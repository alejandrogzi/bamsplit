// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! A BGZF block writer that compresses in parallel **and** reports exact
//! virtual offsets.
//!
//! # Why not use the one in `noodles`
//!
//! `noodles-bgzf` ships two writers and neither does both jobs:
//!
//! | writer | parallel | `virtual_position()` |
//! | --- | --- | --- |
//! | `bgzf::io::Writer` | no | yes |
//! | `bgzf::io::MultithreadedWriter` | yes | **no** |
//!
//! Streaming index construction needs a virtual offset per record, and the
//! flagship `bamsplit chrom` path wants compression spread across the thread
//! budget. Choosing either writer would give up one of the two.
//!
//! # How offsets survive parallelism
//!
//! A record's virtual offset is `(compressed_block_offset << 16) | offset
//! within the block`. The second half is known the instant the record is
//! staged; the first half is only known once every *preceding* block has been
//! compressed, because a block's compressed size is not predictable.
//!
//! So this writer splits the two:
//!
//! ```text
//! write_record()                       flush_batch()
//!   ├─ append to the staging block       ├─ compress N blocks in parallel
//!   ├─ note (block_index, offset)        ├─ write them in order
//!   └─ return a RecordSlot               ├─ now every block's base is known
//!                                        └─ resolve slots -> VirtualRange
//! ```
//!
//! Slots are resolved in write order and handed to the caller through a
//! callback, so an index consumer still sees strictly increasing offsets. The
//! resolution buffer holds at most one batch of records, which bounds memory
//! at roughly `batch_size * 64 KiB` regardless of file size.
//!
//! With one worker the batch is one block and resolution happens immediately,
//! so the single-threaded path pays nothing for the machinery.
//!
//! # Frame layout
//!
//! ```text
//! 1f 8b 08 04                gzip magic, DEFLATE, FEXTRA
//! 00 00 00 00                MTIME = 0, so output is byte-reproducible
//! 00 ff                      XFL, OS = unknown
//! 06 00                      XLEN = 6
//! 42 43 02 00 <BSIZE u16>    the BGZF extra field; BSIZE = block length - 1
//! <deflate data>
//! <CRC32 u32> <ISIZE u32>
//! ```

use std::io::Write;
use std::sync::Arc;

use flate2::Compression;
use flate2::write::DeflateEncoder;

/// The most uncompressed bytes packed into one block.
///
/// `htslib` uses `0xff00`; matching it keeps block boundaries familiar to
/// anyone comparing output with `bgzip`, and leaves ample headroom under the
/// 65 536-byte hard cap for the gzip framing and worst-case DEFLATE expansion.
pub const MAX_BLOCK_UNCOMPRESSED: usize = 0xff00;

/// The 28-byte empty block that marks a well-formed end of file.
pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

const HEADER_LEN: usize = 18;
const TRAILER_LEN: usize = 8;

/// A BGZF virtual offset: `compressed << 16 | uncompressed`.
///
/// Kept as a plain `u64` because that is what both the index builder and the
/// manifest want, and because it makes the arithmetic in this module obvious.
pub type VirtualOffset = u64;

/// Packs a virtual offset, saturating rather than wrapping.
///
/// A compressed offset above 2^48 means a single output has exceeded 256 TiB,
/// which no real run reaches; saturating keeps the function total instead of
/// introducing a panic on a path that can never be exercised.
#[must_use]
pub const fn pack_offset(compressed: u64, uncompressed: u16) -> VirtualOffset {
    let compressed = if compressed > (1 << 48) - 1 {
        (1 << 48) - 1
    } else {
        compressed
    };
    (compressed << 16) | (uncompressed as u64)
}

/// The `[start, end)` virtual-offset range one record occupies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VirtualRange {
    /// Where the record's `block_size` prefix begins.
    pub start: VirtualOffset,
    /// Where the next record begins.
    pub end: VirtualOffset,
}

/// A position inside the not-yet-compressed stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Slot {
    block_index: u64,
    offset_in_block: u32,
}

/// One staged, uncompressed block awaiting compression.
#[derive(Debug)]
struct StagedBlock {
    index: u64,
    data: Vec<u8>,
}

/// A compressed block, ready to write.
#[derive(Debug)]
struct CompressedBlock {
    index: u64,
    frame: Vec<u8>,
    uncompressed_len: usize,
}

/// A record whose virtual range cannot be computed until its batch is written.
///
/// The entry is queued *before* the record's bytes are appended, because
/// appending can seal a block and sealing can compact the block-offset table.
/// If the entry were queued afterwards, a record that straddled a block
/// boundary would find the offset of its own starting block already discarded.
/// `complete` marks the window in which `end` is not yet known; resolution
/// stops at an incomplete entry and picks it up on the next flush.
#[derive(Debug)]
struct PendingRecord<T> {
    start: Slot,
    end: Slot,
    payload: T,
    complete: bool,
}

/// Writes BGZF blocks, resolving virtual offsets as batches complete.
///
/// `T` is whatever the caller wants handed back with each resolved record —
/// the index context, in practice.
pub struct BgzfBlockWriter<W: Write, T> {
    inner: W,
    level: Compression,
    pool: Option<Arc<rayon::ThreadPool>>,
    batch_size: usize,

    staging: Vec<u8>,
    staged: Vec<StagedBlock>,
    spare: Vec<Vec<u8>>,
    pending: Vec<PendingRecord<T>>,

    next_block_index: u64,
    /// Compressed start offsets of the blocks written since
    /// `first_unresolved_index`, in block order.
    ///
    /// The block one past the end starts at `compressed_written`, which is what
    /// makes a record that ends exactly on a block boundary resolvable as soon
    /// as the block before it is written.
    block_offsets: Vec<u64>,
    first_unresolved_index: u64,
    /// Total compressed bytes written to `inner`, excluding the EOF marker.
    compressed_written: u64,
    finished: bool,
}

impl<W: Write, T> std::fmt::Debug for BgzfBlockWriter<W, T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("BgzfBlockWriter")
            .field("batch_size", &self.batch_size)
            .field("staged_blocks", &self.staged.len())
            .field("pending_records", &self.pending.len())
            .field("compressed_written", &self.compressed_written)
            .field("finished", &self.finished)
            .finish_non_exhaustive()
    }
}

impl<W: Write, T> BgzfBlockWriter<W, T> {
    /// Creates a writer.
    ///
    /// `base_compressed_offset` is added to every reported offset, which lets a
    /// parked-and-reopened output continue an index without a discontinuity.
    /// `worker_count` of 1 keeps compression on the calling thread.
    pub fn new(
        inner: W,
        compression_level: u32,
        pool: Option<Arc<rayon::ThreadPool>>,
        worker_count: usize,
        base_compressed_offset: u64,
    ) -> Self {
        let batch_size = worker_count.max(1);
        Self {
            inner,
            level: Compression::new(compression_level.min(9)),
            pool: if batch_size > 1 { pool } else { None },
            batch_size,
            staging: Vec::with_capacity(MAX_BLOCK_UNCOMPRESSED),
            staged: Vec::with_capacity(batch_size),
            spare: Vec::new(),
            pending: Vec::new(),
            next_block_index: 0,
            block_offsets: Vec::new(),
            first_unresolved_index: 0,
            compressed_written: base_compressed_offset,
            finished: false,
        }
    }

    /// Total compressed bytes written, including the base offset.
    #[must_use]
    pub const fn compressed_position(&self) -> u64 {
        self.compressed_written
    }

    /// Borrows the sink.
    pub const fn get_ref(&self) -> &W {
        &self.inner
    }

    fn current_slot(&self) -> Slot {
        Slot {
            block_index: self.next_block_index,
            // Bounded by `MAX_BLOCK_UNCOMPRESSED`, which fits `u32`.
            offset_in_block: self.staging.len() as u32,
        }
    }

    /// Appends `bytes` to the stream, splitting across blocks as needed.
    ///
    /// `resolved` has to be threaded all the way down here: appending can seal
    /// a block, sealing can flush a batch, and flushing is exactly when record
    /// offsets become known. Passing a no-op at any level would silently
    /// discard those resolutions.
    fn append<F>(&mut self, mut bytes: &[u8], resolved: &mut F) -> std::io::Result<()>
    where
        F: FnMut(T, VirtualRange),
    {
        while !bytes.is_empty() {
            let room = MAX_BLOCK_UNCOMPRESSED - self.staging.len();
            let take = room.min(bytes.len());
            self.staging.extend_from_slice(&bytes[..take]);
            bytes = &bytes[take..];
            if self.staging.len() == MAX_BLOCK_UNCOMPRESSED {
                self.seal_block(resolved)?;
            }
        }
        Ok(())
    }

    /// Moves the staging buffer into the batch, flushing if the batch is full.
    fn seal_block<F>(&mut self, resolved: &mut F) -> std::io::Result<()>
    where
        F: FnMut(T, VirtualRange),
    {
        if self.staging.is_empty() {
            return Ok(());
        }
        let mut data = self.spare.pop().unwrap_or_default();
        data.clear();
        std::mem::swap(&mut data, &mut self.staging);
        self.staging.reserve(MAX_BLOCK_UNCOMPRESSED);
        self.staged.push(StagedBlock {
            index: self.next_block_index,
            data,
        });
        self.next_block_index += 1;
        if self.staged.len() >= self.batch_size {
            self.flush_batch(resolved)?;
        }
        Ok(())
    }

    /// Writes one length-prefixed BAM record.
    ///
    /// `payload` is handed back through `resolved` once the record's virtual
    /// range is known, which may be during this call or during a later one.
    ///
    /// # Errors
    ///
    /// Returns any I/O error from the sink.
    pub fn write_record<F>(
        &mut self,
        body: &[u8],
        payload: T,
        resolved: &mut F,
    ) -> std::io::Result<()>
    where
        F: FnMut(T, VirtualRange),
    {
        // `block_size` counts the body only; it is a u32 on disk.
        let prefix = u32::try_from(body.len()).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidInput,
                format!("BAM record body of {} bytes exceeds u32", body.len()),
            )
        })?;

        let start = self.current_slot();
        let reserved = self.pending.len();
        self.pending.push(PendingRecord {
            start,
            end: start,
            payload,
            complete: false,
        });

        self.append(&prefix.to_le_bytes(), resolved)?;
        self.append(body, resolved)?;

        // The reservation is the last entry unless a flush drained earlier
        // ones, which only ever shortens the queue from the front.
        let end = self.current_slot();
        let index = reserved.min(self.pending.len().saturating_sub(1));
        if let Some(record) = self.pending.get_mut(index) {
            record.end = end;
            record.complete = true;
        }
        Ok(())
    }

    /// Writes raw bytes that are not a record, such as the BAM header.
    ///
    /// # Errors
    ///
    /// Returns any I/O error from the sink.
    pub fn write_raw<F>(&mut self, bytes: &[u8], resolved: &mut F) -> std::io::Result<()>
    where
        F: FnMut(T, VirtualRange),
    {
        self.append(bytes, resolved)
    }

    /// Compresses and writes every staged block, in order.
    fn flush_batch<F>(&mut self, resolved: &mut F) -> std::io::Result<()>
    where
        F: FnMut(T, VirtualRange),
    {
        if self.staged.is_empty() {
            return Ok(());
        }

        let level = self.level;
        let staged = std::mem::take(&mut self.staged);

        let mut compressed: Vec<CompressedBlock> = match (&self.pool, staged.len()) {
            (Some(pool), count) if count > 1 => {
                use rayon::prelude::*;
                pool.install(|| {
                    staged
                        .par_iter()
                        .map(|block| CompressedBlock {
                            index: block.index,
                            uncompressed_len: block.data.len(),
                            frame: encode_block(&block.data, level),
                        })
                        .collect()
                })
            }
            _ => staged
                .iter()
                .map(|block| CompressedBlock {
                    index: block.index,
                    uncompressed_len: block.data.len(),
                    frame: encode_block(&block.data, level),
                })
                .collect(),
        };

        // Recycle the uncompressed buffers.
        for block in staged {
            self.spare.push(block.data);
        }

        compressed.sort_unstable_by_key(|block| block.index);
        for block in &compressed {
            debug_assert_eq!(
                block.index,
                self.first_unresolved_index + self.block_offsets.len() as u64,
                "blocks must be recorded in index order"
            );
            self.block_offsets.push(self.compressed_written);
            self.inner.write_all(&block.frame)?;
            self.compressed_written += block.frame.len() as u64;
            let _ = block.uncompressed_len;
        }

        self.try_resolve(resolved);
        Ok(())
    }

    /// Emits every pending record whose block offsets are now known.
    ///
    /// Records are staged in write order and block indices never decrease, so
    /// the resolvable records are always a prefix of `pending`. That turns
    /// resolution into one scan plus one `drain`, rather than a filter that
    /// reallocates the whole list.
    fn try_resolve<F>(&mut self, resolved: &mut F)
    where
        F: FnMut(T, VirtualRange),
    {
        let resolvable = self
            .pending
            .iter()
            .take_while(|record| {
                record.complete
                    && self.base_of(record.start).is_some()
                    && self.base_of(record.end).is_some()
            })
            .count();

        for record in self.pending.drain(..resolvable) {
            // `base_of` returned `Some` for both slots in the scan above, and
            // nothing between then and now changes the block table.
            let start = pack_offset(
                base_or_zero(
                    &self.block_offsets,
                    self.first_unresolved_index,
                    self.compressed_written,
                    record.start,
                ),
                record.start.offset_in_block as u16,
            );
            let end = pack_offset(
                base_or_zero(
                    &self.block_offsets,
                    self.first_unresolved_index,
                    self.compressed_written,
                    record.end,
                ),
                record.end.offset_in_block as u16,
            );
            resolved(record.payload, VirtualRange { start, end });
        }
        self.compact_offsets();
    }

    /// The compressed start offset of `slot`'s block, if it is known.
    ///
    /// The block one past the last written one is known too: it begins exactly
    /// where the written bytes end.
    #[allow(clippy::unnecessary_wraps)]
    fn base_of(&self, slot: Slot) -> Option<u64> {
        let position =
            usize::try_from(slot.block_index.checked_sub(self.first_unresolved_index)?).ok()?;
        match position.cmp(&self.block_offsets.len()) {
            std::cmp::Ordering::Less => self.block_offsets.get(position).copied(),
            std::cmp::Ordering::Equal => Some(self.compressed_written),
            std::cmp::Ordering::Greater => None,
        }
    }

    /// Drops block offsets no pending record can still reference.
    fn compact_offsets(&mut self) {
        let lowest = self.pending.first().map_or(
            self.first_unresolved_index + self.block_offsets.len() as u64,
            |record| record.start.block_index.min(record.end.block_index),
        );
        let drop_count = usize::try_from(lowest.saturating_sub(self.first_unresolved_index))
            .unwrap_or(0)
            .min(self.block_offsets.len());
        if drop_count > 0 {
            self.block_offsets.drain(..drop_count);
            self.first_unresolved_index += drop_count as u64;
        }
    }

    /// Flushes everything staged, without writing the end-of-file marker.
    ///
    /// After this returns, [`compressed_position`](Self::compressed_position)
    /// is a block boundary, which is what makes parking and reopening an output
    /// safe.
    ///
    /// # Errors
    ///
    /// Returns any I/O error from the sink.
    pub fn flush_blocks<F>(&mut self, resolved: &mut F) -> std::io::Result<()>
    where
        F: FnMut(T, VirtualRange),
    {
        self.seal_block(resolved)?;
        self.flush_batch(resolved)?;
        self.inner.flush()
    }

    /// Flushes, writes the end-of-file marker, and returns the sink.
    ///
    /// # Errors
    ///
    /// Returns any I/O error from the sink, or an error if a record's virtual
    /// range could not be resolved — which would indicate a bug in this module
    /// rather than bad input.
    pub fn finish<F>(mut self, resolved: &mut F) -> std::io::Result<W>
    where
        F: FnMut(T, VirtualRange),
    {
        self.flush_blocks(resolved)?;
        if !self.pending.is_empty() {
            return Err(std::io::Error::other(format!(
                "{} BGZF record offsets were never resolved",
                self.pending.len()
            )));
        }
        self.inner.write_all(&BGZF_EOF)?;
        self.compressed_written += BGZF_EOF.len() as u64;
        self.inner.flush()?;
        self.finished = true;
        Ok(self.inner)
    }

    /// Flushes and returns the sink **without** the end-of-file marker.
    ///
    /// Used when an output is parked to stay inside `--max-open-files`; the
    /// marker is written when the output is finally closed.
    ///
    /// # Errors
    ///
    /// Returns any I/O error from the sink.
    pub fn park<F>(mut self, resolved: &mut F) -> std::io::Result<(W, u64)>
    where
        F: FnMut(T, VirtualRange),
    {
        self.flush_blocks(resolved)?;
        if !self.pending.is_empty() {
            return Err(std::io::Error::other(format!(
                "{} BGZF record offsets were never resolved before parking",
                self.pending.len()
            )));
        }
        let position = self.compressed_written;
        self.finished = true;
        Ok((self.inner, position))
    }
}

/// The compressed start offset of `slot`'s block.
///
/// A free function so the borrow checker allows it to be called while
/// `pending` is being drained. Falls back to zero for a slot whose block is not
/// in the table, which the caller has already ruled out via
/// [`BgzfBlockWriter::base_of`].
fn base_or_zero(
    block_offsets: &[u64],
    first_unresolved_index: u64,
    compressed_written: u64,
    slot: Slot,
) -> u64 {
    let Some(position) = slot
        .block_index
        .checked_sub(first_unresolved_index)
        .and_then(|position| usize::try_from(position).ok())
    else {
        return 0;
    };
    match position.cmp(&block_offsets.len()) {
        std::cmp::Ordering::Less => block_offsets[position],
        std::cmp::Ordering::Equal => compressed_written,
        std::cmp::Ordering::Greater => 0,
    }
}

/// Compresses one block into a complete BGZF frame.
fn encode_block(data: &[u8], level: Compression) -> Vec<u8> {
    let mut frame = Vec::with_capacity(HEADER_LEN + data.len() / 2 + TRAILER_LEN + 64);
    frame.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, // magic, DEFLATE, FEXTRA
        0x00, 0x00, 0x00, 0x00, // MTIME = 0 keeps output byte-reproducible
        0x00, 0xff, // XFL, OS = unknown
        0x06, 0x00, // XLEN
        0x42, 0x43, 0x02, 0x00, // "BC", SLEN = 2
        0x00, 0x00, // BSIZE placeholder
    ]);

    let mut encoder = DeflateEncoder::new(frame, level);
    // Writing into a `Vec` cannot fail, and `finish` only surfaces sink errors.
    let _ = encoder.write_all(data);
    let Ok(mut frame) = encoder.finish() else {
        return encode_stored_block(data);
    };

    frame.extend_from_slice(&crc32fast::hash(data).to_le_bytes());
    frame.extend_from_slice(&(data.len() as u32).to_le_bytes());

    // A block that grew past the 65 536-byte frame limit is re-emitted with no
    // compression, which is bounded by `MAX_BLOCK_UNCOMPRESSED + 5 * ceil(n/65535)`
    // and therefore always fits.
    if frame.len() > 0x1_0000 {
        return encode_stored_block(data);
    }

    let bsize = (frame.len() - 1) as u16;
    frame[16..18].copy_from_slice(&bsize.to_le_bytes());
    frame
}

/// Emits a block with DEFLATE's "stored" mode, the worst-case fallback.
fn encode_stored_block(data: &[u8]) -> Vec<u8> {
    let mut frame = Vec::with_capacity(HEADER_LEN + data.len() + 16 + TRAILER_LEN);
    frame.extend_from_slice(&[
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
        0x00, 0x00, 0x00,
    ]);

    let mut remaining = data;
    loop {
        let take = remaining.len().min(0xffff);
        let (chunk, rest) = remaining.split_at(take);
        let final_block = u8::from(rest.is_empty());
        frame.push(final_block);
        frame.extend_from_slice(&(take as u16).to_le_bytes());
        frame.extend_from_slice(&(!(take as u16)).to_le_bytes());
        frame.extend_from_slice(chunk);
        remaining = rest;
        if remaining.is_empty() {
            break;
        }
    }

    frame.extend_from_slice(&crc32fast::hash(data).to_le_bytes());
    frame.extend_from_slice(&(data.len() as u32).to_le_bytes());
    let bsize = (frame.len() - 1) as u16;
    frame[16..18].copy_from_slice(&bsize.to_le_bytes());
    frame
}

#[cfg(test)]
#[allow(clippy::items_after_statements)]
mod tests {
    use super::*;
    use std::io::Read;

    fn inflate(bgzf: &[u8]) -> Vec<u8> {
        let mut reader = noodles_bgzf::io::Reader::new(bgzf);
        let mut out = Vec::new();
        reader.read_to_end(&mut out).expect("valid BGZF");
        out
    }

    fn record(body_len: usize, fill: u8) -> Vec<u8> {
        vec![fill; body_len]
    }

    #[test]
    fn writes_a_stream_noodles_can_read() {
        let mut resolved = Vec::new();
        let mut writer: BgzfBlockWriter<Vec<u8>, usize> =
            BgzfBlockWriter::new(Vec::new(), 6, None, 1, 0);
        writer
            .write_raw(b"BAM\x01", &mut |_, _| {})
            .expect("written");
        for (index, fill) in [b'a', b'b', b'c'].into_iter().enumerate() {
            writer
                .write_record(&record(100, fill), index, &mut |payload, range| {
                    resolved.push((payload, range));
                })
                .expect("written");
        }
        let bytes = writer
            .finish(&mut |payload, range| resolved.push((payload, range)))
            .expect("finished");

        let inflated = inflate(&bytes);
        assert_eq!(&inflated[..4], b"BAM\x01");
        assert_eq!(inflated.len(), 4 + 3 * (4 + 100));
        assert_eq!(resolved.len(), 3);
        assert_eq!(bytes[bytes.len() - 28..], BGZF_EOF);
    }

    #[test]
    fn record_bodies_survive_byte_for_byte() {
        let bodies: Vec<Vec<u8>> = (0..500)
            .map(|index| {
                let length = 32 + (index * 7) % 4096;
                (0..length).map(|byte| (byte % 251) as u8).collect()
            })
            .collect();

        let mut writer: BgzfBlockWriter<Vec<u8>, ()> =
            BgzfBlockWriter::new(Vec::new(), 1, None, 1, 0);
        for body in &bodies {
            writer
                .write_record(body, (), &mut |(), _| {})
                .expect("written");
        }
        let bytes = writer.finish(&mut |(), _| {}).expect("finished");

        let inflated = inflate(&bytes);
        let mut cursor = 0usize;
        for body in &bodies {
            let length =
                u32::from_le_bytes(inflated[cursor..cursor + 4].try_into().expect("four bytes"))
                    as usize;
            assert_eq!(length, body.len());
            cursor += 4;
            assert_eq!(&inflated[cursor..cursor + length], &body[..]);
            cursor += length;
        }
        assert_eq!(cursor, inflated.len());
    }

    #[test]
    fn virtual_offsets_can_be_seeked_to() {
        // A record that begins exactly on a block boundary has two equally
        // valid virtual offsets: `(block, block_len)` and `(block + 1, 0)`.
        // This writer emits the second, which is also what `htslib`'s writer
        // reports after a flush; `noodles`' reader lazily reports the first
        // until it refills. Comparing the two numerically would therefore be
        // testing a representation, not a behaviour. What must hold is that
        // seeking to a reported offset lands on that record — so that is what
        // is asserted.
        use noodles_bgzf::VirtualPosition;

        let bodies: Vec<Vec<u8>> = (0..2_000)
            .map(|index| record(40 + index % 900, (index % 256) as u8))
            .collect();

        let mut resolved: Vec<(usize, VirtualRange)> = Vec::new();
        let mut writer: BgzfBlockWriter<Vec<u8>, usize> =
            BgzfBlockWriter::new(Vec::new(), 6, None, 1, 0);
        for (index, body) in bodies.iter().enumerate() {
            writer
                .write_record(body, index, &mut |payload, range| {
                    resolved.push((payload, range));
                })
                .expect("written");
        }
        let bytes = writer
            .finish(&mut |payload, range| resolved.push((payload, range)))
            .expect("finished");

        assert_eq!(resolved.len(), bodies.len());
        for (position, (index, _)) in resolved.iter().enumerate() {
            assert_eq!(*index, position, "records must resolve in write order");
        }
        // Offsets must be strictly increasing and each record's end must be its
        // successor's start.
        for window in resolved.windows(2) {
            assert!(window[0].1.start < window[0].1.end, "empty range");
            assert!(window[0].1.end <= window[1].1.start, "ranges overlap");
        }

        let mut reader = noodles_bgzf::io::Reader::new(std::io::Cursor::new(bytes));
        let mut prefix = [0u8; 4];
        for (index, body) in bodies.iter().enumerate() {
            reader
                .seek(VirtualPosition::from(resolved[index].1.start))
                .expect("seekable");
            reader.read_exact(&mut prefix).expect("prefix");
            let length = u32::from_le_bytes(prefix) as usize;
            assert_eq!(length, body.len(), "record {index}");
            let mut payload = vec![0u8; length];
            reader.read_exact(&mut payload).expect("body");
            assert_eq!(payload, *body, "record {index}");
        }
    }

    #[test]
    fn parallel_and_serial_writers_agree_on_offsets() {
        let bodies: Vec<Vec<u8>> = (0..3_000)
            .map(|index| record(50 + index % 700, (index % 97) as u8))
            .collect();

        let run = |workers: usize| {
            let pool = Arc::new(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(workers)
                    .build()
                    .expect("pool"),
            );
            let mut resolved: Vec<(usize, VirtualRange)> = Vec::new();
            let mut writer: BgzfBlockWriter<Vec<u8>, usize> =
                BgzfBlockWriter::new(Vec::new(), 6, Some(pool), workers, 0);
            for (index, body) in bodies.iter().enumerate() {
                writer
                    .write_record(body, index, &mut |payload, range| {
                        resolved.push((payload, range));
                    })
                    .expect("written");
            }
            let bytes = writer
                .finish(&mut |payload, range| resolved.push((payload, range)))
                .expect("finished");
            (bytes, resolved)
        };

        let (serial_bytes, serial_offsets) = run(1);
        let (parallel_bytes, parallel_offsets) = run(4);

        // Identical block boundaries and identical compression settings mean
        // the two must be byte-identical, which also pins reproducibility.
        assert_eq!(serial_bytes, parallel_bytes);
        assert_eq!(serial_offsets, parallel_offsets);
        assert_eq!(inflate(&serial_bytes), inflate(&parallel_bytes));
    }

    #[test]
    fn a_record_larger_than_one_block_spans_blocks_correctly() {
        let body = record(MAX_BLOCK_UNCOMPRESSED * 3 + 17, 0x5a);
        let mut resolved = Vec::new();
        let mut writer: BgzfBlockWriter<Vec<u8>, ()> =
            BgzfBlockWriter::new(Vec::new(), 6, None, 1, 0);
        writer
            .write_record(&body, (), &mut |(), range| resolved.push(range))
            .expect("written");
        let bytes = writer
            .finish(&mut |(), range| resolved.push(range))
            .expect("finished");

        assert_eq!(resolved.len(), 1);
        let inflated = inflate(&bytes);
        assert_eq!(inflated.len(), 4 + body.len());
        assert_eq!(&inflated[4..], &body[..]);

        let mut reader = noodles_bgzf::io::Reader::new(&bytes[..]);
        assert_eq!(u64::from(reader.virtual_position()), resolved[0].start);
        let mut sink = Vec::new();
        reader.read_to_end(&mut sink).expect("readable");
    }

    #[test]
    fn incompressible_data_still_fits_a_block() {
        // A deterministic pseudo-random stream that DEFLATE cannot shrink.
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let data: Vec<u8> = (0..MAX_BLOCK_UNCOMPRESSED)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state >> 24) as u8
            })
            .collect();

        let frame = encode_block(&data, Compression::new(9));
        assert!(frame.len() <= 0x1_0000, "frame is {} bytes", frame.len());
        let bsize = u16::from_le_bytes([frame[16], frame[17]]) as usize;
        assert_eq!(bsize + 1, frame.len());

        let mut with_eof = frame.clone();
        with_eof.extend_from_slice(&BGZF_EOF);
        assert_eq!(inflate(&with_eof), data);
    }

    #[test]
    fn the_stored_fallback_round_trips() {
        let data: Vec<u8> = (0..1000).map(|index| (index % 256) as u8).collect();
        let mut frame = encode_stored_block(&data);
        frame.extend_from_slice(&BGZF_EOF);
        assert_eq!(inflate(&frame), data);
    }

    #[test]
    fn an_empty_stream_is_just_the_eof_marker() {
        let writer: BgzfBlockWriter<Vec<u8>, ()> = BgzfBlockWriter::new(Vec::new(), 6, None, 1, 0);
        let bytes = writer.finish(&mut |(), _| {}).expect("finished");
        assert_eq!(bytes, BGZF_EOF.to_vec());
        assert!(inflate(&bytes).is_empty());
    }

    #[test]
    fn parking_omits_the_eof_marker_and_reports_a_block_boundary() {
        let mut writer: BgzfBlockWriter<Vec<u8>, ()> =
            BgzfBlockWriter::new(Vec::new(), 6, None, 1, 0);
        writer
            .write_record(&record(100, b'x'), (), &mut |(), _| {})
            .expect("written");
        let (first_half, position) = writer.park(&mut |(), _| {}).expect("parked");
        assert_eq!(position as usize, first_half.len());
        assert_ne!(
            &first_half[first_half.len().saturating_sub(28)..],
            &BGZF_EOF[..],
            "a parked stream must not be terminated"
        );

        // Resume with the recorded base offset and finish.
        let mut resumed: BgzfBlockWriter<Vec<u8>, ()> =
            BgzfBlockWriter::new(Vec::new(), 6, None, 1, position);
        let mut offsets = Vec::new();
        resumed
            .write_record(&record(100, b'y'), (), &mut |(), range| offsets.push(range))
            .expect("written");
        let second_half = resumed
            .finish(&mut |(), range| offsets.push(range))
            .expect("finished");

        assert_eq!(offsets.len(), 1);
        assert_eq!(offsets[0].start >> 16, position, "offsets must continue");

        let mut whole = first_half;
        whole.extend_from_slice(&second_half);
        let inflated = inflate(&whole);
        assert_eq!(inflated.len(), 2 * (4 + 100));
    }

    #[test]
    fn offsets_are_packed_and_saturate() {
        assert_eq!(pack_offset(0, 0), 0);
        assert_eq!(pack_offset(1, 2), (1 << 16) | 2);
        assert_eq!(pack_offset(u64::MAX, 0) >> 16, (1 << 48) - 1);
    }

    #[test]
    fn the_resolution_buffer_stays_bounded() {
        let mut writer: BgzfBlockWriter<Vec<u8>, ()> =
            BgzfBlockWriter::new(Vec::new(), 1, None, 4, 0);
        let mut peak = 0usize;
        for _ in 0..20_000 {
            writer
                .write_record(&record(64, b'z'), (), &mut |(), _| {})
                .expect("written");
            peak = peak.max(writer.pending.len());
        }
        let _ = writer.finish(&mut |(), _| {}).expect("finished");
        // One batch of 64-byte records is at most 4 * 65280 / 68 ≈ 3 840.
        assert!(peak < 8_000, "pending buffer peaked at {peak}");
    }
}
