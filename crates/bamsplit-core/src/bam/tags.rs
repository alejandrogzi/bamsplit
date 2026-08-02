// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! Bounds-checked decoding of the BAM auxiliary (tag) section.
//!
//! The BAM auxiliary section is a flat sequence of fields:
//!
//! ```text
//! tag[2]  type[1]  value...
//! ```
//!
//! where `type` is one of `A c C s S i I f Z H B`. `B` arrays add a subtype
//! byte and a little-endian `u32` element count.
//!
//! Everything here borrows from the record body: no field is copied, and the
//! decoder performs one linear scan with explicit bounds checks. There is no
//! unchecked indexing and no unaligned pointer dereference — multi-byte scalars
//! are rebuilt from byte slices with [`from_le_bytes`](i32::from_le_bytes).
//!
//! # Example
//!
//! ```
//! use bamsplit_core::bam::tags::{RawTagValue, TagReader};
//!
//! // NM:i:3 followed by RG:Z:rg0
//! let data = b"NMi\x03\x00\x00\x00RGZrg0\x00";
//! let mut reader = TagReader::new(data);
//!
//! let (tag, value) = reader.next().unwrap()?;
//! assert_eq!(&tag, b"NM");
//! assert_eq!(value.as_integer(), Some(3));
//!
//! let (tag, value) = reader.next().unwrap()?;
//! assert_eq!(&tag, b"RG");
//! assert!(matches!(value, RawTagValue::String(b"rg0")));
//! # Ok::<_, bamsplit_core::error::TagError>(())
//! ```

use std::fmt::Write as _;

use crate::error::TagError;

/// A two-character BAM auxiliary tag, e.g. `*b"RG"`.
pub type Tag = [u8; 2];

/// Renders a tag for diagnostics, escaping bytes that are not printable ASCII.
#[must_use]
pub fn render_tag(tag: Tag) -> String {
    let mut rendered = String::with_capacity(2);
    for byte in tag {
        if byte.is_ascii_graphic() {
            rendered.push(byte as char);
        } else {
            let _ = write!(rendered, "\\x{byte:02x}");
        }
    }
    rendered
}

/// The subtype of a `B` (array) auxiliary value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ArraySubtype {
    /// `c`: signed 8-bit.
    Int8,
    /// `C`: unsigned 8-bit.
    UInt8,
    /// `s`: signed 16-bit.
    Int16,
    /// `S`: unsigned 16-bit.
    UInt16,
    /// `i`: signed 32-bit.
    Int32,
    /// `I`: unsigned 32-bit.
    UInt32,
    /// `f`: IEEE-754 binary32.
    Float,
}

impl ArraySubtype {
    /// Decodes a subtype byte.
    fn from_byte(byte: u8) -> Option<Self> {
        Some(match byte {
            b'c' => Self::Int8,
            b'C' => Self::UInt8,
            b's' => Self::Int16,
            b'S' => Self::UInt16,
            b'i' => Self::Int32,
            b'I' => Self::UInt32,
            b'f' => Self::Float,
            _ => return None,
        })
    }

    /// The on-disk width of one element, in bytes.
    #[must_use]
    pub const fn element_size(self) -> usize {
        match self {
            Self::Int8 | Self::UInt8 => 1,
            Self::Int16 | Self::UInt16 => 2,
            Self::Int32 | Self::UInt32 | Self::Float => 4,
        }
    }

    /// The SAM type character for this subtype.
    #[must_use]
    pub const fn type_code(self) -> u8 {
        match self {
            Self::Int8 => b'c',
            Self::UInt8 => b'C',
            Self::Int16 => b's',
            Self::UInt16 => b'S',
            Self::Int32 => b'i',
            Self::UInt32 => b'I',
            Self::Float => b'f',
        }
    }
}

/// A borrowed `B` array value.
///
/// The payload is validated (subtype known, `count * element_size` does not
/// overflow, and the whole payload is in bounds) at construction time, so
/// iteration cannot fail.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RawTagArray<'a> {
    subtype: ArraySubtype,
    count: u32,
    payload: &'a [u8],
}

impl<'a> RawTagArray<'a> {
    /// The element subtype.
    #[must_use]
    pub const fn subtype(&self) -> ArraySubtype {
        self.subtype
    }

    /// The declared element count.
    #[must_use]
    pub const fn len(&self) -> u32 {
        self.count
    }

    /// Whether the array declares zero elements.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.count == 0
    }

    /// The raw little-endian payload, excluding the subtype and count.
    #[must_use]
    pub const fn payload(&self) -> &'a [u8] {
        self.payload
    }

    /// Iterates the elements as `f64`, widening integers losslessly except for
    /// `u32`/`i32` values beyond 2^53, which cannot occur in practice for the
    /// arrays BAM stores but are documented here for completeness.
    pub fn iter_f64(&self) -> impl Iterator<Item = f64> + 'a {
        let subtype = self.subtype;
        let size = subtype.element_size();
        self.payload
            .chunks_exact(size)
            .map(move |chunk| match subtype {
                // Every arm indexes a chunk whose length is exactly `size` by
                // construction of `chunks_exact`, so the slicing cannot panic.
                ArraySubtype::Int8 => f64::from(chunk[0] as i8),
                ArraySubtype::UInt8 => f64::from(chunk[0]),
                ArraySubtype::Int16 => f64::from(i16::from_le_bytes([chunk[0], chunk[1]])),
                ArraySubtype::UInt16 => f64::from(u16::from_le_bytes([chunk[0], chunk[1]])),
                ArraySubtype::Int32 => {
                    f64::from(i32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
                }
                ArraySubtype::UInt32 => {
                    f64::from(u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
                }
                ArraySubtype::Float => {
                    f64::from(f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
                }
            })
    }

    /// Iterates the elements as `i64`, returning [`None`] for a float array.
    #[must_use]
    pub fn iter_i64(&self) -> Option<impl Iterator<Item = i64> + 'a> {
        if self.subtype == ArraySubtype::Float {
            return None;
        }
        let subtype = self.subtype;
        let size = subtype.element_size();
        Some(
            self.payload
                .chunks_exact(size)
                .map(move |chunk| match subtype {
                    ArraySubtype::Int8 => i64::from(chunk[0] as i8),
                    ArraySubtype::UInt8 => i64::from(chunk[0]),
                    ArraySubtype::Int16 => i64::from(i16::from_le_bytes([chunk[0], chunk[1]])),
                    ArraySubtype::UInt16 => i64::from(u16::from_le_bytes([chunk[0], chunk[1]])),
                    ArraySubtype::Int32 => {
                        i64::from(i32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
                    }
                    ArraySubtype::UInt32 => {
                        i64::from(u32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]))
                    }
                    // Unreachable: guarded by the `Float` early return above.
                    ArraySubtype::Float => 0,
                }),
        )
    }
}

/// A borrowed BAM auxiliary value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RawTagValue<'a> {
    /// `A`: a single printable character.
    Character(u8),
    /// `c`: signed 8-bit integer.
    Int8(i8),
    /// `C`: unsigned 8-bit integer.
    UInt8(u8),
    /// `s`: signed 16-bit integer.
    Int16(i16),
    /// `S`: unsigned 16-bit integer.
    UInt16(u16),
    /// `i`: signed 32-bit integer.
    Int32(i32),
    /// `I`: unsigned 32-bit integer.
    UInt32(u32),
    /// `f`: IEEE-754 binary32.
    Float(f32),
    /// `Z`: NUL-terminated string, without the terminator.
    String(&'a [u8]),
    /// `H`: NUL-terminated hex string, without the terminator.
    Hex(&'a [u8]),
    /// `B`: typed array.
    Array(RawTagArray<'a>),
}

impl<'a> RawTagValue<'a> {
    /// The SAM type character (`A`, `c`, …, `B`).
    #[must_use]
    pub const fn type_code(&self) -> u8 {
        match self {
            Self::Character(_) => b'A',
            Self::Int8(_) => b'c',
            Self::UInt8(_) => b'C',
            Self::Int16(_) => b's',
            Self::UInt16(_) => b'S',
            Self::Int32(_) => b'i',
            Self::UInt32(_) => b'I',
            Self::Float(_) => b'f',
            Self::String(_) => b'Z',
            Self::Hex(_) => b'H',
            Self::Array(_) => b'B',
        }
    }

    /// Whether this is a `B` array.
    #[must_use]
    pub const fn is_array(&self) -> bool {
        matches!(self, Self::Array(_))
    }

    /// The value as an `i64`, for the integral types only.
    #[must_use]
    pub const fn as_integer(&self) -> Option<i64> {
        Some(match *self {
            Self::Int8(value) => value as i64,
            Self::UInt8(value) => value as i64,
            Self::Int16(value) => value as i64,
            Self::UInt16(value) => value as i64,
            Self::Int32(value) => value as i64,
            Self::UInt32(value) => value as i64,
            _ => return None,
        })
    }

    /// The value as a byte string, for `Z` and `H` only.
    #[must_use]
    pub const fn as_bytes(&self) -> Option<&'a [u8]> {
        match *self {
            Self::String(value) | Self::Hex(value) => Some(value),
            _ => None,
        }
    }

    /// The array payload, for `B` only.
    #[must_use]
    pub const fn as_array(&self) -> Option<RawTagArray<'a>> {
        match *self {
            Self::Array(array) => Some(array),
            _ => None,
        }
    }

    /// Renders the value as a deterministic, locale-independent routing key.
    ///
    /// The encoding is stable across platforms and Rust versions for every
    /// type except `f` and `B:f`, where it is stable for a given Rust release
    /// (Rust's float `Display` emits the shortest string that round-trips).
    ///
    /// | type | encoding |
    /// | --- | --- |
    /// | `A` | the character itself |
    /// | `c C s S i I` | decimal, `-` for negatives, no separators |
    /// | `f` | shortest round-tripping decimal; `NaN`, `inf`, `-inf` |
    /// | `Z` | the raw bytes, verbatim |
    /// | `H` | the raw hex digits, verbatim |
    /// | `B` | `<subtype>,<v0>,<v1>,…` — only reachable via `--allow-array-tags` |
    ///
    /// Floating-point values never go through locale-sensitive formatting:
    /// Rust's `Display` for `f32`/`f64` is locale-independent by definition.
    #[must_use]
    pub fn to_routing_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        self.write_routing_bytes(&mut out);
        out
    }

    /// Appends [`to_routing_bytes`](Self::to_routing_bytes) to `out`.
    pub fn write_routing_bytes(&self, out: &mut Vec<u8>) {
        match *self {
            Self::Character(value) => out.push(value),
            Self::Int8(value) => append_i64(out, i64::from(value)),
            Self::UInt8(value) => append_i64(out, i64::from(value)),
            Self::Int16(value) => append_i64(out, i64::from(value)),
            Self::UInt16(value) => append_i64(out, i64::from(value)),
            Self::Int32(value) => append_i64(out, i64::from(value)),
            Self::UInt32(value) => append_i64(out, i64::from(value)),
            Self::Float(value) => append_f32(out, value),
            Self::String(value) | Self::Hex(value) => out.extend_from_slice(value),
            Self::Array(array) => {
                out.push(array.subtype().type_code());
                if array.subtype() == ArraySubtype::Float {
                    for element in array.iter_f64() {
                        out.push(b',');
                        append_f64(out, element);
                    }
                } else if let Some(elements) = array.iter_i64() {
                    for element in elements {
                        out.push(b',');
                        append_i64(out, element);
                    }
                }
            }
        }
    }
}

fn append_i64(out: &mut Vec<u8>, value: i64) {
    let mut buffer = itoa_buffer();
    out.extend_from_slice(format_i64(&mut buffer, value));
}

/// A stack buffer wide enough for `-9223372036854775808`.
const fn itoa_buffer() -> [u8; 20] {
    [0u8; 20]
}

/// Formats `value` into `buffer` without allocating, returning the used slice.
fn format_i64(buffer: &mut [u8; 20], value: i64) -> &[u8] {
    let negative = value < 0;
    // `i64::MIN.unsigned_abs()` is exactly `2^63`, which `u64` represents, so
    // this cannot overflow.
    let mut magnitude = value.unsigned_abs();
    let mut cursor = buffer.len();
    loop {
        cursor -= 1;
        buffer[cursor] = b'0' + u8::try_from(magnitude % 10).unwrap_or(0);
        magnitude /= 10;
        if magnitude == 0 {
            break;
        }
    }
    if negative {
        cursor -= 1;
        buffer[cursor] = b'-';
    }
    &buffer[cursor..]
}

fn append_f32(out: &mut Vec<u8>, value: f32) {
    if value.is_nan() {
        out.extend_from_slice(b"NaN");
    } else if value.is_infinite() {
        out.extend_from_slice(if value.is_sign_negative() {
            b"-inf"
        } else {
            b"inf"
        });
    } else {
        // `Display for f32` is locale-independent and emits the shortest
        // decimal that round-trips.
        let _ = write!(ByteWriter(out), "{value}");
    }
}

fn append_f64(out: &mut Vec<u8>, value: f64) {
    if value.is_nan() {
        out.extend_from_slice(b"NaN");
    } else if value.is_infinite() {
        out.extend_from_slice(if value.is_sign_negative() {
            b"-inf"
        } else {
            b"inf"
        });
    } else {
        let _ = write!(ByteWriter(out), "{value}");
    }
}

/// Adapts a `Vec<u8>` to [`std::fmt::Write`] for UTF-8-only formatting.
struct ByteWriter<'a>(&'a mut Vec<u8>);

impl std::fmt::Write for ByteWriter<'_> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        self.0.extend_from_slice(s.as_bytes());
        Ok(())
    }
}

/// A forward-only, bounds-checked scanner over a BAM auxiliary section.
///
/// The scanner never allocates and never panics: every read is preceded by an
/// explicit length check, and every arithmetic step that could overflow uses
/// [`checked_add`](usize::checked_add) or [`checked_mul`](usize::checked_mul).
#[derive(Debug, Clone)]
pub struct TagReader<'a> {
    data: &'a [u8],
    offset: usize,
}

impl<'a> TagReader<'a> {
    /// Creates a scanner over an auxiliary section.
    #[must_use]
    pub const fn new(data: &'a [u8]) -> Self {
        Self { data, offset: 0 }
    }

    /// The byte offset of the next field, relative to the start of the section.
    #[must_use]
    pub const fn offset(&self) -> usize {
        self.offset
    }

    /// Whether the whole section has been consumed.
    #[must_use]
    pub const fn is_exhausted(&self) -> bool {
        self.offset >= self.data.len()
    }

    /// Decodes the next field.
    ///
    /// Returns [`None`] once the section is exhausted.
    ///
    /// # Errors
    ///
    /// Returns [`TagError`] when the section is truncated, a type code is
    /// unknown, an array length overflows, or a string is unterminated.
    #[allow(clippy::should_implement_trait)]
    pub fn next(&mut self) -> Option<Result<(Tag, RawTagValue<'a>), TagError>> {
        if self.is_exhausted() {
            return None;
        }
        Some(self.read_field())
    }

    fn read_field(&mut self) -> Result<(Tag, RawTagValue<'a>), TagError> {
        let start = self.offset;
        let remaining = self.data.len() - start;
        if remaining < 3 {
            self.offset = self.data.len();
            return Err(TagError::TruncatedHeader {
                offset: start,
                available: remaining,
            });
        }

        let tag: Tag = [self.data[start], self.data[start + 1]];
        let type_code = self.data[start + 2];
        let value_start = start + 3;
        let value = self.read_value(tag, type_code, value_start)?;
        Ok((tag, value))
    }

    fn read_value(
        &mut self,
        tag: Tag,
        type_code: u8,
        value_start: usize,
    ) -> Result<RawTagValue<'a>, TagError> {
        let available = self.data.len() - value_start;

        let mut fixed = |width: usize| -> Result<&'a [u8], TagError> {
            if available < width {
                self.offset = self.data.len();
                return Err(TagError::Truncated {
                    tag: render_tag(tag),
                    needed: width,
                    available,
                });
            }
            let slice = &self.data[value_start..value_start + width];
            self.offset = value_start + width;
            Ok(slice)
        };

        let value = match type_code {
            b'A' => RawTagValue::Character(fixed(1)?[0]),
            b'c' => RawTagValue::Int8(fixed(1)?[0] as i8),
            b'C' => RawTagValue::UInt8(fixed(1)?[0]),
            b's' => {
                let bytes = fixed(2)?;
                RawTagValue::Int16(i16::from_le_bytes([bytes[0], bytes[1]]))
            }
            b'S' => {
                let bytes = fixed(2)?;
                RawTagValue::UInt16(u16::from_le_bytes([bytes[0], bytes[1]]))
            }
            b'i' => {
                let bytes = fixed(4)?;
                RawTagValue::Int32(i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
            }
            b'I' => {
                let bytes = fixed(4)?;
                RawTagValue::UInt32(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
            }
            b'f' => {
                let bytes = fixed(4)?;
                RawTagValue::Float(f32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]))
            }
            b'Z' | b'H' => {
                let rest = &self.data[value_start..];
                let Some(nul) = rest.iter().position(|&byte| byte == 0) else {
                    self.offset = self.data.len();
                    return Err(TagError::UnterminatedString {
                        tag: render_tag(tag),
                        type_code: char::from(type_code),
                    });
                };
                let payload = &rest[..nul];
                self.offset = value_start + nul + 1;
                if type_code == b'H' {
                    validate_hex(tag, payload)?;
                    RawTagValue::Hex(payload)
                } else {
                    RawTagValue::String(payload)
                }
            }
            b'B' => {
                if available < 5 {
                    self.offset = self.data.len();
                    return Err(TagError::Truncated {
                        tag: render_tag(tag),
                        needed: 5,
                        available,
                    });
                }
                let subtype_byte = self.data[value_start];
                let Some(subtype) = ArraySubtype::from_byte(subtype_byte) else {
                    self.offset = self.data.len();
                    return Err(TagError::UnknownArraySubtype {
                        tag: render_tag(tag),
                        subtype: char::from(subtype_byte),
                        byte: subtype_byte,
                    });
                };
                let count = u32::from_le_bytes([
                    self.data[value_start + 1],
                    self.data[value_start + 2],
                    self.data[value_start + 3],
                    self.data[value_start + 4],
                ]);
                let payload_start = value_start + 5;
                let payload_available = self.data.len() - payload_start;
                let element_size = subtype.element_size();
                let Some(payload_len) = (count as usize).checked_mul(element_size) else {
                    self.offset = self.data.len();
                    return Err(TagError::ArrayLengthOverflow {
                        tag: render_tag(tag),
                        count,
                        element_size,
                        available: payload_available,
                    });
                };
                if payload_len > payload_available {
                    self.offset = self.data.len();
                    return Err(TagError::ArrayLengthOverflow {
                        tag: render_tag(tag),
                        count,
                        element_size,
                        available: payload_available,
                    });
                }
                self.offset = payload_start + payload_len;
                RawTagValue::Array(RawTagArray {
                    subtype,
                    count,
                    payload: &self.data[payload_start..payload_start + payload_len],
                })
            }
            other => {
                self.offset = self.data.len();
                return Err(TagError::UnknownValueType {
                    tag: render_tag(tag),
                    type_code: char::from(other),
                    byte: other,
                });
            }
        };

        Ok(value)
    }
}

fn validate_hex(tag: Tag, payload: &[u8]) -> Result<(), TagError> {
    if payload.len() % 2 != 0 {
        return Err(TagError::MalformedHex {
            tag: render_tag(tag),
            reason: "an `H` value must contain an even number of hex digits",
        });
    }
    if let Some(bad) = payload.iter().find(|byte| !byte.is_ascii_hexdigit()) {
        let _ = bad;
        return Err(TagError::MalformedHex {
            tag: render_tag(tag),
            reason: "an `H` value must contain only the characters [0-9A-Fa-f]",
        });
    }
    Ok(())
}

/// Finds one tag in an auxiliary section, stopping at the first match.
///
/// This is the hot path for `--key tag:XX` and `--tag XX` routing: it decodes
/// only as far as it must and never allocates.
///
/// # Errors
///
/// Returns [`TagError`] if the section is malformed before the requested tag
/// is reached.
pub fn find_tag(data: &[u8], wanted: Tag) -> Result<Option<RawTagValue<'_>>, TagError> {
    let mut reader = TagReader::new(data);
    while let Some(field) = reader.next() {
        let (tag, value) = field?;
        if tag == wanted {
            return Ok(Some(value));
        }
    }
    Ok(None)
}

/// Validates an entire auxiliary section, optionally rejecting duplicate tags.
///
/// Used by `bamsplit inspect --full` and by the malformed-input test suite.
///
/// # Errors
///
/// Returns [`TagError`] for any structural problem, and for a repeated tag when
/// `reject_duplicates` is set.
pub fn validate_section(data: &[u8], reject_duplicates: bool) -> Result<usize, TagError> {
    let mut reader = TagReader::new(data);
    let mut seen: Vec<Tag> = Vec::new();
    let mut count = 0usize;
    while let Some(field) = reader.next() {
        let (tag, _) = field?;
        count += 1;
        if reject_duplicates {
            if seen.contains(&tag) {
                return Err(TagError::DuplicateTag {
                    tag: render_tag(tag),
                });
            }
            seen.push(tag);
        }
    }
    Ok(count)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn collect(data: &[u8]) -> Result<Vec<(Tag, RawTagValue<'_>)>, TagError> {
        let mut reader = TagReader::new(data);
        let mut out = Vec::new();
        while let Some(field) = reader.next() {
            out.push(field?);
        }
        Ok(out)
    }

    #[test]
    fn decodes_every_scalar_type() {
        let mut data = Vec::new();
        // Each field is `tag[2] type[1] value…`.
        data.extend_from_slice(b"AaA*");
        data.extend_from_slice(b"Ccc");
        data.push(0xff); // -1
        data.extend_from_slice(b"CuC");
        data.push(0xff); // 255
        data.extend_from_slice(b"Sis");
        data.extend_from_slice(&(-2i16).to_le_bytes());
        data.extend_from_slice(b"SuS");
        data.extend_from_slice(&65535u16.to_le_bytes());
        data.extend_from_slice(b"Iii");
        data.extend_from_slice(&(-3i32).to_le_bytes());
        data.extend_from_slice(b"IuI");
        data.extend_from_slice(&4_000_000_000u32.to_le_bytes());
        data.extend_from_slice(b"Flf");
        data.extend_from_slice(&0.5f32.to_le_bytes());
        data.extend_from_slice(b"ZzZhello\x00");
        data.extend_from_slice(b"HxHDEADBEEF\x00");

        let fields = collect(&data).expect("all scalar types decode");
        assert_eq!(fields.len(), 10);
        assert_eq!(fields[0].1, RawTagValue::Character(b'*'));
        assert_eq!(fields[1].1, RawTagValue::Int8(-1));
        assert_eq!(fields[2].1, RawTagValue::UInt8(255));
        assert_eq!(fields[3].1, RawTagValue::Int16(-2));
        assert_eq!(fields[4].1, RawTagValue::UInt16(65535));
        assert_eq!(fields[5].1, RawTagValue::Int32(-3));
        assert_eq!(fields[6].1, RawTagValue::UInt32(4_000_000_000));
        assert_eq!(fields[7].1, RawTagValue::Float(0.5));
        assert_eq!(fields[8].1, RawTagValue::String(b"hello"));
        assert_eq!(fields[9].1, RawTagValue::Hex(b"DEADBEEF"));
    }

    #[test]
    fn decodes_arrays_of_every_subtype() {
        for (subtype, element_size) in [
            (b'c', 1usize),
            (b'C', 1),
            (b's', 2),
            (b'S', 2),
            (b'i', 4),
            (b'I', 4),
            (b'f', 4),
        ] {
            let mut data = Vec::new();
            data.extend_from_slice(b"XXB");
            data.push(subtype);
            data.extend_from_slice(&3u32.to_le_bytes());
            data.extend(std::iter::repeat_n(0u8, 3 * element_size));

            let fields = collect(&data).expect("array decodes");
            let array = fields[0].1.as_array().expect("is an array");
            assert_eq!(array.len(), 3);
            assert_eq!(array.subtype().element_size(), element_size);
            assert_eq!(array.iter_f64().count(), 3);
        }
    }

    #[test]
    fn rejects_array_element_count_overflow() {
        let mut data = Vec::new();
        data.extend_from_slice(b"XXB");
        data.push(b'i');
        data.extend_from_slice(&u32::MAX.to_le_bytes());
        data.extend_from_slice(&[0u8; 4]);

        let error = collect(&data).expect_err("must reject");
        assert!(
            matches!(error, TagError::ArrayLengthOverflow { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_array_payload_out_of_bounds() {
        let mut data = Vec::new();
        data.extend_from_slice(b"XXB");
        data.push(b'S');
        data.extend_from_slice(&10u32.to_le_bytes());
        data.extend_from_slice(&[0u8; 4]); // only 2 of 10 elements present

        let error = collect(&data).expect_err("must reject");
        assert!(
            matches!(error, TagError::ArrayLengthOverflow { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_unterminated_string() {
        let error = collect(b"XXZno-nul").expect_err("must reject");
        assert!(
            matches!(error, TagError::UnterminatedString { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_unknown_type_code() {
        let error = collect(b"XXq\x00").expect_err("must reject");
        assert!(
            matches!(error, TagError::UnknownValueType { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_unknown_array_subtype() {
        let mut data = Vec::new();
        data.extend_from_slice(b"XXB");
        data.push(b'q');
        data.extend_from_slice(&0u32.to_le_bytes());
        let error = collect(&data).expect_err("must reject");
        assert!(
            matches!(error, TagError::UnknownArraySubtype { .. }),
            "{error}"
        );
    }

    #[test]
    fn rejects_truncated_field_header() {
        let error = collect(b"XX").expect_err("must reject");
        assert!(matches!(error, TagError::TruncatedHeader { .. }), "{error}");
    }

    #[test]
    fn rejects_truncated_scalar() {
        let error = collect(b"XXi\x01\x02").expect_err("must reject");
        assert!(matches!(error, TagError::Truncated { .. }), "{error}");
    }

    #[test]
    fn rejects_odd_length_hex() {
        let error = collect(b"XXHABC\x00").expect_err("must reject");
        assert!(matches!(error, TagError::MalformedHex { .. }), "{error}");
    }

    #[test]
    fn rejects_non_hex_digits() {
        let error = collect(b"XXHZZ\x00").expect_err("must reject");
        assert!(matches!(error, TagError::MalformedHex { .. }), "{error}");
    }

    #[test]
    fn finds_a_tag_without_decoding_the_rest() {
        let data = b"NMi\x03\x00\x00\x00RGZrg0\x00";
        let value = find_tag(data, *b"RG").expect("valid").expect("present");
        assert_eq!(value, RawTagValue::String(b"rg0"));
        assert!(find_tag(data, *b"ZZ").expect("valid").is_none());
    }

    #[test]
    fn detects_duplicate_tags_only_when_asked() {
        let data = b"NMi\x01\x00\x00\x00NMi\x02\x00\x00\x00";
        assert_eq!(validate_section(data, false).expect("valid"), 2);
        let error = validate_section(data, true).expect_err("must reject");
        assert!(matches!(error, TagError::DuplicateTag { .. }), "{error}");
    }

    #[test]
    fn routing_bytes_are_deterministic_and_locale_free() {
        assert_eq!(RawTagValue::Int32(-42).to_routing_bytes(), b"-42");
        assert_eq!(
            RawTagValue::UInt32(4_000_000_000).to_routing_bytes(),
            b"4000000000"
        );
        assert_eq!(RawTagValue::Character(b'+').to_routing_bytes(), b"+");
        assert_eq!(RawTagValue::String(b"a b").to_routing_bytes(), b"a b");
        assert_eq!(RawTagValue::Float(0.5).to_routing_bytes(), b"0.5");
        assert_eq!(RawTagValue::Float(1.0).to_routing_bytes(), b"1");
        assert_eq!(RawTagValue::Float(f32::NAN).to_routing_bytes(), b"NaN");
        assert_eq!(RawTagValue::Float(f32::INFINITY).to_routing_bytes(), b"inf");
        assert_eq!(
            RawTagValue::Float(f32::NEG_INFINITY).to_routing_bytes(),
            b"-inf"
        );
    }

    #[test]
    fn routing_bytes_for_i64_extremes() {
        assert_eq!(
            RawTagValue::Int32(i32::MIN).to_routing_bytes(),
            b"-2147483648"
        );
        let mut buffer = itoa_buffer();
        assert_eq!(format_i64(&mut buffer, i64::MIN), b"-9223372036854775808");
        assert_eq!(format_i64(&mut buffer, 0), b"0");
    }

    #[test]
    fn array_routing_bytes_are_canonical() {
        let mut data = Vec::new();
        data.extend_from_slice(b"XXB");
        data.push(b'i');
        data.extend_from_slice(&2u32.to_le_bytes());
        data.extend_from_slice(&1i32.to_le_bytes());
        data.extend_from_slice(&(-2i32).to_le_bytes());

        let value = find_tag(&data, *b"XX").expect("valid").expect("present");
        assert_eq!(value.to_routing_bytes(), b"i,1,-2");
    }

    #[test]
    fn empty_section_yields_no_fields() {
        assert_eq!(collect(b"").expect("valid").len(), 0);
        assert_eq!(validate_section(b"", true).expect("valid"), 0);
    }

    #[test]
    fn render_tag_escapes_non_printable_bytes() {
        assert_eq!(render_tag(*b"RG"), "RG");
        assert_eq!(render_tag([0x00, b'G']), "\\x00G");
    }
}
