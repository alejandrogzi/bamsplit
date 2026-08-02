// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! A reversible, collision-detecting encoder from logical keys to filenames.
//!
//! Routing keys come from untrusted places: reference names in a BAM header,
//! `Z`-tag values, transcript identifiers in a BED file. Any of them can
//! contain `/`, `..`, a NUL byte, a Windows-reserved device name, or 4 000
//! characters of arbitrary UTF-8. Concatenating one into a path is how a split
//! tool writes outside its output directory.
//!
//! # The encoding
//!
//! Bytes in the **safe set** — `A-Z`, `a-z`, `0-9`, `-`, `_`, `.`, `+`, `~` —
//! pass through. Everything else becomes `%XX` with **upper-case** hex, so the
//! encoding is deterministic and the output is always printable ASCII.
//!
//! ```text
//! chr1                 -> chr1
//! chr1/alternate       -> chr1%2Falternate
//! HLA-A*01:01          -> HLA-A%2A01%3A01
//! ..                   -> %2E%2E
//! (empty)              -> %-
//! ```
//!
//! Three cases get extra treatment *after* the byte-wise pass, because they are
//! dangerous as whole names rather than as individual bytes:
//!
//! * `.` and `..` — path traversal. The first byte is re-encoded, giving `%2E`
//!   and, after the trailing-dot rule below also fires, `%2E%2E`.
//! * A Windows reserved device name (`CON`, `NUL`, `COM3`, …), with or without
//!   an extension. The first byte is re-encoded, giving `%43ON`.
//! * A trailing `.`, which Windows silently strips. The last byte is
//!   re-encoded.
//!
//! Each of those substitutions is itself a legal `%XX` escape, so decoding is
//! unchanged and round-tripping still holds.
//!
//! The empty key encodes to the reserved two-byte string `%-`. A `%` in the
//! encoding of a non-empty key is always followed by two hex digits, so `%-`
//! cannot arise any other way and the mapping stays injective.
//!
//! # Long keys
//!
//! Above [`MAX_STEM_LEN`] bytes the encoding is replaced by
//! `<truncated-prefix>--<16 hex digits of xxh3_64(original bytes)>`. That is
//! **not** reversible, and it is the only case that is not: the manifest always
//! carries the complete logical value, so nothing is lost.
//!
//! # Example
//!
//! ```
//! use bamsplit_core::output::filename::{decode, encode};
//!
//! assert_eq!(encode(b"chr1"), "chr1");
//! assert_eq!(encode(b"chr1/alternate"), "chr1%2Falternate");
//! assert_eq!(decode("chr1%2Falternate").as_deref(), Some(&b"chr1/alternate"[..]));
//! ```

use std::collections::HashMap;
use std::path::{Component, Path, PathBuf};

use crate::error::OutputError;

/// The longest encoded stem `bamsplit` will emit, in bytes.
///
/// Most filesystems cap a path component at 255 bytes. `bamsplit` appends
/// `.bam.bai` (8 bytes) and, while writing, a `.part.<pid>.<random>` suffix of
/// up to about 30 bytes, so the stem budget is set well below the limit.
pub const MAX_STEM_LEN: usize = 200;

/// The encoding of the empty key.
pub const EMPTY_KEY_ENCODING: &str = "%-";

/// Whether a byte survives encoding unchanged.
const fn is_safe(byte: u8) -> bool {
    byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'_' | b'.' | b'+' | b'~')
}

const HEX: &[u8; 16] = b"0123456789ABCDEF";

/// Windows device names that cannot be used as a file name, even with an
/// extension. Compared case-insensitively against the stem before the first
/// `.`.
const WINDOWS_RESERVED: [&str; 22] = [
    "CON", "PRN", "AUX", "NUL", "COM1", "COM2", "COM3", "COM4", "COM5", "COM6", "COM7", "COM8",
    "COM9", "LPT1", "LPT2", "LPT3", "LPT4", "LPT5", "LPT6", "LPT7", "LPT8", "LPT9",
];

/// Encodes a logical key into a filesystem-safe stem.
///
/// The result contains only printable ASCII, never a path separator, and is
/// never `.` or `..`. See the module documentation for the exact rules.
#[must_use]
pub fn encode(key: &[u8]) -> String {
    if key.is_empty() {
        return EMPTY_KEY_ENCODING.to_string();
    }

    let mut out = Vec::with_capacity(key.len());
    for &byte in key {
        if is_safe(byte) {
            out.push(byte);
        } else {
            push_escape(&mut out, byte);
        }
    }

    // `out` is ASCII by construction, so this cannot fail.
    let mut encoded = String::from_utf8(out).unwrap_or_default();

    if needs_leading_escape(&encoded) {
        // The first byte of an all-ASCII string is a whole character.
        let first = encoded.as_bytes()[0];
        let mut escaped = Vec::with_capacity(encoded.len() + 2);
        push_escape(&mut escaped, first);
        escaped.extend_from_slice(&encoded.as_bytes()[1..]);
        encoded = String::from_utf8(escaped).unwrap_or_default();
    }

    if encoded.ends_with('.') {
        let mut escaped = encoded.as_bytes()[..encoded.len() - 1].to_vec();
        push_escape(&mut escaped, b'.');
        encoded = String::from_utf8(escaped).unwrap_or_default();
    }

    if encoded.len() > MAX_STEM_LEN {
        encoded = shorten(&encoded, key);
    }

    encoded
}

fn push_escape(out: &mut Vec<u8>, byte: u8) {
    out.push(b'%');
    out.push(HEX[usize::from(byte >> 4)]);
    out.push(HEX[usize::from(byte & 0xf)]);
}

fn needs_leading_escape(encoded: &str) -> bool {
    if encoded == "." || encoded == ".." {
        return true;
    }
    let stem = encoded.split('.').next().unwrap_or(encoded);
    WINDOWS_RESERVED
        .iter()
        .any(|reserved| stem.eq_ignore_ascii_case(reserved))
}

/// Replaces an over-long encoding with a prefix and a stable digest.
///
/// The prefix is cut at an escape boundary so the retained part still decodes,
/// which keeps the name readable; the digest is over the *original* bytes, so
/// two different long keys cannot collide unless xxh3 does.
fn shorten(encoded: &str, original: &[u8]) -> String {
    const DIGEST_LEN: usize = 16;
    const SEPARATOR: &str = "--";

    let digest = format!("{:016x}", xxhash_rust::xxh3::xxh3_64(original));
    let budget = MAX_STEM_LEN - DIGEST_LEN - SEPARATOR.len();

    let bytes = encoded.as_bytes();
    let mut cut = budget.min(bytes.len());
    // Do not cut inside a `%XX` escape: walk back past a trailing partial one.
    for back in 0..3 {
        if cut < back + 1 {
            break;
        }
        if bytes[cut - back - 1] == b'%' {
            cut -= back + 1;
            break;
        }
    }
    let prefix = encoded.get(..cut).unwrap_or("");
    format!("{prefix}{SEPARATOR}{digest}")
}

/// Decodes an encoded stem back to the original bytes.
///
/// Returns [`None`] for input that this encoder could not have produced: a
/// truncated or non-hex escape. Names shortened by the length cap decode to
/// something that is *not* the original key, which is why the manifest carries
/// the logical value verbatim.
#[must_use]
pub fn decode(encoded: &str) -> Option<Vec<u8>> {
    if encoded == EMPTY_KEY_ENCODING {
        return Some(Vec::new());
    }

    let bytes = encoded.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut index = 0;
    while index < bytes.len() {
        if bytes[index] == b'%' {
            let high = bytes.get(index + 1).copied().and_then(from_hex)?;
            let low = bytes.get(index + 2).copied().and_then(from_hex)?;
            out.push((high << 4) | low);
            index += 3;
        } else {
            out.push(bytes[index]);
            index += 1;
        }
    }
    Some(out)
}

const fn from_hex(byte: u8) -> Option<u8> {
    Some(match byte {
        b'0'..=b'9' => byte - b'0',
        b'A'..=b'F' => byte - b'A' + 10,
        b'a'..=b'f' => byte - b'a' + 10,
        _ => return None,
    })
}

/// Whether an encoded stem is safe to join onto an output directory.
///
/// [`encode`] always produces a safe stem; this exists so a stem that arrives
/// from a filename template can be checked with the same rule.
#[must_use]
pub fn is_safe_component(stem: &str) -> bool {
    if stem.is_empty() || stem == "." || stem == ".." {
        return false;
    }
    if stem.ends_with('.') || stem.ends_with(' ') {
        return false;
    }
    if needs_leading_escape(stem) {
        return false;
    }
    !stem
        .bytes()
        .any(|byte| byte < 0x20 || byte == 0x7f || matches!(byte, b'/' | b'\\'))
        && Path::new(stem).components().count() == 1
        && matches!(
            Path::new(stem).components().next(),
            Some(Component::Normal(_))
        )
}

/// A filename template such as `sample1_{key}`.
///
/// Supported placeholders:
///
/// | placeholder | expands to |
/// | --- | --- |
/// | `{key}` | the encoded routing key |
/// | `{index}` | the output's 0-based ordinal |
///
/// Literal braces are written `{{` and `}}`. A template that renders to
/// something unsafe is rejected at render time, not silently sanitized.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FilenameTemplate {
    parts: Vec<TemplatePart>,
    source: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum TemplatePart {
    Literal(String),
    Key,
    Index,
}

impl Default for FilenameTemplate {
    fn default() -> Self {
        Self {
            parts: vec![TemplatePart::Key],
            source: "{key}".to_string(),
        }
    }
}

impl FilenameTemplate {
    /// Parses a template.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError::UnknownTemplatePlaceholder`] for an unrecognized
    /// `{…}` and [`OutputError::InvalidTemplate`] for an unbalanced brace, a
    /// path separator, a control character, or a template with no `{key}`.
    pub fn parse(template: &str) -> Result<Self, OutputError> {
        let invalid = |reason: &str| OutputError::InvalidTemplate {
            template: template.to_string(),
            reason: reason.to_string(),
        };

        if template.is_empty() {
            return Err(invalid("the template is empty"));
        }
        if template.bytes().any(|byte| matches!(byte, b'/' | b'\\')) {
            return Err(invalid(
                "a template names one file, so it cannot contain `/` or `\\`",
            ));
        }
        if template.bytes().any(|byte| byte < 0x20 || byte == 0x7f) {
            return Err(invalid("the template contains a control character"));
        }

        let mut parts = Vec::new();
        let mut literal = String::new();
        let mut rest = template;

        while let Some(open) = rest.find(['{', '}']) {
            let (head, tail) = rest.split_at(open);
            literal.push_str(head);

            if let Some(escaped) = tail.strip_prefix("{{") {
                literal.push('{');
                rest = escaped;
                continue;
            }
            if let Some(escaped) = tail.strip_prefix("}}") {
                literal.push('}');
                rest = escaped;
                continue;
            }
            if tail.starts_with('}') {
                return Err(invalid("unmatched `}`; write `}}` for a literal brace"));
            }

            let Some(close) = tail.find('}') else {
                return Err(invalid("unmatched `{`; write `{{` for a literal brace"));
            };
            let name = &tail[1..close];
            if !literal.is_empty() {
                parts.push(TemplatePart::Literal(std::mem::take(&mut literal)));
            }
            parts.push(match name {
                "key" => TemplatePart::Key,
                "index" => TemplatePart::Index,
                other => {
                    return Err(OutputError::UnknownTemplatePlaceholder {
                        placeholder: other.to_string(),
                        template: template.to_string(),
                    });
                }
            });
            rest = &tail[close + 1..];
        }
        literal.push_str(rest);
        if !literal.is_empty() {
            parts.push(TemplatePart::Literal(literal));
        }

        if !parts.contains(&TemplatePart::Key) {
            return Err(invalid(
                "the template must contain `{key}`, or every output would share one name",
            ));
        }

        Ok(Self {
            parts,
            source: template.to_string(),
        })
    }

    /// The template as the user wrote it.
    #[must_use]
    pub fn source(&self) -> &str {
        &self.source
    }

    /// Renders the template for one output.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError::InvalidTemplate`] if the rendered stem is not a
    /// safe path component.
    pub fn render(&self, encoded_key: &str, index: usize) -> Result<String, OutputError> {
        let mut out = String::with_capacity(encoded_key.len() + self.source.len());
        for part in &self.parts {
            match part {
                TemplatePart::Literal(text) => out.push_str(text),
                TemplatePart::Key => out.push_str(encoded_key),
                TemplatePart::Index => {
                    use std::fmt::Write as _;
                    let _ = write!(out, "{index}");
                }
            }
        }
        if !is_safe_component(&out) {
            return Err(OutputError::InvalidTemplate {
                template: self.source.clone(),
                reason: format!("rendered to the unsafe filename {out:?}"),
            });
        }
        Ok(out)
    }
}

/// The paths one logical output occupies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OutputPaths {
    /// The logical key, verbatim.
    pub logical_key: Vec<u8>,
    /// The encoded key, before any template is applied.
    pub encoded_key: String,
    /// The final stem, after the template is applied.
    pub stem: String,
    /// The BAM path.
    pub bam: PathBuf,
    /// The index path, when one will be produced.
    pub index: Option<PathBuf>,
}

/// Assigns filenames, refusing to let two logical keys share one.
///
/// A collision is possible in principle — two keys that differ only in bytes a
/// shortened name discards — so it is checked rather than assumed away.
#[derive(Debug)]
pub struct FilenameEncoder {
    directory: PathBuf,
    template: FilenameTemplate,
    assigned: HashMap<String, Vec<u8>>,
    next_index: usize,
}

impl FilenameEncoder {
    /// Creates an encoder writing into `directory`.
    pub fn new(directory: impl Into<PathBuf>, template: FilenameTemplate) -> Self {
        Self {
            directory: directory.into(),
            template,
            assigned: HashMap::new(),
            next_index: 0,
        }
    }

    /// The output directory.
    #[must_use]
    pub fn directory(&self) -> &Path {
        &self.directory
    }

    /// How many names have been assigned.
    #[must_use]
    pub fn len(&self) -> usize {
        self.assigned.len()
    }

    /// Whether no name has been assigned yet.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.assigned.is_empty()
    }

    /// Assigns paths for one logical key.
    ///
    /// Calling this twice with the same key returns the same paths.
    ///
    /// # Errors
    ///
    /// Returns [`OutputError::EncodedKeyCollision`] when a different logical
    /// key already owns the rendered name, and
    /// [`OutputError::InvalidTemplate`] when the template renders to something
    /// unsafe.
    pub fn assign(
        &mut self,
        logical_key: &[u8],
        index_extension: Option<&str>,
    ) -> Result<OutputPaths, OutputError> {
        let encoded_key = encode(logical_key);
        let stem = self.template.render(&encoded_key, self.next_index)?;

        match self.assigned.get(&stem) {
            Some(existing) if existing.as_slice() == logical_key => {}
            Some(existing) => {
                return Err(OutputError::EncodedKeyCollision {
                    first: String::from_utf8_lossy(existing).into_owned(),
                    second: String::from_utf8_lossy(logical_key).into_owned(),
                    encoded: stem,
                });
            }
            None => {
                self.assigned.insert(stem.clone(), logical_key.to_vec());
                self.next_index += 1;
            }
        }

        let bam = self.directory.join(format!("{stem}.bam"));
        let index =
            index_extension.map(|extension| self.directory.join(format!("{stem}.bam.{extension}")));

        Ok(OutputPaths {
            logical_key: logical_key.to_vec(),
            encoded_key,
            stem,
            bam,
            index,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn safe_names_pass_through() {
        for name in [
            "chr1",
            "chrX",
            "scaffold_12",
            "GL000191.1",
            "HLA-A",
            "a+b~c",
        ] {
            assert_eq!(encode(name.as_bytes()), name, "{name}");
        }
    }

    #[test]
    fn path_separators_are_escaped() {
        assert_eq!(encode(b"chr1/alternate"), "chr1%2Falternate");
        assert_eq!(encode(b"a\\b"), "a%5Cb");
        assert_eq!(encode(b"../../etc/passwd"), "..%2F..%2Fetc%2Fpasswd");
    }

    #[test]
    fn dot_names_cannot_traverse() {
        assert_eq!(encode(b"."), "%2E");
        // The leading escape leaves a trailing `.`, which is then escaped too.
        assert_eq!(encode(b".."), "%2E%2E");
        for key in [b".".as_slice(), b"..".as_slice()] {
            let encoded = encode(key);
            assert!(is_safe_component(&encoded), "{encoded}");
            assert_eq!(decode(&encoded).as_deref(), Some(key));
        }
    }

    #[test]
    fn percent_is_escaped_so_the_encoding_is_injective() {
        assert_eq!(encode(b"%2F"), "%252F");
        assert_eq!(decode("%252F").as_deref(), Some(&b"%2F"[..]));
        assert_ne!(encode(b"%2F"), encode(b"/"));
    }

    #[test]
    fn control_characters_and_spaces_are_escaped() {
        assert_eq!(encode(b"a b"), "a%20b");
        assert_eq!(encode(b"a\tb"), "a%09b");
        assert_eq!(encode(b"a\nb"), "a%0Ab");
        assert_eq!(encode(&[0u8, 1, 0x7f]), "%00%01%7F");
    }

    #[test]
    fn windows_reserved_characters_and_names_are_escaped() {
        assert_eq!(encode(b"a<b>c:d\"e|f?g*h"), "a%3Cb%3Ec%3Ad%22e%7Cf%3Fg%2Ah");
        assert_eq!(encode(b"CON"), "%43ON");
        assert_eq!(encode(b"nul"), "%6Eul");
        assert_eq!(encode(b"COM3.bam"), "%43OM3.bam");
        // A name that merely starts with a reserved word is fine.
        assert_eq!(encode(b"CONTIG"), "CONTIG");
        for key in [&b"CON"[..], b"nul", b"COM3.bam"] {
            assert_eq!(decode(&encode(key)).as_deref(), Some(key));
        }
    }

    #[test]
    fn a_trailing_dot_is_escaped() {
        assert_eq!(encode(b"chr1."), "chr1%2E");
        assert_eq!(decode("chr1%2E").as_deref(), Some(&b"chr1."[..]));
    }

    #[test]
    fn non_ascii_bytes_are_escaped() {
        // "chré" in UTF-8.
        assert_eq!(encode("chré".as_bytes()), "chr%C3%A9");
        assert_eq!(decode("chr%C3%A9").as_deref(), Some("chré".as_bytes()));
    }

    #[test]
    fn the_empty_key_has_a_reserved_encoding() {
        assert_eq!(encode(b""), EMPTY_KEY_ENCODING);
        assert_eq!(decode(EMPTY_KEY_ENCODING).as_deref(), Some(&b""[..]));
        assert!(is_safe_component(EMPTY_KEY_ENCODING));
        // No non-empty key can produce it.
        assert_ne!(encode(b"-"), EMPTY_KEY_ENCODING);
    }

    #[test]
    fn long_keys_are_shortened_with_a_stable_digest() {
        let key = vec![b'x'; 4096];
        let encoded = encode(&key);
        assert!(encoded.len() <= MAX_STEM_LEN, "{}", encoded.len());
        assert_eq!(encoded, encode(&key), "the digest must be stable");
        assert!(encoded.contains("--"));
        assert!(is_safe_component(&encoded));

        let mut other = key.clone();
        other.push(b'y');
        assert_ne!(encoded, encode(&other));
    }

    #[test]
    fn shortening_does_not_cut_inside_an_escape() {
        let key: Vec<u8> = std::iter::repeat_n(b'/', 4096).collect();
        let encoded = encode(&key);
        assert!(encoded.len() <= MAX_STEM_LEN);
        let (prefix, _) = encoded.rsplit_once("--").expect("has a digest");
        assert!(
            decode(prefix).is_some(),
            "prefix {prefix:?} must still decode"
        );
    }

    #[test]
    fn decoding_rejects_malformed_escapes() {
        assert!(decode("%").is_none());
        assert!(decode("%2").is_none());
        assert!(decode("%zz").is_none());
        assert!(decode("abc%").is_none());
    }

    #[test]
    fn encoded_names_are_always_single_safe_components() {
        for key in [
            &b"chr1"[..],
            b"chr1/alt",
            b"..",
            b".",
            b"",
            b"CON",
            b"a b",
            &[0u8, 0xff],
            "染色体".as_bytes(),
        ] {
            let encoded = encode(key);
            assert!(is_safe_component(&encoded), "{key:?} -> {encoded:?}");
            let joined = Path::new("/out").join(&encoded);
            assert!(
                joined.starts_with("/out") && joined.components().count() == 3,
                "{joined:?} escaped the output directory"
            );
        }
    }

    #[test]
    fn round_trips_every_single_byte() {
        for byte in 0..=255u8 {
            let key = [byte];
            let encoded = encode(&key);
            assert_eq!(decode(&encoded).as_deref(), Some(&key[..]), "byte {byte}");
        }
    }

    #[test]
    fn round_trips_multi_byte_keys() {
        let keys: Vec<Vec<u8>> = vec![
            b"chr1".to_vec(),
            b"chr1/alt:1-2".to_vec(),
            b"%%%".to_vec(),
            (0..=255u8).collect(),
            b"...".to_vec(),
        ];
        for key in keys {
            let encoded = encode(&key);
            if encoded.len() >= MAX_STEM_LEN {
                // Shortened names are documented as lossy; the manifest keeps
                // the logical value instead.
                continue;
            }
            assert_eq!(decode(&encoded).as_deref(), Some(&key[..]), "{key:?}");
        }
    }

    #[test]
    fn templates_render_and_validate() {
        let template = FilenameTemplate::parse("sample1_{key}").expect("valid");
        assert_eq!(template.render("chr1", 0).expect("safe"), "sample1_chr1");

        let indexed = FilenameTemplate::parse("{index}-{key}").expect("valid");
        assert_eq!(indexed.render("chrX", 7).expect("safe"), "7-chrX");

        let braces = FilenameTemplate::parse("{{{key}}}").expect("valid");
        assert_eq!(braces.render("chr1", 0).expect("safe"), "{chr1}");
    }

    #[test]
    fn templates_reject_dangerous_and_malformed_input() {
        assert!(matches!(
            FilenameTemplate::parse("../{key}"),
            Err(OutputError::InvalidTemplate { .. })
        ));
        assert!(matches!(
            FilenameTemplate::parse("{sample}"),
            Err(OutputError::UnknownTemplatePlaceholder { .. })
        ));
        assert!(matches!(
            FilenameTemplate::parse("prefix-only"),
            Err(OutputError::InvalidTemplate { .. })
        ));
        assert!(matches!(
            FilenameTemplate::parse("{key"),
            Err(OutputError::InvalidTemplate { .. })
        ));
        assert!(matches!(
            FilenameTemplate::parse("key}"),
            Err(OutputError::InvalidTemplate { .. })
        ));
        assert!(matches!(
            FilenameTemplate::parse(""),
            Err(OutputError::InvalidTemplate { .. })
        ));
    }

    #[test]
    fn a_template_that_renders_to_a_reserved_name_is_rejected() {
        let template = FilenameTemplate::parse("{key}").expect("valid");
        // `encode` would never produce this, but a template could.
        let error = template.render("CON", 0).expect_err("must reject");
        assert!(
            matches!(error, OutputError::InvalidTemplate { .. }),
            "{error}"
        );
    }

    #[test]
    fn the_encoder_assigns_stable_paths() {
        let mut encoder = FilenameEncoder::new("/out", FilenameTemplate::default());
        let first = encoder.assign(b"chr1", Some("bai")).expect("assigned");
        assert_eq!(first.bam, Path::new("/out/chr1.bam"));
        assert_eq!(first.index.as_deref(), Some(Path::new("/out/chr1.bam.bai")));
        assert_eq!(first.encoded_key, "chr1");

        let again = encoder.assign(b"chr1", Some("bai")).expect("assigned");
        assert_eq!(first, again);
        assert_eq!(encoder.len(), 1);
    }

    #[test]
    fn the_encoder_detects_a_collision() {
        // Two keys whose encodings collide only because a template discards the
        // difference; simulated here with a template that ignores the index.
        let mut encoder = FilenameEncoder::new("/out", FilenameTemplate::default());
        encoder.assign(b"chr1", None).expect("assigned");
        // Force a collision by re-registering the same stem for a different key.
        let clash = encoder.assign(b"chr1", None);
        assert!(clash.is_ok(), "the same key must be idempotent");

        let mut manual = FilenameEncoder::new("/out", FilenameTemplate::default());
        manual
            .assigned
            .insert("chr1".to_string(), b"other".to_vec());
        let error = manual.assign(b"chr1", None).expect_err("must reject");
        assert!(
            matches!(error, OutputError::EncodedKeyCollision { .. }),
            "{error}"
        );
    }

    #[test]
    fn the_index_extension_is_optional() {
        let mut encoder = FilenameEncoder::new("/out", FilenameTemplate::default());
        let paths = encoder.assign(b"chr1", None).expect("assigned");
        assert!(paths.index.is_none());

        let csi = encoder.assign(b"chr2", Some("csi")).expect("assigned");
        assert_eq!(csi.index.as_deref(), Some(Path::new("/out/chr2.bam.csi")));
    }

    #[test]
    fn unsafe_components_are_recognized() {
        assert!(is_safe_component("chr1"));
        assert!(!is_safe_component(""));
        assert!(!is_safe_component("."));
        assert!(!is_safe_component(".."));
        assert!(!is_safe_component("a/b"));
        assert!(!is_safe_component("a\\b"));
        assert!(!is_safe_component("a\0b"));
        assert!(!is_safe_component("trailing."));
        assert!(!is_safe_component("trailing "));
        assert!(!is_safe_component("CON"));
    }
}
