// Copyright (c) 2025 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
// Distributed under the terms of the Apache License, Version 2.0.

//! `bamsplit shard`.

use bamsplit_core::bam::header::ReadGroupField;
use bamsplit_core::error::ConfigError;
use bamsplit_core::routing::shard::{DEFAULT_SEED, ShardKeySource};
use bamsplit_core::{Result, ShardOptions};

use crate::cli::{GlobalOptions, ShardArgs};

/// Parses a `--key` value.
///
/// `tag:XX` is handled here rather than as a `clap` enum because its payload is
/// part of the value.
///
/// # Errors
///
/// Returns [`ConfigError`] for an unrecognized key or a malformed tag.
pub fn parse_key(value: &str) -> std::result::Result<ShardKeySource, ConfigError> {
    if let Some(tag) = value.strip_prefix("tag:") {
        return Ok(ShardKeySource::Tag(bamsplit_core::routing::parse_tag(tag)?));
    }
    Ok(match value {
        "qname" => ShardKeySource::QName,
        "record" => ShardKeySource::Record,
        "read-group" => ShardKeySource::ReadGroupField(ReadGroupField::ReadGroup),
        "sample" => ShardKeySource::ReadGroupField(ReadGroupField::Sample),
        "library" => ShardKeySource::ReadGroupField(ReadGroupField::Library),
        "platform-unit" => ShardKeySource::ReadGroupField(ReadGroupField::PlatformUnit),
        other => {
            return Err(ConfigError::OutOfRange {
                option: "--key",
                constraint: "one of qname, record, tag:XX, read-group, sample, library, \
                             platform-unit",
                value: other.to_string(),
            });
        }
    })
}

/// Runs a shard split.
///
/// # Errors
///
/// Propagates every [`bamsplit_core::Error`].
pub fn run(globals: &GlobalOptions, args: ShardArgs) -> Result<()> {
    let run = globals.to_run_options(crate::cli::command_line(), args.index, None);
    let options = ShardOptions {
        out_dir: args.out_dir,
        shards: args.shards,
        key: parse_key(&args.key)?,
        seed: args.seed.unwrap_or(DEFAULT_SEED),
        missing: args.missing.into(),
        missing_name: args.missing_name,
        shard_width: args.shard_width,
    };
    let report = bamsplit_core::split_by_shard(&args.input, &options, &run)?;
    super::report_summary(&report);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_documented_key_parses() {
        assert_eq!(parse_key("qname").expect("valid"), ShardKeySource::QName);
        assert_eq!(parse_key("record").expect("valid"), ShardKeySource::Record);
        assert_eq!(
            parse_key("tag:CB").expect("valid"),
            ShardKeySource::Tag(*b"CB")
        );
        assert_eq!(
            parse_key("platform-unit").expect("valid"),
            ShardKeySource::ReadGroupField(ReadGroupField::PlatformUnit)
        );
    }

    #[test]
    fn an_unknown_key_is_rejected_with_the_valid_set() {
        let error = parse_key("barcode").expect_err("must reject");
        assert!(error.to_string().contains("qname"), "{error}");
    }

    #[test]
    fn a_malformed_tag_key_is_rejected() {
        assert!(parse_key("tag:TOOLONG").is_err());
        assert!(parse_key("tag:").is_err());
    }
}
