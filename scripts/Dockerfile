# ---------- Build Stage ----------
FROM rust:1.93.0-bookworm AS builder

WORKDIR /app

COPY Cargo.toml Cargo.lock ./
COPY src ./src
COPY scripts ./scripts
COPY crates ./crates

RUN cargo build --release --all-features --locked --workspace && \
    strip target/release/bamsplit

# ---------- Runtime Stage ----------
FROM debian:bookworm-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    ca-certificates \
    procps \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/bamsplit /usr/local/bin/bamsplit

# Set up non-root user
RUN useradd -m -u 1000 cuser && \
    chmod +x /usr/local/bin/bamsplit

USER cuser
WORKDIR /data

RUN bamsplit --help
