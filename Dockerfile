FROM rust:1.66-slim-buster AS builder

WORKDIR /build

COPY ./ .

RUN cargo build --release

FROM debian:buster-slim

# Import from builder
# Copy our build
COPY --from=builder /build/target/release/pnps-utils /bin/

CMD ["/bin/pnps-utils"]
