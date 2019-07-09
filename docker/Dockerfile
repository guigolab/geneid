# Dockerfile for geneid
#
FROM debian:stretch AS builder

ARG GENEID_VER=1.4.5

RUN apt-get update && \
    apt-get install --no-install-recommends -y \
        curl \
        ca-certificates \
        build-essential

RUN curl -L https://github.com/guigolab/geneid/archive/v${GENEID_VER}.tar.gz | \
    tar xz && \
    cd geneid-${GENEID_VER} && \
    make BIN=/build

FROM debian:stretch

LABEL author.name="Emilio Palumbo"
LABEL author.email="emiliopalumbo@gmail.com"

COPY --from=builder /build/geneid /usr/local/bin/

ENTRYPOINT [ "/usr/local/bin/geneid" ]