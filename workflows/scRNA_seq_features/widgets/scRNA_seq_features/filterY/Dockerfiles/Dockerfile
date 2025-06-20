ARG samVersion=1.18
FROM debian:bookworm-slim as builder
ARG DEBIAN_FRONTEND=noninteractive
ARG TARGETPLATFORM
ARG samVersion

RUN apt-get update && apt-get -y install build-essential curl bzip2 zlib1g-dev 	libncurses-dev libbz2-dev liblzma-dev  libcurl4-openssl-dev librtmp-dev
RUN curl -L https://github.com/samtools/samtools/releases/download/$samVersion/samtools-$samVersion.tar.bz2  | tar -jvxf - 
WORKDIR samtools-$samVersion
RUN make -j 8
#Make the amd and arm64 directories so that they always exist for the next stage
RUN mkdir -p /lib/x86_64-linux-gnu
RUN mkdir -p /lib/aarch64-linux-gnu


FROM debian:bookworm-slim
ARG samVersion
ARG TARGETPLATFORM
COPY --from=builder /samtools-$samVersion/samtools /usr/local/bin/samtools
COPY --from=builder /samtools-$samVersion/misc/ /usr/local/bin/
ENV SOURCEPATH="/lib/x86_64-linux-gnu"
COPY --from=builder ${SOURCEPATH}/libcurl.so* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libncurses* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libnghttp* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libssh* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libps* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libssl* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/librtmp* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libcrypto* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libgssapi* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libldap* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/liblber* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libbrotl* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libkrb5* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libk5crypto* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libsas* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libkeyutils* ${SOURCEPATH}/
ENV SOURCEPATH="/lib/aarch64-linux-gnu"
COPY --from=builder ${SOURCEPATH}/libcurl.so* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libncurses* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libnghttp* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libssh* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libps* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libssl* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/librtmp* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libcrypto* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libgssapi* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libldap* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/liblber* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libbrotl* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libkrb5* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libk5crypto* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libsas* ${SOURCEPATH}/
COPY --from=builder ${SOURCEPATH}/libkeyutils* ${SOURCEPATH}/
COPY filterChrY.sh /usr/local/bin/filterChrY.sh