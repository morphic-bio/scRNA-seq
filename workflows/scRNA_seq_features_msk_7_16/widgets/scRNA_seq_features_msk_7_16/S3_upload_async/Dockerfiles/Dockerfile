FROM debian:bookworm-slim
ARG DEBIAN_FRONTEND=noninteractive
ARG VERSION="2.0.30"
ARG ARCH  # Allow ARCH to be passed at build time, but it could be empty

# Run apt-get and architecture detection
RUN apt-get update && apt-get install -y curl zip unzip \
    # Check if ARCH is defined, if not detect using uname
    && if [ -z "$ARCH" ]; then \
        ARCH=$(uname -m); \
        if [ "$ARCH" = "aarch64" ]; then ARCH="arm64"; fi; \
        if [ "$ARCH" = "x86_64" ]; then ARCH="x86_64"; fi; \
    fi \
    # Download and install AWS CLI based on detected or passed ARCH
    && curl "https://awscli.amazonaws.com/awscli-exe-linux-$ARCH-$VERSION.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf awscliv2.zip aws \
    && apt-get remove -y curl zip unzip \
    && apt-get autoremove -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Copy the upload script
COPY upload.sh /usr/local/bin/upload.sh
RUN chmod +x /usr/local/bin/upload.sh
