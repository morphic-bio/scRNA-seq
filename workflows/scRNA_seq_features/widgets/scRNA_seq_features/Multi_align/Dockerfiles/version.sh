#!/bin/bash

if [ $( uname -m ) = 'x86_64' ]; then
    ln -sf /usr/local/bin/"${VERSION}"/Linux_x86_64_static/STAR /usr/local/bin/STAR
else
    ln -sf /usr/local/bin/"${VERSION}"/arm64/STAR /usr/local/bin/STAR
fi
exec "$@"