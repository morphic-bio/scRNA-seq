FROM alpine:3.17.1
RUN apk update upgrade --no-cache && apk add bash curl wget unzip tar bzip2\
 && rm -rf /var/cache/apk*
COPY gather_files.sh /usr/local/bin/gather_files.sh
RUN chmod +x /usr/local/bin/gather_files.sh
