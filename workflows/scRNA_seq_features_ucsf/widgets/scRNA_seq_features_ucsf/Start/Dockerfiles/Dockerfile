FROM alpine:3.12.1
RUN apk update upgrade --no-cache && apk add bash && rm -rf /var/cache/apk* 
COPY start.sh /usr/local/bin/.
