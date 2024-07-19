## Build arg for all software versions
ARG VERSION_UBUNTU=22.04

## Start FROM ubuntu image with version specified
FROM ubuntu:${VERSION_UBUNTU}

## Disable interactive questions
ARG DEBIAN_FRONTEND=noninteractive

## Update Ubuntu
RUN apt update && \
	apt-get update && \
	apt -y install python3-pip

RUN apt-get -y install fastp trimmomatic

## Add scripts and to callable bin
ADD scripts/ /scripts/
RUN chmod -R +x /scripts/

## Add adapter information
ADD adapters/ /adapters/

CMD ["/scripts/main.py"]
