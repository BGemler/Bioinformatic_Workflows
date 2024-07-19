## Build arg for all software versions
ARG VERSION_UBUNTU=22.04

## Start FROM ubuntu image with version specified
FROM ubuntu:${VERSION_UBUNTU}

## Disable interactive questions
ARG DEBIAN_FRONTEND=noninteractive

## Update Ubuntu
RUN apt update && \
	apt -y install python3-pip

## Seperate out to take advantage of cacheing
RUN apt-get -y install bowtie2 kraken2
    
## Add fastp and move to callable bin
ADD scripts/ /scripts/
RUN chmod -R +x /scripts/

CMD ["/scripts/main.py"]
