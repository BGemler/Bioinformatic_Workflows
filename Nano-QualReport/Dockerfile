## Build arg for all software versions
ARG VERSION_UBUNTU=22.04

## Start FROM ubuntu image with version specified
FROM ubuntu:${VERSION_UBUNTU}

## Disable interactive questions
ARG DEBIAN_FRONTEND=noninteractive

## Update Apt & Install FastQC
RUN apt update && \
	apt -y install fastqc

## Install pip & Multiqc
RUN apt -y install python3-pip && \
	pip install multiqc

## Add scripts to container, make executable
ADD scripts /scripts
RUN chmod +x /scripts/main.py

CMD ["/scripts/main.py"]