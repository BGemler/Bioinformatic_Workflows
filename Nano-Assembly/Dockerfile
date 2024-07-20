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

RUN pip install matplotlib

## Add scripts and to callable bin
ADD src/ /src/

# Unzip SPAdes executable from /src/
RUN cd /src/
RUN tar -xvf /src/SPAdes-4.0.0-Linux.tar.gz

ADD scripts/ /scripts/

# Make executable permissions
RUN chmod -R +x /scripts/
RUN chmod +x /SPAdes-4.0.0-Linux/bin/spades.py

CMD ["/scripts/main.py"]
