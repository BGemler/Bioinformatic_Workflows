## Build arg for all software versions
ARG VERSION_UBUNTU=22.04

## Start FROM ubuntu image with version specified
FROM ubuntu:${VERSION_UBUNTU}

## Disable interactive questions
ARG DEBIAN_FRONTEND=noninteractive

## Update Ubuntu
RUN apt update && \
	apt-get update && \
	apt -y install python3-pip zlib1g-dev

RUN pip install matplotlib

## Add scripts and to callable bin
ADD src/ /src/

# Unzip SPAdes and QUAST executable from /src/
RUN tar -xvf /src/SPAdes-4.0.0-Linux.tar.gz
RUN tar -xvf /src/quast-5.2.0.tar.gz

# Install QUAST
WORKDIR /quast-5.2.0/
RUN python3 setup.py install_full

ADD scripts/ /scripts/

# Make executable permissions
RUN chmod -R +x /scripts/
RUN chmod +x /SPAdes-4.0.0-Linux/bin/spades.py
RUN chmod +x /quast-5.2.0/quast.py

CMD ["/scripts/main.py"]
