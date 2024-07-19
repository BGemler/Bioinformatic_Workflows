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
RUN apt-get -y install kraken2 python-is-python3
    
## Add bracken, download it
ADD src/ /src/
RUN cd "/src/" && \
	tar -xvf "v2.9.tar.gz"

RUN cd "/src/Bracken-2.9/" && \
		bash "./install_bracken.sh"

## Install Python packages
RUN pip3 install numpy scipy

ADD KrakenTools/ /KrakenTools/
RUN chmod -R +x /KrakenTools/

ADD scripts/ /scripts/
RUN chmod -R +x /scripts/

CMD ["/scripts/utils.py"]