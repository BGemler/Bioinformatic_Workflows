#!/bin/bash


DockerContainers=("Nano-HostRemoval" "Nano-K2TaxClass" "Nano-QualReport" "Nano-ReadTrim" "Nano-SummarizeUtils")

for container in ${DockerContainers[@]}; do
    # Enter container
    cd $container

    # Build the container
    ./build

    # Enter root directory
    cd ..
done
