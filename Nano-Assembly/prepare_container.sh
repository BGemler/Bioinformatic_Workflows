#!/bin/bash

# Run this script prior to building container. Pulls SPAdes executable file for linux from their GitHub

SRC_Loc="src/"
if ! [[ -d $SRC_Loc ]]; then
    # redirect to stderr so function only returns prinft output 
    echo -e "\nCreating src directory: ${SRC_Loc}\n" >&2 
    mkdir -p "${SRC_Loc}"
fi

# Download github SPAdes executable package
cd $SRC_Loc
wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz
