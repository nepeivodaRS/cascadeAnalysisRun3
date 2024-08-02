#! /usr/bin/env bash

directory="results"
# Check if the directory exists
if [ -d "$directory" ]; then
    echo "Directory '$directory' exists."
else
    echo "Directory '$directory' does not exist."
    mkdir $directory
fi

workingDir=\"${PWD}\"

root -l "../topoStudy1D.C(2,
             \"DATA.txt\",
             \"MC.txt\",
             \"R_{v0}, cm\",
             0.4, 2.5,
             0.80, 1.12,
             0.95, 1.05,
             true,
             $workingDir,
             \"_v0radius\")"