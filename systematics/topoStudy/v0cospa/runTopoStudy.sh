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
             \"v0 cos(PA)\",
             0.969, 0.999,
             -0.35, 1.55,
             0.9, 1.1,
             true,
             $workingDir,
             \"_v0cospa\")"