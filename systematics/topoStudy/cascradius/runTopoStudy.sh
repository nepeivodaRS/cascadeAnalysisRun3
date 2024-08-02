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
             \"R_{casc.}, cm\",
             0.4, 1.3,
             0.55, 1.35,
             0.85, 1.15,
             true,
             $workingDir,
             \"_cascradius\")"