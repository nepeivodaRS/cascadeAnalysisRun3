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
             \"c#tau\",
             2, 8,
             0.5, 1.25,
             0.9, 1.1,
             false,
             $workingDir,
             \"_proplifetime\")"