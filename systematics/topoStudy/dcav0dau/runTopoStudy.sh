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
             \"DCA_{V0 dau.}, cm\",
             0.1, 1.0,
             0.55, 1.35,
             0.85, 1.15,
             false,
             $workingDir,
             \"_dcav0dau\")"