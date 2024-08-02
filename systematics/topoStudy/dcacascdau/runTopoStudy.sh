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
             \"DCA_{cascdau}, cm\",
             0.1, 1.0,
             0.45, 1.25,
             0.8, 1.2,
             false,
             $workingDir,
             \"_dcacascdau\")"