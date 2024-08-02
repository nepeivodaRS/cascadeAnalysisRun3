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
             \"n_{#sigma, TPC}\",
             2, 6,
             0., 1.45,
             0.46, 1.14,
             false,
             $workingDir,
             \"_nsigmatpc\")"