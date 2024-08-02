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
             \"N_{tpccrrows}\",
             50, 80,
             0.15, 1.35,
             0.9, 1.1,
             true,
             $workingDir,
             \"_tpccrrows\")"