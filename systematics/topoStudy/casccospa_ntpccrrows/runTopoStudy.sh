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

root -l "../topoStudy2D.C(2,
             \"DATA.txt\",
             \"MC.txt\",
             \"casc cos(PA)\", \"N_{tpccrrows}\",
             0.969, 0.999,
             50, 80,
             16, 16,
             1,
             $workingDir,
             \"_casccospa-ntpccrrows\")"