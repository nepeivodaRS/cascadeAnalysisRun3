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
             \"#Delta#it{M}_{#Lambda}, GeV/#it{c}^{2}\",
             0.001, 0.020,
             -0.4, 1.6,
             0.8, 1.2,
             false,
             $workingDir,
             \"_lambdamasswin\")"