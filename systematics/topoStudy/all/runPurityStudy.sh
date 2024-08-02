#! /usr/bin/env bash

directory="results"
# Check if the directory exists
if [ -d "$directory" ]; then
    echo "Directory '$directory' exists."
else
    echo "Directory '$directory' does not exist."
    mkdir $directory
fi

ulimit -n 10000

root -l "../topoStudy.C(2,
             \"DATA.txt\",
             \"MC.txt\",
             \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial\",
             \"_purityStudy\")"