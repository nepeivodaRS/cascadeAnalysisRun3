#! /usr/bin/env bash

directory="results"
# Check if the directory exists
if [ -d "$directory" ]; then
    echo "Directory '$directory' exists."
else
    echo "Directory '$directory' does not exist."
    mkdir $directory
fi

#DEF
o2-analysis-lf-cascpostprocessing -b --configuration json://config/configDEF.json
mv AnalysisResults.root $directory/DATA-LHC22o_apass6-DEF.root

rm dpl-config.json