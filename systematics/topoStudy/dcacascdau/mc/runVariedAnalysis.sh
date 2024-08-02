#! /usr/bin/env bash

directory="results"
# Check if the directory exists
if [ -d "$directory" ]; then
    echo "Directory '$directory' exists."
else
    echo "Directory '$directory' does not exist."
    mkdir $directory
fi

#default
o2-analysis-lf-cascpostprocessing -b --configuration json://config/configDEF.json
mv AnalysisResults.root results/MC-LHC24b1-DEF.root

#systematics
for i in {1..20}
do
    o2-analysis-lf-cascpostprocessing -b  --configuration json://config/systjson$i.json
    mv AnalysisResults.root results/MC-LHC24b1-$i.root
done

rm dpl-config.json