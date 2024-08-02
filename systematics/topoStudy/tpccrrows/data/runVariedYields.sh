#! /usr/bin/env bash

directory="yieldsOut"
# Check if the directory exists
if [ -d "$directory" ]; then
    echo "Directory '$directory' exists."
else
    echo "Directory '$directory' does not exist."
    mkdir $directory
fi

cut="tpccrrows"

workingDir=\"${PWD}\"

pathToDATAfileList="/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/fileLists/input_data-AR.txt"
# Read the first line of the file into a variable using head and command substitution
pathToDATA=\"$(head -n 1 "$pathToDATAfileList")\"

#default
root -l -b -q "../../../../yield.C(2, 0, 0, 0,
               $pathToDATA,
               \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/topoStudy/$cut/data/results/DATA-LHC22o_apass6-DEF.root\",
               $workingDir,
               \"_DEF\")"

#systematics
for i in {1..31}
do
    root -l -b -q "../../../../yield.C(2, 0, 0, 0,
                   $pathToDATA,
                   \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/topoStudy/$cut/data/results/DATA-LHC22o_apass6-$i.root\",
                   $workingDir,
                   \"_$i\")"        
done
