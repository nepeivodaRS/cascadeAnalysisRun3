#! /usr/bin/env bash

directory="yieldsOut"
# Check if the directory exists
if [ -d "$directory" ]; then
    echo "Directory '$directory' exists."
else
    echo "Directory '$directory' does not exist."
    mkdir $directory
fi

workingDir=\"${PWD}\"

pathToMCfileList="/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/fileLists/input_mc-AR.txt"
# Read the first line of the file into a variable using head and command substitution
pathToMC=\"$(head -n 1 "$pathToMCfileList")\"

#default
root -l -b -q "../../../yield.C(2, 0, 0, 1,
               $pathToMC,
               \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/mc/results/MC-LHC24b1-DEF.root\",
               $workingDir,
               \"_DEF\")"

#loosest
root -l -b -q "../../../yield.C(2, 0, 0, 1,
               $pathToMC,
               \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/mc/results/MC-LHC24b1-LOOSEST.root\",
               $workingDir,
               \"_LOOSEST\")"

#systematics
for i in {1..1000}
do
    root -l -b -q "../../../yield.C(2, 0, 0, 1,
                   $pathToMC,
                   \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/mc/results/MC-LHC24b1-$i.root\",
                   $workingDir,
                   \"_$i\")"
done
