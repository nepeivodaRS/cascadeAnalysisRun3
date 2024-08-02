#! /usr/bin/env bash

directory="/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/data/yieldsOutEffCorr"
# Check if the directory exists
if [ -d "$directory" ]; then
    echo "Directory '$directory' exists."
else
    echo "Directory '$directory' does not exist."
    mkdir $directory
fi

directory="/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/data/efficiencies"
# Check if the directory exists
if [ -d "$directory" ]; then
    echo "Directory '$directory' exists."
else
    echo "Directory '$directory' does not exist."
    mkdir $directory
fi

workingDir="\"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/data\""

pathToMCfileList="/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/fileLists/input_mc-AR.txt"
# Read the first line of the file into a variable using head and command substitution
pathToMC=\"$(head -n 1 "$pathToMCfileList")\"

pathToDATAfileList="/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/fileLists/input_data-AR.txt"
# Read the first line of the file into a variable using head and command substitution
pathToDATA=\"$(head -n 1 "$pathToDATAfileList")\"

pathToEfficiencyCalib="\"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/data/efficiencies/effCalibration_LHC24b1b.root\""

root -l -b -q "../../../calibrationMC.C(
                $pathToMC,
                $pathToDATA,
                $workingDir,
                \"_LHC24b1b\")"

root -l -b -q "../../../effCorr.C(2, 0,
                $pathToMC,
                \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/mc/results/MC-LHC24b1-DEF.root\",
                $pathToEfficiencyCalib,
                $workingDir,
                \"_DEF\",
                \"_LHC24b1b\")"

root -l -b -q "../../../effCorr.C(2, 0,
                $pathToMC,
                \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/mc/results/MC-LHC24b1-LOOSEST.root\",
                $pathToEfficiencyCalib,
                $workingDir,
                \"_LOOSEST\",
                \"_LHC24b1b\")"

for i in {1..1000}
    do
        root -l -b -q "../../../effCorr.C(2, 0,
                        $pathToMC,
                        \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/mc/results/MC-LHC24b1-$i.root\",
                        $pathToEfficiencyCalib,
                        $workingDir,
                        \"_$i\",
                        \"_LHC24b1b\")"
    done

