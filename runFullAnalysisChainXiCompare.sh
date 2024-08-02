#!/bin/bash
# Xi analysis
workingDir="\"${PWD}\""

date="2024-07-24"

pathToDATA="\"${PWD}/data/24jul-lhc22o-pass6/AnalysisResults.root\""
pathToDATAPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc22o_pass6/$date/AnalysisResults.root\""

pathToMC="\"${PWD}/data/24jul-lhc24b1b/AnalysisResults.root\""
pathToMCPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc24b1b/$date/AnalysisResults.root\""

partType=2

for i in `seq 0 10`
  do
    root -l -b -q "yield.C($partType, $i, 0, 0,
                   $pathToDATA,
                   $pathToDATAPP,
                   $workingDir,
                   \"_run2\")"
  done

root -l "effCorr.C($partType, 0,
                $pathToMC,
                $pathToMCPP,
                $pathToDATA,
                $workingDir,
                \"_run2\",
                \"_LHC24b1b_run2\")"

# root -l "yieldInMult.C($partType, 0, $workingDir, $pathToDATA, \"_run2\",  \"_LHC24b1b\")"
root -l "compPublXiMB.C($partType, 0, $workingDir, \"_run2\")"
