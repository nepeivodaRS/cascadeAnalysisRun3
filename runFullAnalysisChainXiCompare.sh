#!/bin/bash
# Xi analysis
workingDir="\"${PWD}\""

date="2024-08-07"

pathToDATA="\"${PWD}/data/6aug-lhc22o-pass6-tight/AnalysisResults.root\""
pathToDATAPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc22o-pass6-medium-tight/$date/AnalysisResults.root\""

pathToMC="\"${PWD}/data/4aug-lhc24b1b-tight/AnalysisResults.root\""
pathToMCPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc24b1b-tight/$date/AnalysisResults.root\""

partType=2
inel=0

# for i in `seq 0 10`
#   do
#     root -l -b -q "yield.C($partType, $i, $inel, 0,
#                    $pathToDATA,
#                    $pathToDATAPP,
#                    $workingDir,
#                    \"_run2-medium-tight\")"
#   done

# root -l "effCorr.C($partType, $inel,
#                 $pathToMC,
#                 $pathToMCPP,
#                 $pathToDATA,
#                 $workingDir,
#                 \"_run2-medium-tight\",
#                 \"_LHC24b1b\")"

# root -l "yieldInMult.C($partType, 0, $workingDir, $pathToDATA, \"_run2-medium-tight\",  \"_LHC24b1b\")"
root -l "compPublXiMB.C($partType, $inel, $workingDir, \"_run2-medium-tight\")"
