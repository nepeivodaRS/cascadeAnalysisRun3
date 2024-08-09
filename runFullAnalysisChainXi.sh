#!/bin/bash
# Xi analysis
workingDir="\"${PWD}\""

date="2024-08-07"

pathToDATA="\"${PWD}/data/6aug-lhc22o-pass6-tight/AnalysisResults.root\""
pathToDATAPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc22o-pass6-medium-tight/$date/AnalysisResults.root\""

pathToMC="\"${PWD}/data/4aug-lhc24b1b-tight/AnalysisResults.root\""
pathToMCPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc24b1b-tight/$date/AnalysisResults.root\""

# for i in `seq 0 10`
#   do
#     root -l -b -q "yield.C(2, $i, 0, 0,
#                    $pathToDATA,
#                    $pathToDATAPP,
#                    $workingDir,
#                    \"-medium-tight\")"
#   done

# root -l "effCorr.C(2, 0,
#                 $pathToMC,
#                 $pathToMCPP,
#                 $pathToDATA,
#                 $workingDir,
#                 \"-medium-tight\",
#                 \"_LHC24b1b\")"

# root -l "meanInMult.C(2, 0, $workingDir)"
# root -l "sigmaInMult.C(2, 0, $workingDir)"
# root -l "purityInMult.C(2, 0, $workingDir)"
# root -l "yieldInMult.C(2, 0, $workingDir, $pathToDATA,  \"medium-tight\",  \"_LHC24b1b\")"
root -l "yieldInMultFitted.C(2, 0, 3, $workingDir, \"-medium-tight\")"
