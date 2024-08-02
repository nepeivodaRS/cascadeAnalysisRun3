#!/bin/bash
# Xi MC Closure
workingDir="\"${PWD}\""

date="2024-07-24"

pathToDATA="\"${PWD}/data/24jul-lhc22o-pass6/AnalysisResults.root\""
pathToDATAPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc22o_pass6/$date/AnalysisResults.root\""

pathToMC="\"${PWD}/data/24jul-lhc24b1b/AnalysisResults.root\""
pathToMCPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc24b1b/$date/AnalysisResults.root\""

# for i in `seq 0 1`
#   do
#     root -l -b -q "yield.C(2, $i, 0, 1,
#                    $pathToMC,
#                    $pathToMCPP,
#                    $workingDir,
#                    \"_mc-closure\")"
#   done

# root -l "effCorr.C(2, 0,
#                 $pathToMC,
#                 $pathToMCPP,
#                 $pathToMC,
#                 $workingDir,
#                 \"_mc-closure\",
#                 \"_LHC24b1b\")"

root -l "mcClosure.C(2, 0,
         $workingDir,
         $pathToMC,
         $pathToMCPP)"

# root -l "yieldInMultFitted.C(2, 0, 3, $workingDir, \"_mc-closure\")"
