#!/bin/bash
# Xi MC Closure
workingDir="\"${PWD}\""

date="2024-08-07"

pathToDATA="\"${PWD}/data/6aug-lhc22o-pass6-tight/AnalysisResults.root\""
pathToDATAPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc22o-pass6-medium-tight/$date/AnalysisResults.root\""

pathToMC="\"${PWD}/data/4aug-lhc24b1b-tight/AnalysisResults.root\""
pathToMCPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc24b1b-tight/$date/AnalysisResults.root\""

# for i in `seq 0 10`
#   do
#     root -l -b -q "yield.C(2, $i, 0, 1,
#                    $pathToMC,
#                    $pathToMCPP,
#                    $workingDir,
#                    \"-mc-closure-tight\")"
#   done

root -l "effCorr.C(2, 0,
                $pathToMC,
                $pathToMCPP,
                $pathToMC,
                $workingDir,
                \"-mc-closure-tight\",
                \"-LHC24b1b-mc-closure\")"

root -l "mcClosure.C(2, 0,
         $workingDir,
         $pathToMC,
         $pathToMCPP,
         \"-mc-closure-tight\")"

# root -l "yieldInMultFitted.C(2, 0, 3, $workingDir, \"-mc-closure-tight\")"
