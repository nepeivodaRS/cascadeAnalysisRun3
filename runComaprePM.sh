#!/bin/bash
# Xi analysis
workingDir="\"${PWD}\""

date="2024-08-07"

pathToDATA="\"${PWD}/data/6aug-lhc22o-pass6-tight/AnalysisResults.root\""
pathToDATAPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc22o-pass6-medium-tight/$date/AnalysisResults.root\""

pathToMC="\"${PWD}/data/4aug-lhc24b1b-tight/AnalysisResults.root\""
pathToMCPP="\"${PWD}/O2PhysicsRuns/results/run3_13tev/xi/lhc24b1b-tight/$date/AnalysisResults.root\""

for i in `seq 0 10`
  do
    root -l -b -q "yield.C(0, $i, 0, 0,
                   $pathToDATA,
                   $pathToDATAPP,
                   $workingDir,
                   \"-pm\")"
  done

for i in `seq 0 10`
  do
    root -l -b -q "yield.C(1, $i, 0, 0,
                   $pathToDATA,
                   $pathToDATAPP,
                   $workingDir,
                   \"-pm\")"
  done

root -l "effCorr.C(0, 0,
                $pathToMC,
                $pathToMCPP,
                $pathToDATA,
                $workingDir,
                \"_pm\",
                \"_LHC24b1b\")"

root -l "effCorr.C(1, 0,
                $pathToMC,
                $pathToMCPP,
                $pathToDATA,
                $workingDir,
                \"_pm\",
                \"_LHC24b1b\")"

# root -l "meanInMult.C(0, 0, $workingDir)"
# root -l "sigmaInMult.C(0, 0, $workingDir)"
# root -l "purityInMult.C(0, 0, $workingDir)"
# root -l "yieldInMult.C(0, 0, $workingDir, $pathToDATA,  \"\",  \"_LHC24b1b\")"
# root -l "yieldInMult.C(1, 0, $workingDir, $pathToDATA,  \"\",  \"_LHC24b1b\")"
root -l "pmCompare.C(0, 0, $workingDir,  \"_pm\")" # 0 for Xi, 1 for Omega
