#!/bin/bash
# log file where the terminal output will be saved
LOGFILE="cascpostprocessingLog.txt"

formatted_date=$(date +"%Y-%m-%d")  # Format date as YYYY-MM-DD

# directory of this script
DIR_THIS="$(dirname "$(realpath "$0")")"

OPTION="-b --configuration json:/${PWD}/jsonConfigs/ppRun3xi13tevMC.json"

o2-analysis-lf-cascpostprocessing ${OPTION} \
> "./logFiles/$LOGFILE" 2>&1

#--aod-file /Users/rnepeiv/workLund/PhD_work/run3omega/data/used_in_analysis_note/AO2D_500khz_without_inelgt0.root

# report status
rc=$?
if [ $rc -eq 0 ]; then
  echo "No problems!"
  resultsDir="${PWD}/results/run3_13tev/xi/lhc24b1b/$formatted_date"
  # Check if the directory exists
  if [ -d "$resultsDir" ]; then
      echo "Directory '$resultsDir' exists."
  else
      mkdir -p $resultsDir
      echo "Directory '$resultsDir' is created."
  fi
  mv "AnalysisResults.root" "$resultsDir/AnalysisResults.root"
  mv "dpl-config.json" "${PWD}/jsonConfigs/latestUsed.json"
else
  echo "Error: Exit code $rc"
  echo "Check the log file $LOGFILE"
  exit $rc
fi