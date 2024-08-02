#!/bin/bash
# log file where the terminal output will be saved
LOGFILE="fullLog.txt"

formatted_date=$(date +"%Y-%m-%d")  # Format date as YYYY-MM-DD

# directory of this script
DIR_THIS="$(dirname "$(realpath "$0")")"

OPTION="-b --configuration json:/${PWD}/jsonConfigs/ppRun3DataFull.json"

o2-analysis-timestamp ${OPTION} |
o2-analysis-ft0-corrected-table ${OPTION} |
o2-analysis-track-propagation ${OPTION} |
o2-analysis-event-selection ${OPTION} |
o2-analysis-multiplicity-table ${OPTION} |
o2-analysis-centrality-table ${OPTION} |
o2-analysis-pid-tof ${OPTION} |
o2-analysis-pid-tof-base ${OPTION} |
o2-analysis-pid-tpc-base ${OPTION} |
o2-analysis-pid-tpc ${OPTION} |
o2-analysis-lf-lambdakzerobuilder ${OPTION} |
o2-analysis-lf-cascadebuilder ${OPTION} |
o2-analysis-trackselection ${OPTION} |
o2-analysis-lf-cascpostprocessing ${OPTION} |
o2-analysis-lf-cascqaanalysis ${OPTION} \
> "./logFiles/$LOGFILE" 2>&1

#--aod-file /Users/rnepeiv/workLund/PhD_work/run3omega/data/used_in_analysis_note/AO2D_500khz_without_inelgt0.root

# report status
rc=$?
if [ $rc -eq 0 ]; then
  echo "No problems!"
  resultsDir="${PWD}/results/run3_13tev/full/lhc22o_pass6/$formatted_date"
  # Check if the directory exists
  if [ -d "$resultsDir" ]; then
      echo "Directory '$resultsDir' exists."
  else
      echo "Directory '$resultsDir' does not exist."
      mkdir -p $resultsDir
  fi
  mv "AnalysisResults.root" "$resultsDir/AnalysisResults.root"
  mv "dpl-config.json" "${PWD}/jsonConfigs/latestUsed.json"
else
  echo "Error: Exit code $rc"
  echo "Check the log file $LOGFILE"
  exit $rc
fi