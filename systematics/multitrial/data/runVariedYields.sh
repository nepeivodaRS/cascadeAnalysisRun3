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

pathToDATAfileList="/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/fileLists/input_data-AR.txt"
# Read the first line of the file into a variable using head and command substitution
pathToDATA=\"$(head -n 1 "$pathToDATAfileList")\"

# #default
# for i in `seq 0 10`
#     do
#         root -l -b -q "../../../yield.C(2, $i, 0, 0,
#                        $pathToDATA,
#                        \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/data/results/DATA-LHC22o_apass6-DEF.root\",
#                        $workingDir,
#                        \"_DEF\")"
#     done

# #loosest
# for i in `seq 0 10`
#     do
#         root -l -b -q "../../../yield.C(2, $i, 0, 0,
#                        $pathToDATA,
#                        \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/data/results/DATA-LHC22o_apass6-LOOSEST.root\",
#                        $workingDir,
#                        \"_LOOSEST\")"
#     done

#systematics
for i in {40..1000}
do
    for j in `seq 0 10`
        do
            root -l -b -q "../../../yield.C(2, $j, 0, 0,
                           $pathToDATA,
                           \"/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisSQM/systematics/multitrial/data/results/DATA-LHC22o_apass6-$i.root\",
                           $workingDir,
                           \"_$i\")"
        done
done
