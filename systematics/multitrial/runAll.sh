#! /usr/bin/env bash
# rm -r ./data/config/*
# rm -r ./data/results/*
# rm -r ./data/yieldsOut/*
# rm -r ./data/efficiencies/*
# rm -r ./data/yieldsOutEffCorr/*
# rm -r ./mc/results/*
# rm -r ./mc/config/*
# rm -r ./mc/yieldsOut/*
# rm -r ./mc/yieldsOutEffCorr/*
# rm -r ./results/*
# rm DATA.txt
# rm MC.txt

# root -l -q "generatejson_systematics_multitrial.cpp"

cd data
# ./runVariedAnalysis.sh
./runVariedYields.sh
cd ../mc
./runEffCorr.sh
# cd ../../
