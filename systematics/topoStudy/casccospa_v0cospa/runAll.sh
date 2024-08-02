#! /usr/bin/env bash
rm -r ./data/config/*
rm -r ./data/results/*
rm -r ./data/yieldsOut/*
rm -r ./mc/results/*
rm -r ./mc/config/*
rm -r ./mc/yieldsOut/*
rm -r ./results/*
rm DATA.txt
rm MC.txt

root -l -q "generatejson_systematics.cpp"

cd data
./runVariedAnalysis.sh
./runVariedYields.sh
cd ../mc
./runVariedAnalysis.sh
./runVariedYields.sh
cd yieldsOut/
find "$(pwd)" -type f -name "*.root" > ../../MC.txt
cd ../../data/yieldsOut
find "$(pwd)" -type f -name "*.root" > ../../DATA.txt
cd ../../
# ./runTopoStudy.sh
