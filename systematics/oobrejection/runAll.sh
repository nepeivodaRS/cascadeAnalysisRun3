#! /usr/bin/env bash
rm -r ./data/results/*
rm -r ./data/yieldsOut/*
rm -r ./data/efficiencies/*
rm -r ./data/yieldsOutEffCorr/*
rm -r ./mc/results/*

cd data
./runVariedAnalysis.sh
./runVariedYields.sh
cd ../mc
./runVariedAnalysis.sh
./runEffCorr.sh
