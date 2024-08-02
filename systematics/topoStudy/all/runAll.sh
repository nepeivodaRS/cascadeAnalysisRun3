#! /usr/bin/env bash
rm DATA.txt
rm MC.txt

cd ../../multitrial/mcPurity/yieldsOut
find "$(pwd)" -type f -name "*.root" > ../../../topoStudy/all/MC.txt
cd ../../dataPurity/yieldsOut
find "$(pwd)" -type f -name "*.root" > ../../../topoStudy/all/DATA.txt
cd ../../
#./runTopoStudy.sh
