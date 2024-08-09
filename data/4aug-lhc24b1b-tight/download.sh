#!/bin/bash

path="/Users/rnepeiv/workLund/PhD_work/run3omega/cascadeAnalysisPostSQM/data/4aug-lhc24b1b-tight"

touch ${path}/tmp/ddd.sh
while read p; do
  echo "alien_cp $p file::${path}${p#alien:///alice/cern.ch*}" >> ${path}/tmp/ddd.sh
done < ${path}/files.txt

chmod u+x ${path}/tmp/ddd.sh
time parallel --progress -j 10 -a ${path}/tmp/ddd.sh
#source ${path}/tmp/ddd.sh
mv ${path}/tmp/ddd.sh ${path}/tmp/prevDownload.sh