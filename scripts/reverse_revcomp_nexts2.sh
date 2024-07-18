#!/bin/bash
data_dir=$1
cd $data_dir
mv RunInfo.xml backup_changed_RunInfo.xml
cat backup_original_RunInfo.xml | sed '/Read Number/,/3/ s/IsReverseComplement="N"/IsReverseComplement="Y"/' >RunInfo.xml
# first find a line with read number, then in the next line search for a 3, in that line with a 3 make the conversion from True to false 
#cat backup_original_RunInfo.xml | sed 's/IsReverseComplement="Y"/IsReverseComplement="N"/g' >RunInfo.xml # first trial, set everything the same
# sed '/Read Number/,/3/ s/IsReverseComplement="Y"/IsReverseComplement="N"/'
# made according to : https://superuser.com/questions/1689644/find-and-replace-text-in-a-file-after-match-of-pattern-only-for-first-occurrence/1689959#1689959
