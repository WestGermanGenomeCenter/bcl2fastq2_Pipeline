#!/bin/bash
data_dir=$1
cd $data_dir
mv RunInfo.xml backup_changed_RunInfo.xml
mv backup_original_RunInfo.xml RunInfo.xml
