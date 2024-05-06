cat *.yaml | grep - | grep -v http | grep -v nodefaults | grep -v name >complete_software_env.tsv
