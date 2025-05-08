cho "software_versions:" >software_mqc_versions.yml
cat *.yaml | grep = | grep -v name |grep -v "#" | sed 's/- //g' | sed 's/=/:"/g' | awk '{print "\t",$1."\""}' >>software_mqc_versions.yml