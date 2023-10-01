cat $1 | awk -F'\t' '{split($1,a,"_"); split($6,b,"_");  if (a[1]=="bahahaHasm" && b[1]=="bahahascf" && a[2]==b[2] ) print $0; }'
