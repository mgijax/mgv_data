BEGIN { print "["; }
{
  sep="";
  if ( NR > 1 ) sep= ",";
  print sep"[\"" $1 "\",\"" $2 "\",\"" substr($3,1,1)substr($4,1,1) "\"]"
}
END { print "]"; }
