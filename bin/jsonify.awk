BEGIN { print "["; }
{
  sep="";
  if ( NR > 1 ) sep= ",";
  openb = "[\""
  csep = "\",\""
  closeb = "\"]"
  sub(/NCBITaxon:/,"",$2)
  sub(/NCBITaxon:/,"",$4)
  print sep openb $1 csep $2 csep $3 csep $4 csep substr($5,1,1)substr($6,1,1) closeb
}
END { print "]"; }
