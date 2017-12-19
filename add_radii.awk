BEGIN{
  while(( getline line<"newcell.lmp") > 0 ) {
     print line
  }
}
