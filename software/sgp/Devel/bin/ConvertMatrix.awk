BEGIN{
  getline<ARGV[1];
  for (i=1;i<=NF;i++)
    AA[i+1]=$i;
  
  while(getline<ARGV[1]>0)
    for (i=2;i<=NF;i++)
      MAT[$1,AA[i]]=$i;


  for (i in MAT) {
    split(i,AAs,SUBSEP);
    print AAs[1],AAs[2],MAT[i];
  }

}
