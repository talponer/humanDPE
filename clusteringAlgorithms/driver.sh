# generation of numerical, R-friendly sequence library files: 

  ./seq2mono.pl dm6_epd05.seq >  dm6_epd05.dat
  ./seq2mono.pl hg19_epd05.seq > hg19_epd05.dat

# clustering of promoter regions: 

  Rscript --vanilla dm6_epd05_em1.R
  Rscript --vanilla dm6_epd05_em2.R
  Rscript --vanilla hg19_epd05_em1.R
  Rscript --vanilla hg19_epd05_em2.R

# visualiztation of the motifs via an HTML table 

  ./mk_logo_table.pl 3 1 cons_dm6_K3_s1  cons_dm6_K3_s1
  ./mk_logo_table.pl 6 2 cons_dm6_K6_s1  cons_dm6_K6_s1
  ./mk_logo_table.pl 3 1 cons_hg19_K3_s1 cons_hg19_K3_s1
  ./mk_logo_table.pl 6 2 cons_hg19_K6_s1 cons_hg19_K6_s1

