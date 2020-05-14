#!/usr/bin/perl

$sq = "";
while(<STDIN>) {
   if(/^> *(\S+)/) {if($sq ne "") {process_seq(); $sq = ""}; $name=$1}
   else {chomp; $sq .= $_}
   }

process_seq();

sub process_seq {print "$name";
   $sq =~ tr/a-z/A-Z/; $sq =~ s/[^A-Z]/0/g; $sq =~ tr/ACGT/1234/;
   $l = length($sq); @sq = split //, $sq;
   for($i=0; $i<$l; $i++) {printf "%2i", $sq[$i]} 
   print "\n"; 
   } 
