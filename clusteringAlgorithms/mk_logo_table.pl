#!/usr/bin/perl
#
$n=$ARGV[0];
$m=$ARGV[1];
$prefix=$ARGV[2];
$w=int(850/$m);

print "<table border=1 cellpadding=2 width=900>\n";
for($i=1; $i<=$n; $i += $m) {print "<tr>";
   for($j=$i; $j<$i+$m; $j++) {print "<td><font size=2>";
      $name=$prefix . sprintf("_%03i",$j);
      $file=`ls img/$name*`; $name=$file; $name =~ s/img\///; $name =~ s/.jpg//; 
      print "$name <br><img src=$file width=$w></td>\n"}
   print "</tr>"}
print "</table>";

