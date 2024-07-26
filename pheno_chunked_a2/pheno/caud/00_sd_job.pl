use warnings;

mkdir "out";
open(OUT, ">sd_jobs.txt");
foreach $i (1..22){
	print OUT "Rscript sd.R $i\n";
}
close(OUT);