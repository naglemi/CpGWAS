use warnings;

@cells = ("Astrocyte","Microglia","Neuronal","Oligo");
foreach $cell (@cells){
	$bed1 = $cell.".enhancer.bed";
	open(OUT, ">temp.bed");
	foreach $cell2 (@cells){
		if($cell2 ne $cell){
			$bed2 = $cell2.".enhancer.bed";
			open(IN, "../$bed2");
			while(<IN>){
				chomp;
				print OUT "$_\n";
			}
			close(IN);
		}
	}
	close(OUT);	
	$outfile = 	$cell.".bed";
	`bedtools subtract -a ../$bed1 -b temp.bed -A > $outfile`;
}
`rm temp.bed`;