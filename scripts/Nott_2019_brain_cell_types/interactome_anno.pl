use warnings;

# read gene anno
open(IN, "Neuronal_interactome1_anno.bed");
while(<IN>){
	chomp;
	@tokens=split(' ',$_);
	$bin=$tokens[0]."_".$tokens[1]."_".$tokens[2];
	$record=$tokens[3]."_".$tokens[4]."_".$tokens[5]."_".$tokens[6]."_".$tokens[7];
	push @{$anno1{$bin}},$record;
}
close(IN);

open(IN, "Neuronal_interactome2_anno.bed");
while(<IN>){
	chomp;
	@tokens=split(' ',$_);
	$bin=$tokens[0]."_".$tokens[1]."_".$tokens[2];
	$record=$tokens[3]."_".$tokens[4]."_".$tokens[5]."_".$tokens[6]."_".$tokens[7];
	push @{$anno2{$bin}},$record;
}
close(IN);

# output contacts linked to promoters
open(IN, "Neuronal_interactome_way1.bed");
open(OUT, ">Neuronal_interactome_way1_anno.bed");
while(<IN>){
	chomp;
	@tokens=split(' ',$_);
	$bin=$tokens[3]."_".$tokens[4]."_".$tokens[5];
	if(defined $anno2{$bin}){
		foreach $record (@{$anno2{$bin}}){
			print OUT "$_\t";
			@temp=split(/_/,$record);
			foreach $rec(@temp){
				print OUT "$rec\t";
			}
			print OUT "\n";
		}		
	}
}
close(IN);
close(OUT);

open(IN, "Neuronal_interactome_way2.bed");
open(OUT, ">Neuronal_interactome_way2_anno.bed");
while(<IN>){
	chomp;
	@tokens=split(' ',$_);
	$bin=$tokens[3]."_".$tokens[4]."_".$tokens[5];
	if(defined $anno1{$bin}){
		foreach $record (@{$anno1{$bin}}){
			print OUT "$_\t";
			@temp=split(/_/,$record);
			foreach $rec(@temp){
				print OUT "$rec\t";
			}
			print OUT "\n";
		}		
	}
}
close(IN);
close(OUT);
