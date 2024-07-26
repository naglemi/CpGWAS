use warnings;

open(IN, "/dcs04/lieber/statsgen/shizhong/wgbs/finemap/dlpfc/pheno/chr1/1/p1");
while(<IN>){
	chomp;
	@tokens=split(' ',$_);
	$id = $tokens[0];
	$id =~ s/Br0/Br/;
	$tag{$id} = 1;
}
close(IN);

open(OUT, ">dlpfc_brnum");
foreach $id (keys %tag){
	print OUT "$id\n";
}
close(OUT);
	