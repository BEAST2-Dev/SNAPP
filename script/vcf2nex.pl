#!/bin/perl

# Attempts to convert vcf file to nexus that can be read by BEAUti with SNAPP template
# YMMV
# Usage: perl vcf2nex.pl < file.cvs > file.nex


# process header
while ($header !~ /CHROM/) {
	$header = <>;
}

@s = split('\s+', $header);
for ($i = 15; $i <= $#s; $i++) {
	$taxon[$i-15] = $s[$i];
}

# process body
$k = 0;
while ($s = <>) {
	@s = split("[\\s:]+", $s);
	for ($i = 15; $i < $#s; $i +=3) {
		$s = $s[$i];
		if ($s =~ /(.*)\/(.*)/) {
			$data[$i/3-5][0] = $data[$i/3-5][0].$1;
			$data[$i/3-5][1] = $data[$i/3-5][1].$2;
		} else {
			$data[$i/3-5][0] = $data[$i/3-5][0].'?';
			$data[$i/3-5][1] = $data[$i/3-5][1].'?';
		}	
	}
	$k++;
}

print "#NEXUS\n";
print "Begin data;\n";
print "        Dimensions ntax=".(2*$#taxon+2)." nchar=$k;\n";
print "        Format datatype=binary symbols=\"01\" gap=-;\n";
print "        Matrix\n";
for ($i = 0; $i <= $#taxon; $i++) {
	print "$taxon[$i]_1  $data[$i][0]\n";
	print "$taxon[$i]_2  $data[$i][1]\n";
}
print "End;\n"
