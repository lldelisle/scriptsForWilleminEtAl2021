use strict;
use Getopt::Long;

sub gen_split {
	my($l1, $l2, $l3, $l4, $outf, $re) = @_;
	my($l_re) = length($re);
	my($i) = -1;
	my($base) = 0;
	my($pos) = index($l2, $re, $base+1);
	my($seg) = 0;

	my(@head) = split(" ", $l1);
	my($head_suf) = join(" ", @head[1..$#head]);

	while($pos > $base) {
		printf $outf "$head[0]"."SEG$seg:%u:%u $head_suf\n",$base+1,$pos+$l_re;
		print $outf substr($l2, $base, $l_re + $pos - $base)."\n";
		print $outf "$l3\n";
		print $outf  substr($l4, $base, $l_re + $pos - $base)."\n";
		$seg++; 
		$base = $pos;	
		$pos = index($l2, $re, $base+1);
	}
	$pos = length($l2);
	printf $outf "$head[0]"."SEG$seg:%u:%u $head_suf\n",$base+1,length($l2);
	print $outf substr($l2, $base, $l_re + $pos - $base)."\n";
	print $outf "$l3\n";
	print $outf substr($l4, $base, $l_re + $pos - $base)."\n";
}

################code####################

die "usage: perl splitOneFastqOnRe.pl --input inputFile [can be gzip] --output outputFile [will be gzip] --re_seq GATC\n" if scalar(@ARGV) < 1;
my ($inputFile, $outputFile, $re_seq);

GetOptions("input=s" => \$inputFile,
		"output=s" => \$outputFile,
		"re_seq=s" => \$re_seq)
	or die "not all args supplied\n";


open(my $out_fh, "| gzip > $outputFile") or die "cannot write to $outputFile";

my $in_fh;
if ( $inputFile =~ /.*\.gz$/ ) {
	open( $in_fh, "zcat $inputFile |" ) or die "cannot open $inputFile";

} else {
	open( $in_fh, "<" ,$inputFile ) or die "cannot open $inputFile";
}
print STDERR "proccessing $inputFile\n";

my $read_cnt;

while (<$in_fh>) {
	$read_cnt++;
	my $l11 = $_; chomp $l11;
	my $l12 = <$in_fh>; chomp $l12;
	my $l13 = <$in_fh>; chomp $l13;
	my $l14 = <$in_fh>; chomp $l14;
	gen_split($l11,$l12,$l13,$l14, $out_fh, $re_seq);
	if ($read_cnt % 100000 == 0) {
		printf STDERR "proccessed %uK reads... ", $read_cnt/1000;	
	}
}

print "$read_cnt\n";
