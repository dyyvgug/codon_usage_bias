use strict;
use warnings;
use Statistics::Descriptive; 
use Number::Format qw/round/;

my $stat = Statistics::Descriptive::Full->new();

my $win_size = 30;

my $slide = 1;

my %cai;
	open (hand2, "CBI_ref_hg38.txt") or die $!;

	while (<hand2>)  {
		chomp;
		next if /^AA/;
		my ($aa, $codon, $CAI, $rand) = split /\t/;
		$cai{$codon} = $CAI;
		}

	close hand2;

die "Usage: $0 <dir> <extion>\n" unless @ARGV == 2;
my $Dir = $ARGV[0] ;
my $Ext = $ARGV[1] ;
opendir(DH, "$Dir") or die "Can't open: $!\n" ;

my @list = grep {/$Ext$/ && -f "$Dir/$_" } readdir(DH) ;
closedir(DH) ;
chdir($Dir) or die "Can't cd dir: $!\n" ;
foreach my $file (@list){
	open(FH, "$file") or die "Can't open: $!\n" ;
	print "$file:\n" ;
	open (OUTFILE, ">$Dir/$file.csv");
	while(<FH>){
		my $seq;
		print "$_\n";
		$seq .= $_;
		$seq =~ s/\s+//;  
	

	for (my $i = 0; $i < length($seq)- $win_size * 3; $i += $slide * 3)  {
		#calculate sliding window of 10aa
		my @val = ();
		for (my $j = $i; $j < $i + $win_size*3; $j += 3)  {
			my $cod = substr($seq, $j, 3);
			push @val, $cai{$cod};
			}
		print join(",",@val)."\n";
		$stat->add_data(@val);
		my $cai_avg = round($stat->geometric_mean(), 3);
		$stat->clear();
		#my $position = int($i/3)+1;
		print OUTFILE "$i,$cai_avg\n";  }
		$seq = '';

	}
close(FH) ; close OUTFILE;
}




	
