#!/usr/bin/perl -w
#Alexander Bowles, University of Essex 2018
#Beginning Perl for Bioinformatics
use strict;
use Term::ANSIColor;

#Script to parse the output of MCL, checking for presence/absence of each species in each group of homology
#Each line of the MCL output is an homology group, each line will be checked for the presence of species labels/tags

#Program USAGE
my $USAGE = "Script to parse the output of MCL, checking for presence/absence of each species in each group of homology
USAGE: $0 MCL_output\n";
unless (@ARGV) {
	print color ("red"), $USAGE, color ("reset") and exit;
}

#Check if the file exists, open it
my $filename = $ARGV[0];
unless (check_file ($filename)) {
	open(FILE, $filename) or die "Could not open file $filename\n";
}

#Open output file
my $outfile = $filename;
$outfile =~ s/...(.*)\....$/$1/g;
$outfile = "02_".$outfile."_gene_numbers_parsed.out";
print "Outfile name: $outfile\n";
unless ( open(OUTPUT, ">$outfile") ) {
    print "Cannot open file \"$outfile\" to write to!!\n\n";
    exit;
}

#Count each taxa ocurrence in each line,
my @taxa = array_taxa ();
print "Taxa: @taxa\n";
my $line_number = 0;
my @output = ();

while (<FILE>) {
	$line_number++;
	@output = ();
	push (@output, $line_number);
	foreach my $taxon (@taxa){
		my $count = 0;
		while ($_ =~ /$taxon/g) {
			$count++;
		}
		push (@output, $count);
	}
	print OUTPUT "@output\n"; #Print results to output file, one line at a time.
}
print "\nResults printed in file $outfile\n";

#Close file handles
close(OUTPUT);
close FILE;

#End of program
print color ("red"), "\nend\n\n", color ("reset");
exit;


################################################################################
#                                Subroutines                                   #
################################################################################

#CHECK_FILE: checks file properties, checks if file exists, is flat, is empty, or can be opened
#-e exists, -f flat file, -s empty
sub check_file {
    my ($filename) = @_;
    unless (-e $filename || (undef $filename)) {print "File does NOT exists\n" and exit;}
    unless (-f $filename) {print "File is NOT a flat file\n" and exit;}
    unless (-s $filename) {print "File is empty\n" and exit;}
    unless (open (FH, $filename)) {print "Can not be opened\n" and exit;}
    close FH;
    return;
}

#ARRAY_TAXA: just an array definition to not clutter the main program.
sub array_taxa {	
	my  (@taxa)  = qw/acas ppal ddis ttra falb spun amac tmel scer cfra sros mlei lgig ctel cele dmel drer xtro ggal acar hsap ngru bnat ehux gthe ptet tgon pine ptro tpse cpar gsul gphl cmer ppur bpra oluc otau mcom mpus apro cvar psp. csub czof gpec vcar crei ceus kfla mpol ppat smoe gbil pabi gmon atri spol zmar zmue drot ashe dcat pequ aoff pdac egui cnus acom obra opun oglu oruf omer obar ogla oniv osai osaj lper pedu bdis hvul atau tura taes ttur lpee ecru sita zmay sbic etef zjap zmat zpac otho miti mbal macu ecal mcor nnuc kfed rcre vvin lang adur aipa ccaj pang pvul vrad gmax gsoj cari mtru tpra fves pavi pmum pper pbre pcom mdom zjuj mnot cmax cmos clan lsic cmel csat ptri ppru lusi jcur mesc hbra rcom cfol pgra egra dlon abux ccle csin cich cmed cpap thas esal tpar bnap bole brap siri bvug crub cgra atha alyr tcac ccap coli dzib garb ghir grai fesc dcay bvul sole cqui ahyp cacu csie achi lsat ebre hann pgin dcar ccan cgig fexc oeur bhyg mgut sind himp gaur ugib inil itri pinf paxi nobt nsyl ntab ntom slyc spen spim stub cann cbac cchi/;
return @taxa;
}
