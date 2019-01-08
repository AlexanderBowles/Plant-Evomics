#! /usr/bin/perl -w
#Alexander Bowles, University of Essex 2018; with help of Jordi Paps (Oxford, Essex & Bristol) and Patrick Gemmell (Oxford)
#Beginning Perl for Bioinformatics
use strict;
use AnyDBM_File;
use Fcntl;
use Term::ANSIColor;

#Script to extract groups of homology from a parsed MCL output, based on the taxonomic distribution

#Introduce file names
my $MCL_out = "Input/mcl_output_all_2"; 
my $MCL_columns_parsed = "Input/02_mcl_output_all_2_gene_numbers_parsed.out";

#Check if files exists, open them
check_file ($MCL_out);
check_file ($MCL_columns_parsed);

#Create DBM files from the MCL output file for fast lookups.
my %MCL_out;
tie (%MCL_out,"AnyDBM_File", "$MCL_out.db", O_CREAT|O_RDWR, 0777) or die + "Can't open database: $!\n";
#%MCL_out = MCL_output_to_hash ($MCL_out);
print "\nMCL out (showing some random lines out of ", scalar (keys %MCL_out), "):\n";
my @MCL_out_keys = keys %MCL_out;
foreach (@MCL_out_keys[0..4]){
	my $first_elements = substr ($MCL_out{$_}, 0, 80);
	print "$_ => $first_elements...\n";
} print "...\n";

#Create DBM files from the MCL parsed file for fast lookups.
my %MCL_columns_parsed;
tie (%MCL_columns_parsed,"AnyDBM_File", "$MCL_columns_parsed.db", O_CREAT|O_RDWR, 0777) or die + "Can't open database: $!\n";
#%MCL_columns_parsed = MCL_columns_parsed_to_hash ($MCL_columns_parsed);
print "\nMCL parsed (each column is one species, showing some random lines out of ", scalar (keys %MCL_columns_parsed), "):\n\n";
print "Species: acas ppal ddis ttra falb spun amac tmel scer cfra sros mlei lgig ctel cele dmel drer xtro ggal acar hsap ngru bnat ehux gthe ptet tgon pine ptro tpse cpar gsul gphl cmer ppur bpra oluc otau mcom mpus apro cvar psp. csub czof gpec vcar crei ceus kfla mpol ppat smoe gbil pabi gmon atri spol zmar zmue drot ashe dcat pequ aoff pdac egui cnus acom obra opun oglu oruf omer obar ogla oniv osai osaj lper pedu bdis hvul atau tura taes ttur lpee ecru sita zmay sbic etef zjap zmat zpac otho miti mbal macu ecal mcor nnuc kfed rcre vvin lang adur aipa ccaj pang pvul vrad gmax gsoj cari mtru tpra fves pavi pmum pper pbre pcom mdom zjuj mnot cmax cmos clan lsic cmel csat ptri ppru lusi jcur mesc hbra rcom cfol pgra egra dlon abux ccle csin cich cmed cpap thas esal tpar bnap bole brap siri bvug crub cgra atha alyr tcac ccap coli dzib garb ghir grai fesc dcay bvul sole cqui ahyp cacu csie achi lsat ebre hann pgin dcar ccan cgig fexc oeur bhyg mgut sind himp gaur ugib inil itri pinf paxi nobt nsyl ntab ntom slyc spen spim stub cann cbac cchi\n\n";
my @MCL_columns_parsed_keys = keys %MCL_columns_parsed;
foreach (@MCL_columns_parsed_keys[0...4]){
	print "$_ => $MCL_columns_parsed{$_}\n";
} print "...\n";

#Create hash from the clade definition in subrout CLADES
my %spp;
my @value_list = ();
%spp = hash_spp();

#Calculate how many columns (species/terminal tip) are in the hash
walk_hash (\%spp, \@value_list);
my @total_columns = scalar @value_list;
print "\nNumber of species: ", (scalar @value_list), "\n";

#Ask user for clade/spp to check genes taxonomic distribution, perform search
my $user_input = '';	#To store the searching criteria fom user
my @arguments = ();	#To store the split of the searching criteria
my @search= ();  	#To store the taxa and options of the searching criteria
my $taxa  = '';		#To store the taxa from @search
my $option = '';	#To store options from @search (present, absent, minus, atleast,...)
my @columns = ();	#To store columns to search
my $final_search;	#To store columns to search plus the options (present, absent, minus, atleast,...)
my %final_searches = ();#To store ALL searches columns and options (present, absent, minus, atleast,...)
my @true_flags = (); 	#To store the flags that indicate if a homology group/line passes all the queries checks
my @outgroup = ();	#To store the columns left over at the end after extracting all the ingroups columns
my @good_homology_groups_spp_names = ();	#To store the groups of homology fullfilling the search criteria
my @good_homology_groups_spp_names_copy = ();	#Backup
my @good_homology_groups_columns_parsed = ();	#To store the columns of the groups of homology fullfilling the search criteria

OUTER: do {
	print "\nPlease, enter name of the clade/species to search (\"example\" for some samples, \"tree\" to print the evolutionary tree, Enter to exit):\n";
	$user_input = <STDIN>;
	chomp $user_input;
	unless ($user_input =~ /^\s*$|^exit|^quit/i) {			#Unless user wants to exit...
		if ($user_input =~ /example|help/i) {
			print_examples();			#Print examples of commands if user requests it
		} elsif ($user_input =~ /tree/i) {
			print_hash_colors(%spp);		#Print examples of commands if user requests it
		} else {
			#Here the real search starts, emptying variables for next loop and parsing the user input
			%final_searches = ();			#Empty the hash containing all the search conditions
			@good_homology_groups_spp_names = ();	#Empty the hash containing the results from previous search
			@good_homology_groups_columns_parsed = (); #Empty the hash containing the columns of the groups of homology fullfilling the search criteria
			@arguments = split (" ", $user_input);	#Decompose the user input in different arguments entered by user, each taxa query delimited by a space (" ")
			#print "\@arguments: @arguments\n";			
			foreach (@arguments) {
				#print "\$_ in \@arguments: $_\n";
				@search = split ("-", $_);	#Decompose each argument into taxa (first item) and options (second item)
				#print "\@search: @search\n";
				$taxa = $search[0];
				$option = $search[1];
				unless (defined $option && $option =~ /pre|all|abs|none|min|but|atl|onl|jus/gi){
						print "Taxa $taxa is missing valid options (present, absent, minus, atleast, only,...)\n";
						goto OUTER;
				}
				if ($taxa =~ /^out|^rest|^other/i) {		#We store outgroup conditions in the hash, as "outgroup" does not exists in the species hash of hashes
					$final_searches {$taxa."_".$option} = "outgroup_".$option;       #Store ALL the searches in a hash, keys are taxa and values the columns and options
				} else {
					@columns = obtain_taxa_columns (\$taxa, \%spp, \%MCL_columns_parsed); #Obtain the columns belonging to each taxa
					#print "Columns that will be inspected: @columns\n";
					if (scalar @columns == 0 ){
						print "Taxa $taxa not found in taxa list\n";
						goto OUTER;
					}
					unless (defined $search[1] && $search[1] =~ /pre|all|abs|none|min|but|atl|onl|jus/gi){
						print "Taxa $taxa is missing valid options (present, absent, minus, atleast, only)\n";
						goto OUTER;
					}
					$final_search = join ("_", @columns, $option); 		#Join the columns and the options to store them later in %final_search
					$final_searches {$taxa."_".$option} = $final_search;	#Store ALL the searches in a hash, keys are taxa and values the columns and options
				}	
			}
			
			#Read the hash containing the file with the MCL output parsed to columns, check the search criteria
			my @keys_MCL_columns_parsed = keys %MCL_columns_parsed;			
			#print "Colums to be inspected: @keys_MCL_columns_parsed\n";
			my $HG_counter = 0;
			foreach my $homology_group (@keys_MCL_columns_parsed){ 				#For each line of the MCL parsed file...
				#print "\n", '$homology_group: ', $homology_group, "\n";
				@true_flags = ();							#Empty the flags that indicate if a homology group/line passes all the queries checks
				my @MCL_columns = split (" ", $MCL_columns_parsed{$homology_group});	#Explode the string of the homology group to columns				
				@outgroup = @MCL_columns;						#Storing all columns in @outgroups, later the ingroups columns will be spliced out o this array
				my @keys = keys %final_searches;
				#print 'keys in %final_searches: ', "@keys" , "\n";
				#print "MCL columns_parsed:\n@MCL_columns\n";
				foreach my $query (keys %final_searches) {				#Now go query from query stored in %final_searches
					#print '$query ', $query,"\n";
					unless ($query =~ /outg/i) {					#Leave outgroup check for the end
						my @query = split ("_", $final_searches{$query});	#Explode the string of the columns to check
						#print "\@query: @query ";
						my $condition = pop @query;				#Save the options (present, absent, minus, atleast,...)
						#print "\$condition: $condition \n";
						@query = @query;				
						#print "\@query columns:  @query\n";
						my $total_taxa_in_clade = scalar @query;		#Stores total number of taxa in queried clade, used to later check options
						my $sum = 0;
						#print "Columns contens: ";
						foreach my $query_column (@query) {			#For each column corresponding to the taxa in the query, we extract the value of the column in MCL file
							#print "$MCL_columns[$query_column] ";							
							if ($MCL_columns[$query_column] != 0) {
								$sum++;
							}
							splice (@outgroup, $query_column, 1, "X"); 	#Remove column from outgroup array, to check later the remaining columns
						}
						#print "\nTotal taxa in clade:$total_taxa_in_clade\n";
						#print "Taxa of clade in this MCL: $sum\n";
						#Check if the group of homology matches the conditions requested by user for this clade
						my $check = check_conditions($condition, $total_taxa_in_clade, $sum);
						if ($check eq "wrong") {
							goto OUTER;
						} else {
							push (@true_flags, $check);
						}
					}
				}
				#OUTGROUP CHECK
				foreach (keys %final_searches) {				#Now we check if user asked for outgroup conditions and is stored in %final_searches
					if ($_ =~ /outg/i) {
						#print "OUTGROUP\n";
						#print "\@outgroup:\n@outgroup\n";
						#print "\@MCL_columns:\n@MCL_columns\n";
						my @query = split ("_", $final_searches{$_});	#Explode the string of the columns to check
						my $condition = pop @query;				#Save the options (present, absent, minus, atleast,...)
						#print '$condition in outgroup: ', $condition,"\n";
						@outgroup = grep { $_ ne "X" } @outgroup;
						#print 'MCL contens in outgroup: ', "\n@outgroup","\n";
						my $total_taxa_in_clade = scalar @outgroup;		#Stores total number of taxa in queried clade, used to later check options
						#print 'Total taxa in outgroup: ', $total_taxa_in_clade, "\n";
						my $sum = 0;
						#print "\$MCL_columns[\$query_column] in outgroup:\n";
						foreach (@outgroup) {			#For each column corresponding to the taxa in the query, we extract the value of the column in MCL file
							#print "$_ ";
							if ($_ != 0) {
								$sum++;
							}
						}
						#print 'Outgroup in this MCL: ', $sum, "\n";
						my $check = check_conditions($condition, $total_taxa_in_clade, $sum);
							if ($check eq "wrong") {
								goto OUTER;
							} else {
								push (@true_flags, $check);
							}
					}
				}
				#print "\@true_flags: @true_flags\n";
				#Check if all flags are true, then store the group of homology in %results
				my $flags = 0;
				
				foreach (@true_flags) {
					if ($_ !~ 'true' ) {
						$flags++;
					}					
				}
				if ($flags == 0) {
					push (@good_homology_groups_spp_names, "$homology_group\t$MCL_out{$homology_group}\n");
					push (@good_homology_groups_columns_parsed, "$homology_group\t$MCL_columns_parsed{$homology_group}\n");
					$HG_counter++;
				}				
			}
			@good_homology_groups_spp_names_copy = @good_homology_groups_spp_names;
			print "\nNumber of groups of homology found: $HG_counter \n";
			if ($HG_counter == 0) { goto OUTER; };
			
			#Save 4 different output files
			print "\nDo you want to see results (gene names, MCL groups, etc) and save them in files? Yes/No\n";
			my $save_files = <STDIN>;
			chomp $save_files;
			if ($save_files =~ /y/gi) {
				my $output_filename = join ("_", @arguments);
				$output_filename = "Output/".$output_filename."_$HG_counter\_HGs";
				
				# 1) Save the groups of homology that match the query, with spp names
				#print "\nShowing first sequence names for first groups of homology:\n";
				#foreach (@good_homology_groups_spp_names) {
				#	print substr($_, 0, 160), "\n\n";
				#}
				unless (open (OUTPUT1, ">$output_filename"."_MCL_genes_IDs.out")) {
					print "Can't open file to save";
				}
				print OUTPUT1 @good_homology_groups_spp_names;
				close OUTPUT1;
				
				# 2) Save the columns from MCL parsed file for groups of homology that match the query
				#print "\nShowing columns for the few first groups of homology:\n";
				#print "Species: Crei Ppat Smoe Atri Atha Ehux Bnat Rfil Ttra Falb Spun Amac Scer Sarc Cfra Cow_ Mbre Sros Aque Ocar Mley Pbac Tadh Nvec Adig Hmag Gsal Sjap Sman Egra Emul Hmic Avag Cgig Pfuc Lgig Ctel Hrob Tspi Rcul Cele Bmal Smar Isca Smim Mmar Dpul Znev Tcas Dmel Skow Spur Bflo Cint Csav Bsch Odio Drer Xtro Ggal Acar Hsap\n\n";
				#foreach (@good_homology_groups_columns_parsed) {
				#	print $_;
				#}
				unless (open (OUTPUT2, ">$output_filename"."_MCL_columns_parsed.out")) {
					print "Can't open file to save";
				}
				print OUTPUT2 "HG\tacas\tppal\tddis\tttra\tfalb\tspun\tamac\ttmel\tscer\tcfra\tsros\tmlei\tlgig\tctel\tcele\tdmel\tdrer\txtro\tggal\tacar\thsap\tngru\tbnat\tehux\tgthe\tptet\ttgon\tpine\tptro\ttpse\tcpar\tgsul\tgphl\tcmer\tppur\tbpra\toluc\totau\tmcom\tmpus\tapro\tcvar\tpsp.\tcsub\tczof\tgpec\tvcar\tcrei\tceus\tkfla\tmpol\tppat\tsmoe\tgbil\tpabi\tgmon\tatri\tspol\tzmar\tzmue\tdrot\tashe\tdcat\tpequ\taoff\tpdac\tegui\tcnus\tacom\tobra\topun\toglu\toruf\tomer\tobar\togla\toniv\tosai\tosaj\tlper\tpedu\tbdis\thvul\tatau\ttura\ttaes\tttur\tlpee\tecru\tsita\tzmay\tsbic\tetef\tzjap\tzmat\tzpac\totho\tmiti\tmbal\tmacu\tecal\tmcor\tnnuc\tkfed\trcre\tvvin\tlang\tadur\taipa\tccaj\tpang\tpvul\tvrad\tgmax\tgsoj\tcari\tmtru\ttpra\tfves\tpavi\tpmum\tpper\tpbre\tpcom\tmdom\tzjuj\tmnot\tcmax\tcmos\tclan\tlsic\tcmel\tcsat\tptri\tppru\tlusi\tjcur\tmesc\thbra\trcom\tcfol\tpgra\tegra\tdlon\tabux\tccle\tcsin\tcich\tcemd\tcpap\tthas\tesal\ttpar\tbnap\tbole\tbrap\tsiri\tbvug\tcrub\tcgra\tatha\talyr\ttcac\tccap\tcoli\tdzib\tgarb\tghir\tgrai\tfesc\tdcay\tbvul\tsole\tcqui\tahyp\tcacu\tcsie\tachi\tlsat\tebre\thann\tpgin\tdcar\tccan\tcgig\tfexc\toeur\tbhyg\tmgut\tsind\thimp\tgaur\tugib\tinil\titri\tpinf\tpaxi\tnobt\tnsyl\tntab\tntom\tslyc\tspen\tspim\tstub\tcann\tcbac\tcchi\n";
				print OUTPUT2 @good_homology_groups_columns_parsed;
				close OUTPUT2;
				
				# 3) Now save in return format, one taxa per line
				#@good_homology_groups_spp_names = @good_homology_groups_spp_names_copy;
				my @gene_names_return = ();
				print "\nShowing first sequence names for groups of homology:\n";
				foreach (@good_homology_groups_spp_names) {
					chomp $_;
					my @homology_group = split (" ", $_);
					my $gene_names = parse_gene_names_return_format(@homology_group);
					push (@gene_names_return, $gene_names);
					print "\n",substr($gene_names, 0, 480),"...\n";
				}
				unless (open (OUTPUT3, ">$output_filename"."_MCL_annotated_genes.out")) {
					print "Can't open file to save";
				}
				print OUTPUT3 @gene_names_return;
				close OUTPUT3;
				
				# 4) Save the names of the taxa present in the groups of homology 
				my @taxa_names = ();
				foreach (@good_homology_groups_spp_names) {
					chomp $_;
					my @homology_group = split (" ", $_);
					push (@taxa_names, parse_taxa_names(@homology_group));
				}
				#print "\nShowing the labels of all the taxa for each homology group:\n";
				#foreach (@taxa_names) {
				#	print "$_";
				#	}
				unless (open (OUTPUT4, ">$output_filename"."_taxa_names.out")) {
					print "Can't open file to save";
				}
				print OUTPUT4 @taxa_names;
				close OUTPUT4;
				print "\nResults saved to files $output_filename\n";
			}	
		}
	}
} until ($user_input =~ /^\s*$|^exit|^quit/i);

#Untie hash files
untie %MCL_out;
untie %MCL_columns_parsed;

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
    unless (-e $filename) {print "File $filename does NOT exists\n" and exit;}
    unless (-f $filename) {print "File $filename is NOT a flat file\n" and exit;}
    unless (-s $filename) {print "File $filename is empty\n" and exit;}
    unless (open (FH, $filename)) {print "File $filename can not be opened\n" and exit;}
    close FH;
    return;
}

#MCL_COLUMNS_PARSED_TO_HASH: subrout to parse a MCL output file and store it in a hash, in which the
#keys are the group homology and the values the genes found in that group
sub MCL_columns_parsed_to_hash {
	my ($filename) = @_;
	my %hash;
	my @line = '';
	my $line_counter = 0;
	
	open(FH, $filename);
	my @file = <FH>;
	foreach my $line (@file) {   	#Parse file one line at a time
		$line_counter++;
		chomp $line;
		$line =~ s/^\d*\s//;	#To remove first digits, which indicate the group of homology number/ID
		#$line =~ s/\s/\t/g;
		#$line =~ s/\t$//g;
		my $key = "HG_".$line_counter;
		$hash{$key} = "$line";	
	}
	close FH;
	return %hash;
}

#MCL_TO_HASH: subrout to parse a MCL output file and store it in a hash, in which the
#keys are the group homology and the values the genes found in that group
sub MCL_output_to_hash {
	my ($filename) = @_;
	my %hash;
	my $line_counter = 0;
	
	open(FH, $filename);
	my @file = <FH>;
	foreach my $line (@file) {   			#Parse file one line at a time
		$line_counter++;
		chomp $line;
		#$line =~ s/\s/\t/g;
		#$line =~ s/\t$//g;
		my $key = "HG_".$line_counter;
		$hash{$key} = "$line";	
	}
	close FH;
	return %hash;
}

#WALK_HASH: subroutine to traverse the hash of hashes, modified after in http://www.perlmonks.org/?node_id=116162
sub walk_hash { 
my ($hash, $value_list, $key_list) = @_;
while (my ($key, $value) = each %$hash) {
	push @$key_list, $key;
	if (ref($value) eq 'HASH') {
		walk_hash ($value, $value_list, $key_list);
	} else {
		push @$value_list, $value;
	}
	pop @$key_list;
	}
}

#QUERY: intermediate subroutine that sends the taxa search to the to recoursive subroutines PG_WALK and PRINT_EVERYTHING
sub obtain_taxa_columns { 
	my ($taxa, $spp, $MCL_parsed) = @_; 
	my @results = ();
	
	#Send query and hash      

	pg_walk (\%$spp, [], $$taxa, \@results);
	@results = sort {$a <=> $b} @results;
	return @results; 
}

#PG_WALK: subroutine to traverse the hash of hashes till finding the queried label
#Hat tip to Patrick Gemmell (Oxford) for the help, he modified the subroutine
#found in http://www.perlmonks.org/?node_id=116162
sub pg_walk {
	my ($hash, $key_list, $query_text, $localresults) = @_;
	while (my ($key, $value ) = each %$hash) {
		push @$key_list, $key;
		if ($key =~ /^$query_text/gi) {
			print "Taxa that will be searched: $key\n";
			print_everything($value , $localresults);
		} else {
			if (ref($value ) eq 'HASH') {
				pg_walk($value , $key_list, $query_text, $localresults);
			}
		}
		pop @$key_list;
    }
}

#PRINT_EVERYTHING: subroutine to traverse the hash of hashes to print the values corresponding to the query in sub PG_WALK (see above).
#Hat tip to Patrick Gemmell (Oxford) for the help, he modified the subroutine
#found in http://www.perlmonks.org/?node_id=116162
sub print_everything {
	my ($hash, $localresults, $key_list) = @_;
	while (my ($key, $value) = each %$hash) {
		push @$key_list, $key;
		if (ref($value) eq 'HASH') {
			print_everything($value, $localresults, $key_list);
		} else {
			push @$localresults, $value;
		}
		pop @$key_list;
	}
}

#CHECK_CONDITIONS: subroutine to check the conditions of presence/absence specified by the user, returns TRUE or FALSE
sub check_conditions {
	my ($condition, $total_taxa_in_clade, $sum) = @_;
	#Perform the check according to the different options (present, absent, minus, atleast, only...
	if ($condition =~ /pre|all/gi) {
		if ($sum == $total_taxa_in_clade) {
			#print "True all\n";
			return "true";
		} else {
			#print "False all\n";
			return "false";
		}
	} elsif ($condition =~ /abs|none/gi) {
		if ($sum == 0) {
			#print "True none\n";
			return "true";
		} else {
			#print "False none\n";
			return "false";
		}
	} elsif ($condition =~ /atl\D*/gi ) {
		$condition =~ s/$&//g;
		#print "\$condition: $condition\n";
		if ($condition > $total_taxa_in_clade || $condition == 0) {
			print "\nSpecified \"at least\" value ($condition) is greater than the number of taxa in the clade ($total_taxa_in_clade), or is not not a number, or is 0\n";
			return "wrong";
		}								
		if ($sum >= $condition) {
			#print "True atleast $condition\n";
			return "true";
		} else {
			#print "Atleast $condition false\n";
			return "false";
		}
	} elsif ($condition =~ /min\D*|but\D*/gi ) {
		$condition =~ s/$&//g;							
		#print "\$condition: $condition \n";
		if ($condition > $total_taxa_in_clade || $condition == 0) {
			print "\nSpecified \"minus\" value ($condition) is greater than the number of taxa in the clade ($total_taxa_in_clade), or is not not a number, or is 0\n";
			return "wrong";
		}								
		if ($sum == ($total_taxa_in_clade - $condition)) {
			#print "True minus $condition\n";
			return "true";
		} else {
			#print "Minus $condition false\n";
			return "false";
		}
	} elsif ($condition =~ /onl\D*|jus\D*/gi ) {
		$condition =~ s/$&//g;								
		#print "\$condition: $condition \n";
		if ($condition > $total_taxa_in_clade || $condition == 0) {
			print "\nSpecified \"only\" value ($condition) is greater than the number of taxa in the clade ($total_taxa_in_clade), or is not not a number, or is 0\n";
			return "wrong";
		}								
		if ($sum == $condition) {
			#print "True only $condition\n";
			return "true";
		} else {
			#print "Only $condition false\n";
			return "false";
		}
	}
}

#PARSE_GENE_NAMES_TEXT_FORMAT: subrout to extract the gene names from the annotated genomes for the list of groups of homology
# producing an output in text format
sub parse_gene_names_text_format {
	my (@homology_group) = @_;
	my @list_annotated_genomes = qw /HG_ hsap dmel cele drer cpar kfla vcar ppat smoe atri osai taes mtru atha cann/;
	#@list_annotated_genomes = reverse @list_annotated_genomes;
	my $results;

	foreach my $annotated_taxa (@list_annotated_genomes) {
		#print "\$gene: $gene\n";
		foreach my $gene (@homology_group) {
			#print "\$annotated_taxa: $annotated_taxa ";
			if ($gene =~ /$annotated_taxa/i) {
				#print color ("red"), "Match!\n", color ("reset");
				$results .= "$gene ";
			}
		}
	}
	#print "\@results in sub: $results\n";
	if ($results =~ /^HG_\d*\s*$/) {
		$results .= "This group of homology does not contain any annotated genome";
	}
	#print "Results @results\n";
	#@results = sort @results;
	$results .= "\n";
	return $results;
}

#PARSE_GENE_NAMES_RETURN_FORMAT: subrout to extract the gene names from the annotated genomes for the list of groups of homology
# producing an output in a format with returns
sub parse_gene_names_return_format {
	my (@homology_group) = @_;
	my @list_annotated_genomes = qw /HG_ hsap dmel cele drer cpar kfla vcar ppat smoe atri osai taes mtru atha cann/;
	#@list_annotated_genomes = reverse @list_annotated_genomes;
	my $results;
	my $flag = 0;

	foreach my $annotated_taxa (@list_annotated_genomes) {
		$flag = 0;
		#print "\$gene: $gene\n";
		foreach my $gene (@homology_group) {
			#print "\$annotated_taxa: $annotated_taxa ";
			if ($gene =~ /$annotated_taxa/i) {
				#print color ("red"), "Match!\n", color ("reset");
				$results .= "$gene\t";
				$flag = 1;
			}
		}
		if ($flag == 1) { $results .= "\n"; }
	}
	
	#print "\@results in sub: $results\n";
	if ($results =~ /^HG_\d*\s*$/) {
		$results .= "This group of homology does not contain any annotated genome\n";
	}
	#print "Results @results\n";
	#@results = sort @results;
	#$results .= "\n";
	return $results;
}

#PARSE_TAXA_NAMES: subrout to extract the taxa names from the list of groups of homology
sub parse_taxa_names {
	my (@homology_group) = @_;
	my %taxons;
	my $results;
	my $group_ID = splice (@homology_group, 0 ,1);

	foreach my $taxon (@homology_group) {
		$taxon =~ /^.{4}/;
		$taxons{$&} = '';
	}
	my @keys = sort keys %taxons;
	$results = join ("\t", @keys);
	$results = "$group_ID "."\t$results"."\n";
	return $results;
}

#PRINT_EXAMPLES: subroutine to print some examples to user
sub print_examples {
	print '
Clade/species names can be truncated, but the start of the clade name should match the table printed above.
Search is case insensitive.

Some search examples (first 4 digits in examples stand for rest of taxa, the other 4 for ingroup):	
  "Vertebrata-present" => genes found in ALL vertebrate species, present or absent in other clades/rest of taxa
	Rest of taxa ???? Ingroup 1111
  "Vertebrata-present Rest-absent" => genes found in ALL vertebrate species, absent in other clades/rest of taxa
	Rest of taxa 0000 Ingroup 1111
  "Vertebrata-present Rest-present" => genes found in ALL vertebrate species, present in other clades/rest of taxa
	Rest of taxa 1111 Ingroup 1111
  "Vertebrata-absent Rest-present" => genes found in rest of taxa species, absent in Vertebrata
	Rest of taxa 1111 Ingroup 0000
  "Homo-present Mus-present Rest-absent" => genes only found in humans and mice. Species can be specified one by one.

The number of species presenting/missing for a gene can be fine-tuned with minus#, atleast#, only# for both ingroup and rest of taxa:
  "Vertebrata-minus1" => found in ALL vertebrate species but one, present or absent in other clades/rest of taxa
	Rest of taxa ???? Ingroup 1110 / 1101 / 1011 / 0111
  "Vertebrata-minus2 Rest-minus1" => genes found in ALL vertebrate species but one, absent in other clades/rest of taxa
	Rest of taxa 1110 / 1101 / 1011 / 0111 Ingroup 1100 / 1010 / 1001 / 0110 / 0101 / 0011
  "Vertebrata-atleast1 Rest-atleast1" => genes found in at least 1 vertebrate species and 1 rest of taxa species
	Rest of taxa 1000 / 1100 / 1110 / 1111 / 1010 / 1011 / 1001 / 1101 / 0110 / 0111 etc.
	Ingroup  1000 / 1100 / 1110 / 1111 / 1010 / 1011 / 1001 / 1101 / 0110 / 0111 etc.
  "Vertebrata-only3" => return genes found in just 3 vertebrate species, present or absent in other clades/rest of taxa
	Rest of taxa ???? Ingroup 1110 / 1101 / 1011 / 0111
  
Different criteria can be combined in a single search:
  "Vertebrata-minus1 Echinodermata-atleast2" => genes found in ALL vertebrate species but one, AND present in at least two echinoderms, absent/present in other clades/rest of taxa
	Rest of taxa ???? Vertebrata 1110 / 1101 / 1011 / 0111 Echinodermata 1100 / 1010 / 1001 / 0110 / 0101 / 0011
  "Vertebrata-atleast2 Urochordata-atleast2" => genes found in 2 or more vertebrate species OR 2 or more urochordates, independently if they are present/absent in other clades/rest of taxa
	Rest of taxa ???? Vertebrata 1100 / 1010 / 1001 / 0110 / 0101 / 0011 Urochordata 1100 / 1010 / 1001 / 0110 / 0101 / 0011
  "Nematoda-absent Platyhelminthes-absent Rest-present" => genes found in clades/rest of taxa, absent (convergently lost) in round worms and flatworms
	Rest of taxa 1111 Nematoda 0000 Platyhelminthes 0000
	
Carefull with nested taxa!!! Start with the greater group taking into account the conditions for the smaller group:
  To find genes in ALL chordates but missing only in humans => "Chordata-minus1 Hsap-absent"
  To find genes in ALL chordates but missing only in vertebrates => "Chordata-minus5 Vertebrata-absent"
  To find genes in at least one clade of chordates, but missing only in vertebrates => "Cephalocordata-atleast1 Urochordata-atleast1 Vertebrata-absent"
  ';
	return;
}

#HASH_SPP: a subroutine to define define the hash of hashes containing of all the clades and spp included, put here to not clutter the main program.
#Each species is assigned a numeric value, same as the column they occupy in the parsed MCL output, so Crei is the first column and Hsap occupies the last column.
#Thus, when user asks for a group, these values can be used as index to lookup lines/arrays.
#Each element shoud have the same number of levels, or the subrout "PRINT_HASH_COLORS" won't work.
sub hash_spp {
	my (%spp) = ();

##   All spp, one by one
##     DOMAIN  ## SUBDOMAIN ## SUPERGROUP     ##  KINGDOM      ## SUBKINGDOM1 ## SUBKINGDOM2 ##  SUBKINGDOM3 ## CLADE1      ##  SUBCLADE1      ##  CLADE2        ##  SUBCLADE2_1 ## SUBCLADE2_2 ## CLADE3 ## SUBCLADE3 ## ORDER     ## FAMILY       ##  SPECIES
##   {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Malpighiales'}{'Euphorbiaceae'}{'Hevea_brasiliensis(hbra)'}

# Amorphea
## Apusozoa
## Unikonta
### Non-holozoan_Amorphea
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Amoebozoa'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Longamoebia'}{'Acanthamoebidae'}{'Acanthamoeba_castellanii_(acas)'}{''} = 0;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Amoebozoa'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Dictyosteliidae'}{'Dictyosteliidae'}{'Polysphondylium_pallidum_(ppal)'}{''} = 1;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Amoebozoa'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Dictyosteliidae'}{'Dictyosteliidae'}{'Dictyostelium_discoideum_(ddis)'}{''} = 2;
$spp {'Eukaryota'}{'Amorphea'}{'Apusozoa'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Thecomonadea'}{'Thecomonadea'}{'Apusomonadidae'}{'Thecamonas_trahens_(ttra)'}{''} = 3;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Holomycota/Nucletmycea'}{''}{''}{''}{''}{''}{''}{''}{'Holomycota/Nucletmycea'}{'Fonticulida'}{'Fonticulida'}{'Fonticulida'}{'Fonticula_alba_(falb)'}{''} = 4;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Holomycota/Nucletmycea'}{''}{''}{''}{''}{''}{''}{''}{'Chytridiomycota'}{'Chytridiomycetes'}{'Spizellomycetales'}{'Spizellomycetaceae'}{'Spizellomyces_punctatus_(spun)'}{''} = 5;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Holomycota/Nucletmycea'}{''}{''}{''}{''}{''}{''}{''}{'Blastocladiomycota'}{'Blastocladiomycetes'}{'Blastocladiales'}{'Blastocladiaceae'}{'Allomyces_macrogynus_(amac)'}{''} = 6;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Holomycota/Nucletmycea'}{''}{'Dikarya'}{'Ascomycota'}{''}{''}{''}{'Pezizomycotina'}{'Pezizomycotina'}{'Pezizomycetes'}{'Pezizales'}{'Tuberaceae'}{'Tuber_melanosporum_(tmel)'}{''} = 7;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Holomycota/Nucletmycea'}{''}{'Dikarya'}{'Ascomycota'}{''}{''}{''}{'Saccharomycotina'}{'Saccharomycotina'}{'Saccharomycetes'}{'Saccharomycetales'}{'Saccharomycetaceae'}{'Saccharomyces_cerevisiae_(scer)'}{''} = 8;

#####  Non-metazoan_Holozoa
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Teretosporea'}{''}{''}{''}{''}{''}{''}{''}{''}{'Ichthyosporea'}{'Creolimacidae'}{'Creolimacidae'}{'Creolimax_fragrantissima_(cfra)'}{''} = 9;

### Non-holozoan metazoa
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Choanoflagellata'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Salpingoecidae'}{'Salpingoecidae'}{'Salpingoeca_rosetta_(sros)'}{''} = 10;

### Non-bilaterian metazoa
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Basal2'}{'Ctenophora'}{'Tentaculata'}{''}{''}{''}{''}{''}{''}{''}{'Lobata'}{'Bolinopsidae'}{'Mnemiopsis_leidyi_(mlei)'}{''} = 11;

##### Lophotrochozoa
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Basal2'}{'Bilateria'}{'Protostomia'}{'Lophotrochozoa'}{'Lopho1'}{'Lopho2'}{'Mollusca'}{''}{''}{'Gastropoda'}{'Lottiidae'}{'Lottiidae'}{'Lottia_gigantea_(lgig)'}{''} = 12;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Basal2'}{'Bilateria'}{'Protostomia'}{'Lophotrochozoa'}{'Lopho1'}{'Lopho2'}{'Lopho3'}{'Annelida'}{''}{'Polychaeta'}{'Capitellida'}{'Capitellidae'}{'Capitella_teleta_(ctel)'}{''} = 13;

#####  Ecdysozoa
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Basal2'}{'Bilateria'}{'Protostomia'}{'Ecdysozoa'}{'Nematoida'}{'Nematoda'}{''}{''}{''}{'Secernentea '}{'Rhabditidae'}{'Rhabditidae'}{'Caenorhabditis_elegans_(cele)'}{''} = 14;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Basal2'}{'Bilateria'}{'Protostomia'}{'Ecdysozoa'}{'Panarthropoda'}{'Arthropoda'}{'Mandibulata'}{'Mandib1'}{'Hexapoda'}{'Insecta'}{'Diptera'}{'Drosophilidae'}{'Drosophila_melanogaster_(dmel)'}{''} = 15;

#####  Deuterostomes
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Basal2'}{'Bilateria'}{'Deuterostomia'}{'Chordata'}{'Olfactores'}{'Vertebrata'}{'Teleost'}{''}{''}{'Actinopterygii'}{'Cypriniformes'}{'Cyprinidae'}{'Danio_rerio_(drer)'}{''} = 16;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Basal2'}{'Bilateria'}{'Deuterostomia'}{'Chordata'}{'Olfactores'}{'Vertebrata'}{'Tetrapoda'}{''}{''}{'Batrachia'}{'Anura'}{'Pipidae'}{'Xenopus_tropicalis_(xtro)'}{''} = 17;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Basal2'}{'Bilateria'}{'Deuterostomia'}{'Chordata'}{'Olfactores'}{'Vertebrata'}{'Tetrapoda'}{'Amniota'}{'Sauropsida'}{'Tetrapoda'}{'Galliformes '}{'Phasianidae'}{'Gallus_gallus_(ggal)'}{''} = 18;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Basal2'}{'Bilateria'}{'Deuterostomia'}{'Chordata'}{'Olfactores'}{'Vertebrata'}{'Tetrapoda'}{'Amniota'}{'Sauropsida'}{'Reptilia'}{'Squamata'}{'Polychrotidae'}{'Anolis_carolinensis_(acar)'}{''} = 19;
$spp {'Eukaryota'}{'Amorphea'}{'Unikonta'}{'Opisthokonta'}{'Basal2'}{'Bilateria'}{'Deuterostomia'}{'Chordata'}{'Olfactores'}{'Vertebrata'}{'Tetrapoda'}{'Amniota'}{'Mammalia'}{'Mammalia'}{'Primates'}{'Hominidae'}{'Homo_sapiens_(hsap)'}{''} = 20;

# Bikonta
## Non-Archaeplastida Bikonta (SAR and Excavata)
$spp {'Eukaryota'}{'Bikonta'}{'Excavata'}{'Amoebida'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Vahlkampfidae'}{'Vahlkampfidae'}{'Naegleria_gruberi_(ngru)'}{''} = 21;
$spp {'Eukaryota'}{'Bikonta'}{'SAR'}{'Rhizaria'}{'Cercozoa'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Chlorarachniophyceae'}{'Bigelowiella_natans_(bnat)'}{''} = 22;
$spp {'Eukaryota'}{'Bikonta'}{'SAR'}{'Chromoalveolata'}{'Haptophyta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Isochrysidales'}{'Noelaerhabdaceae'}{'Emiliana_huxleyii_(ehux)'}{''} = 23;
$spp {'Eukaryota'}{'Bikonta'}{'SAR'}{'Chromoalveolata'}{'Cryptopyhta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Pyrenomonadales'}{'Geminigeraceae'}{'Guillardia_theta_(gthe)'}{''} = 24;
$spp {'Eukaryota'}{'Bikonta'}{'SAR'}{'Chromoalveolata'}{'Alveolata'}{''}{''}{''}{''}{''}{''}{''}{''}{'Oligohymenophorea'}{'Peniculida'}{'Parameciidae'}{'Paramecium_tetraurelia_(ptet)'}{''} = 25;
$spp {'Eukaryota'}{'Bikonta'}{'SAR'}{'Chromoalveolata'}{'Alveolata'}{''}{''}{''}{''}{''}{''}{''}{''}{'Conoidasida'}{'Eucoccidiorida'}{'Sarcocystidae'}{'Toxoplasma_gondii_(tgon)'}{''} = 26;
$spp {'Eukaryota'}{'Bikonta'}{'SAR'}{'Stramenopiles'}{'Oomycota'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Peronosporales'}{'Peronosporaceae'}{'Phytophthora_infestans_(pine)'}{''} = 27;
$spp {'Eukaryota'}{'Bikonta'}{'SAR'}{'Stramenopiles'}{'Bacillariophyta'}{''}{''}{''}{''}{''}{''}{''}{''}{'Bacillariophyceae'}{'Naviculales'}{'Phaeodactylaceae'}{'Phaeodactylum_tricornutum_(ptro)'}{''} = 28;
$spp {'Eukaryota'}{'Bikonta'}{'SAR'}{'Stramenopiles'}{'Bacillariophyta'}{''}{''}{''}{''}{''}{''}{''}{''}{'Coscinodiscophyceae'}{'Thalassiosirales'}{'Thalassiosiraceae'}{'Thalassiosira_pseudonana_(tpse)'}{''} = 29;

## Archaeplastida
### Non-Viridiplante Archaeplastida
##### Glaucophyta
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{''}{'Glaucophyta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Glaucophyceae'}{'Glaucocystales'}{'Glaucocystaceae'}{'Cyanophora_paradoxa_(cpar)'}{''} = 30;

##### Rhodophyta
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Rhodophyta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Cyanidiophyceae'}{'Cyanidiales'}{'Cyanidiaceae'}{'Galdieria_sulphuraria_(gsul)'}{''} = 31;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Rhodophyta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Cyanidiophyceae'}{'Cyanidiales'}{'Cyanidiaceae'}{'Galdieria_phlegrea_(gphl)'}{''} = 32;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Rhodophyta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Cyanidiophyceae'}{'Cyanidiales'}{'Cyanidiaceae'}{'Cyanidioschyzon_merolae_(cmer)'}{''} = 33;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Rhodophyta'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Bangiophyceae'}{'Porphyridiales'}{'Porphyridiaceae'}{'Porphyridium_purpureum_(ppur)'}{''} = 34;

### Viridiplantae
#####  Chlorophyta
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Prasinophytes'}{'Mamiellophyceae'}{''}{''}{''}{''}{''}{''}{''}{'Mamiellales'}{'Bathyococcaceae'}{'Bathycoccus_prasinos_(bpra)'}{''} = 35;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Prasinophytes'}{'Mamiellophyceae'}{''}{''}{''}{''}{''}{''}{''}{'Mamiellales'}{'Bathyococcaceae'}{'Ostreococcus_lucimarinus_(oluc)'}{''} = 36;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Prasinophytes'}{'Mamiellophyceae'}{''}{''}{''}{''}{''}{''}{''}{'Mamiellales'}{'Bathyococcaceae'}{'Ostreococcus_tauri_(otau)'}{''} = 37;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Prasinophytes'}{'Mamiellophyceae'}{''}{''}{''}{''}{''}{''}{''}{'Mamiellales'}{'Mamiellaceae'}{'Micromonas_commoda_(mcom)'}{''} = 38;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Prasinophytes'}{'Mamiellophyceae'}{''}{''}{''}{''}{''}{''}{''}{'Mamiellales'}{'Mamiellaceae'}{'Micromonas_pusilla_(mpus)'}{''} = 39;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Trebouxiophyceae'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chlorellales'}{'Chlorellaceae'}{'Auxenochlorella_protothecoides_(apro)'}{''} = 40;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Trebouxiophyceae'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chlorellales'}{'Chlorellaceae'}{'Chlorella_variabilis_(cvar)'}{''} = 41;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Trebouxiophyceae'}{''}{''}{''}{''}{''}{''}{''}{''}{'Trebouxiophyceae_incertae_sedis'}{'Coccomyxaceae'}{'Picochlorum_sp._(psp.)'}{''} = 42;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Trebouxiophyceae'}{''}{''}{''}{''}{''}{''}{''}{''}{'Trebouxiophyceae_incertae_sedis'}{'Coccomyxaceae'}{'Coccomyxa_subellipsoidea_(csub)'}{''} = 43;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Chlorophyceae'}{''}{''}{''}{''}{''}{''}{''}{''}{'Sphaeropleales'}{'Chromochloridaceae'}{'Chromochloris_zofingiensis_(czof)'}{''} = 44;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Chlorophyceae'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chlamydomonadales'}{'Volvocaceae'}{'Gonium_pectorale_(gpec)'}{''} = 45;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Chlorophyceae'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chlamydomonadales'}{'Volvocaceae'}{'Volvox_carteri_(vcar)'}{''} = 46;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Chlorophyceae'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chlamydomonadales'}{'Chlamydomonadaceae'}{'Chlamydomonas_reinhardtii_(crei)'}{''} = 47;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Chlorophyta'}{'Chlorophyceae'}{''}{''}{''}{''}{''}{''}{''}{''}{'Chlamydomonadales'}{'Chlamydomonadaceae'}{'Chlamydomonas_eustigma_(ceus)'}{''} = 48;

##### Charophyta
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Klebsormidiophyceae'}{''}{''}{''}{''}{''}{''}{''}{''}{''}{'Klebsormidiales'}{'Klebsormidiaceae'}{'Klebsormidium_flaccidum_(kfla)'}{''} = 49;

#####  Bryophytes (Liverworts/ mosses/ hornworts) 
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Bryophytes'}{'Marchantiophyta'}{'Marchantiopsida'}{'Marchantiidae'}{''}{''}{''}{''}{''}{'Marchantiales'}{'Marchantiaceae'}{'Marchantia_polymorpha_(mpol)'}{''} = 50; 
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Bryophytes'}{'Bryophyta'}{'Bryophytina'}{'Bryopsida'}{'Funariidae'}{''}{''}{''}{''}{'Funariales'}{'Funariaceae'}{'Physcomitrella_patens_(ppat)'}{''} = 51;

#####  Lycophytes
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Lycopodiopsida'}{''}{''}{''}{''}{''}{''}{''}{'Selaginellales'}{'Selaginellaceae'}{'Selaginella_moellendorffii_(smoe)'}{''} = 52;

#####  Gymnosperms
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Acrogymnospermae'}{'Ginkgoidae'}{''}{''}{''}{''}{''}{'Ginkgoales'}{'Ginkgoaceae'}{'Ginkgo_biloba_(gbil)'}{''} = 53;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Acrogymnospermae'}{'Pinidae'}{''}{''}{''}{''}{''}{'Pinales'}{'Pinaceae'}{'Picea_abies_(pabi)'}{''} = 54;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Acrogymnospermae'}{'Gnetidae'}{''}{''}{''}{''}{''}{'Gnetales'}{'Gnetaceae'}{'Gnetum_monatum_(gmon)'}{''} = 55;

##### Angiosperms
#######  Basal angiosperms
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Basal_Magnoliophyta'}{''}{''}{''}{''}{''}{'Amborellales'}{'Amborellaceae'}{'Amborella_trichopoda_(atri)'}{''} = 56;

#######  Mesangiospermae
######### Monocots (Liliopsida)
########### Non-PetrosaviidaeLiliopsida
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{''}{''}{''}{''}{'Alismatales'}{'Araceae'}{'Spirodela_polyrhiza_(spol)'}{''} = 57;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{''}{''}{''}{''}{'Alismatales'}{'Zosteraceae'}{'Zostera_marina_(zmar)'}{''} = 58;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{''}{''}{''}{''}{'Alismatales'}{'Zosteraceae'}{'Zostera_muelleri_(zmue)'}{''} = 59;

########### Non-Commelinids Petrosaviidae
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{''}{''}{''}{'Dioscoreales'}{'Dioscoreaceae'}{'Dioscorea_rotundata_(drot)'}{''} = 60;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{''}{''}{''}{'Asparagales'}{'Orchidaceae'}{'Apostasia_shenzhenica_(ashe)'}{''} = 61;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{''}{''}{''}{'Asparagales'}{'Orchidaceae'}{'Dendrobium_catenatum_(dcat)'}{''} = 62;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{''}{''}{''}{'Asparagales'}{'Orchidaceae'}{'Phalaenopsis_equestris_(pequ)'}{''} = 63;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{''}{''}{''}{'Asparagales'}{'Asparagaceae'}{'Asparagus_officinalis_(aoff)'}{''} = 64;

########### Arecales
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Arecales'}{'Arecaceae'}{'Phoenix_dactylifera_(pdac)'}{''} = 65;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Arecales'}{'Arecaceae'}{'Elaeis_guineensis_(egui)'}{''} = 66;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Arecales'}{'Arecaceae'}{'Cocos_nucifera_(cnus)'}{''} = 67;

########### Poales
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Bromeliaceae'}{'Ananas_comosus_(acom)'}{''} = 68;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oryza_brachyantha_(obra)'}{''} = 69;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oryza_punctata_(opun)'}{''} = 70;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oryza_glumipatula_(oglu)'}{''} = 71;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oryza_rufipogon_(oruf)'}{''} = 72;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oryza_meridionalis_(omer)'}{''} = 73;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oryza_barthii_(obar)'}{''} = 74;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oryza_glaberrima_(ogla)'}{''} = 75;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oryza_nivara_(oniv)'}{''} = 76;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oryza_sativa_Indica_(osai)'}{''} = 77;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oryza_sativa_Japonica_(osaj)'}{''} = 78;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Leersia_perrieri_(lper)'}{''} = 79;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Phyllostachys_edulis_(pedu)'}{''} = 80;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Brachypodium_distachyon_(bdis)'}{''} = 81;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Hordeum_vulgare_(hvul)'}{''} = 82;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Aegilops_tauschii_(atau)'}{''} = 83;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Triticum_urartu_(tura)'}{''} = 84;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Triticum_aestivum_(taes)'}{''} = 85;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Triticum_turgidium_(ttur)'}{''} = 86;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Lolium_perenne_(lpee)'}{''} = 87;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Echinochloa_crus-galli_(ecru)'}{''} = 88;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Setaria_italica_(sita)'}{''} = 89;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Zea_mays_(zmay)'}{''} = 90;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Sorghum_bicolor_(sbic)'}{''} = 91;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Eragrostis_tef_(etef)'}{''} = 92;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Zoysia_japonica_(zjap)'}{''} = 93;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Zoysia_matrella_(zmat)'}{''} = 94;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Zoysia_pacifica_(zpac)'}{''} = 95;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Poales'}{'Poaceae'}{'Oropetium_thomaeum_(otho)'}{''} = 96;

########## Zingiberales
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Zingiberales'}{'Musaceae'}{'Musa_itinerans_(miti)'}{''} = 97;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Zingiberales'}{'Musaceae'}{'Musa_balbisiana_(mbal)'}{''} = 98;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Liliopsida'}{'Petrosaviidae'}{'Commelinids'}{''}{''}{'Zingiberales'}{'Musaceae'}{'Musa_acuminata_(macu)'}{''} = 99;

########## Early diverging eudicots
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{''}{''}{''}{'BaseAng'}{'Proteales'}{'Nelumbonaceae'}{'Nelumbo_nucifera_(nnuc)'}{''} = 100;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{''}{''}{''}{'BaseAng'}{'Ranunculales'}{'Papaveraceae'}{'Eschscholzia_californica_(ecal)'}{''} = 101;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{''}{''}{''}{'BaseAng'}{'Ranunculales'}{'Papaveraceae'}{'Macleaya_cordata_(mcor)'}{''} = 102;

########## Gunneridae
#################  Non Rosids/Asterids Pentapetalae
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{''}{''}{'Saxifragales'}{'Crassulaceae'}{'Kalanchoe_fedtschenkoi_(kfed)'}{''} = 103;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{''}{''}{'Saxifragales'}{'Crassulaceae'}{'Rhodiola_crenulata_(rcre)'}{''} = 104;

##########  Rosids
#################  Rosids incertae sedis
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Rosids_incertae_sedis'}{'Vitales'}{'Vitaceae'}{'Vitis_vinifera_(vvin)'}{''} = 105;

#################  Fabids 
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Lupinus_angustifolius_(lang)'}{''} = 106;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Arachis_duranensis_(adur)'}{''} = 107;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Arachis_ipaensis_(aipa)'}{''} = 108;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Cajanus_cajan_(ccaj)'}{''} = 109;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Phaseolus_angularis_(pang)'}{''} = 110;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Phaseolus_vulgaris_(pvul)'}{''} = 111;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Vigna_radiata_(vrad)'}{''} = 112;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Glycine_max_(gmax)'}{''} = 113;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Glycine_soja_(gsoj)'}{''} = 114;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Cicer_arietinum_(cari)'}{''} = 115;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Medicago_truncatula_(mtru)'}{''} = 116;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Fabales'}{'Fabaceae'}{'Trifolium_pratense_(tpra)'}{''} = 117;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Rosales'}{'Rosaceae'}{'Fragaria_vesca_(fves)'}{''} = 118;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Rosales'}{'Rosaceae'}{'Prunus_avium_(pavi)'}{''} = 119;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Rosales'}{'Rosaceae'}{'Prunus_mume_(pmum)'}{''} = 120;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Rosales'}{'Rosaceae'}{'Prunus_persica_(pper)'}{''} = 121;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Rosales'}{'Rosaceae'}{'Pyrus_bretschneideri_(pbre)'}{''} = 122;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Rosales'}{'Rosaceae'}{'Pyrus_communis_(pcom)'}{''} = 123;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Rosales'}{'Rosaceae'}{'Malus_domestica_(mdom)'}{''} = 124;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Rosales'}{'Rhamnaceae'}{'Ziziphus_jujuba_(zjuj)'}{''} = 125;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Rosales'}{'Moraceae'}{'Morus_notabilis_(mnot)'}{''} = 126;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Cucurbitales'}{'Cucurbitaceae'}{'Cucurbita_maxima_(cmax)'}{''} = 127;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Cucurbitales'}{'Cucurbitaceae'}{'Cucurbita_moschata_(cmos)'}{''} = 128;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Cucurbitales'}{'Cucurbitaceae'}{'Citrullus_lanatus_(clan)'}{''} = 129;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Cucurbitales'}{'Cucurbitaceae'}{'Lagenaria_siceraria_(lsic)'}{''} = 130;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Cucurbitales'}{'Cucurbitaceae'}{'Cucumis_melo_(cmel)'}{''} = 131;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Cucurbitales'}{'Cucurbitaceae'}{'Cucumis_sativus_(csat)'}{''} = 132;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Malpighiales'}{'Saliaceae'}{'Populus_trichocarpa_(ptri)'}{''} = 133;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Malpighiales'}{'Saliaceae'}{'Populus_pruinosa_(ppru)'}{''} = 134;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Malpighiales'}{'Linaceae'}{'Linum_usitatissimum_(lusi)'}{''} = 135;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Malpighiales'}{'Euphorbiaceae'}{'Jatropha_curcas_(jcur)'}{''} = 136;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Malpighiales'}{'Euphorbiaceae'}{'Manihot_esculenta_(mesc)'}{''} = 137;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Malpighiales'}{'Euphorbiaceae'}{'Hevea_brasiliensis_(hbra)'}{''} = 138;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Malpighiales'}{'Euphorbiaceae'}{'Ricinus_communis_(rcom)'}{''} = 139;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Fabids'}{'Oxalidales'}{'Cephalotaceae'}{'Cephalotus_follicularis_(cfol)'}{''} = 140;

#################  Malvids 
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Myrtales'}{'Lythraceae'}{'Punica_granatum_(pgra)'}{''} = 141;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Myrtales'}{'Myrtaceae'}{'Eucalyptus_grandis_(egra)'}{''} = 142;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Sapindales'}{'Sapindaceae'}{'Dimocarpus_longan_(dlon)'}{''} = 143;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Sapindales'}{'Rutaceae'}{'Atalantia_buxifolia_(abux)'}{''} = 144;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Sapindales'}{'Rutaceae'}{'Citrus_clementina_(ccle)'}{''} = 145;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Sapindales'}{'Rutaceae'}{'Citrus_sinensis_(csin)'}{''} = 146;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Sapindales'}{'Rutaceae'}{'Citrus_ichangensis_(cich)'}{''} = 147;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Sapindales'}{'Rutaceae'}{'Citrus_medica_(cmed)'}{''} = 148;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Caricaceae'}{'Carica_papaya_(cpap)'}{''} = 149;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Cleomaceae'}{'Tarenaya_hassleriana_(thas)'}{''} = 150;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Eutrema_salsugineum_(esal)'}{''} = 151;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Thellungiella_parvula_(tpar)'}{''} = 152;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Brassica_napus_(bnap)'}{''} = 153;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Brassica_oleracea_(bole)'}{''} = 154;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Brassica_rapa_(brap)'}{''} = 155;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Sisymbrium_irio_(siri)'}{''} = 156;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Barbarea_vulgaris_(bvug)'}{''} = 157;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Capsella_rubella_(crub)'}{''} = 158;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Capsella_grandiflora_(cgra)'}{''} = 159;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Arabidopsis_thaliana_(atha)'}{''} = 160;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Brassicales'}{'Brassicaceae'}{'Arabidopsis_lyrata_(alyr)'}{''} = 161;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Malvales'}{'Malvaceae'}{'Theobroma_cacao_(tcac)'}{''} = 162;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Malvales'}{'Malvaceae'}{'Corchorus_capsularis_(ccap)'}{''} = 163;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Malvales'}{'Malvaceae'}{'Corchorus_olitorius_(coli)'}{''} = 164;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Malvales'}{'Malvaceae'}{'Durio_zibethinus_(dzib)'}{''} = 165;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Malvales'}{'Malvaceae'}{'Gossypium_arboreum_(garb)'}{''} = 166;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Malvales'}{'Malvaceae'}{'Gossypium_hirsutum_(ghir)'}{''} = 167;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Rosids'}{'Malvids'}{'Malvales'}{'Malvaceae'}{'Gossypium_raimondii_(grai)'}{''} = 168;

#################  Non-Rosid/Asterid Pentapetalae
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{''}{''}{'Caryophylalles'}{'Polygonaceae'}{'Fagopyrum_esculentum_(fesc)'}{''} = 169;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{''}{''}{'Caryophylalles'}{'Caryophyllaceae'}{'Dianthus_caryophyllus_(dcay)'}{''} = 170;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{''}{''}{'Caryophylalles'}{'Chenopodiaceae'}{'Beta_vulgaris_(bvul)'}{''} = 171;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{''}{''}{'Caryophylalles'}{'Chenopodiaceae'}{'Spinacia_oleracea_(sole)'}{''} = 172;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{''}{''}{'Caryophylalles'}{'Chenopodiaceae'}{'Chenopodium_quinoa_(cqui)'}{''} = 173;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{''}{''}{'Caryophylalles'}{'Amaranthaceae'}{'Amaranthus_hypochondriacus_(ahyp)'}{''} = 174;

#################  Non-Campulids/Lamiids Asterids 
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{''}{'Cornales'}{'Nyssaceae'}{'Campthotheca_acuminata_(cacu)'}{''} = 175;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{''}{'Ericales'}{'Theaceae'}{'Camellia_sinensis_(csie)'}{''} = 176;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{''}{'Ericales'}{'Actinidiaceae'}{'Actinidia_chinensis_(achi)'}{''} = 177;

#################  Campanulids 
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Campanulids'}{'Asterales'}{'Asteraceae'}{'Lactuca_sativa_(lsat)'}{''} = 178;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Campanulids'}{'Asterales'}{'Asteraceae'}{'Erigeron_breviscapus_(ebre)'}{''} = 179;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Campanulids'}{'Asterales'}{'Asteraceae'}{'Helianthus_annuus_(hann)'}{''} = 180;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Campanulids'}{'Apiales'}{'Araliaceae'}{'Panax_ginseng_(pgin)'}{''} = 181;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Campanulids'}{'Apiales'}{'Apiaceae'}{'Daucus_carota_(dcar)'}{''} = 182;

#################  Lamiids
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Gentiales'}{'Rubiaceae'}{'Coffea_canephora_(ccan)'}{''} = 183;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Gentiales'}{'Apocynaceae'}{'Calotropis_gigantea_(cgig)'}{''} = 184;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Lamiales'}{'Oleaceae'}{'Fraxinus_excelsior_(fexc)'}{''} = 185;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Lamiales'}{'Oleaceae'}{'Olea_europaea_(oeur)'}{''} = 186;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Lamiales'}{'Gesneriaceae'}{'Boea_hygrometrica_(bhyg)'}{''} = 187;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Lamiales'}{'Phyrmaceae'}{'Mimulus_guttatus_(mgut)'}{''} = 188;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Lamiales'}{'Pedaliaceae'}{'Sesamum_indicum_(sind)'}{''} = 189;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Lamiales'}{'Bignoniaceae'}{'Handroanthis_impetiginosus_(himp)'}{''} = 190;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Lamiales'}{'Lentibulariaceae'}{'Genlisea_aurea_(gaur)'}{''} = 191;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Lamiales'}{'Lentibulariaceae'}{'Utricularia_gibba_(ugib)'}{''} = 192;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Convolvulaceae'}{'Ipomoea_nil_(inil)'}{''} = 193;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Convolvulaceae'}{'Ipomoea_trifida_(itri)'}{''} = 194;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Petunia_inflata_(pinf)'}{''} = 195;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Petunia_axillaris_(paxi)'}{''} = 196;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Nicotiana_obtusifolia_(nobt)'}{''} = 197;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Nicotiana_sylvestris_(nsyl)'}{''} = 198;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Nicotiana_tabacum_(ntab)'}{''} = 199;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Nicotiana_tomentosiformis_(ntom)'}{''} = 200;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Solanum_lycopersicum_(slyc)'}{''} = 201;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Solanum_pennellii_(spen)'}{''} = 202;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Solanum_pimpinenllifolium_(spim)'}{''} = 203;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Solanum_tuberosum_(stub)'}{''} = 204;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Capsicum_annuum_(cann)'}{''} = 205;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Capsicum_baccatum_(cbac)'}{''} = 206;
$spp {'Eukaryota'}{'Bikonta'}{'Archaeplastida'}{'RhodoVirid'}{'Viridiplantae'}{'Streptophyta'}{'Embryophyta'}{'Tracheophyta'}{'Euphyllophyta'}{'Angiosperm'}{'Mesangiospermae'}{'Eudicotyledons'}{'Gunneridae'}{'Pentapetalae'}{'Asterids'}{'Lamiids'}{'Solanales'}{'Solanaceae'}{'Capsicum_chinense_(cchi)'}{''} = 207;
##     DOMAIN  ## SUBDOMAIN ## SUPERGROUP     ##  KINGDOM      ## SUBKINGDOM1 ## SUBKINGDOM2 ##  SUBKINGDOM3 ## CLADE1      ##  SUBCLADE1      ##  CLADE2        ##  SUBCLADE2_1 ## SUBCLADE2_2 ## CLADE3 ## SUBCLADE3 ## ORDER     ## FAMILY       ##  SPECIES ##

print_hash_colors (%spp);
return %spp;
}

#PRINT_HASH_COLORS: subroutine to traverse the hash of hashes and print the contents, after http://www.perlmonks.org/?node_id=116162
sub print_hash_colors {
    my (%spp) = @_; 
#Define variables (taxonomic ranks)for the hash of hashes
my $domain;
my $subdomain;
my $supergroup;
my $kingdom;
my $subkingdom1;
my $subkingdom2;
my $subkingdom3;
my $clade1;
my $subclade1;
my $clade2;
my $subclade2_1;
my $subclade2_2;
my $clade3;
my $subclade3;
my $order;
my $family;
my $species;

#Print all the keys of the hash of hashes with a nested FOR loop
##     DOMAIN  ## SUBDOMAIN ## SUPERGROUP     ##  KINGDOM      ## SUBKINGDOM1 ## SUBKINGDOM2 ##  SUBKINGDOM3 ## CLADE1      ##  SUBCLADE1      ##  CLADE2        ##  SUBCLADE2_1 ## SUBCLADE2_2 ## CLADE3 ## SUBCLADE3 ## ORDER     ## FAMILY       ##  SPECIES ##
print "\nTree:";
for $domain (keys %spp) {
 print "\n$domain\n";
  for $subdomain (keys %{ $spp{$domain} }) {
   unless ($subdomain eq '') {print color ("blue"),"      $subdomain\n", color ("reset");}
    for $supergroup (keys %{ $spp{$domain}{$subdomain} }){
     unless ($supergroup eq '') {print color ("yellow"),"      $supergroup\n", color ("reset");}
      for $kingdom (keys %{ $spp{$domain}{$subdomain}{$supergroup} }){
       unless ($kingdom eq '') {print color ("green"),"      $kingdom\n", color ("reset");}
        for $subkingdom1 (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom} }){
         unless ($subkingdom1 eq '') {print color ("magenta"),"      $subkingdom1\n", color ("reset");}
          for $subkingdom2 (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1} }){
           unless ($subkingdom2 eq '') {print color ("cyan"),"      $subkingdom2\n", color ("reset");}
            for $subkingdom3 (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2} }){
             unless ($subkingdom3 eq '') {print color ("red"),"      $subkingdom3\n", color ("reset");}
              for $clade1 (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3} }){
               unless ($clade1 eq '') {print color ("blue"),"      $clade1\n", color ("reset");}
                for $subclade1 (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3}{$clade1} }){
                 unless ($subclade1 eq '') {print color ("yellow"),"      $subclade1\n", color ("reset");}
                  for $clade2 (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3}{$clade1}{$subclade1} }){
                   unless ($clade2 eq '') {print color ("green"),"      $clade2\n", color ("reset");}
                    for $subclade2_1 (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3}{$clade1}{$subclade1}{$clade2} }){
                     unless ($subclade2_1 eq '') {print color ("magenta"),"      $subclade2_1\n", color ("reset");}
                      for $subclade2_2 (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3}{$clade1}{$subclade1}{$clade2}{$subclade2_1} }){
                       unless ($subclade2_2 eq '') {print color ("cyan"),"      $subclade2_2\n", color ("reset");}
                        for $clade3 (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3}{$clade1}{$subclade1}{$clade2}{$subclade2_1}{$subclade2_2} }){
                         unless ($clade3 eq '') {print color ("red"),"     $clade3\n", color ("reset");}
                          for $subclade3 (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3}{$clade1}{$subclade1}{$clade2}{$subclade2_1}{$subclade2_2}{$clade3} }){
                           unless ($subclade3 eq '') {print color ("blue"),"     $subclade3\n", color ("reset");}
                            for $order (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3}{$clade1}{$subclade1}{$clade2}{$subclade2_1}{$subclade2_2}{$clade3}{$subclade3} }){
                             unless ($order eq '') {print color ("yellow"),"    $order\n", color ("reset");}
                              for $family (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3}{$clade1}{$subclade1}{$clade2}{$subclade2_1}{$subclade2_2}{$clade3}{$subclade3}{$order} }){
                               unless ($family eq '') {print color ("green"),"   $family\n", color ("reset");}
                                for $species (keys %{ $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3}{$clade1}{$subclade1}{$clade2}{$subclade2_1}{$subclade2_2}{$clade3}{$subclade3}{$order}{$family} }){
                                 print color ("white"), " $species => $spp{$domain}{$subdomain}{$supergroup}{$kingdom}{$subkingdom1}{$subkingdom2}{$subkingdom3}{$clade1}{$subclade1}{$clade2}{$subclade2_1}{$subclade2_2}{$clade3}{$subclade3}{$order}{$family}{$species}{''}\n", color ("reset");
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
return;
}
