use warnings;
use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use List::Util;
use FindBin qw($Bin);
use lib "$Bin";
#Collect Inputs
=head1 NAME

strain_approximation.pl -  A script to approximate ST's using Multiple Sequence
                           Alignments and the Kimura distance between these
                           alignments.
=head1 SYNOPSIS

  USAGE: strain_approximation.pl
   --typer_input_file <Sequence File List (-i option from typer.pl)>
	 --st_results <ST_all.out file from a finished typer run>
	 --help

=head1 OPTIONS
B<--typer_input_file, i> : Sequence File List

B<--st_results, s>   : ST_all.out file from typer.pl run.

B<--type_alleles, t> : Number of alleles in schema.

B<--attributes, a> : Number of attributes in the input.

B<--help, h>         : Display this help message.

=head1  DESCRIPTION

This program uses assigned ST's to approximate an ST for those assigned "NEW".

=head1  CONTACT

    Chris Greco
    cgreco@jcvi.org

=cut

my %opts;
GetOptions(\%opts,
	'help|h',
	'typer_input_file|i=s',
	'st_results|s=s',
	'type_alleles|t=i',
	'attributes|a=i',
) || die "Error getting options! $!";

pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

my $stfile = $opts{st_results};
my $typer_input_file = $opts{typer_input_file};

my (@header, %st_designations, %st_calls, %new_genomes);
open(my $st, '<', $stfile) or die "Couldn't open $stfile\n";
my $header_seen = 0;
while (<$st>){
	my $line = $_;
	chomp $line;
	if ($header_seen == 0){
		my @pre_header = split("\t",$line);
		@header = @pre_header[0 .. $opts{type_alleles} + 1];
		$header_seen = 1;
} else {
		my @line_values = split("\t", $line);
		chomp @line_values;
        if ("NEW" ~~ @line_values){
            $new_genomes{$line_values[0]} = "APPROXIMATE";
        }
		$st_designations{$line_values[0]} = join("\t", @line_values[1 .. $opts{type_alleles} + 1]);
		$st_calls{$line_values[0]} = $line_values[1];
	}
}

my $type_strains = 0;

my %attributeHash;
open(my $th, '<', $typer_input_file) or die "Couldn't open $typer_input_file\n";
while (<$th>){
	my $line = $_;
	chomp $line;
	my @line_values = split(/\t/,$line,3);
	my $genome = $line_values[0];
	my $path = $line_values[1];
	$attributeHash{$genome} = "";
	if ($line_values[2]){
		$type_strains++;
		$line_values[2] =~ s/^\s+|\s+$//g;
		$attributeHash{$genome} = $line_values[2];
	}
}

my @typeStrains;
while (my ($key, $value) = each(%attributeHash)) {
    if ($value ne ""){
			push (@typeStrains, $key);
		}
}

&generate_dist_mat($type_strains);
my $values_hash = &parse_dist_mat();
my $type_strain_hashes = &sort_hashes(\%$values_hash, \@typeStrains);
my %type_strain_hashes = %$type_strain_hashes;


my @st_out_header = @header;
my @add_to_st_out_header = ("Best Hit", "Second Best Hit");
push @st_out_header, @add_to_st_out_header;

my @st_approx_header = @header;
my @add_to_st_approx_header = ("Best Hit ID", "Best Hit Attribute", "Best Hit ST", "Best Hit Distance", "Second Best Hit ID", "Second Best Hit Attribute", "Second Best Hit ST", "Second Best Hit Distance");
push @st_approx_header, @add_to_st_approx_header;


open(my $st_all, '>', "ST_all_test.out") or die "Couldn't open ST_all.out\n";
open(my $st_approx, '>', "ST_approx.details") or die "Couldn't open ST_approx.details\n";

print $st_all join("\t", @st_out_header) . "\n";
print $st_approx join("\t", @st_approx_header) . "\n";

foreach my $key (keys %st_designations){

    my @st_approx_list = $key;
	push (@st_approx_list,  split("\t", $st_designations{$key}));
	my @st_out_list = @st_approx_list;
	#Add best two Type Strains to ST_OUT
	#Looking at a defined type strain-- only report itself
	if ($attributeHash{$key}){
			#Add itself to ST_all.out
			push (@st_out_list, $attributeHash{$type_strain_hashes->{$key}->{1}->{"Sample"}});
			#Add info to ST_approx
			push (@st_approx_list, $type_strain_hashes->{$key}->{1}->{"Sample"});
			push (@st_approx_list, $attributeHash{$type_strain_hashes->{$key}->{1}->{"Sample"}});
			push (@st_approx_list, $st_calls{$type_strain_hashes->{$key}->{1}->{"Sample"}});
			push (@st_approx_list, $type_strain_hashes->{$key}->{1}->{"Identity"});
	} else {
        if (exists $new_genomes{$key}){
    		push (@st_out_list, $attributeHash{$type_strain_hashes->{$key}->{1}->{"Sample"}});
    		push (@st_out_list, $attributeHash{$type_strain_hashes->{$key}->{2}->{"Sample"}});

    		#Add info to ST_APPROX
    		push (@st_approx_list, $type_strain_hashes->{$key}->{1}->{"Sample"});
    		push (@st_approx_list, $attributeHash{$type_strain_hashes->{$key}->{1}->{"Sample"}});
    		push (@st_approx_list, $st_calls{$type_strain_hashes->{$key}->{1}->{"Sample"}});
    		push (@st_approx_list, $type_strain_hashes->{$key}->{1}->{"Identity"});

    		push (@st_approx_list, $type_strain_hashes->{$key}->{2}->{"Sample"});
    		push (@st_approx_list, $attributeHash{$type_strain_hashes->{$key}->{2}->{"Sample"}});
    		push (@st_approx_list, $st_calls{$type_strain_hashes->{$key}->{2}->{"Sample"}});
    		push (@st_approx_list, $type_strain_hashes->{$key}->{2}->{"Identity"});
    	 }
     }
		print $st_all join( "\t", map { defined $_ ? $_ : '' } @st_out_list ) . "\n";
		print $st_approx join( "\t", map { defined $_ ? $_ : '' } @st_approx_list ) . "\n";
}


sub generate_dist_mat{
	my $type_strains = shift;
	if ($type_strains < 2){
		print "ERROR: There are too few defined type strains from the input. Can't approximate type strain.\n";
		exit();
	} else {
	#Run muscle to create aligned
	my $aligned_file = "allGenomesJoinedAlleles.fasta";
	if (-e $aligned_file){
	my $distance_matrix_command = "Rscript";
	$distance_matrix_command .= " $Bin/distance_matrix_generation.R $aligned_file";
	system($distance_matrix_command) == 0 || die "\n";
} else {
	exit("Couldn't find the MUSCLE generated alignment file.\n");
}
	}
}

sub parse_dist_mat{
	open(my $fh, '<', "allGenomesJoinedAlleles.dist") or die "Couldn't open allGenomesJoinedAlleles.dist\n";
	my $header = 0;
	my %header_hash;

	my %values_hash;
	while(<$fh>){
		my $line = $_;
		my @split_line = split(",", $line);
		if ($header == 0){
			# $i = 0 is the index column
			for (my $i=0; $i < scalar @split_line; $i++){
				my $genome = $split_line[$i];
				chomp $genome;
				$header_hash{$i} = $genome;
				$header = 1;
				}
		} else {
			my $sample = $split_line[0];
			for (my $i=1; $i < scalar @split_line; $i++){
				my $val = $split_line[$i];
				chomp $val;
				my $current_genome = $header_hash{$i};
				$values_hash{$sample}{$current_genome} = $val;
				}
		}
	}
	return (\%values_hash);
}

sub sort_hashes{
	my $hash_of_hashes = shift;
	my $typeStrains = shift;
	my %out_hash;
	foreach my $sample (keys %$hash_of_hashes){
		my $outfile = $sample . "/" . $sample . "_sorted_distances.txt";
		open (my $of, ">", $outfile) or die "Couldn't open file $outfile.\n";
		my $current_hash = $hash_of_hashes->{$sample};
		my @hit_names = sort {$current_hash->{$a} <=> $current_hash->{$b} } keys (%$current_hash);
		my @hit_values = @{$current_hash}{@hit_names};
		my (@out_hit_names, @out_hit_values);
		for (my $i=0; $i < scalar @hit_names; $i++){
			push (@out_hit_names, $hit_names[$i]);
		}
		for (my $i=0; $i < scalar @hit_values; $i++){
			push (@out_hit_values, $hit_values[$i])
		}
		print $of join("\t", @out_hit_names) . "\n";
		print $of join("\t", @out_hit_values) . "\n";
		close $of;
		&run_ckmeans($outfile);
		my $clusterFile = $sample . "/" . $sample . "_cluster_identity.txt";
		my $type_strain_hash = &generate_strain_approx($clusterFile, $typeStrains);
		$out_hash{$sample} = $type_strain_hash;
		}
		return (\%out_hash);
	}

sub run_ckmeans{
	my $input_file = shift;
	my $cmd = "Rscript $Bin/1d_clustering.R $input_file";
	system($cmd) == 0 || die "\n";
}

sub generate_strain_approx{
		my ($clusterFile, $typeStrains) = @_;
		my %out_hash;
		open (my $fh, "<", $clusterFile) or die "Couldn't open file $clusterFile.\n";
		my $genome = (split "/", $clusterFile)[0];
		my (@samples, @identities, @cluster_id);
		my $counter = 0;
		while (<$fh>){
			my $line = $_;
			chomp $line;
			my @split_line = split("\t", $line);
			$counter++;
			if ($counter == 1){
				@samples = @split_line;
			} elsif ($counter == 2){
				@identities = @split_line;
			} elsif ($counter == 3){
				@cluster_id = @split_line;
			}
		}
		my @type_strain_indices;
		my @typeStrainsPresent = @$typeStrains;
		for (my $i = 0; $i < scalar @samples; $i++){
				if ($samples[$i] ~~ @typeStrainsPresent){
					push(@type_strain_indices, $i);
				}
		}
		my $num_of_type_strains = 0;
		for (my $index_of_interest = 0; $index_of_interest < scalar @type_strain_indices; $index_of_interest++){
			my $index = $type_strain_indices[$index_of_interest];
			if ($cluster_id[$index] == 1){
				$num_of_type_strains++;
				$out_hash{$num_of_type_strains}{"Sample"} = $samples[$index];
				$out_hash{$num_of_type_strains}{"Identity"} = $identities[$index];
			}
		}

		if (not exists($out_hash{1})){
			$out_hash{1}{"Sample"} = "";
			$out_hash{1}{"Identity"} = "";
		}

		if (not exists($out_hash{2})){
			$out_hash{2}{"Sample"} = "";
			$out_hash{2}{"Identity"} = "";
		}

		return (\%out_hash);
}
