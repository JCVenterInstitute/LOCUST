use warnings;
use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use List::Util;
use Bio::TreeIO;
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
) || die "Error getting options! $!";

pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

# Main Function
#Step 1 -- Get a list of all "type" strains or other defined genomes (column 3 of genomes.list)

my $attributeHash = &parse_type_file($opts{typer_input_file}); #%$attributeHash
my %attributeHash = %$attributeHash;
my $type_strains = &find_type_strains(%attributeHash); #@$type_strains
my @type_strains = @$type_strains;

if (scalar @type_strains < 2) {
    print "Can't do strain approximation with less than 2 defined genomes.\n";
    exit();
}

my $st_hash = &get_st_vals;
my %st_hash = %$st_hash;

#Step 2-- Generate Unsupervised Clustering Data and Append Results to distance matrix
#This is all taken care of via the r script (distance matrix generation, clustering, etc.)

#Output generates a file "tsne_mclust_results.txt"
&generate_dist_mat_and_cluster;
my $clusters_with_distance = "mclust_results.txt";

#Step 3 -- Find best type strain hit per genome
#read in cluster results
#write per genome results
my $main_info = &read_clusters(\@$type_strains, $clusters_with_distance);
my %main_info = %$main_info;

#Step 4 -- Add info to the ST_all.out file
# ST_approx.details contains 2 * (Best Hit ID, Best Hit Attribute, Best Hit ST, Best Hit Distance, Best Hit Cluster)
# Must also add Genome's Cluster Membership

#best hit id = @type_strain
#best hit attribute = $attributeHash{$hit id}
#best hit st = ST_all.out file

#ST_all_approx.out contains Best Hit and Second Best Hit

my $st_file = "ST_all.out";
my $details_file = "ST_approx.details";
my $approx_file = "ST_all_approx.out";

open(my $st_fh, "<", $st_file) or die "Couldn't open $st_file\n";
open(my $dt_fh, ">", $details_file) or die "Couldn't open $details_file\n";
open(my $ap_fh, ">", $approx_file) or die "Couldn't open $approx_file\n";
my @st_approx_details_header;
my @st_all_approx_header;
my $attrib_cols;
my $firstline = 1;
while(<$st_fh>){
    my $line = $_;
    chomp $line;
    #Add to Headers
    if ($firstline){
        $firstline = 0;
        my @split_line = split("\t", $line);

        $attrib_cols = scalar @split_line - ($opts{type_alleles} + 2); #Sample and ST columns

        push (@st_all_approx_header, @split_line);
        push (@st_all_approx_header, ("Best Hit", "Second Best Hit"));
        print $ap_fh join("\t", @st_all_approx_header) . "\n";
        push (@st_approx_details_header, @split_line);
        push (@st_approx_details_header, ("Genome Cluster", "Best Hit ID", "Best Hit Attribute", "Best Hit ST", "Best Hit Distance", "Best Hit Cluster", "Second Best Hit ID", "Second Best Hit Attribute", "Second Best Hit ST", "Second Best Hit Distance", "Second Best Hit Cluster"));
        print $dt_fh join("\t", @st_approx_details_header) . "\n";
    } else {
        my @st_approx_details_body;
        my @st_all_approx_body;

        my @split_line = split("\t", $line);
        @split_line = @split_line[0..$opts{type_alleles} + 1];
        my $genome = $split_line[0];
        push(@st_approx_details_body, @split_line);
        push(@st_all_approx_body, @split_line);
        my $best_hit_file = $genome . "/" . $genome . "_sorted_distances.txt";
        open(my $bh, "<", $best_hit_file) or die "Couldn't open $best_hit_file\n";
        my @filelines;
        while(<$bh>){
            my $line = $_;
            chomp $line;
            push(@filelines, $line);
        }
        my @type_strains = split("\t", $filelines[0]);
        my @type_strain_distances = split("\t", $filelines[1]);
        my @type_strain_clusters = split("\t", $filelines[2]);
        #Current Genome Cluster
        my $genome_cluster = $main_info{$genome}->{"Cluster"};
        #Best Hit
        my $best_hit_id = $type_strains[0];
        my $best_hit_attribute = $attributeHash{$best_hit_id};
        my $best_hit_st = $st_hash{$best_hit_id};
        my $best_hit_distance = $type_strain_distances[0];
        my $best_hit_cluster = $type_strain_clusters[0];
        #Second Best Hit
        my $second_best_hit_id = $type_strains[1];
        my $second_best_hit_attribute = $attributeHash{$second_best_hit_id};
        my $second_best_hit_st = $st_hash{$second_best_hit_id};
        my $second_best_hit_distance = $type_strain_distances[1];
        my $second_best_hit_cluster = $type_strain_clusters[1];
        close($bh);

        for (my $i=0;$i < $attrib_cols;$i++){
            push(@st_all_approx_body, "");
            push(@st_approx_details_body, "");
        }

        #ST_all_approx.out
        push(@st_all_approx_body, ($best_hit_attribute, $second_best_hit_attribute));
        print $ap_fh join("\t", @st_all_approx_body) . "\n";
        #ST_approx.details
        push(@st_approx_details_body, ($genome_cluster, $best_hit_id, $best_hit_attribute, $best_hit_st, $best_hit_distance, $best_hit_cluster));
        push(@st_approx_details_body, ($second_best_hit_id, $second_best_hit_attribute, $second_best_hit_st, $second_best_hit_distance, $second_best_hit_cluster));
        print $dt_fh join("\t", @st_approx_details_body) . "\n";


    }
}

close($dt_fh);
close($ap_fh);

sub get_st_vals{
    my $in_file = "ST_all.out";
    my %st_hash;
    open(my $if, "<", $in_file) or die "Couldn't open $in_file\n";
    while (<$if>){
        my $line = $_;
        chomp $line;
        my @split_line = split("\t", $line);
        $st_hash{$split_line[0]} = $split_line[1];
    }
    return (\%st_hash);
}


sub read_clusters{
    my $type_strains = shift;
    my $cluster_file = shift;
    my $firstline = 1;
    my %cluster_cols;
    my %header_hash;
    my %out_hash;
    my %cluster_hash;
    my @sorted_clusters;
    my $cluster_identity;
    open (my $ch, "<", $cluster_file) or die "Couldn't open $cluster_file\n";
    while(<$ch>){
        my $line = $_;
        chomp $line;
        if ($firstline){
            $firstline = 0;
            my @split_line = split("\t", $line);
            my %typestrains = map { $_ => 1} @$type_strains;
            for (my $i=0; $i < scalar @split_line; $i++){
                my $genome = $split_line[$i];
                chomp $genome;
                if ($genome eq "Cluster"){
                    $cluster_identity = $i + 1;
                } elsif ($genome =~ /^Cluster/){
                    $cluster_cols{$genome} = $i + 1;
                    push (@sorted_clusters, $genome);
                }
                if(exists($typestrains{$genome})){
                    $header_hash{$i + 1} = $genome;
                }
            }
        } else {
            my @split_line = split("\t", $line);
            my $genome_row = $split_line[0];
            while ( my ($key, $value) = each (%header_hash)){
                $out_hash{$genome_row}{$value} = $split_line[$key];
            }
            $out_hash{$genome_row}{"Cluster"} = $split_line[$cluster_identity];
            my $ci_file = $genome_row . "/" . $genome_row . "_cluster_identity.txt"; #Cluster Identity
            open (my $ci, ">", $ci_file) or die "Couldn't open $ci_file";
            #Header of Cluster Identity Columns
            print $ci "Cluster" . "\t" . join("\t", @sorted_clusters) . "\n";
            #Get Cluster and Cluster Uncertainty in Same Order
            my @genome_specific_cluster_row;
            push (@genome_specific_cluster_row, $split_line[$cluster_identity]);
            foreach (@sorted_clusters){
                my $index = $cluster_cols{$_};
                push(@genome_specific_cluster_row, $split_line[$index - 1]);
            }
            print $ci join("\t", @genome_specific_cluster_row) . "\n";
            close $ci;


            }
        }
        for my $genome (sort keys %out_hash){
            #Write sorted distances (lowest to greatest) for genome to each type strain
            my $sorted_distances_file = $genome . "/" . $genome . "_sorted_distances.txt";
            open (my $sd, ">", $sorted_distances_file) or die "Couldn't open $sorted_distances_file";

            my $h2 = $out_hash{$genome};
            my @type_strains_ordered = sort { $h2->{$a} <=> $h2->{$b}} keys %$h2;
            @type_strains_ordered = grep { $_ ne "Cluster"} @type_strains_ordered;
            my @type_strain_distances = map { $h2->{$_} } @type_strains_ordered;
            my @cluster_membership = map { $out_hash{$_}->{"Cluster"} } @type_strains_ordered;
            print $sd join("\t", @type_strains_ordered) . "\n";
            print $sd join("\t", @type_strain_distances) . "\n";
            print $sd join("\t", @cluster_membership) . "\n";
        }
        return(\%out_hash);
    }

sub parse_type_file{
    my $typer_input_file = shift;
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
    		$line_values[2] =~ s/^\s+|\s+$//g;
    		$attributeHash{$genome} = $line_values[2];
    	}
    }
    return(\%attributeHash);
}

sub find_type_strains{
    my $attribute_hash = shift;
    my @typeStrains;
    while (my ($key, $value) = each(%$attributeHash)) {
        if ($value ne ""){
    			push (@typeStrains, $key);
    		}
    }
    return(\@typeStrains);
}

sub generate_dist_mat_and_cluster{
	#Run muscle to create aligned
	my $aligned_file = "allGenomesJoinedAlleles.fasta";
	if (-e $aligned_file){
	       my $distance_matrix_command = "Rscript";
	       $distance_matrix_command .= " $Bin/mclust_cluster.R $aligned_file";
	       system($distance_matrix_command) == 0 || die "\n";
    } else {
	exit("Couldn't find the MUSCLE generated alignment file.\n");
        }
	}
