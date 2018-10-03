use warnings;
use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use List::Util;
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
) || die "Error getting options! $!";

pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

my $stfile = $opts{st_results};
my $typer_input_file = $opts{typer_input_file};

my (@header, %st_designations, %st_calls);
open(my $st, '<', $stfile) or die "Couldn't open $stfile\n";
my $header_seen = 0;
while (<$st>){
	my $line = $_;
	if ($header_seen == 0){
		@header = split("\t",$line);
		$header_seen = 1;
} else {
		my @line_values = split("\t", $line);
		chomp @line_values;
		$st_designations{$line_values[0]} = join("\t", @line_values[1 .. $#line_values]);
		$st_calls{$line_values[0]} = $line_values[1];
	}
}

my %attributeHash;
open(my $th, '<', $typer_input_file) or die "Couldn't open $typer_input_file\n";
while (<$th>){
	my $line = $_;
	chomp $line;
	my @line_values = split(/\t/,$line,3);
	my $genome = $line_values[0];
	my $path = $line_values[1];
	$attributeHash{$genome} = "";
	my $genomeAttributes = $line_values[2];
	$genomeAttributes =~ s/^\s+|\s+$//g;
	$attributeHash{$genome} = $genomeAttributes;
}

#Run muscle to create aligned
my $aligned_file = "allGenomesJoinedAlleles.fasta";
my $distmat_out = "allGenomesJoinedAlleles.distmat";
my $distmat_command = "distmat -nucmethod 2 -outfile $distmat_out $aligned_file";
system($distmat_command) == 0 || die "\n";

my @n_by_n_distance_array;
open(my $fh, '<', $distmat_out) or die "Couldn't open $distmat_out\n";

my $row = 0;
my $num_of_genomes;
my %genome_idx;
while (<$fh>){
  my $line = $_;
  if ($line =~ /^\t/){
    #distmat output = 1st line of interest begins with tab
    my @split_line = split("\t", $line);
    #Header Row -- Used to initialize Matrix
    if ($row == 0){
      $num_of_genomes = $#split_line;
      push @n_by_n_distance_array, [(0) x $num_of_genomes] for 0 .. $num_of_genomes; #initialize n x n array of size # of genomes
      for (my $col =0; $col < $num_of_genomes; $col++){
        my $val = $split_line[$col +1];  #First line is empy tab
        $val =~ s/^\s+|\s+$//g;
        $n_by_n_distance_array[$row][$col] = $val;
      }
      $row = $row + 1;
    } else {
      my ($genome, $idx) = split(" ", $split_line[-1]);
      $genome_idx{$idx - 1} = $genome;
      for (my $col =0; $col < $num_of_genomes; $col++){
        my $val = $split_line[$col +1];  #First line is empy tab
        $val =~ s/^\s+|\s+$//g;
        $n_by_n_distance_array[$row][$col] = $val;
      }
        $row++;
    }
  }
}

for (my $row_idx = 1; $row_idx <= $num_of_genomes; $row_idx++){ #row_idx 0 is the header
  for (my $col_idx = 0; $col_idx < $num_of_genomes; $col_idx++){
    if ($n_by_n_distance_array[$row_idx][$col_idx] eq ""){
      #print "$row_idx\t$col_idx\n";
      #print "ROW:" . ($row_idx - 1) . "\n";
      #print "COL:" . ($col_idx + 1) . "\n";
      $n_by_n_distance_array[$row_idx][$col_idx] = $n_by_n_distance_array[$col_idx + 1][$row_idx - 1];
    }
  }
}

#Now iterate over each row to get closest
my @row_hashes;
for (my $matrix_row = 1; $matrix_row <= $num_of_genomes; $matrix_row ++){
  my %row_values;
  my $genome_being_examined = $matrix_row - 1; #Key for %genome_idx
  my $genome = $genome_idx{$genome_being_examined};
  #Iterate over values in row
  for (my $col = 0; $col < $num_of_genomes; $col++){
      my $row_value = abs($n_by_n_distance_array[$matrix_row][$col]);
      $row_values{$col} = $row_value;
    
  }
  push (@row_hashes, \%row_values);
}

#Hash of Genome to Best Hit(s)
my %out_st_line;
my %st_all_out;
#Generate Sorted Arrays of Hits and Values (smallest -> largest)
for (my $i = 0; $i < $#row_hashes; $i++){
  my %row_hash = %{$row_hashes[$i]};
  my $genome = $genome_idx{$i};
	my @hit_names;
  my @hit_values;
  my $best_hit = 10000; #Arbitrarily High
  foreach my $hit_name (sort {$row_hash{$a} <=> $row_hash{$b} } keys %row_hash){
    push @hit_names, $genome_idx{$hit_name};
    push @hit_values, $row_hash{$hit_name};
  }
  my $out_file = $genome . "/" . $genome . "_sorted_distances.txt";
  open(my $of, ">", $out_file) or die "Couldn't open file $out_file.\n";
    print $of join("\t", @hit_names) . "\n";
    print  $of join("\t", @hit_values) . "\n";
	my @allele_calls = split("\t", $st_designations{$genome});
	my @st_array = @allele_calls;
	my $st_call = $allele_calls[0];
	my @possible_indices;
	for (my $i = 0; $i < $#hit_values; $i++){
		 	my $hit_name = $hit_names[$i];
			if ($attributeHash{$hit_name}){ #If current hit name has a defined 3rd column of metadata
				push @possible_indices, $i;
		}
 }
 my %possible_hit_info;
 for my $val (@possible_indices){
	 			my $genome = $hit_names[$val];
				my $dist = $hit_values[$val];
				push (@{$possible_hit_info{$dist}}, $genome);
 }
	my %best_hit;
	my %second_best_hit;
	my $count = 0;
	foreach my $val (sort {$a <=> $b} keys %possible_hit_info){
		$count = $count + 1;
		if ($count == 1){
			$best_hit{$val} = $possible_hit_info{$val};
	} elsif ($count == 2){
		$second_best_hit{$val} = $possible_hit_info{$val};
	}
}
	my $best_hit_val = (keys %best_hit)[0];
	my @best_hits = (values %best_hit)[0];
	my $second_best_hit_val = (keys %second_best_hit)[0];
	my @second_best_hits = (values %second_best_hit)[0];
	my @best_hit_names;
	my @best_hits_sts;
	my @best_hit_attribute;
	my @second_best_hit_names;
	my @second_best_hits_sts;
	my @second_best_hit_attribute;
	foreach my $hit_list (@best_hits){
		foreach my $hit (@{$hit_list}){
			push(@best_hit_names, $hit);
			push(@best_hits_sts, $st_calls{$hit});
			push(@best_hit_attribute, $attributeHash{$hit});
		}
	}
	foreach my $hit_list (@second_best_hits){
		foreach my $hit (@{$hit_list}){
			push (@second_best_hit_names, $hit);
			push(@second_best_hits_sts, $st_calls{$hit});
			push(@second_best_hit_attribute, $attributeHash{$hit});
		}
	}
	push (@allele_calls, join("/", @best_hit_names));
	push (@allele_calls, $best_hit_val);
	push (@allele_calls, join("/", @best_hits_sts));
	push (@allele_calls, join("/", @best_hit_attribute));
	push (@st_array, join("/", @best_hit_attribute));
	push (@allele_calls, join("/", @second_best_hit_names));
	push (@allele_calls, $second_best_hit_val);
	push (@allele_calls, join("/", @second_best_hits_sts));
	push (@allele_calls, join("/", @second_best_hit_attribute));
	my $joined_allele_calls = join("\t", @allele_calls);
	my $joined_st_array = join("\t", @st_array);
	$out_st_line{$genome} = $joined_allele_calls;
	$st_all_out{$genome} = $joined_st_array;
}

#Write results
my $out_approx_file = "ST_approx.out";
open(my $of, ">", $out_approx_file) or die "Couldn't open file $out_approx_file.\n";
my $header = join("\t", @header);
chomp $header;
print $of $header . "\tBest Hit\tBest Kimura Distance\tBest Sequence Type\tBest Proxy/Type Strain\tSecond Best Hit\tSecond Best Kimura Distance\tSecond Best Sequence Type\tSecond Best Proxy/Type Strain\n";
while (my ($key, $value) = each %out_st_line){
	print $of "$key\t$value\n";
}

my $out_st_file = "ST_all.out";
open(my $of, ">", $out_st_file) or die "Couldn't open file $out_st_file.\n";
my $header = join("\t", @header);
chomp $header;
print $of $header . "\tBest Hit\n";
while (my ($key, $value) = each %st_all_out){
	print $of "$key\t$value\n";
}
