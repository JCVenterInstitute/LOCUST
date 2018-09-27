use warnings;
use strict;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

#Collect Inputs
=head1 NAME

strain_approximation.pl -  A script to approximate ST's using Multiple Sequence
                           Alignments and the Kimura distance between these
                           alignments.

=head1 SYNOPSIS

  USAGE: strain_approximation.pl
   --input_file <ST Output File from Locust>
                          --help

=head1 OPTIONS

B<--input_file, i>   : ST_all.out file from typer.pl run.

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
	'input_file|i=s',
) || die "Error getting options! $!";

pod2usage( { -exitval => 1, -verbose => 2 } ) if $opts{help};

my $stfile = "ST_all.out"







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

#Look at all genomes for matches? Or just those that are Type Strain representatives?
#Look at only genomes that have length within N bases of longest sequence? Avoid those that are SHORT or MISSING?

#Now iterate over each row to get closest
my @row_hashes;
for (my $matrix_row = 1; $matrix_row <= $num_of_genomes; $matrix_row ++){
  my %row_values;
  my $genome_being_examined = $matrix_row - 1; #Key for %genome_idx
  my $genome = $genome_idx{$genome_being_examined};
  #Iterate over values in row
  for (my $col = 0; $col < $num_of_genomes; $col++){
    if ($col != $genome_being_examined){
      my $row_value = abs($n_by_n_distance_array[$matrix_row][$col]);
      $row_values{$col} = $row_value;
    }
  }
  push (@row_hashes, \%row_values);
}

#Hash of Genome to Best Hit(s)
my %best_hits;
#Generate Sorted Arrays of Hits and Values (smallest -> largest)
for (my $i = 0; $i < $#row_hashes; $i++){
  my %row_hash = %{$row_hashes[$i]};
  my $genome = $genome_idx{$i};
  my @hit_names;
  my @hit_values;
  my $best_hit = 10000; #Arbitrarily High
  my @best_hits;
  foreach my $hit_name (sort {$row_hash{$a} <=> $row_hash{$b} } keys %row_hash){
    push @hit_names, $genome_idx{$hit_name};
    push @hit_values, $row_hash{$hit_name};
    if ($row_hash{$hit_name} < $best_hit){
      $best_hit = $row_hash{$hit_name};
      push @best_hits, $genome_idx{$hit_name};
    } elsif ($row_hash{$hit_name} == $best_hit){
      push @best_hits, $genome_idx{$hit_name};
    }
  }
  $best_hits{$genome} = \@best_hits;
  my $out_file = $genome . "/" . $genome . "_sorted_distances.txt";
  open(my $of, ">", $out_file) or die "Couldn't open file $out_file.\n";
    print $of join("\t", @hit_names) . "\n";
    print  $of join("\t", @hit_values) . "\n";
}
