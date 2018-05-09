use File::Slurp;
use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;

my %opts;
GetOptions(\%opts,
	'help|h',
	'input_file|i=s',
) || die "Error getting options! $!";

if ( $opts{help} ) {
	die "\nInput an ST_Out.txt (-i/--input_file) file from a Locust run to generate an ITOL Colored Strip Annotation File.\nThere are 220 available colors. If you have more ST's than this, colors will be repeated.\n\n"
}

my @color_palette = ("rgba(102,194,165,1)", "rgba(252,141,98,1)", "rgba(141,160,203,1)", "rgba(231,138,195,1)", "rgba(166,216,84,1)", "rgba(255,217,47,1)", "rgba(229,196,148,1)", "rgba(179,179,179,1)", "rgba(228,26,28,1)", "rgba(55,126,184,1)", "rgba(77,175,74,1)", "rgba(152,78,163,1)", "rgba(255,127,0,1)", "rgba(255,255,51,1)", "rgba(166,86,40,1)", "rgba(247,129,191,1)", "rgba(153,153,153,1)", "rgba(141,211,199,1)", "rgba(255,255,179,1)", "rgba(190,186,218,1)", "rgba(251,128,114,1)", "rgba(128,177,211,1)", "rgba(253,180,98,1)", "rgba(179,222,105,1)", "rgba(252,205,229,1)", "rgba(217,217,217,1)", "rgba(188,128,189,1)", "rgba(204,235,197,1)", "rgba(255,237,111,1)", "rgba(137,197,218,1)", "rgba(218,87,36,1)", "rgba(116,217,68,1)", "rgba(206,80,202,1)", "rgba(63,73,33,1)", "rgba(192,113,124,1)", "rgba(203,213,136,1)", "rgba(95,127,199,1)", "rgba(103,55,112,1)", "rgba(211,217,62,1)", "rgba(56,51,62,1)", "rgba(80,133,120,1)", "rgba(215,193,177,1)", "rgba(104,144,48,1)", "rgba(173,111,59,1)", "rgba(205,155,205,1)", "rgba(209,66,133,1)", "rgba(109,222,136,1)", "rgba(101,41,38,1)", "rgba(127,220,192,1)", "rgba(200,66,72,1)", "rgba(133,105,213,1)", "rgba(94,115,143,1)", "rgba(209,163,61,1)", "rgba(138,124,100,1)", "rgba(89,152,97,1)", "rgba(102,194,165,0.5)", "rgba(252,141,98,0.5)", "rgba(141,160,203,0.5)", "rgba(231,138,195,0.5)", "rgba(166,216,84,0.5)", "rgba(255,217,47,0.5)", "rgba(229,196,148,0.5)", "rgba(179,179,179,0.5)", "rgba(228,26,28,0.5)", "rgba(55,126,184,0.5)", "rgba(77,175,74,0.5)", "rgba(152,78,163,0.5)", "rgba(255,127,0,0.5)", "rgba(255,255,51,0.5)", "rgba(166,86,40,0.5)", "rgba(247,129,191,0.5)", "rgba(153,153,153,0.5)", "rgba(141,211,199,0.5)", "rgba(255,255,179,0.5)", "rgba(190,186,218,0.5)", "rgba(251,128,114,0.5)", "rgba(128,177,211,0.5)", "rgba(253,180,98,0.5)", "rgba(179,222,105,0.5)", "rgba(252,205,229,0.5)", "rgba(217,217,217,0.5)", "rgba(188,128,189,0.5)", "rgba(204,235,197,0.5)", "rgba(255,237,111,0.5)", "rgba(137,197,218,0.5)", "rgba(218,87,36,0.5)", "rgba(116,217,68,0.5)", "rgba(206,80,202,0.5)", "rgba(63,73,33,0.5)", "rgba(192,113,124,0.5)", "rgba(203,213,136,0.5)", "rgba(95,127,199,0.5)", "rgba(103,55,112,0.5)", "rgba(211,217,62,0.5)", "rgba(56,51,62,0.5)", "rgba(80,133,120,0.5)", "rgba(215,193,177,0.5)", "rgba(104,144,48,0.5)", "rgba(173,111,59,0.5)", "rgba(205,155,205,0.5)", "rgba(209,66,133,0.5)", "rgba(109,222,136,0.5)", "rgba(101,41,38,0.5)", "rgba(127,220,192,0.5)", "rgba(200,66,72,0.5)", "rgba(133,105,213,0.5)", "rgba(94,115,143,0.5)", "rgba(209,163,61,0.5)", "rgba(138,124,100,0.5)", "rgba(89,152,97,0.5)", "rgba(102,194,165,0.75)", "rgba(252,141,98,0.75)", "rgba(141,160,203,0.75)", "rgba(231,138,195,0.75)", "rgba(166,216,84,0.75)", "rgba(255,217,47,0.75)", "rgba(229,196,148,0.75)", "rgba(179,179,179,0.75)", "rgba(228,26,28,0.75)", "rgba(55,126,184,0.75)", "rgba(77,175,74,0.75)", "rgba(152,78,163,0.75)", "rgba(255,127,0,0.75)", "rgba(255,255,51,0.75)", "rgba(166,86,40,0.75)", "rgba(247,129,191,0.75)", "rgba(153,153,153,0.75)", "rgba(141,211,199,0.75)", "rgba(255,255,179,0.75)", "rgba(190,186,218,0.75)", "rgba(251,128,114,0.75)", "rgba(128,177,211,0.75)", "rgba(253,180,98,0.75)", "rgba(179,222,105,0.75)", "rgba(252,205,229,0.75)", "rgba(217,217,217,0.75)", "rgba(188,128,189,0.75)", "rgba(204,235,197,0.75)", "rgba(255,237,111,0.75)", "rgba(137,197,218,0.75)", "rgba(218,87,36,0.75)", "rgba(116,217,68,0.75)", "rgba(206,80,202,0.75)", "rgba(63,73,33,0.75)", "rgba(192,113,124,0.75)", "rgba(203,213,136,0.75)", "rgba(95,127,199,0.75)", "rgba(103,55,112,0.75)", "rgba(211,217,62,0.75)", "rgba(56,51,62,0.75)", "rgba(80,133,120,0.75)", "rgba(215,193,177,0.75)", "rgba(104,144,48,0.75)", "rgba(173,111,59,0.75)", "rgba(205,155,205,0.75)", "rgba(209,66,133,0.75)", "rgba(109,222,136,0.75)", "rgba(101,41,38,0.75)", "rgba(127,220,192,0.75)", "rgba(200,66,72,0.75)", "rgba(133,105,213,0.75)", "rgba(94,115,143,0.75)", "rgba(209,163,61,0.75)", "rgba(138,124,100,0.75)", "rgba(89,152,97,0.75)", "rgba(102,194,165,0.25)", "rgba(252,141,98,0.25)", "rgba(141,160,203,0.25)", "rgba(231,138,195,0.25)", "rgba(166,216,84,0.25)", "rgba(255,217,47,0.25)", "rgba(229,196,148,0.25)", "rgba(179,179,179,0.25)", "rgba(228,26,28,0.25)", "rgba(55,126,184,0.25)", "rgba(77,175,74,0.25)", "rgba(152,78,163,0.25)", "rgba(255,127,0,0.25)", "rgba(255,255,51,0.25)", "rgba(166,86,40,0.25)", "rgba(247,129,191,0.25)", "rgba(153,153,153,0.25)", "rgba(141,211,199,0.25)", "rgba(255,255,179,0.25)", "rgba(190,186,218,0.25)", "rgba(251,128,114,0.25)", "rgba(128,177,211,0.25)", "rgba(253,180,98,0.25)", "rgba(179,222,105,0.25)", "rgba(252,205,229,0.25)", "rgba(217,217,217,0.25)", "rgba(188,128,189,0.25)", "rgba(204,235,197,0.25)", "rgba(255,237,111,0.25)", "rgba(137,197,218,0.25)", "rgba(218,87,36,0.25)", "rgba(116,217,68,0.25)", "rgba(206,80,202,0.25)", "rgba(63,73,33,0.25)", "rgba(192,113,124,0.25)", "rgba(203,213,136,0.25)", "rgba(95,127,199,0.25)", "rgba(103,55,112,0.25)", "rgba(211,217,62,0.25)", "rgba(56,51,62,0.25)", "rgba(80,133,120,0.25)", "rgba(215,193,177,0.25)", "rgba(104,144,48,0.25)", "rgba(173,111,59,0.25)", "rgba(205,155,205,0.25)", "rgba(209,66,133,0.25)", "rgba(109,222,136,0.25)", "rgba(101,41,38,0.25)", "rgba(127,220,192,0.25)", "rgba(200,66,72,0.25)", "rgba(133,105,213,0.25)", "rgba(94,115,143,0.25)", "rgba(209,163,61,0.25)", "rgba(138,124,100,0.25)", "rgba(89,152,97,0.25)");
my %st_dict;
my %color_dict;
my $input_file = $opts{input_file} || die "Couldn't open the input file.\n";

my @full_st_array;
open(my $fh, "<", $input_file)
    or die "Failed to open file: $!\n";
while(<$fh>) {
    chomp;
    push @full_st_array, $_;
}
close $fh;

my @headerArray;
foreach (split '\t', $full_st_array[0]){
  push @headerArray, $_;
}


for (my $i = 1; $i <  scalar @full_st_array; $i++){
	my @row = split "\t", $full_st_array[$i];
	my $genome = $row[0];
	my $ST_CALL = $row[1];
	my $row_length = scalar @row;
	my @allele_designations = @row[2 .. $row_length];
	#Captures any NEW allele calls-- immediately designates ST as new
	if ( "NEW" ~~ @allele_designations) {
		$ST_CALL = "new";
	}
	#Captures Short or Missing allele calls-- immediately designates ST as unknown
	if ("SHORT" ~~ @allele_designations || "MISSING" ~~ @allele_designations){
		$ST_CALL = "unknown";
	}
	#Captures those with all alleles called but unknown to ST schema file. All UNKNOWN ST's should have been removed already except for these.
	if ((not "NEW" ~~ @allele_designations) && $ST_CALL eq "UNKNOWN") {
		$ST_CALL = "new";
	}

	$st_dict{$genome} = $ST_CALL;
}

my $count = 0;
for my $st (values %st_dict){
	if ($color_dict{$st}){
		next;
	}
	else{
		if ($count < scalar @color_palette){
			$color_dict{$st} = $color_palette[$count];
		}
		else {
			$color_dict{$st} = $color_palette[$count % scalar @color_palette];
		}
	}
	$count ++;
}

my @legend_shapes;
my @legend_colors;
my @legend_labels;
#print Dumper \%color_dict;
while(my($k, $v) = each %color_dict){
	push(@legend_shapes, 1);
	push(@legend_colors, $v);
	push(@legend_labels, $k);
}

my $legend_shape_str = join("\t", @legend_shapes);
my $legend_colors_str = join("\t", @legend_colors);
my $legend_labels_str = join("\t", @legend_labels);

my $outfile = "ITOL_ST_ANNOTATION.txt";
open (my $oh, '>', $outfile) or die "Couldn't open the outfile: $outfile";
print $oh "DATASET_COLORSTRIP\n";
print $oh "SEPERATOR TAB\n";
print $oh "DATASET_LABEL\tSequence Type\n";
print $oh "LEGEND_TITLE\tSequence Type\n";
print $oh "LEGEND_SHAPES\t$legend_shape_str\n";
print $oh "LEGEND_COLORS\t$legend_colors_str\n";
print $oh "LEGEND_LABELS\t$legend_labels_str\n";
print $oh "COLOR\t#ff0000\n";
print $oh "DATA\n";
while (my ($key, $value) = each (%st_dict)){
	print $oh "$key\t$color_dict{$value}\t$value\n";
}
close $oh;
if (scalar @color_palette <= keys %color_dict){
	print "Warning: There are less available colors than there are Sequence Types. Colors will be repeated.\n"
}
