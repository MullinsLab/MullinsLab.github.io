#!/usr/bin/perl

# August 16, 2010
# run initial mosaik alignment with given input parameters, split the .ace files into fwd and rev files,  run refinement on the fwd and reverse files
use strict;
use Getopt::Long;
use File::Basename;



my %option = (
	      'seq' => '',
	      'qual' => '',
	      'ref' => '',
	      'out_seq' => '',
	      'rur' => '',
	      'hs' => 15,
	      'mmp' => 0.30,
	      'minp' => 0.7,
	      'st' => '454',
	      'script_loc' => '.'
	      );

my $usage = "usage: align_splitace_refine.pl [-option value]

options:  
-seq		input reads fasta file
-qual		input quality file
-ref		input reference file to align
-out_seq        output ace file from MOSAIK
-rur		output file for the information of unaligned reads
-hs		hash size for MosaikAligner (default: 15)
-mmp		maximum mismatch percentage (default: 30% )
-minp		minimum percentage of the aligned read length (default: 0.7)
-st             sequencing technology [454, helicos, illumina, sanger, solid]
 
";

GetOptions (\%option, 'seq=s', 'qual=s', 'ref=s', 'out_seq=s', 'rur=s', 'hs=i', 'mmp=f', 'minp=f', 'st=s');

my $seq = $option{'seq'} or die $usage;
my $qual = $option{'qual'};
my $ref = $option{'ref'} or die $usage;
my $rur = $option{'rur'};
my $out_seq = $option{'out_seq'};
my $hs = $option{'hs'};
my $mmp = $option{'mmp'};
my $minp = $option{'minp'};
my $st = $option{'st'};
my $header = $option{'script_loc'};

my $ref_name = "";

#### Print parameters on screen ###########

print "\n============== parameters ==============\n";
print "Sequence file: $seq\n";
if ($qual) {
  print "Quality file: $qual\n";
}
print "Reference file: $ref\n";
print "rur: $rur\n";
print "hs: $hs\n";
print "mmp: $mmp\n";
print "minp: $minp\n";
print "st: $st\n";
$header =~ s/\s//g;
$header =~ s/\t//g;
$header =~ s/\r//g;
$header =~ s/\n//g;

if(length $header > 0){
  print "Script location: $header\n";
}
else{
  print "Script location: current working directory\n";
  $header = ".";
}
print "========================================\n\n";

##### Build the ref sequence #########
system("/Applications/cli_apps/MosaikBuild -fr $ref -oa ref.dat | tee MOSAIK_run.log");

###### Do the same for reads #########
if ($qual){
  system ("/Applications/cli_apps/MosaikBuild -fr $seq -fq $qual -st 454 -out seq.dat | tee MOSAIK_run.log");
}

#### Mosaik Align parameters ###########
system ("/Applications/cli_apps/MosaikAligner -in seq.dat -out seq_aligned.dat -ia ref.dat -rur $rur -hs $hs -mmp $mmp -minp $minp -p 8 | tee -a MOSAIK_run.log");

##### Mosaik sort parameters ##########
system ("/Applications/cli_apps/MosaikSort -in seq_aligned.dat -out seq_resolved.dat | tee -a MOSAIK_run.log");

###### Mosaik assembler #################
system ("/Applications/cli_apps/MosaikAssembler -in seq_resolved.dat -ia ref.dat -out $out_seq -f ace | tee -a MOSAIK_run.log");

exit;




