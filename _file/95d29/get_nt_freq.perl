#!/usr/bin/perl

#########################
# Script to pick out SNP from a given alignment
# The alignment should be a nt alignment

## 
#########################

my $usage = "perl get_nt_freq.pl in_file out_file sample_id";

my $in_file = shift or die $usage;
my $sample_id = shift or die $usage;


my @alignment = ();         # store input file
our %ref_pos = ();
our %all_read_start = ();
our %all_read_end = ();
our %all_read_seq = ();

our $ref_seq = "";
our $ref_len = 0;
our %ref_pos = (); # Index for storing positions

our %all_pos_char = ();
our %consensus_char_for_pos = ();
our %var_pos = ();
our %all_pos_freq = ();

our %ref_wo_gap = ();

our %var_per_pos = ();
our %indel_pos = ();  # store positions that have INDELs

my $freq_file = "$sample_id" . "_nt_variant_freq_table.out";
my $var_file = "$sample_id" . "_nt_coverage_variation_table.out";


open(INPUT, $in_file);
@alignment = <INPUT>;
close INPUT;

for( my $i = 0; $i < scalar @alignment; $i++){
    chomp $alignment[$i];
    $alignment[$i] =~ s/\s//g; $alignment[$i] =~ s/\r//g; $alignment[$i] =~ s/\n//g;
    
    if($alignment[$i] =~ />(\S+)/){
	my $header = $1;
	$header =~ s/\s//g; $header =~ s/\t//g;
	$header =~ s/\r//g; $header =~ s/\n//g;
	
	if( ($header =~ /reference/i)||  ($header =~ /hxb2/i ) || ($header =~ /consensus/i) ){
	    # deal with the reference sequence
	    for(my $k=$i+1; $k < scalar @alignment; $k++){
		chomp $alignment[$k];
		$alignment[$k] =~ s/\s//g; $alignment[$k] =~ s/\r//g; $alignment[$k] =~ s/\n//g;
		if($alignment[$k] =~ />(\S+)/){
		    #next header
		    $k = scalar @alignment;
		}
		else{
		    $ref_seq = "$ref_seq" . "$alignment[$k]";
		}
	    }
	    my @ref = split(//, $ref_seq);
	    $ref_len = length $ref_seq;
	    %ref_pos = get_nt_pos(\@ref); # 1 based index for positions
	}
	## All other sequences besides the reference sequence
	else{
	    my $curr_read = "";
	    my $read_name = $header;
	    
	    
	    for(my $k=$i+1; $k < scalar @alignment; $k++){
		chomp $alignment[$k];
		$alignment[$k] =~ s/\s//g; $alignment[$k] =~ s/\r//g; $alignment[$k] =~ s/\n//g;
		if($alignment[$k] =~ />(\S+)/){
		    #next header
		    $k = scalar @alignment;
		}
		else{
		    $curr_read = "$curr_read" . "$alignment[$k]";
		}
	    }
	    my @read = split(//, $curr_read);
	    my ($start, $end) = get_start_end(\@read);
	    
	    $all_read_start{$read_name} = $start;
	    $all_read_end{$read_name} = $end;
	    if(exists $all_read_seq{$read_name}){
		print "THIS READ $read_name SEEMS TO BE PRESENT IN THIS ALIGNMENT MORE THAN ONCE.\n";
		print "Please fix this error and then run the program again. Program will exit.\n";
		exit;
	    }
	    else{
		$all_read_seq{$read_name} = $curr_read;
	    }
	}
    }
}

%ref_wo_gap = get_ref(\%ref_pos);   ## STORE THE REFERENCE BASED POSITION

my $num_reads = scalar (keys %all_read_seq);
print "Number of unique sequences in this file: $num_reads\n";
print "calculating frequencies ....\n";



#### GET THE TYPE OF CHARACTERS AT EACH POSTION #######
foreach my $p(sort {$a <=> $b} keys %ref_pos){
  my %char_ct = ();
  
  foreach my $rd(keys %all_read_seq){
    
    if( ($p >= $all_read_start{$rd})  && ($p <= $all_read_end{$rd})){
    
      my $rdch = substr($all_read_seq{$rd}, ($p-1), 1);
   
      if(exists $all_pos_char{$p}{$rdch}){
	my $ct = $all_pos_char{$p}{$rdch};
	$all_pos_char{$p}{$rdch} = $ct + 1;
	$char_ct{$rdch} = $all_pos_char{$p}{$rdch};
      }
      else{
	$all_pos_char{$p}{$rdch} = 1;
	$char_ct{$rdch} = $all_pos_char{$p}{$rdch};
      }
    }
    
  } # done going through all reads for this position
  ## Get the consensus character for this position
  my $conchar = get_key_with_highest_val(\%char_ct);
  $consensus_char_for_pos{$p} = $conchar; ## even if consensus char is a gap at that position
}





############# CALCULATE FREQUENCIES ##########
foreach my $p(sort {$a <=> $b} keys %ref_pos){
  
  my $coverage = get_read_coverage($p);
  if($coverage > 0){
    
    foreach my $ch(keys %{$all_pos_char{$p}}){

      $all_pos_freq{$p}{$ch} = $all_pos_char{$p}{$ch}/$coverage;
      
      if( ($all_pos_freq{$p}{$ch} > 0)  && ($all_pos_freq{$p}{$ch} < 1) ){
	$var_pos{$p}{$ch} = $all_pos_freq{$p}{$ch};
      }
    }
  }
}

## find the number of INDELs and SNPs at each position
foreach my $p(sort {$a <=> $b} keys %ref_pos){
  my $num_var = 0;

  foreach my $ch(keys %{$all_pos_char{$p}}){
    
    if( ($all_pos_freq{$p}{$ch} > 0) && ($all_pos_freq{$p}{$ch} < 1) ){
      $num_var++;
    }
    else{
      $num_var = 0;
    }
   
    if($ch eq '-'){
      $indel_pos{$p} = 1;
    }
  }
  $var_per_pos{$p} = $num_var;
}

# print frequencies to file
open(FREQFILE, ">$freq_file");

print FREQFILE "------------ ALL VARIANT  FREQUENCIES --------------\n";
print FREQFILE "\tA\tG\tT\tC\t-\n";

## print frequencies of positions with variations
foreach my $varpos(sort {$a <=> $b} keys %var_pos){
  my $ref_based_pos = $ref_wo_gap{$varpos};

  print FREQFILE "$ref_based_pos\t$var_pos{$varpos}{'A'}\t$var_pos{$varpos}{'G'}\t$var_pos{$varpos}{'T'}\t$var_pos{$varpos}{'C'}\t$var_pos{$varpos}{'-'}\n";
}


print FREQFILE "------- ALL FREQUENCIES --------\n";
print FREQFILE "\tA\tG\tT\tC\t-\n";

## print frequencies of positions with variations                                                                                   
foreach my $varpos(sort {$a <=> $b} keys %all_pos_freq){
    my $ref_based_pos = $ref_wo_gap{$varpos};

  print FREQFILE "$ref_based_pos\t$all_pos_freq{$varpos}{'A'}\t$all_pos_freq{$varpos}{'G'}\t$all_pos_freq{$varpos}{'T'}\t$all_pos_freq{$varpos}{'C'}\t$all_pos_freq{$varpos}{'-'}\n";
}


close FREQFILE;


#### print coverage and number of variations to file ####
open(OUTPUT_FILE, ">$var_file");
print OUTPUT_FILE "position\t type\t readcoverage\tnum of variants\n";
foreach my $p(sort {$a <=> $b} keys %var_per_pos){
  my $type = "";
  
  if( (exists $indel_pos{$p}) && ($var_per_pos{$p} > 0) ){
    $type = "INDEL";
  }
 
  my $ref_based_pos = $ref_wo_gap{$p};
  my $coverage = get_read_coverage($p);
  
  print OUTPUT_FILE "$ref_based_pos\t $type\t$coverage\t$var_per_pos{$p}\n";
}

close OUTPUT_FILE;


exit;

##########SUBROUTINES ############

###################
# Get read consensus based on given alignment
###################
sub get_key_with_highest_val{
  my ($href1) = @_;
  my %vals = %$href1;

  my %curr_ct = (); # store all freq for this position and their

  foreach my $ch(keys %vals){
  
    my $ct = $vals{$ch};
    $curr_ct{$ct} = $ch;
  }

  my @sorted_cts = sort ({$b <=> $a} keys %curr_ct);
  my $high_ct = shift @sorted_cts;
  my $high_char = $curr_ct{$high_ct};
  return $high_char;
}


##############
# for a given ref sequence position, find the number of reads that cover that position,
##############
sub get_read_coverage{
  my($position) = @_;
  my $out = 0;
  
  # cycle through all reads
  foreach my $rd(keys %all_read_seq){
    
    #get start and end positions for this read
    my $s = $all_read_start{$rd};
    my $e = $all_read_end{$rd};
    
    if( ($position >= $s) && ($position <= $e) ){
      ## Since the red names includes frequencies, include freq information
    
      #print "$rd read freq: $rdfreq\n";
      $out = $out + 1;
    }
  }
  return $out;
}


#############
# Given an array, split it into hash, with key being 1 based index
# and value being the array element
#############
sub get_nt_pos{
    my($aref) = @_;
    my %out = ();
    my $count = 1;
    foreach my $rec(@$aref){
	$rec =~ s/\s//g;
	$out{$count} = $rec;
	$count++;
    }
    return %out;
}


#############
# Given an array, get the start and end position
#############
sub get_start_end{
    my($aref) = @_;

    my @arr = @$aref;

    # Both start and end positions are based based on 1-based indeces
    my $st = 0;
    my $end = 0;
    # start position
    for(my $i =0; $i < scalar @arr; $i++){
      $arr[$i] =~ s/\s//g; $arr[$i] =~ s/\t//g; $arr[$i] =~ s/\r//g; $arr[$i] =~ s/\n//g;
      if($arr[$i] ne '-'){
	$st = $i+1;
	$i = scalar @arr;
      }
    }
    # end position
    for(my $j = (scalar @arr -1); $j > 0; $j--){
      $arr[$j] =~ s/\s//g; $arr[$j] =~ s/\t//g; $arr[$j] =~ s/\r//g; $arr[$j] =~ s/\n//g;
      if($arr[$j] ne '-'){
	$end = $j+1;
	$j = 0;
      }
    }
    return($st, $end);
}


###################
# Get the nt characters in ref sequence 
####################
sub get_ref{
  my($href) = @_;
  my %align_pos = %$href;
  my %out = ();
  
  my @sorted_pos = sort {$a <=> $b} keys %align_pos;
  my $align_st = shift @sorted_pos;
  my $align_end = pop @sorted_pos;
  
  print "align start $align_st  align end $align_end\n"; 
  my $ct = 1; # count for orig pos

  #####################################
  my $insertion_count = 0;
  for(my $k = $align_st; $k <= $align_end; $k++){

    my $char = $align_pos{$k};
    $char =~ s/\s//g;
    if($char ne '-'){
      # not a gap, process
      $insertion_count = 0; # reset the insertion counter;
      if($ct <= 1){
	$out{$k} = "1";
	$ct = $ct+1;
      }
      else{
	$out{$k} = "$ct";
	$ct = $ct+1;
      }
    }
    else{
      # place in alignment with possible insertions
      my $st_pos = "";
      my $end_pos = "";

      # increment the insertion counter for each position where there is a gap
      $insertion_count = $insertion_count + 1;

      # find start of insertion
      for(my $i = $k; $i >= $align_st; $i--){
	if($i == 1){
	  $st_pos = "1";
	}
	else{
	  if($align_pos{$i} ne '-'){
	    if(exists $out{$i}){
	      $st_pos = $out{$i};
	      last;
	    }
	  }
	}
      }
      for(my $i = $k; $i <= $align_end; $i++){
	if($align_pos{$i} ne '-'){
	  # found the end point
	  $end_pos = "$ct";
	  last;	
	}
      }
      $out{$k} = "$st_pos" . "-" . "$end_pos" . "." . "$insertion_count";
    }
  }
  return %out;
}
