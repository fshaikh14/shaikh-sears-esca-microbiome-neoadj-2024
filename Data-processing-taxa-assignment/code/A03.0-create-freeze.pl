#!/usr/bin/perl
use Data::Dumper;
use Getopt::Std;
use Scalar::Util qw/looks_like_number/;
use List::Util qw/shuffle sum min max/;
use POSIX qw/ceil floor/;
use File::Basename;
use warnings;
use strict;
# --------------------------------------------------------------------
#  Author: James Robert White, PhD
#  Email: jwhite@respherabio.com
# --------------------------------------------------------------------
#  Create freeze 16S rRNA dataset
# --------------------------------------------------------------------
my $outDir = "../analysis/".basename($0); $outDir =~ s/\.pl$//g;
if (-e $outDir){
  `rm -r $outDir`;
}
`mkdir $outDir`;
# --------------------------------------------------------------------
my $PREFIX    = "jhmi-esca-microbiota";
# --------------------------------------------------------------------
my $freezeDir = "$outDir";
# --------------------------------------------------------------------
my $extDir = "../external/RI.2022-08-08/insight_results";
# --------------------------------------------------------------------
# full
`mkdir $freezeDir/total-counts`;
# otu.biom
`cp $extDir/otu_table.biom.gz $freezeDir/total-counts/$PREFIX.biom.gz`;
# picrust.biom
`cp $extDir/picrust.predictions.biom.gz $freezeDir/total-counts/$PREFIX.picrust.biom.gz`;
# reps.fna
`cp $extDir/insight.reps.fna.gz $freezeDir/total-counts/$PREFIX.reps.fna.gz`;
# summary
`cp $extDir/otu_table.summary.txt $freezeDir/total-counts/$PREFIX.stat-summary.txt`;
# otu table
`cp $extDir/otu_table.txt.gz $freezeDir/total-counts/$PREFIX.counts.species-otu.txt.gz`;
# summarize taxa L2, L3, L4, L5, L6
my %Ls = ("L2" => "phylum",
          "L3" => "class",
          "L4" => "order",
          "L5" => "family",
          "L6" => "genus" );
foreach my $k1 (sort keys %Ls){
  `cp $extDir/summarize-taxa/otu_table_$k1.txt $freezeDir/total-counts/$PREFIX.counts.$Ls{$k1}.txt`;
}

# normalized 5000 ----------------------------------------------------------------------------
`mkdir $freezeDir/normalized.5000`;
$extDir = "../external/RI.2022-08-08/insight_results/rd5000";
# biom
`cp $extDir/otu_table_even5000.biom.gz $freezeDir/normalized.5000/$PREFIX.biom.gz`;
# reps
`cp $extDir/reps_even5000.fna.gz $freezeDir/normalized.5000/$PREFIX.reps.fna.gz`;
# tree
`cp $extDir/reps_even5000.tre.gz $freezeDir/normalized.5000/$PREFIX.tre.gz`;
# beta
`cp $extDir/beta_diversity/distance_matrices/bray_curtis_dm.txt $freezeDir/normalized.5000/$PREFIX.beta-diversity.bray_curtis.txt`;
`cp $extDir/beta_diversity/distance_matrices/unweighted_unifrac_dm.txt $freezeDir/normalized.5000/$PREFIX.beta-diversity.unweighted_unifrac.txt`;
`cp $extDir/beta_diversity/distance_matrices/weighted_unifrac_dm.txt $freezeDir/normalized.5000/$PREFIX.beta-diversity.weighted_unifrac.txt`;

# alpha
format_to_columns("$extDir/alpha_diversity/alphadiv_even5000.metrics_x_samples.txt", "$freezeDir/normalized.5000/$PREFIX.alpha-diversity.txt", "no");

# picrust tables
format_to_columns("$extDir/summarize_picrust/ko_kegglevel3.txt", "$freezeDir/normalized.5000/$PREFIX.picrust.txt", "yes");

# taxa tables
foreach my $k1 (sort keys %Ls){
  format_to_columns("$extDir/summarize_taxa/qiime_out/otu_table_even5000_$k1.txt", "$freezeDir/normalized.5000/$PREFIX.$Ls{$k1}.txt", "yes");
}
# otu table
format_to_columns("$extDir/summarize_taxa/qiime_out/otu_table_even5000.txt", "$freezeDir/normalized.5000/$PREFIX.species-otu.txt", "yes");

`gzip $freezeDir/normalized.5000/$PREFIX.species-otu.txt`;

# --------------------------------------------------------------------------------
sub format_to_columns
{
  my ($input, $output, $prct) = @_;
  # convert rows to columns ---------
  # format names for easier eval in R
  # optionally convert final columns into % abundance
  # values yes or no
  if ($prct ne "yes" and $prct ne "no"){
    die "Error format_to_columns needs a yes or a no on % values.\n";
  }
  open IN, "$input" or die;
  my %data     = ();
  my %total    = ();
  my @header   = ();
  my @features = ();
  while(<IN>){
    chomp($_);
    next if ($_ =~ /\#\ Constructed\ from\ biom\ file/);
    my @A = split "\t", $_;
    if (!defined($header[1])){
      @header    = @A;
      if ($header[$#header] eq "taxonomy"){
        pop @header;
      }
    }else{
      my $feature = $A[0];
      $feature    =~ s/\ //g;
      $feature    =~ s/\;/\./g;
      $feature    =~ s/\:/\./g;
      $feature    =~ s/\_\_/\_/g;
      $feature    =~ s/\-/\_/g;
      $feature    =~ s/\///g;
      $feature    =~ s/\.\./\./g;
      $feature    =~ s/\_\_/\_/g;
      $feature    =~ s/\_\_/\_/g;

      # print "$feature\n";
      push @features, $feature;
      for my $i (1 .. $#header){
        $data{$header[$i]}{$feature} = $A[$i];
        $total{$header[$i]}         += $A[$i];
      }
    }
  }
  close IN;

  if ($prct eq "yes"){
    foreach my $k1 (sort keys %data){
      foreach my $f (@features){
        $data{$k1}{$f} = 100*$data{$k1}{$f}/$total{$k1};
      }
    }
  }

  my @colorder = (@features);
  open OUT, ">$output" or die;
  print OUT "SampleID";
  foreach my $k2 (@colorder){
    print OUT "\t$k2";
  }
  print OUT "\n";

  foreach my $k1 (sort keys %data){
    print OUT "$k1";
    foreach my $k2 (@colorder){
      print OUT "\t$data{$k1}{$k2}";
    }
    print OUT "\n";
  }
  close OUT;
} # end format_to_columns
