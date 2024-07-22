#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use File::Basename;

# load 5,000x subsampled results
my $inDir  = "../analysis/A03.0-create-freeze/normalized.5000";
my $outDir = "../analysis/A03.1-add-metadata";
my @files = qw/jhmi-esca-microbiota.alpha-diversity.txt
jhmi-esca-microbiota.phylum.txt
jhmi-esca-microbiota.class.txt
jhmi-esca-microbiota.order.txt
jhmi-esca-microbiota.family.txt
jhmi-esca-microbiota.genus.txt
jhmi-esca-microbiota.species-otu.txt.gz
jhmi-esca-microbiota.picrust.txt/;

my $inMeta = "../data/meta.2022-11-07.txt";
my %meta = ();
my @mhead = ();
open IN, "$inMeta" or die;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if (!defined($mhead[1])){
    @mhead = @A;
  }else{
    # only feces to be included --
    next if ($A[4] ne "Feces");
    for my $i (1 .. $#mhead){
      if (!defined($A[$i]) or $A[$i] eq ""){
        $A[$i] = "n/a";
      }
      $meta{$A[0]}{$mhead[$i]} = $A[$i];
    }
  }
}
close IN;

if (-e "$outDir"){
  `rm -r $outDir`;
}
`mkdir -p $outDir`;
my %spPrct = ();
foreach my $f (@files){
  print "$f...\n";
  my @header = ();
  if ($f !~ /gz/){
    open IN, "$inDir/$f" or die;
  }else{
    open IN, "gunzip -c $inDir/$f |" or die;
  }
  open OUT, ">$outDir/$f" or die;
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
    if (!defined($header[1])){
      @header = @A;
      print OUT "SampleID";
      for my $mi (1 .. $#mhead){
        print OUT "\t$mhead[$mi]";
      }
      for my $i (1 .. $#A){
        print OUT "\t$header[$i]";
      }
      print OUT "\n";
    }else{
      if (!defined($meta{$A[0]})){
        next;
      }
      print OUT "$A[0]";
      for my $mi (1 .. $#mhead){
        print OUT "\t$meta{$A[0]}{$mhead[$mi]}";
      }
      for my $i (1 .. $#A){
        print OUT "\t$A[$i]";
        if ($f =~ /species/){
          push @{$spPrct{$header[$i]}}, $A[$i];
        }
      }
      print OUT "\n";
    }
  }
  close IN;
  close OUT;

  if ($f =~ /species/){
    my %gte0pt1Prct = ();
    foreach my $k1 (sort keys %spPrct){
      my @k1 = @{$spPrct{$k1}};
      if (average(@k1) >= 0.01){
        $gte0pt1Prct{$k1} = 1;
      }
    }
    open OUT, ">$outDir/jhmi-esca-microbiota.species-otu.avgGTE0.01Prct.txt" or die;
    open IN, "$outDir/$f" or die;
    my @spheader = ();
    while(<IN>){
      chomp($_);
      my @A = split "\t", $_;
      if (!defined($spheader[1])){
        @spheader = @A;
      }
      print OUT "$A[0]";
      for my $i (1 .. $#A){
        if ($i <= $#mhead or defined($gte0pt1Prct{$spheader[$i]})){
          print OUT "\t$A[$i]";
        }
      }
      print OUT "\n";
    }
    close IN;
    close OUT;
  }
}

# add beta-diversity distance matrices with this sample set
my @bfiles = qw/jhmi-esca-microbiota.beta-diversity.bray_curtis.txt
jhmi-esca-microbiota.beta-diversity.unweighted_unifrac.txt
jhmi-esca-microbiota.beta-diversity.weighted_unifrac.txt/;
foreach my $f (@bfiles){
  print "$f...\n";
  my @header = ();
  open IN, "$inDir/$f" or die;
  open OUT, ">$outDir/$f" or die;
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
    if (!defined($header[1])){
      @header = @A;
      for my $i (1 .. $#A){
        if (defined($meta{$header[$i]})){
          print OUT "\t$header[$i]";
        }
      }
      print OUT "\n";
    }else{
      if (defined($meta{$A[0]})){
        print OUT "$A[0]";
        for my $i (1 .. $#A){
          if (defined($meta{$header[$i]})){
            print OUT "\t$A[$i]";
          }
        }
        print OUT "\n";
      }
    }
  }
  close IN;
  close OUT;
} # end of beta-diversity files

# compile into a single Excel file
my $dex = 1;
foreach my $f (@files){
  if ($f =~ /species/){
    $f = "jhmi-esca-microbiota.species-otu.avgGTE0.01Prct.txt";
  }
  my $str = $f;
  $str =~ s/jhmi-esca-microbiota.//g;
  $str =~ s/\.txt//g;
  $str = "S$dex.0 | $str";
  print "Excel file $f...$str\n";
  `./lib/excel2csv -p -S \"$str\" $outDir/$f $outDir/jhmi-esca-microbiota.merged.xlsx`;
  $dex++;
}

sub average
{
  my (@B)       = @_;
  if (!defined($B[0])){
    return("na");
  }
  my $sum       = 0;
  foreach my $s (@B){
    $sum       += $s;
  }
  my $favg      = $sum/($#B+1);
  return($favg);
}
