#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_f $opt_a $opt_b $opt_d $opt_n $opt_w);

# Usage
my $usage = "
FstFromJoinedFreqTablesWindow.pl - reads an joined allele frequency table ala multiple outputs of MakeFreqTableFromPooledPileup.pl and calculates mean pairwise Fst per window
Copyright (C) 2020 by Jacob A Tennessen 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: perl FstFromJoinedFreqTablesWindow.pl options
 required:
  -f	an allele frequency table
  -a  comma-delimited list of sample freq positions [zero-based numbering], first group
  -b  comma-delimited list of sample freq positions [zero-based numbering], second group
 optional:
  -d  mimimum depth per sample [default = 5]
  -n  minimum number of samples per group [default = 5]
  -w  window size [default = 10000]

";

#############

# command line processing.
getopts('f:a:b:d:n:w:');
die $usage unless ($opt_f);
die $usage unless ($opt_a);
die $usage unless ($opt_b);

my ($freqfile, $grouppos1, $grouppos2, $mindepth, $minsamp, $window);

$freqfile = $opt_f if $opt_f;
$grouppos1 = $opt_a if $opt_a;
$grouppos2 = $opt_b if $opt_b;
if (defined $opt_d) {
  $mindepth = $opt_d;
} else {
  $mindepth = 5;
}
if (defined $opt_n) {
  $minsamp = $opt_n;
} else {
  $minsamp = 5;
}
if (defined $opt_w) {
  $window = $opt_w;
} else {
  $window = 10000;
}

my @pileupfiledata = split /\//, $freqfile;

my $genotypes = pop @pileupfiledata;

my $outputfile = "Mean_Fsts_Windows_$window"."_$genotypes";

if (defined $pileupfiledata[0]) {
    my $dir = join "/", @pileupfiledata;
    $genotypes = "$dir/$genotypes";
    $outputfile = "$dir/$outputfile";
}

my @g1f = split ",", $grouppos1;

my @g2f = split ",", $grouppos2;

my $total;

my %HoHoAfsts;

my %HoHoAheterozygosity;

open(IN, "$genotypes") || die "can't open $genotypes\n";

while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split /\s+/, $line;
    unless (defined $total) {
        $total = scalar (@data);
    }
    my @switch;
    my $switch = 0;
    foreach my $d (@data) {
        if ($d =~ /,/) {
            if ($d =~ /$data[1]/) {
                $switch = 0;
            } else {
               $switch = 1; 
            }
        }
        if ($switch == 1) {
            push @switch, 1;
        } else {
            push @switch, 0;
        }
    }
    my $group1alleletotal = 0;
    my $good1samples = 0;
    my $group1depth = 0;
    foreach my $g1 (@g1f) {
        if (($data[$g1+1] =~ /\d/)&&($data[$g1+1] >= $mindepth)) {
            my $allelecount1 = $data[$g1]*$data[$g1+1];
            if ($switch[$g1] == 1) {
                $allelecount1 = (1-$data[$g1])*$data[$g1+1];
            }
            $group1alleletotal += $allelecount1;
            $good1samples +=1;
            $group1depth += $data[$g1+1];
        }
    }
    my $group2alleletotal = 0;
    my $good2samples = 0;
    my $group2depth = 0;
    foreach my $g2 (@g2f) {
        if (($data[$g2+1] =~ /\d/)&&($data[$g2+1] >= $mindepth)) {
            my $allelecount2 = $data[$g2]*$data[$g2+1];
            if ($switch[$g2] == 1) {
                $allelecount2 = (1-$data[$g2])*$data[$g2+1];;
            }
            $group2alleletotal += $allelecount2;
            $good2samples +=1;
            $group2depth += $data[$g2+1];
        }
    }
    if (($good1samples >= $minsamp)&&($good2samples >= $minsamp)) {
        my $depthratio = $group1depth/$group2depth;
        if (($depthratio > (2/3))&&($depthratio < (3/2))) {
          my $group1freq = $group1alleletotal/$group1depth;
          my $group2freq = $group2alleletotal/$group2depth;
          my $fstat = Fst ($group1depth,$group2depth,$group1freq,$group2freq);
          my $meanfreq = ($group1freq+$group2freq)/2;
          my $het = 2*$meanfreq*(1-$meanfreq);
          my @chromdata = split "_", $data[0];
          my @shortchrom = split "F", $chromdata[0];
          my $rounded1 = (int($chromdata[1]/$window))*$window;
          my $rounded2 = (int(($chromdata[1]-($window/2))/$window))*$window + ($window/2);
          push @{$HoHoAfsts{$shortchrom[0]}{$rounded1}}, $fstat;
          push @{$HoHoAfsts{$shortchrom[0]}{$rounded2}}, $fstat;
          push @{$HoHoAheterozygosity{$shortchrom[0]}{$rounded1}}, $het;
          push @{$HoHoAheterozygosity{$shortchrom[0]}{$rounded2}}, $het;
        }
    }
}

close (IN);

my @out;

foreach my $chrom (sort (keys %HoHoAfsts)) {

  foreach my $window (sort (keys %{$HoHoAfsts{$chrom}})) {
    my $wtotalfst = 0;
    my $wcountfst = 0;
    foreach my $f (@{$HoHoAfsts{$chrom}{$window}}) {
      $wtotalfst += $f;
      $wcountfst +=1;
    }
    my $wtotalhet = 0;
    my $wcounthet = 0;
    foreach my $h (@{$HoHoAheterozygosity{$chrom}{$window}}) {
      $wtotalhet += $h;
      $wcounthet +=1;
    }
    if ($wcountfst > 0) {
      my $meanfst = $wtotalfst/$wcountfst;
      my $meanhet = $wtotalhet/$wcounthet;
      push @out, "$chrom\t$window\t$meanfst\t$wcountfst\t$meanhet";
    }
  }
}

if (defined $out[0]) {
    my $result = join "\n", @out;
    unless ( open(OUT, ">$outputfile") ) {
        print "Cannot open file \"$outputfile\" to write to!!\n\n";
        exit;
    }
    print OUT "Scaffold\tWindow\tFst\tSites\tHetE\n$result";
    close (OUT);
}


#########################

sub Fst {
    my ($size1, $size2, $freq1, $freq2) = @_;
    
    my @freqs = ($freq1, $freq2);
    my @samplesizes = ($size1, $size2);
    
    my $num_sub_pops = 2;

    my $Fst;
    my ($TS_sub1,$TS_sub2); # numerator and denominator for theta (?)
    my @alleles = (0,1);

    my $avg_samp_size         = 0; # n-bar
    my $avg_allele_freq       = 0; # p-tilda-A-dot
    my $total_samples_squared = 0; # s-squared?

    foreach (my $c = 0; $c < 2; $c ++) {
	my $s = $samplesizes[$c];
	$avg_samp_size += $s;
	$total_samples_squared += $s**2;
	my $all_freq = $freqs[$c];
	$avg_allele_freq += $s * $all_freq;
    }
    
    my $total_samples =  $avg_samp_size;	# sum of n over i sub-populations
    $avg_samp_size /= $num_sub_pops;
    $avg_allele_freq /= $total_samples;

    my $adj_samp_size = ( 1/ ($num_sub_pops - 1)) * ( $total_samples - ( $total_samples_squared/$total_samples));   # n-sub-c

    my $variance              = 0; # s-squared-sub-A
    my $sum_variance          = 0;
    my $i = 0;		# we have cached the marker info
    for (my $d = 0; $d <2; $d ++) {
	my $s = $samplesizes[$d];
	$sum_variance += $s * (($freqs[$d] - $avg_allele_freq)**2);
    }
    $variance = ( 1 / (( $num_sub_pops-1)*$avg_samp_size))*$sum_variance;

	       $TS_sub1 = $variance - 
		   ( ( 1/($avg_samp_size-1))*
		     ( ($avg_allele_freq*(1-$avg_allele_freq))-
		       ( (($num_sub_pops-1)/$num_sub_pops)*$variance)));
	       $TS_sub2 = ( (($adj_samp_size-1)/($avg_samp_size-1))*
			      $avg_allele_freq*(1-$avg_allele_freq) ) +
			      ( 1 + ( (($num_sub_pops-1)*
				       ($avg_samp_size-$adj_samp_size))/ 
				      ($avg_samp_size - 1))) * 
				      ($variance/$num_sub_pops);
    
    unless (($freq1 == $freq2)&&(($freq1 == 0)||($freq1 == 1))) {
        $Fst = $TS_sub1 / $TS_sub2;
    } else {
        $Fst = 0;
    }
   
    return $Fst;
}
