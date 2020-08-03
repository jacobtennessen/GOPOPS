#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_s $opt_l $opt_c $opt_m );

# Usage
my $usage = "
SimulatePools.pl - simulates Fst differences due to chance in pooled sequencing
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

Usage: perl SimulatePools.pl options
 optional:
  -s  sample size for each pool to simulate [default = 600]
  -l  number of loci to simulate [default = 100,000]
  -c  sequencing coverage to simulate [default = 500]
  -m  minimum allele frequency to simulate [default = 0.2]

";

#############

# command line processing.
getopts('s:l:c:m:');

my ($samplesize, $locicount, $coverage, $minfreq);

if (defined $opt_s) {
  $samplesize = $opt_s;
} else {
  $samplesize = 600;
}
if (defined $opt_l) {
  $locicount = $opt_l;
} else {
  $locicount = 100000;
}
if (defined $opt_c) {
  $coverage = $opt_c;
} else {
  $coverage = 500;
}
if (defined $opt_m) {
  $minfreq = $opt_m;
} else {
  $minfreq = 0.2;
}

my %fsts;

my @notablespots = (0.5);

for (my $p = 10; $p <= $locicount; $p*= 10) {
    my $n = 1 - 1/$p;
    push @notablespots, $n;
}

my $rf = 0.5 - $minfreq;

for (my $l = 0; $l < $locicount; $l++) {
    my $freq = (rand($rf)) + $minfreq;
    my $pop1 = 0;
    my $pop2 = 0;
    for (my $p1 = 0; $p1 < (2*$samplesize); $p1++) {
        my $allele1 = rand(1/$freq);
        if ($allele1 <= 1) {
            $pop1 += 1;
        }
        my $allele2 = rand(1/$freq);
        if ($allele2 <= 1) {
            $pop2 += 1;
        }
    }
    my $freq1 = $pop1/(2*$samplesize);
    my $freq2 = $pop2/(2*$samplesize);
    my $reads1 = 0;
    my $reads2 = 0;
    for (my $c = 0; $c < $coverage; $c++) {
        my $r1 = rand(1/$freq1);
        if ($r1 <= 1) {
            $reads1 += 1;
        }
        my $r2 = rand(1/$freq2);
        if ($r2 <= 1) {
            $reads2 += 1;
        }
    }
    my $finalfreq1 = $reads1/$coverage;
    my $finalfreq2 = $reads2/$coverage;
    my $div = sprintf "%.3f", (Fst($coverage,$coverage,$finalfreq1,$finalfreq2));
    if (defined $fsts{$div}) {
        $fsts{$div} += 1;
    } else {
        $fsts{$div} = 1;
    }
}

my @fsts = sort by_number (keys %fsts);

my $tail = 0;
my %seenspots;

foreach my $fk (@fsts) {
  my $prop = $fsts{$fk}/$locicount;
  $tail += $prop;
	foreach my $n (@notablespots) {
	    if (($tail >= $n)&&(!defined $seenspots{$n})) {
            print "$n lower bound for neutral at $fk.\n";
            $seenspots{$n} = 1;
	    }
	}
}

#################################

sub by_number {
    if ($a < $b) {-1} elsif ($a > $b) {1} else {0}
}

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
