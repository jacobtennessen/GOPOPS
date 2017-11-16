#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_g $opt_d $opt_m );

# Usage
my $usage = "
AlleleFreqsFromPileup.pl - reads a pileup file and returns within-sample frequencies of putative SNPs
Copyright (C) 2017 by Jacob A Tennessen

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

Usage: perl AlleleFreqsFromPileup.pl options
 required:
  -g	a pileup file
 optional:
  -d	minimum depth to call genotypes [default = 8]
  -m	maximum number of samples with missing data [default = 1]

";

#############

# command line processing.
getopts('g:d:m:');
die $usage unless ($opt_g);

my ($pileup, $mindepth, $missingthresh);

$pileup	= $opt_g if $opt_g;
if (defined $opt_d) {
    $mindepth = $opt_d;
} else {
    $mindepth = 8; 
}
if (defined $opt_m) {
    $missingthresh = $opt_m;
} else {
    $missingthresh = 1; 
}

my @pileupfiledata = split /\//, $pileup;

my $genotypes = pop @pileupfiledata;

my $outputfile = "Allele_Frequencies_$genotypes";

if (defined $pileupfiledata[0]) {
    my $dir = join "/", @pileupfiledata;
    $genotypes = "$dir/$genotypes";
    $outputfile = "$dir/$outputfile";
}

my $total; #number of columns in pileup file

my @indnumbers; #columns with genotype info for each individual

my $indcount; #number of individuals

my @out; #allele frequency output

open(IN, "$genotypes") || die "can't open $genotypes\n";

while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    unless (defined $total) {
        $total = scalar (@data);
        $indcount = ($total - 3)/3;
        my @title = ("Allele");
        my $titleindcount = 1;
        for (my $x = 3; $x < $total; $x +=3) {
            push @indnumbers, $x;
            push @title, "Freq$titleindcount\tSize$titleindcount";
            $titleindcount +=1;
        }
        my $title = join "\t", @title;
        push @out, $title;
    }
    my %HoHbasesperind;
    my $indtrack = 1;
    my @localbases;
    my %seenlocalbases;
    my %allgood;
    push @localbases, $data[2];
    $seenlocalbases{$data[2]} = 1;
    foreach my $id (@indnumbers) {
        $allgood{$indtrack} = 0;
        if ($data[$id] >= $mindepth) {
            my @seqdata = split "", $data[$id+1];
            my $extra = 0;
            my $indel = 0;
            my $postindel = 0;
            foreach my $s (@seqdata) {
                if ($indel == 1) {
                    if ($s =~ /\d/) {
                        $extra = $s;
                    } else {
                        next;
                    }
                    $indel = 0;
                    $postindel = 1;
                } elsif ($extra > 0) {
                    if ($postindel == 1) {
                        if ($s =~ /\d/) { #indel multiple digits
                            $extra = ($extra*10)+$s;
                        } else {
                            $postindel = 0;
                            $extra -=1;
                        }
                    } else {
                        $extra -=1;
                    }
                } elsif (($s =~ /,/)||($s =~ /\./)) {
                    if (defined $HoHbasesperind{$data[2]}{$indtrack}) {
                        $HoHbasesperind{$data[2]}{$indtrack} +=1;
                    } else {
                       $HoHbasesperind{$data[2]}{$indtrack} = 1;
                    }
                    $allgood{$indtrack} +=1;
                } elsif (($s =~ /A/gi)||($s =~ /C/gi)||($s =~ /G/gi)||($s =~ /T/gi)) {
                    my $uc = uc $s;
                    unless (defined $seenlocalbases{$uc}) {
                        push @localbases, $uc;
                        $seenlocalbases{$uc} = 1;
                    }
                    if (defined $HoHbasesperind{$uc}{$indtrack}) {
                        $HoHbasesperind{$uc}{$indtrack} +=1;
                    } else {
                        $HoHbasesperind{$uc}{$indtrack} = 1;
                    }
                    $allgood{$indtrack} +=1;
                } elsif ($s =~ /\^/) {
                    $extra = 1;
                } elsif (($s =~ /-/)||($s =~ /\+/)) {
                    $indel = 1;
                } elsif (($s =~ /\$/)||($s =~ /\*/gi)||($s =~ /</gi)||($s =~ />/gi)||($s =~ /N/gi)) {
                    next;
                } else {
                    next;
                }
            }
        }
        $indtrack +=1;
    }
    my @goodfreqs;
    foreach my $base (@localbases) {
        my @freqs;
        my $freqsum = 0;
        my $goodtally = 0;
        for (my $i = 1; $i <= $indcount; $i++) {
            unless (defined $HoHbasesperind{$base}{$i}) {
                $HoHbasesperind{$base}{$i} = 0;
            }
            my $freq = "-";
            if ($allgood{$i} >= $mindepth) {
                $freq = sprintf "%.4f", ($HoHbasesperind{$base}{$i}/$allgood{$i});
                $freqsum += $freq;
                $goodtally +=1;
            }
            push @freqs, "$freq\t$allgood{$i}";
        }
        if ($goodtally >= ($indcount - $missingthresh)) {
            my $meanfreq = $freqsum/$goodtally;
            if (($meanfreq > 0)&&($meanfreq < 1)) {
                my $freqlist = join "\t", @freqs;
                my $outline = "$data[0]_$data[1]_$base\t$freqlist";
                push @goodfreqs, $outline;
            }
        }
    }
    if ((scalar(@goodfreqs)) > 1) {
        my $reject = shift @goodfreqs;
    }
    if (defined $goodfreqs[0]) {
        my $retained = join "\n", @goodfreqs;
        push @out, $retained;
    }
}

close (IN);

my $result = join "\n", @out;

unless ( open(OUT, ">$outputfile") ) {
    print "Cannot open file \"$outputfile\" to write to!!\n\n";
    exit;
}

print OUT $result;

close (OUT);
