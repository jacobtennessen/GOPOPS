#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw( $opt_g $opt_a $opt_d $opt_n $opt_o );

# Usage
my $usage = "
MakeFreqTableFromPooledPileup.pl - reads a pileup file of pooled samples and calculates allele frequencies
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

Usage: perl MakeFreqTableFromPooledPileup.pl options
 required:
  -g	a pileup file PILEUPNAME. All samples are considered to be pools
 optional:
  -a	minimum mean allele frequency across all samples [default = 0.01]
  -d  mimimum mean depth across all samples [default = 2]
  -n  comma-delimited list of sample names in order
  -o  output file name [default is 'Freqs_PILEUPNAME in same directory as PILEUPNAME ]

";

#############

# command line processing.
getopts('g:a:d:n:o:');
die $usage unless ($opt_g);

my ($pileup, $minfreq, $mindepth, $names, $outputfile);

$pileup	= $opt_g if $opt_g;
my @pileupfiledata = split /\//, $pileup;
my $genotypes = pop @pileupfiledata;

if (defined $opt_a) {
  $minfreq = $opt_a;
} else {
  $minfreq = 0.01;
}
if (defined $opt_d) {
  $mindepth = $opt_d;
} else {
  $mindepth = 2;
}
$names = $opt_n if $opt_n;
$outputfile = $opt_o if $opt_o;

if (defined $pileupfiledata[0]) {
    my $dir = join "/", @pileupfiledata;
    unless (defined $outputfile) {
        $outputfile = "$dir/Freqs_$genotypes";
    }
    $genotypes = "$dir/$genotypes";
}

unless (defined $outputfile) {
    $outputfile = "Freqs_$genotypes";
}

my $highestallelenumber = 2; #maximum number of alleles allowed per SNP.

my $dumplimit = 100000;
    
my $total; #number of columns in pileup file

my @indnumbers; #columns with genotype info for each individual

my $goodline = 0;

my $indcount; #number of individuals

my @out;

my $title;

open(IN, "$genotypes") || die "can't open $genotypes\n";

while (<IN>) {
    my $line = $_;
    $line =~ s/\r|\n//g;
    my @data = split "\t", $line;
    unless (defined $total) {
        $total = scalar (@data);
        $indcount = int( (($total - 3)/3) + 0.99);
        for (my $x = 3; $x < $total; $x +=3) {
            push @indnumbers, $x;
        }
    }
    my %allseenbases;
    $allseenbases{$data[2]} = 1;
    my $refbase = $data[2];
    my @freqs;
    my @depths;
    my $chromsite = "$data[0]_$data[1]";
    foreach my $id (@indnumbers) {
        my $goodcount = 0;
        my $freq = 0;
        if (defined $data[$id+1]) {
            my @seqdata = split "", $data[$id+1];
            my $extra = 0;
            my $indel = 0;
            my $postindel = 0;
            my %bases;
            my $allgood = 0;
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
                    if (defined $bases{$data[2]}) {
                        $bases{$data[2]} +=1;
                    } else {
                        $bases{$data[2]} = 1;
                    }
                    $allgood +=1;
                } elsif (($s =~ /A/gi)||($s =~ /C/gi)||($s =~ /G/gi)||($s =~ /T/gi)) {
                    my $uc = uc $s;
                    if (defined $bases{$uc}) {
                        $bases{$uc} +=1;
                    } else {
                        $bases{$uc} = 1;
                    }
                    $allgood +=1;
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
            my $altcount = 0;
            foreach my $b (keys %bases) {
                unless (defined $allseenbases{$b}) {
                    my $bcount = (scalar(keys %allseenbases))+1;
                    $allseenbases{$b} = $bcount;
                }
                unless ($b =~ /^$refbase$/) {
                    $altcount += $bases{$b};
                }
                $goodcount += $bases{$b};
            }
            if ($goodcount > 0) {
                $freq = sprintf "%.3f", $altcount/$goodcount;
            }
        }
        push @freqs, $freq;
        push @depths, $goodcount;
    }
    if ((scalar(keys %allseenbases) <= 1)||(scalar(keys %allseenbases) > $highestallelenumber)) {
        next;
    }
    my $depthtotal = 0;
    foreach my $d (@depths) {
        $depthtotal +=$d;
    }
    my $meandepth = $depthtotal/$indcount;
    if ($meandepth < $mindepth) {
        next;
    }
    my $freqtotal = 0;
    foreach my $f (@freqs) {
        $freqtotal +=$f;
    }
    my $meanfreq = $freqtotal/$indcount;
    if (($meanfreq < $minfreq)||($meanfreq > (1-$minfreq))) {
        next;
    }
    my @outline;
    for (my $i = 0; $i < $indcount; $i++) {
        if ($meanfreq > 0.5) {
            $freqs[$i] = 1 - $freqs[$i];
        }
        push @outline, "$freqs[$i]\t$depths[$i]";
    }
    my $alt;
    foreach my $a (keys %allseenbases) {
        unless ($a =~ /^$refbase$/) {
            $alt = $a;
        }
    }
    my $alleles = "$refbase,$alt";
    if ($meanfreq > 0.5) {
        $alleles = "$alt,$refbase";
    }
    if ($meanfreq > 0.5) {
        $meanfreq = 1 - $meanfreq;
    }
    $meanfreq = sprintf "%.3f", $meanfreq;
    my $outline = join "\t", @outline;
    push @out, "$chromsite\t$alleles\t$outline\t$meanfreq";
    $goodline +=1;
    if ($goodline >= $dumplimit) {
        if (defined $title) {
            my $result = join "\n", @out;
            unless ( open(OUT, ">>$outputfile") ) {
                print "Cannot open file \"$outputfile\" to write to!!\n\n";
                exit;
            }
            print OUT "\n$result";
            close (OUT);
        } else {
            my @title = ("Marker","Alleles");
            if (defined $names) {
                my @names = split ",", $names;
                foreach my $n (@names) {
                    push @title, "Freq_$n\tDepth_$n";
                }
            } else {
                for (my $i = 1; $i <= $indcount; $i++) {
                    push @title, "Freq_Sample$i\tDepth_Sample$i";
                }
            }
            push @title, "MeanFreq";
            $title = join "\t", @title;
            unshift @out, $title;
            
            my $result = join "\n", @out;
            unless ( open(OUT, ">$outputfile") ) {
                print "Cannot open file \"$outputfile\" to write to!!\n\n";
                exit;
            }
            print OUT $result;
            close (OUT);
        }
        @out = ();
        $goodline = 0;
    }
}

close (IN);

if (defined $out[0]) {
    if (defined $title) {
        my $result = join "\n", @out;
        unless ( open(OUT, ">>$outputfile") ) {
            print "Cannot open file \"$outputfile\" to write to!!\n\n";
            exit;
        }
        print OUT "\n$result";
        close (OUT);
    } else {
        my @title = ("Marker","Alleles");
        if (defined $names) {
            my @names = split ",", $names;
            foreach my $n (@names) {
                push @title, "Freq_$n\tDepth_$n";
            }
        } else {
            for (my $i = 1; $i <= $indcount; $i++) {
                push @title, "Freq_Sample$i\tDepth_Sample$i";
            }
        }
        push @title, "MeanFreq";
        $title = join "\t", @title;
        unshift @out, $title;
        
        my $result = join "\n", @out;
        unless ( open(OUT, ">$outputfile") ) {
            print "Cannot open file \"$outputfile\" to write to!!\n\n";
            exit;
        }
        print OUT $result;
        close (OUT);
    }
    @out = ();
    $goodline = 0;
}
