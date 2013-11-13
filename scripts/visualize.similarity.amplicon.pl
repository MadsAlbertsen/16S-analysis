#!/usr/bin/env perl
###############################################################################
#
#    visualize.similarity.amplicon.pl
#
#	 Takes an input fastafile and calculates the number of 1 
#    nucleotide mismatches to the X most abundant sequences.
#    Outputs a network connection overview.
#    
#    Copyright (C) 2013 Mads Albertsen
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use strict;
use warnings;

#core Perl modules
use Getopt::Long;

#locally-written modules
BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# get input params
my $global_options = checkParams();

my $inseqs;
my $networkstats;
my $nons;
my $ncatch;

$inseqs = &overrideDefault("inseqs.txt",'inseqs');
$nons = &overrideDefault(0,'nons');
$networkstats = &overrideDefault("networkstats.txt",'networkstats');
$ncatch = &overrideDefault("ncatch.txt",'ncatch');

my %seqs;
my $header = "";
my $seq;
my %printed;
my $countn = 0;
my $count = 0;
my %ctype;
my %stype;
my %missing;
my $remove = 0;

######################################################################
# CODE HERE
######################################################################

open(OUT, ">$networkstats") or die("Cannot create file: $networkstats\n");
open(OUT2, ">$ncatch") or die("Cannot create file: $ncatch\n");

open(INseqs, $inseqs) or die("Cannot read file: $inseqs\n");                                       #Add the sequences to a hash table
while ( my $line = <INseqs> ) {
	$line =~ s/\r\n//g;
	if ($line =~ m/>/) { 	
		$count++;
		if ($header ne ""){
			if ($seq !~ m/N/ or $nons == 0){
				$seqs{$seq} = $header;
			}
			else{
				$countn++;
			}
		}
		$line =~ s/>//; 	
		$header = $line;
		$seq = "";
	}
	else{
		$seq = $seq.$line;
	}
}

if ($seq !~ m/N/ or $nons == 0){ #to catch the last sequence
	$seqs{$seq} = $header;
}
else{
	$countn++;
}

if ($nons == 1){
	print "$countn sequences with Ns out of $count sequences in total\n";
}

close INseqs;

print OUT "node1\tnode2\ttype\tdetail\tsensus\tposition\n";
print OUT2 "OTU\tN\n";

#Pretty stats
$ctype{"AT"} = "A-T";
$ctype{"TA"} = "A-T";
$ctype{"AC"} = "A-C";
$ctype{"CA"} = "A-C";
$ctype{"AG"} = "A-G";
$ctype{"GA"} = "A-G";
$ctype{"AN"} = "A-N";
$ctype{"NA"} = "A-N";
$ctype{"TC"} = "T-C";
$ctype{"CT"} = "T-C";
$ctype{"TG"} = "T-G";
$ctype{"GT"} = "T-G";
$ctype{"TN"} = "T-N";
$ctype{"NT"} = "T-N";
$ctype{"CG"} = "C-G";
$ctype{"GC"} = "C-G";
$ctype{"CN"} = "C-N";
$ctype{"NC"} = "C-N";
$ctype{"GN"} = "G-N";
$ctype{"NG"} = "G-N";

$stype{"AT"} = "Mismatch";
$stype{"TA"} = "Mismatch";
$stype{"AC"} = "Mismatch";
$stype{"CA"} = "Mismatch";
$stype{"AG"} = "Mismatch";
$stype{"GA"} = "Mismatch";
$stype{"AN"} = "Nmatch";
$stype{"NA"} = "Nmatch";
$stype{"TC"} = "Mismatch";
$stype{"CT"} = "Mismatch";
$stype{"TG"} = "Mismatch";
$stype{"GT"} = "Mismatch";
$stype{"TN"} = "Nmatch";
$stype{"NT"} = "Nmatch";
$stype{"CG"} = "Mismatch";
$stype{"GC"} = "Mismatch";
$stype{"CN"} = "Nmatch";
$stype{"NC"} = "Nmatch";
$stype{"GN"} = "Nmatch";
$stype{"NG"} = "Nmatch";

foreach my $sequence (keys %seqs){
	$remove = 0;
	if ($sequence !~ m/N/){
		print OUT2 "$seqs{$sequence}\t0\n";
	}
	else{
		print OUT2 "$seqs{$sequence}\t1\n";
	}
	my @nucl =  split(//, $sequence); 
	for (my $count = 0; $count < length($sequence); $count++) {
		if ("A" ne $nucl[$count]){                                    #Should just make a loop that jumps though the 4 nucleotides...
			my @temp = @nucl;
			$temp[$count] = "A";
			my $newotu = join("",@temp);			
			if (exists($seqs{$newotu})){
				$missing{$newotu} = 1;
				$missing{$sequence} = 1;
				my $rev = $seqs{$newotu}."m".$seqs{$sequence};
				if (!exists($printed{$rev})){
					my $pretty = $nucl[$count]."A";
					print OUT "$seqs{$sequence}\t$seqs{$newotu}\t$stype{$pretty}\t$nucl[$count]-A\t$ctype{$pretty}\t$count\n";
					my $hit = $seqs{$sequence}."m".$seqs{$newotu};
					$printed{$hit} = 1;
				}
			}
		}
		if ("T" ne $nucl[$count]){
			my @temp = @nucl;
			$temp[$count] = "T";
			my $newotu = join("",@temp);			
			if (exists($seqs{$newotu})){
				$missing{$newotu} = 1;
				$missing{$sequence} = 1;
				my $rev = $seqs{$newotu}."m".$seqs{$sequence};
				if (!exists($printed{$rev})){
					my $pretty = $nucl[$count]."T";
					print OUT "$seqs{$sequence}\t$seqs{$newotu}\t$stype{$pretty}\t$nucl[$count]-T\t$ctype{$pretty}\t$count\n";
					my $hit = $seqs{$sequence}."m".$seqs{$newotu};
					$printed{$hit} = 1;
				}
			}
		}
		if ("C" ne $nucl[$count]){
			my @temp = @nucl;
			$temp[$count] = "C";
			my $newotu = join("",@temp);			
			if (exists($seqs{$newotu})){
				$missing{$newotu} = 1;
				$missing{$sequence} = 1;
				my $rev = $seqs{$newotu}."m".$seqs{$sequence};
				if (!exists($printed{$rev})){
					my $pretty = $nucl[$count]."C";
					print OUT "$seqs{$sequence}\t$seqs{$newotu}\t$stype{$pretty}\t$nucl[$count]-C\t$ctype{$pretty}\t$count\n";
					my $hit = $seqs{$sequence}."m".$seqs{$newotu};
					$printed{$hit} = 1;
				}
			}
		}	
		if ("G" ne $nucl[$count]){
			my @temp = @nucl;
			$temp[$count] = "G";
			my $newotu = join("",@temp);			
			if (exists($seqs{$newotu})){
				$missing{$newotu} = 1;
				$missing{$sequence} = 1;
				my $rev = $seqs{$newotu}."m".$seqs{$sequence};
				if (!exists($printed{$rev})){
					my $pretty = $nucl[$count]."G";
					print OUT "$seqs{$sequence}\t$seqs{$newotu}\t$stype{$pretty}\t$nucl[$count]-G\t$ctype{$pretty}\t$count\n";
					my $hit = $seqs{$sequence}."m".$seqs{$newotu};
					$printed{$hit} = 1;
				}
			}
		}
		if ("N" ne $nucl[$count]){
			my @temp = @nucl;
			$temp[$count] = "N";
			my $newotu = join("",@temp);			
			if (exists($seqs{$newotu})){
				$missing{$newotu} = 1;
				$missing{$sequence} = 1;
				my $rev = $seqs{$newotu}."m".$seqs{$sequence};
				if (!exists($printed{$rev})){
					my $pretty = $nucl[$count]."N";					
					print OUT "$seqs{$sequence}\t$seqs{$newotu}\t$stype{$pretty}\t$nucl[$count]-N\t$ctype{$pretty}\t$count\n";
					my $hit = $seqs{$sequence}."m".$seqs{$newotu};
					$printed{$hit} = 1;
				}
			}
		}			
		#InDels
		my @tempD = @nucl;
		$tempD[$count] = "";
		my $newotuD = join("",@tempD);	
		if (exists($seqs{$newotuD})){
			$missing{$newotuD} = 1;
			$missing{$sequence} = 1;
			my $rev = $seqs{$newotuD}."m".$seqs{$sequence};
			if (!exists($printed{$rev})){
				print OUT "$seqs{$sequence}\t$seqs{$newotuD}\tIndel\t$nucl[$count]-D\t$nucl[$count]-D\t$count\n";
				my $hit = $seqs{$sequence}."m".$seqs{$newotuD};
				$printed{$hit} = 1;
			}
		}
	}	
}

foreach my $sequence (keys %seqs){
	if (!exists($missing{$sequence})){
		print OUT "$seqs{$sequence}\t$seqs{$sequence}\tnone\tnone\tnone\tnone\n"
	}
}

close OUT;
close OUT2;

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "inseqs|i:s", "networkstats|o:s", "nons|n+", "ncatch|c:s");
    my %options;

    # Add any other command line options, and the code to handle them
    # 
    GetOptions( \%options, @standard_options );
    
	#if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsosy items
    #if(!exists $options{'infile'} ) { print "**ERROR: $0 : \n"; exec("pod2usage $0"); }

    return \%options;
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

__DATA__

=head1 NAME

    visualize.error.amplicon.pl

=head1 COPYRIGHT

   copyright (C) 2013 Mads Albertsen

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

=head1 DESCRIPTION



=head1 SYNOPSIS

visualize.error.amplicon.pl  -i [-o -p -h]

 [-help -h]           Displays this basic usage information
 [-inseqs -i]         Sequence file with all sequences
 [-networkstats -o]   Network stats of 1bp simlilarity
 [-nons -n]           Discard sequences with Ns (flag, default no)
 [-ncatch -c]         Outputs overview of sequences with Ns (default: ncatch.txt).
 