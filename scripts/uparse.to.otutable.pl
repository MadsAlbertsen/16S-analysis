#!/usr/bin/env perl
###############################################################################
#
#    uparse.to.otutable.pl
#
#	 Takes a .uc file and a sampleid file and makes an OTU table..
#    
#    Copyright (C) 2012 Mads Albertsen
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

my $in_uc;
my $in_sampleid;
my $out;

$in_uc = &overrideDefault("in.uc",'in_uc');
$in_sampleid = &overrideDefault("sampleid.txt",'in_sampleid');
$out = &overrideDefault("otutable.txt",'out');
 
my %count;
my %sampleid;

######################################################################
# CODE HERE
######################################################################


### Read in the sample id
open(in_sampleid_fh, $in_sampleid) or die("Cannot read file: $in_sampleid\n");

while ( my $line = <in_sampleid_fh> ) {
	chomp $line;   	
	$line =~ s/>//;
	#$line =~ s/\r\n//g;
	my @id = split(/\t/,$line);
	my $sample = $id[1];
	my @id_split = split(/:/,$id[0]);
	my $read = $id_split[0] . $id_split[1] . $id_split[2] . $id_split[9];
	$sampleid{$read} = $sample;
}

close in_sampleid_fh;

### Go through the read mapping file and extract OTU counts per sample
open(in_uc_fh, $in_uc) or die("Cannot read file: $in_uc\n");

while ( my $line = <in_uc_fh> ) {
	chomp $line;   	
	#$line =~ s/\r\n//g;
	my @id = split(/\t/,$line);
	next if ($id[0] eq "N");
	my @id_split = split(/:/,$id[8]);
	my $sample = $id_split[0] . $id_split[1] . $id_split[2] . $id_split[9];
	my $OTU = $id[9];
		
	if(exists($count{$OTU}{$sampleid{$sample}})) {
		$count{$OTU}{$sampleid{$sample}}++
	} else {
		$count{$OTU}{$sampleid{$sample}} = 1;
	}
}
close in_uc_fh;

#=cut
### Save the data to a file
open(out_fh, ">$out") or die("Cannot create file: $out\n");

my $header = "OTU";
foreach my $sample (keys %sampleid){
		$header = "$header\t$sampleid{$sample}";
	}
print out_fh "$header\n";

foreach my $OTU (keys %count){
	my $OTU_count = $OTU;
	foreach my $sample (keys %sampleid){
		if (exists($count{$OTU}{$sampleid{$sample}})){
			$OTU_count = "$OTU_count\t$count{$OTU}{$sampleid{$sample}}";
		} else{
			$OTU_count = "$OTU_count\t0";
		}
	}
	print out_fh "$OTU_count\n";
}

close out_fh;
#=cut

######################################################################
# TEMPLATE SUBS
######################################################################
sub checkParams {
    #-----
    # Do any and all options checking here...
    #
    my @standard_options = ( "help|h+", "in_uc|u:s", "in_sampleid|s:s", "out|o:s");
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

    vprobes.generateprobes.pl

=head1 COPYRIGHT

   copyright (C) 2012 Mads Albertsen

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

uparse.to.otutable.pl -u -s [-h -o]

 [-help -h]           Displays this basic usage information
 [-in_uc -u]          Input .uc file.
 [-in_sampleid -s]    Input sampleid file.
 [-out -o]            Output OTU table.
 
=cut
