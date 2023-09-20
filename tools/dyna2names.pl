#!/usr/bin/perl

## Copyright (2000) Sandia Corporation. Under the terms of Contract
## DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains 
## certain rights in this software.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in the
##     documentation and/or other materials provided with the
##     distribution.
##
##   * Neither the name of Sandia nor the names of any contributors may
##     be used to endorse or promote products derived from this software
##     without specific prior written permission.
##
##   * Modified source versions must be plainly marked as such, and must
##     not be misrepresented as being the original software.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
## ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR 
## ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
## LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
## OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
## DAMAGE.

## $Source$
## $Revision$
## $Date$
##
## Create a .names file from an LsDyna input file. It will be named 
## d3plot.names.
##
## Modifications
##   - 5/15/2006 WSD Created.
##
## Fix? should we check for *KEYWORD to make sure it is proper type of file?
##
## Fix? could use *INCLUDE to decide other .k files to troll through?
##
## Fix? support multiple files on command line?
##
## Fix? does LSDyna have element sets?

use strict;

my $output_file = "d3plot.names";

# Test input
die "Usage: dyna2names <lsdyna_input_file>\n" if (@ARGV < 1);
die "ERROR: Currently you can only specify 1 input file. If this is a problem\n"
    . "please contact the FCLib developers.\n" if (@ARGV > 1);
my $input_deck = $ARGV[0];
die "Could not find input file $input_deck\n" if (! -e $input_deck); 

# Scan input, collect part names & ids, collect parameter defs
# assuming that any *<KEYWORD> will have no spaces before it on a line
# assuming all *PARAMETER definitions will be encountered before they are used
my $num;
my @names;
my @ids;
my %parameters;
open INPUT_DECK, $input_deck or die "Could not open input file $input_deck\n";
$num = 0;
while (<INPUT_DECK>) {
    chomp;
    if (/^\*PART($|\s|_)/i) {
	# The first line after PART keyword is the name of the part
	while (<INPUT_DECK>) {
	    chomp;
	    if (m/^\$/) {
		; # do nothing - skipping comments
	    }
	    else {
		m/^\s*(.*\S)\s*$/;
		$names[$num] = $1;
		last;
	    }
	}
	# The next line, the first int, is part identifier
	# If it starts with '&' we have to look up id, otherwise is id
	while (<INPUT_DECK>) {
	    chomp;
	    if (m/^\$/) {
		; # do nothing - skipping comments
	    }
	    else {
		m/^\s*(&?)(\w+)/;
		if ($1) {
		    $ids[$num] = $parameters{$2};
		}
		else {
		    $ids[$num] = $2;
		}
		last;
	    }
	}
	$num++;
    }
    elsif (/^\*PARAMETER($|\s)/i) {
	# Each line is a parameter spec
	while (<INPUT_DECK>) {
	    chomp;
	    if (m/^\$/) {
		; # do nothing - skipping commnets
	    }
	    else { # assuming 1 parameter per line
		m/(I|R|i|r)\s*(\w+)\s+([\w\.\-]+)/;
		#print "$_\n";
		#print "** $2 = $3\n";
		$parameters{$2} = $3;
		last;
	    }
	}
    }
}
close INPUT_DECK;

#debug - print the parameters
#print "%parameters =\n";
#my $key;
#my $value;
#while ( ($key, $value) = each %parameters ) {
#    print "$key => $value\n";
#}

# print the names file
open OUTPUT_FILE, ">$output_file" or die "Could not open output file $output_file\n";
my $i;
print OUTPUT_FILE "blocks\n";
for ($i = 0; $i < $num; $i++) {
    #printf "$ids[$i]      $names[$i]\n";
    printf OUTPUT_FILE "%-5d %s\n", $ids[$i], $names[$i];
}
close OUTPUT_DECK;
