#!/home/ivan/bin/perl

use warnings;
use strict;
use blib;
use lib '../../Mol/blib/lib';
use Chemistry::File::Mopac;
use Chemistry::Mol;

my $mol = Chemistry::Mol->read("test.mop", type => "mop");
#my $mol = Chemistry::Mol->read("test.mop");

print $mol->print;
