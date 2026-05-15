#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw($Bin);
use lib "$Bin/lib";
use ADAF::Parameters qw(convert_legacy_file write_toml_file);

usage() unless @ARGV == 2;

my ($legacy_path, $toml_path) = @ARGV;
my $data = convert_legacy_file($legacy_path);
write_toml_file($toml_path, $data);

print "Wrote $toml_path\n";

sub usage {
    die "Usage: $0 legacy.dat model.toml\n";
}
