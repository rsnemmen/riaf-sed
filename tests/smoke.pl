#!/usr/bin/perl

use strict;
use warnings;

use Cwd qw(abs_path);
use File::Spec;
use FindBin qw($Bin);
use lib File::Spec->catdir($Bin, File::Spec->updir(), 'perl', 'lib');
use ADAF::Diagnostics qw(classify_solution);

my $repo_root = abs_path(File::Spec->catdir($Bin, File::Spec->updir()));
my $fortran_dir = File::Spec->catdir($repo_root, 'fortran');

run_command('make', '-C', $fortran_dir);

my $fixtures_dir = File::Spec->catdir($repo_root, 'tests', 'testcases_for_code');
my @good_fixtures = qw(out_nice.dat out_nice03.dat out_nice04.dat);
my @bad_fixtures = qw(
  out_bad.dat
  out_bad02.dat
  out_bad03.dat
  out_bad04.dat
  out_bad05.dat
  out_discont.dat
  out_discont02.dat
  out_spiky.dat
);

for my $fixture (@good_fixtures) {
    my $path = File::Spec->catfile($fixtures_dir, $fixture);
    my $result = classify_solution($path);
    die "$fixture should pass the dyn.pl-style diagnostics.\n" unless $result->{is_nice};
}

for my $fixture (@bad_fixtures) {
    my $path = File::Spec->catfile($fixtures_dir, $fixture);
    my $result = classify_solution($path);
    die "$fixture should fail the dyn.pl-style diagnostics.\n" if $result->{is_nice};
}

my @good_examples = (
    File::Spec->catfile($repo_root, 'tests', 'n1097', 'out_std'),
    File::Spec->catfile($repo_root, 'tests', 'm81', 'out_03'),
);

for my $path (@good_examples) {
    my $result = classify_solution($path);
    die "$path should pass the dyn.pl-style diagnostics.\n" unless $result->{is_nice};
}

my @spectra = (
    File::Spec->catfile($repo_root, 'tests', 'n1097', 'spec_std'),
    File::Spec->catfile($repo_root, 'tests', 'm81', 'spec_03'),
);

for my $path (@spectra) {
    validate_spectrum($path);
}

print "Smoke checks passed.\n";

sub run_command {
    my (@command) = @_;
    system(@command) == 0 or die "Command failed: @command\n";
}

sub validate_spectrum {
    my ($path) = @_;

    open(my $fh, '<', $path) or die "Can't open $path: $!\n";

    my $count = 0;
    my $previous_x;
    while (my $line = <$fh>) {
        next if $line =~ /^\s*#/;
        next if $line =~ /^\s*$/;

        my @fields = split ' ', $line;
        next unless @fields >= 2;

        my $x = to_number($fields[0]);
        my $y = to_number($fields[1]);
        next unless defined $x && defined $y;

        die "$path has a non-increasing frequency grid.\n" if defined $previous_x && $x <= $previous_x;
        die "$path has a non-positive luminosity value.\n" if $y <= 0.0;

        $previous_x = $x;
        $count++;
    }

    close($fh);

    die "$path does not contain enough spectrum samples.\n" if $count < 20;
}

sub to_number {
    my ($value) = @_;
    return undef unless defined $value;

    $value =~ s/[dD]/e/g;
    return undef if $value !~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[-+]?\d+)?$/i;

    return 0.0 + $value;
}
