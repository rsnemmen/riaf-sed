#!/usr/bin/perl

use strict;
use warnings;

use Cwd qw(abs_path);
use File::Spec;
use FindBin qw($Bin);

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

sub classify_solution {
    my ($path) = @_;

    my ($x_ref, $y_ref, $flags_ref) = parse_solution($path);
    my @x = @{$x_ref};
    my @y = @{$y_ref};
    my %flags = %{$flags_ref};

    my $result = {
        has_output => scalar(@x) > 0,
        is_nice => 0,
    };

    return $result unless $result->{has_output};

    my @dydx = first_derivative(\@x, \@y);
    my @d2ydx2 = second_derivative(\@x, \@y);

    my $increase = 1;
    my $discont = 0;
    my $sonic = 'Problem!';
    my $testsonic = 1;

    for my $i (0 .. $#x) {
        if ($dydx[$i] > 0.0 && abs($dydx[$i]) > 0.01 && $x[$i] <= 90.0) {
            $increase = 0;
        }

        if ($d2ydx2[$i] <= -0.09 && abs($d2ydx2[$i]) > 0.01 && $x[$i] > 3.0) {
            $discont = 1;
        }

        if ($y[$i] > 1.0 && $testsonic) {
            $sonic = $i > 0 ? $x[$i - 1] : $x[$i];
            $testsonic = 0;
        }
    }

    $result->{is_nice} =
        $increase == 1
        && $discont == 0
        && $flags{weirdam} == 0
        && $sonic ne 'Problem!'
        && $flags{nan} == 0
        && $flags{failed} == 0;

    return $result;
}

sub parse_solution {
    my ($path) = @_;

    open(my $fh, '<', $path) or die "Can't open $path: $!\n";

    my (@x, @y);
    my %flags = (
        nan => 0,
        failed => 0,
        weirdam => 0,
    );

    while (my $line = <$fh>) {
        $flags{nan} = 1 if $line =~ /NaN/i;
        $flags{failed} = 1 if $line =~ /FAILED/i;
        $flags{weirdam} = 1 if $line =~ /ssll < 0/;

        next if $line =~ /^\s*#/;
        next if $line =~ /^\s*$/;

        my @fields = split ' ', $line;
        next unless @fields >= 2;

        my $x = to_number($fields[0]);
        my $y = to_number($fields[1]);
        next unless defined $x && defined $y;

        push @x, $x;
        push @y, $y;
    }

    close($fh);

    return (\@x, \@y, \%flags);
}

sub first_derivative {
    my ($x_ref, $y_ref) = @_;
    my @x = @{$x_ref};
    my @y = @{$y_ref};
    my $n = scalar(@x);
    return () if $n == 0;
    return (0.0) if $n == 1;

    my @dydx;
    for my $i (0 .. $n - 1) {
        my ($left, $right) = $i == 0 ? (0, 1) : $i == $n - 1 ? ($n - 2, $n - 1) : ($i - 1, $i + 1);
        my $dx = $x[$right] - $x[$left];
        push @dydx, $dx == 0.0 ? 0.0 : ($y[$right] - $y[$left]) / $dx;
    }

    return @dydx;
}

sub second_derivative {
    my ($x_ref, $y_ref) = @_;
    my @x = @{$x_ref};
    my @y = @{$y_ref};
    my $n = scalar(@x);
    return (0.0) x $n if $n < 3;

    my @d2ydx2 = (0.0) x $n;
    for my $i (1 .. $n - 2) {
        my $h1 = $x[$i] - $x[$i - 1];
        my $h2 = $x[$i + 1] - $x[$i];
        next if $h1 == 0.0 || $h2 == 0.0;

        my $s1 = ($y[$i] - $y[$i - 1]) / $h1;
        my $s2 = ($y[$i + 1] - $y[$i]) / $h2;
        $d2ydx2[$i] = 2.0 * ($s2 - $s1) / ($h1 + $h2);
    }

    $d2ydx2[0] = $d2ydx2[1];
    $d2ydx2[$n - 1] = $d2ydx2[$n - 2];

    return @d2ydx2;
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
