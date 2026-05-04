#!/usr/bin/perl

use strict;
use warnings;

use Cwd qw(abs_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use lib File::Spec->catdir($Bin, File::Spec->updir(), 'perl', 'lib');
use ADAF::Diagnostics qw(classify_solution);
use ADAF::Paths qw(parameter_file_from_args);

my $repo_root = abs_path(File::Spec->catdir($Bin, File::Spec->updir()));
my $fortran_dir = File::Spec->catdir($repo_root, 'fortran');

test_parameter_file_selector();

run_command('make', '-C', $fortran_dir);

test_largeR_dynamics_regression();
test_largeR_spectrum_regression();

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

sub test_largeR_dynamics_regression {
    my $workdir = tempdir('adaf-dyn-XXXX', TMPDIR => 1, CLEANUP => 1);
    my $dyn_script = File::Spec->catfile($repo_root, 'perl', 'dyn.pl');
    my $input = File::Spec->catfile($repo_root, 'examples', 'largeR.dat');
    my $candidate = File::Spec->catfile($workdir, 'out');
    my $reference = File::Spec->catfile($repo_root, 'tests', 'reference', 'largeR_dyn.out');
    my $comparator = File::Spec->catfile($repo_root, 'tests', 'compare_dynamics.py');
    my $plotter = File::Spec->catfile($repo_root, 'tests', 'plot_dynamics.py');
    my $plot = File::Spec->catfile($repo_root, 'tests', 'artifacts', 'largeR_dynamics.png');

    run_command_in_dir($workdir, 'perl', $dyn_script, $input);

    my $result = classify_solution($candidate);
    die "$candidate should pass the dyn.pl-style diagnostics.\n" unless $result->{is_nice};
    die "$candidate should have 91 profile shells.\n" unless $result->{linesout} == 91;

    run_command('python3', $comparator, $reference, $candidate);
    run_command('python3', $plotter, $reference, $candidate, $plot);
}

sub test_largeR_spectrum_regression {
    my $workdir = tempdir('adaf-spectrum-XXXX', TMPDIR => 1, CLEANUP => 1);
    my $dyn_script = File::Spec->catfile($repo_root, 'perl', 'dyn.pl');
    my $spectrum_script = File::Spec->catfile($repo_root, 'perl', 'spectrum.pl');
    my $input = File::Spec->catfile($repo_root, 'examples', 'largeR.dat');
    my $dynamics_candidate = File::Spec->catfile($workdir, 'out');
    my $candidate = File::Spec->catfile($workdir, 'spectrum');
    my $reference = File::Spec->catfile($repo_root, 'tests', 'reference', 'largeR_spectrum.out');
    my $comparator = File::Spec->catfile($repo_root, 'tests', 'compare_spectrum.py');
    my $plotter = File::Spec->catfile($repo_root, 'tests', 'plot_spectrum.py');
    my $plot = File::Spec->catfile($repo_root, 'tests', 'artifacts', 'largeR_spectrum.png');

    run_command_in_dir($workdir, 'perl', $dyn_script, $input);

    my $result = classify_solution($dynamics_candidate);
    die "$dynamics_candidate should pass the dyn.pl-style diagnostics.\n"
        unless $result->{is_nice};
    die "$dynamics_candidate should have 91 profile shells.\n"
        unless $result->{linesout} == 91;

    run_command_in_dir($workdir, 'perl', $spectrum_script, $input);
    run_command('python3', $comparator, $reference, $candidate);
    run_command('python3', $plotter, $reference, $candidate, $plot);
}

sub run_command_in_dir {
    my ($workdir, @command) = @_;
    my $oldcwd = Cwd::getcwd();

    chdir $workdir or die "Can't chdir to $workdir: $!\n";
    my $ok = system(@command) == 0;
    chdir $oldcwd or die "Can't chdir back to $oldcwd: $!\n";

    die "Command failed: @command\n" unless $ok;
}

sub test_parameter_file_selector {
    my $default = parameter_file_from_args('dyn.pl');
    die "No-argument parameter selection should use in.dat.\n" unless $default eq 'in.dat';

    my $custom = parameter_file_from_args('dyn.pl', 'model.dat');
    die "One-argument parameter selection should use that file.\n" unless $custom eq 'model.dat';

    my $error;
    eval { parameter_file_from_args('dyn.pl', 'model.dat', 'extra.dat'); };
    $error = $@;
    die "Multiple parameter arguments should fail with usage.\n"
        unless $error =~ /^Usage: dyn\.pl \[parameter-file\]/;
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
