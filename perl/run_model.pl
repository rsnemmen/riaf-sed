#!/usr/bin/env perl

use strict;
use warnings;

use Cwd qw(abs_path getcwd);
use File::Basename qw(basename);
use File::Spec;
use FindBin qw($Bin);
use Term::ANSIColor qw(colored);

use lib "$Bin/lib";
use ADAF::Diagnostics qw(classify_solution);

usage() if @ARGV != 1;

my $input = abs_path($ARGV[0]);
die colored("error", "red") . ": can't find parameter file $ARGV[0]\n"
    unless defined $input && -f $input;

my $workdir = getcwd();
my $repo_root = abs_path(File::Spec->catdir($Bin, File::Spec->updir()));
die colored("error", "red") . ": can't determine repository root from $Bin\n"
    unless defined $repo_root;

my $dyn_script = File::Spec->catfile($Bin, 'dyn.pl');
my $spectrum_script = File::Spec->catfile($Bin, 'spectrum.pl');
my $plot_script = File::Spec->catfile($Bin, 'plot_sed.py');
my $diag = read_parameter($input, 'diag');
my $spectrum = read_parameter($input, 'spec');
my $png = output_png_for($input);

die colored("error", "red") . ": missing diag= in $input\n" unless defined $diag;
die colored("error", "red") . ": missing spec= in $input\n" unless defined $spectrum;

print_stage("ADAF model run");
print_info("Parameter file", $input);
print_info("Working directory", $workdir);

print_stage("Computing dynamics");
run_command($^X, $dyn_script, $input);

print_stage("Checking dynamics solution");
my $diag_path = File::Spec->rel2abs($diag, $workdir);
my $result = classify_solution($diag_path);
if (!$result->{is_nice}) {
    print colored("No physical global solution found.", "bold red") . "\n";
    print_info("Dynamics output", $diag_path);
    print_info("Shells", $result->{linesout});
    print_info("Sonic point", $result->{sonic});
    print_info("Mach max", $result->{largest});
    exit 1;
}
print colored("Physical solution found.", "bold green") . "\n";
print_info("Dynamics output", $diag_path);
print_info("Shells", $result->{linesout});
print_info("Sonic point", $result->{sonic});

print_stage("Computing spectrum");
run_command($^X, $spectrum_script, $input);

print_stage("Plotting SED");
my $spectrum_path = File::Spec->rel2abs($spectrum, $workdir);
my $png_path = File::Spec->rel2abs($png, $workdir);
run_command('python3', $plot_script, $spectrum_path, $png_path);

print_stage("Done");
print_info("Spectrum", $spectrum_path);
print_info("SED plot", $png_path);

sub usage {
    die "Usage: $0 parameter-file\n";
}

sub read_parameter {
    my ($path, $key) = @_;

    open(my $fh, '<', $path) or die "Can't open $path: $!\n";
    while (my $line = <$fh>) {
        $line =~ s/#.*$//;
        next if $line =~ /^\s*$/;
        my ($name, $value) = split /=/, $line, 2;
        next unless defined $name && defined $value;
        trim($name);
        trim($value);
        if ($name eq $key) {
            close($fh);
            return $value;
        }
    }
    close($fh);

    return undef;
}

sub trim {
    $_[0] =~ s/^\s+//;
    $_[0] =~ s/\s+$//;
}

sub output_png_for {
    my ($path) = @_;
    my $name = basename($path);
    $name =~ s/\.[^.]*$//;
    return "$name.png";
}

sub print_stage {
    my ($message) = @_;
    print "\n" . colored("==>", "bold cyan") . " " . colored($message, "bold") . "\n";
}

sub print_info {
    my ($label, $value) = @_;
    print colored(sprintf("%-20s", "$label:"), "cyan") . "$value\n";
}

sub run_command {
    my (@command) = @_;
    print colored('$ ', "yellow") . join(' ', @command) . "\n";
    system(@command);
    my $status = $?;

    if ($status == -1) {
        die colored("error", "red") . ": failed to run @command: $!\n";
    }
    if ($status & 127) {
        die colored("error", "red") . ": @command died with signal " . ($status & 127) . "\n";
    }
    if ($status != 0) {
        die colored("error", "red") . ": @command exited with status " . ($status >> 8) . "\n";
    }
}
