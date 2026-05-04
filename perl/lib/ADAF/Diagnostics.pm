package ADAF::Diagnostics;

use strict;
use warnings;
use Exporter qw(import);
use Math::Derivative qw(Derivative1 Derivative2);

our @EXPORT_OK = qw(classify_solution);

sub classify_solution {
    my ($path) = @_;

    my ($x_ref, $y_ref, $flags_ref) = _parse_solution($path);
    my @x = @{$x_ref};
    my @y = @{$y_ref};
    my %flags = %{$flags_ref};

    my $result = {
        has_output => scalar(@x) > 0,
        is_nice    => 0,
        nooutput   => scalar(@x) > 0 ? 0 : 1,
        linesout   => scalar(@x),
        increase   => 1,
        discont    => 0,
        sonic      => 'Problem!',
        largest    => scalar(@y) ? $y[0] : 'Problem!',
        largestR   => scalar(@x) ? $x[0] : 'Problem!',
        weirdam    => $flags{weirdam},
        nan        => $flags{nan},
        failed     => $flags{failed},
    };

    return $result unless $result->{has_output};

    my @dydx = _first_derivative(\@x, \@y);
    my @d2ydx2 = _second_derivative(\@x, \@y);
    my $testsonic = 1;

    for my $i (0 .. $#x) {
        if ($dydx[$i] > 0.0 && abs($dydx[$i]) > 0.01 && $x[$i] <= 90.0) {
            $result->{increase} = 0;
        }

        if ($d2ydx2[$i] <= -0.09 && abs($d2ydx2[$i]) > 0.01 && $x[$i] > 3.0) {
            $result->{discont} = 1;
        }

        if ($y[$i] > 1.0 && $testsonic) {
            $result->{sonic} = $i > 0 ? $x[$i - 1] : $x[$i];
            $testsonic = 0;
        }

        if ($y[$i] > $result->{largest}) {
            $result->{largest} = $y[$i];
            $result->{largestR} = $x[$i];
        }
    }

    $result->{is_nice} =
        $result->{increase} == 1
        && $result->{discont} == 0
        && $result->{weirdam} == 0
        && $result->{sonic} ne 'Problem!'
        && $result->{nan} == 0
        && $result->{failed} == 0;

    return $result;
}

sub _parse_solution {
    my ($path) = @_;

    open(my $fh, '<', $path) or die "Can't open $path: $!\n";

    my (@x, @y);
    my %flags = (
        nan     => 0,
        failed  => 0,
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

        my $x = _to_number($fields[0]);
        my $y = _to_number($fields[1]);
        next unless defined $x && defined $y;

        push @x, $x;
        push @y, $y;
    }

    close($fh);

    return (\@x, \@y, \%flags);
}

sub _first_derivative {
    my ($x_ref, $y_ref) = @_;
    my $n = scalar(@{$x_ref});
    return () if $n == 0;
    return (0.0) if $n == 1;

    return Derivative1($x_ref, $y_ref);
}

sub _second_derivative {
    my ($x_ref, $y_ref) = @_;
    my $n = scalar(@{$x_ref});
    return (0.0) x $n if $n < 3;

    return Derivative2($x_ref, $y_ref);
}

sub _to_number {
    my ($value) = @_;
    return undef unless defined $value;

    $value =~ s/[dD]/e/g;
    return undef if $value !~ /^[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:e[-+]?\d+)?$/i;

    return 0.0 + $value;
}

1;
