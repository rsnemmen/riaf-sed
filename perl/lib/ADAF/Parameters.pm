package ADAF::Parameters;

use strict;
use warnings;

use Exporter qw(import);

our @EXPORT_OK = qw(read_parameters read_parameter convert_legacy_file write_toml_file);

my %SCHEMA = (
    dynamics => {
        gamma    => 'gamai',
        mass     => 'm',
        beta     => 'beta',
        alpha    => 'alfa',
        delta    => 'delta',
        mdot_out => 'dotm0',
        r_out    => 'rout',
        p_wind   => 'pp0',
    },
    shooting => {
        sl0_initial => 'sl0i',
        sl0_final   => 'sl0f',
        n_models    => 'nmodels',
    },
    boundary => {
        ti   => 'ti',
        te   => 'te',
        mach => 'vcs',
    },
    runtime => {
        diag          => 'diag',
        search_method => 'search_method',
        eig_tol       => 'eig_tol',
        dyn_timeout   => 'dyn_timeout',
        max_workers   => 'max_workers',
    },
    spectrum => {
        distance_pc => 'distance',
        theta_deg   => 'theta',
        filename    => 'spec',
    },
);

my %REQUIRED = (
    dynamics => [qw(gamma mass beta alpha delta mdot_out r_out p_wind)],
    shooting => [qw(sl0_initial sl0_final n_models)],
    boundary => [qw(ti te mach)],
    runtime  => [qw(diag)],
    spectrum => [qw(distance_pc theta_deg filename)],
);

my %LEGACY = (
    gamai  => [qw(dynamics gamma)],
    m      => [qw(dynamics mass)],
    beta   => [qw(dynamics beta)],
    alfa   => [qw(dynamics alpha)],
    delta  => [qw(dynamics delta)],
    dotm0  => [qw(dynamics mdot_out)],
    rout   => [qw(dynamics r_out)],
    pp0    => [qw(dynamics p_wind)],
    sl0i   => [qw(shooting sl0_initial)],
    sl0f   => [qw(shooting sl0_final)],
    nmodels => [qw(shooting n_models)],
    ti     => [qw(boundary ti)],
    te     => [qw(boundary te)],
    vcs    => [qw(boundary mach)],
    diag   => [qw(runtime diag)],
    search_method => [qw(runtime search_method)],
    eig_tol       => [qw(runtime eig_tol)],
    dyn_timeout   => [qw(runtime dyn_timeout)],
    max_workers   => [qw(runtime max_workers)],
    distance => [qw(spectrum distance_pc)],
    theta    => [qw(spectrum theta_deg)],
    spec     => [qw(spectrum filename)],
);

my %STRING_KEYS = map { $_ => 1 } qw(search_method diag filename);

sub read_parameters {
    my ($path) = @_;
    my $data = _read_toml($path);
    _validate_toml($path, $data);
    return _flatten($data);
}

sub read_parameter {
    my ($path, $key) = @_;
    my $params = read_parameters($path);
    return $params->{$key};
}

sub convert_legacy_file {
    my ($path) = @_;

    open(my $fh, '<', $path) or die "Can't open $path: $!\n";
    my %toml;
    my $line_number = 0;
    while (my $line = <$fh>) {
        $line_number++;
        $line =~ s/#.*$//;
        _trim($line);
        next if $line eq '';

        my ($legacy_key, $value) = split /=/, $line, 2;
        die "$path:$line_number: expected legacy key=value syntax.\n"
            unless defined $legacy_key && defined $value;
        _trim($legacy_key);
        _trim($value);

        my $mapping = $LEGACY{$legacy_key};
        die "$path:$line_number: unknown legacy parameter '$legacy_key'.\n"
            unless defined $mapping;
        my ($table, $key) = @$mapping;
        $toml{$table}{$key} = _legacy_value($key, $value);
    }
    close($fh);

    _validate_toml($path, \%toml);
    return \%toml;
}

sub write_toml_file {
    my ($path, $data) = @_;

    open(my $fh, '>', $path) or die "Can't write $path: $!\n";
    print {$fh} _format_toml($data);
    close($fh);
}

sub _read_toml {
    my ($path) = @_;

    open(my $fh, '<', $path) or die "Can't open $path: $!\n";
    local $/;
    my $content = <$fh>;
    close($fh);

    my $data = _read_with_toml_tiny($content);
    return $data if defined $data;
    return _parse_basic_toml($path, $content);
}

sub _read_with_toml_tiny {
    my ($content) = @_;

    return undef unless eval { require TOML::Tiny; 1; };
    return undef unless TOML::Tiny->can('from_toml');

    my $data = TOML::Tiny::from_toml($content);
    return $data;
}

sub _parse_basic_toml {
    my ($path, $content) = @_;

    my %data;
    my $table;
    my $line_number = 0;
    for my $line (split /\n/, $content) {
        $line_number++;
        $line =~ s/#.*$//;
        _trim($line);
        next if $line eq '';

        if ($line =~ /^\[([A-Za-z0-9_]+)\]$/) {
            $table = $1;
            $data{$table} ||= {};
            next;
        }

        die "$path:$line_number: parameter must be inside a TOML table.\n"
            unless defined $table;
        my ($key, $value) = split /=/, $line, 2;
        die "$path:$line_number: expected TOML key = value syntax.\n"
            unless defined $key && defined $value;
        _trim($key);
        _trim($value);
        $data{$table}{$key} = _parse_basic_value($path, $line_number, $value);
    }

    return \%data;
}

sub _parse_basic_value {
    my ($path, $line_number, $value) = @_;

    if ($value =~ /^"(.*)"$/) {
        my $string = $1;
        $string =~ s/\\"/"/g;
        $string =~ s/\\\\/\\/g;
        return $string;
    }

    die "$path:$line_number: unsupported TOML value '$value'.\n"
        unless $value =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
    return $value;
}

sub _validate_toml {
    my ($path, $data) = @_;

    for my $table (keys %$data) {
        die "$path: unknown TOML table [$table].\n" unless exists $SCHEMA{$table};
        for my $key (keys %{ $data->{$table} }) {
            die "$path: unknown parameter [$table].$key.\n"
                unless exists $SCHEMA{$table}{$key};
            _validate_value($path, $table, $key, $data->{$table}{$key});
        }
    }

    for my $table (keys %REQUIRED) {
        die "$path: missing TOML table [$table].\n" unless exists $data->{$table};
        for my $key (@{ $REQUIRED{$table} }) {
            die "$path: missing required parameter [$table].$key.\n"
                unless exists $data->{$table}{$key};
        }
    }
}

sub _validate_value {
    my ($path, $table, $key, $value) = @_;

    die "$path: missing value for [$table].$key.\n"
        unless defined $value && "$value" ne '';

    return if $STRING_KEYS{$key};

    die "$path: [$table].$key must be numeric, got '$value'.\n"
        unless "$value" =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
}

sub _flatten {
    my ($data) = @_;
    my %params;

    for my $table (keys %SCHEMA) {
        for my $key (keys %{ $SCHEMA{$table} }) {
            next unless exists $data->{$table}{$key};
            $params{ $SCHEMA{$table}{$key} } = $data->{$table}{$key};
        }
    }

    return \%params;
}

sub _legacy_value {
    my ($key, $value) = @_;

    return $value if $STRING_KEYS{$key};
    $value =~ s/[dD]/e/g;
    $value =~ s/\.([eE])/.0$1/;
    $value =~ s/\.\z/.0/;
    return $value;
}

sub _format_toml {
    my ($data) = @_;
    my @tables = qw(dynamics shooting boundary runtime spectrum);
    my @parts;

    for my $table (@tables) {
        next unless exists $data->{$table};
        push @parts, "[$table]\n";
        for my $key (_ordered_keys($table)) {
            next unless exists $data->{$table}{$key};
            push @parts, "$key = " . _format_value($key, $data->{$table}{$key}) . "\n";
        }
        push @parts, "\n";
    }

    return join('', @parts);
}

sub _ordered_keys {
    my ($table) = @_;

    return @{ $REQUIRED{$table} }, grep {
        my $key = $_;
        !grep { $_ eq $key } @{ $REQUIRED{$table} }
    } sort keys %{ $SCHEMA{$table} };
}

sub _format_value {
    my ($key, $value) = @_;

    if ($STRING_KEYS{$key}) {
        $value =~ s/\\/\\\\/g;
        $value =~ s/"/\\"/g;
        return qq("$value");
    }

    return $value;
}

sub _trim {
    $_[0] =~ s/^\s+//;
    $_[0] =~ s/\s+$//;
}

1;
