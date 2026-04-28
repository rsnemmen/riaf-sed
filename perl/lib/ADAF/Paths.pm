package ADAF::Paths;

use strict;
use warnings;

use Cwd qw(abs_path);
use Exporter qw(import);
use File::Spec;

our @EXPORT_OK = qw(fortran_binary);

sub fortran_binary {
    my ($script_dir, $binary_name) = @_;

    my $repo_root = abs_path(File::Spec->catdir($script_dir, File::Spec->updir()));
    die "Can't determine repository root from $script_dir.\n" unless defined $repo_root;

    my $binary = File::Spec->catfile($repo_root, 'fortran', $binary_name);
    die "Can't find $binary. Run `cd fortran && make` first.\n" unless -e $binary;
    die "Can't execute $binary. Run `cd fortran && make` first.\n" unless -x $binary;

    return $binary;
}

1;
