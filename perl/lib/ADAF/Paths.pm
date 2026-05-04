package ADAF::Paths;

use strict;
use warnings;

use Cwd qw(abs_path);
use Exporter qw(import);
use File::Copy qw(copy);
use File::Spec;

our @EXPORT_OK = qw(fortran_binary stage_data_files);

sub fortran_binary {
    my ($script_dir, $binary_name) = @_;

    my $repo_root = abs_path(File::Spec->catdir($script_dir, File::Spec->updir()));
    die "Can't determine repository root from $script_dir.\n" unless defined $repo_root;

    my $binary = File::Spec->catfile($repo_root, 'bin', $binary_name);
    die "Can't find $binary. Run `make build` first.\n" unless -e $binary;
    die "Can't execute $binary. Run `make build` first.\n" unless -x $binary;

    return $binary;
}

sub stage_data_files {
    my ($script_dir, $workdir, @files) = @_;

    my $repo_root = abs_path(File::Spec->catdir($script_dir, File::Spec->updir()));
    die "Can't determine repository root from $script_dir.\n" unless defined $repo_root;

    my $data_dir = File::Spec->catdir($repo_root, 'data');

    for my $f (@files) {
        my $src = File::Spec->catfile($data_dir, $f);
        die "Missing data file $src — is the data/ directory present?\n" unless -e $src;
        my $dst = File::Spec->catfile($workdir, $f);
        next if -e $dst;
        symlink($src, $dst) || copy($src, $dst)
            || die "Can't stage $src in $workdir: $!\n";
    }
}

1;
