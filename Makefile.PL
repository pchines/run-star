use ExtUtils::MakeMaker;
# See lib/ExtUtils/MakeMaker.pm for details of how to influence
# the contents of the Makefile that is written.
WriteMakefile(
    'NAME'	    => 'Run::STAR',
    'VERSION'   => '0.01',
    'EXE_FILES' => exe_files("scripts"),
    'AUTHOR'    => 'Peter Chines <pchines@mail.nih.gov>',
    'LICENSE'   => 'U.S. Government Work; no restrictions on use',
    'ABSTRACT'  => 'Scripts to efficiently submit STAR aligner jobs to a cluster',
);

sub exe_files {
    my ($dir) = @_;
    opendir EXEDIR, $dir or die "Can't open '$dir', $!\n";
    my @files = grep { -x $_ && !-d $_ }
        map { "$dir/$_" } readdir EXEDIR;
    return \@files;
}
