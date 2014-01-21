#!/usr/bin/perl -w
# $Id: check_star_run.pl 6356 2013-04-10 19:53:07Z pchines $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::NFSLock qw(uncache);
use GTB::File qw(Open wait_for_file);
our %Opt;
our $Lock;  # NFS-compatible lock on output file is retained until exit

=head1 NAME

check_star_run.pl - checks that STAR ran successfully

=head1 SYNOPSIS

Check that STAR ran and wrote a Log.final.out in the specified directory:

  check_star_run.pl -dir RUN_DIR

Also check that the number of reads processed matches input sequence count:

  check_star_run.pl -dir RUN_DIR -counts counts.txt

Also append STAR run statistics to report.txt:

  check_star_run.pl -dir RUN_DIR -output report.txt

For complete documentation, run C<check_star_run.pl -man>

=head1 DESCRIPTION

This program checks for the Log.final.out that STAR writes upon a successful
conclusion.  If it doesn't find one, or it doesn't match the other criteria,
this program will die with a message, e.g. "No Log.final.out in ..."

If a C<-counts> file is provided, the "Number of input reads" in the
Log.final.out must match.  Of not, this program dies with a message, e.g.
"Wrong number of input reads in ..."

If an C<-output> file is provided, a tab-delimited summary of the
Log.final.out file will be written to that file.  If the file exists, output
will be appended to that file, and only column names that exist will be
written.  Missing data will be represented by "NA".

If the C<-output> file does not exist, all columns present in the
Log.final.out file (plus any that may be specified on the command line with
the C<-field> option) will be included, and the first line of the
file will be a header line with the field names.  Note that many field names
are shortened from what is in Log.final.out, and made R-compatible.

Note that an NFS-compatible lock is obtained on the output file before
writing/appending, so it is safe to have multiple jobs writing to the same
output file.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();
my $expected_in = read_counts_file();
if (!-f "$Opt{dir}/Log.final.out") {
    die "No Log.final.out in $Opt{dir}\n";
}
my ($ofh, $ra_cols) = open_output_file();
if ($expected_in || $ofh) {
    my ($rh, $ra_order) = read_final_log("$Opt{dir}/Log.final.out");
    if ($expected_in && $rh->{input_reads} != $expected_in) {
        die "Wrong number of input reads in $Opt{dir}: "
            . "got $rh->{input_reads}, expected $expected_in\n";
    }
    # TODO: add other (optional) tests for % reads mapped, etc.?
    if ($ofh) {
        if (!@$ra_cols) {
            $ra_cols = ['dir', (sort keys %{ $Opt{field} }), @$ra_order];
            print $ofh join("\t", @$ra_cols), "\n";
        }
        $rh->{dir} = $Opt{dir};
        for my $f (keys %{ $Opt{field} }) {
            $rh->{$f} ||= $Opt{field}{$f};
        }
        print $ofh join("\t", map { defined $rh->{$_} ? $rh->{$_} : "NA" }
                @$ra_cols), "\n";
    }
}
print "OK - Log.final.out is present\n";

#------------
# End MAIN
#------------

sub read_final_log {
    my ($file) = @_;
    my $ifh = Open($file);
    my %x;
    my @cols;
    while (<$ifh>) {
        chomp;
        my ($k, $v) = split /\s*\|\s*/;
        if (!defined $v) {
            next;
        }
        $k = lc $k;
        $k =~ s/^\s+//;
        $k =~ s/\W+/_/g;
        $k =~ s/^_//;
        $k =~ s/_$//;
        $k =~ s/^(number_of_|of_reads_|mapping_speed_)//;
        $v =~ s/%$//;
        $x{$k} = $v;
        push @cols, $k;
    }
    return (\%x, \@cols);
}

sub open_output_file {
    my $ofh;
    my @cols;
    if ($Opt{output}) {
        if ($Opt{output} ne '-') {
            $Lock = File::NFSLock->new($Opt{output},'EX',5*60,30*60);
            if (!$Lock) {
                warn "Could not obtain lock on $Opt{output}, "
                        . "[$File::NFSLock::errstr]";
                return ($ofh, \@cols);
            }
            if (-f $Opt{output}) {
                my $ifh = Open($Opt{output});
                my $head = <$ifh>;
                chomp $head;
                @cols = split /\t/, $head;
                close $ifh;
            }
        }
        $ofh = Open($Opt{output}, 'a');
    }
    return ($ofh, \@cols);
}

sub read_counts_file {
    my $in;
    if ($Opt{counts}) {
        if ($Opt{counts} =~ /^\d+$/) {
            return $Opt{counts};
        }
        uncache($Opt{counts});  # NFS-lock specific to ensure update of file
        wait_for_file($Opt{counts});
        my $cfh = Open($Opt{counts});
        my $line = <$cfh>;
        if ($line =~ /^(\d+)(?:\s.+,\s*(\d+)\s+excluded)?/) {
            $in = $1;
            if ($2) {
                $Opt{field}{filtered_out} ||= $2;
            }
        }
    }
    return $in;
}

sub process_commandline {
    # Set defaults here
    %Opt = (field   => {},
            );
    GetOptions(\%Opt, qw(counts=s dir=s field=s output=s
                manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "check_star_run.pl, ", q$Revision: 6356 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    pod2usage("Not expecting non-option arguments") if @ARGV;
    if (!$Opt{dir}) {
        pod2usage("-dir is required");
    }
    if (!-d $Opt{dir}) {
        die "-dir $Opt{dir} does not exist\n";
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--counts> 1234567

=item B<--counts> /path/to/file

Count of input sequences, either as an integer or in a file in the format
produced by filter_one_fastq.pl: one line, "## sequences written, ## excluded".

=item B<--dir> /path/to/run/dir

REQUIRED.  Directory where STAR was run.

=item B<--field> FieldName=Value

Additional fields to include in output file; typically used to identify the
particular library/lane on which STAR was run.  Multiple fields may be
defined.  Ignored if no output file is given.

=item B<--output> /path/to/file

File to write tab-delimited summary from Log.final.out to; suitable for
loading into R (or Excel, if you must).  If not provided, no summary will be
written.

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=back

=head1 AUTHOR

 Peter Chines - pchines@mail.nih.gov

=head1 LEGAL

This software/database is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software/database is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S.  Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut
