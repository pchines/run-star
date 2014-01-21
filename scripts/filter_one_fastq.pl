#!/usr/bin/perl -w
# $Id: filter_one_fastq.pl 6352 2013-04-10 15:02:47Z pchines $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
our %Opt;

=head1 NAME

filter_one_fastq.pl - filter NPF reads from Casava 1.8-style FASTQ file

=head1 SYNOPSIS

Filter out reads that fail the chastity filter (and have deflines in the
style of CASAVAv1.8), writing to stdout:

  filter_one_fastq.pl in_R1_001.fq.gz [in_R1_002.fq.gz...] > Read1

For complete documentation, run C<filter_one_fastq.pl -man>

=head1 DESCRIPTION

CASAVA 1.8 writes standard FASTQ format files.  The defline (first line) of
each sequence includes additional information about the sequence, including
the read number, a flag to indicate whether the read failed the chastity
filter, and the actual barcode sequence read.

This program reads a FASTQ file, assuming that each sequence is represented
by exactly four lines (defline, sequence, plus-line, qualities).  If a
CASAVA-style defline is found, it will filter out reads that fail the
chastity filter.  If one or more reads without such a defline is found, a
single warning will be output to stderr, but all of these reads will be
retained in the output.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();
my $ofh = Open($Opt{output}, 'w');
my $warned;
my $incl = 0;
my $excl = 0;
for my $infile (@ARGV) {
    my $ifh = Open($infile);
    while (<$ifh>) {
        if (/^\@\S+\s\d+:([YN]):/) {
            if ($1 eq 'Y') {
                for (my $i = 0; $i < 3; ++$i) {
                    $_ = <$ifh>;
                }
                ++$excl;
                next;   # failed chastity filter
            }
        }
        elsif (/^\@/) {
            if (!$warned) {
                $warned = warn "FASTQ file '$infile' is not CASAVAv1.8-style; "
                    . "will pass all reads through\n";
            }
        }
        else {
            chomp;
            die "FASTQ file is corrupt or does not have 4 lines per sequence.\n"
                . "invalid defline: '$_' ";
        }
        # output record
        ++$incl;
        print $ofh $_;
        if ($Opt{'length'}) {
            for (my $i = 0; $i < 3; ++$i) {
                $_ = <$ifh>;
                chomp;
                print $ofh substr($_,0,$Opt{'length'}), "\n";
            }
        }
        else {
            for (my $i = 0; $i < 3; ++$i) {
                $_ = <$ifh>;
                print $ofh $_;
            }
        }
    }
} # loop @ARGV
if ($Opt{counts}) {
    my $cfh;
    if (uc($Opt{counts}) eq "STDERR") {
        print STDERR "$incl sequences written, $excl excluded\n";
    }
    else {
        $cfh = Open($Opt{counts}, 'w');
        print $cfh "$incl sequences written, $excl excluded\n";
    }
}

#------------
# End MAIN
#------------

sub process_commandline {
    # Set defaults here
    %Opt = (
            output  => '-',
            );
    GetOptions(\%Opt, qw(counts=s length=i output=s
                manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "filter_one_fastq.pl, ", q$Revision: 6352 $, "\n"; }
    # If non-option arguments are required, uncomment next line
    if (!@ARGV) {
        pod2usage("Expecting one or more FASTQ file(s).\n"
              . "All must be for the same read, e.g. all read 1 or all read2");
    }
    if ($Opt{counts} && uc($Opt{counts}) eq "STDOUT") {
        $Opt{counts} = '-';
    }
    if ($Opt{counts} && $Opt{counts} eq $Opt{output}) {
        die "Will not write FASTQ and counts to the same location\n";
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--counts> stderr

=item B<--counts> /path/to/file

Destination for counts of sequences: one line, in the format "## sequences
written, ## excluded".

=item B<--length> N

Maximum length of sequences to output.  If sequences are longer, they will be
trimmed at the 3' end.

=item B<--output> /path/to/file

File to write filtered output to; defaults to STDOUT.  May end with .gz to
write gzip-compressed file.

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
