#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use GTB::File qw(Open);
use GTB::Run qw(as_number);
our %Opt;
our $OutfileFormat; # sprintf format for generating filenames
our @Outfiles;  # filenames for all potential output files
our @Ofh;       # filehandles for actual output files
our @Count;     # counts of reads for each filehandle
our %Mate;      # potentially large hash of unmated (yet) paired reads
                # keys are read names (less /1 or /2); values are batch#
                # concatenated with FASTQ for known mate or
                # -1, signifying mate was seen, but won't be written, or
                # -2, signifying that mate pair has already been written (see
                # unique option)

use constant {
    SAM_QNAME   => 0,
    SAM_FLAG    => 1,
    SAM_CHROM   => 2,
    SAM_POS     => 3,
    SAM_SEQ     => 9,
    SAM_QUAL    => 10,
};
use constant {
    FLAG_PAIRED     => 0x001,
    FLAG_UNMAPPED   => 0x004,
    FLAG_REVERSE    => 0x010,
    FLAG_FIRST      => 0x040,
    FLAG_SECOND     => 0x080,
};

=head1 NAME

bam2fastq.pl - create FASTQ files from BAM files

=head1 SYNOPSIS

Create gzipped FASTQ files from a BAM file.  Output files will have same
prefix as BAM file, plus ".0.fq.gz", ".1.fq.gz" or ".2.fq.gz" extension:

  bam2fastq.pl [options] <bam_file>

See bam2fastq.pl -man for detailed description of options.

=head1 DESCRIPTION

Create FASTQ files from BAM files.  The output files are gzip-compressed by
default, and will be written to the current directory, with the same name as
the BAM file, plus an extension.  The extension is ".1.fq.gz" and ".2.fq.gz"
ifor paired reads, and ".0.fq.gz" for unpaired reads, whether because from a
single-end run, or because orphaned by filtering options.  See C<--prefix>
and C<--suffix> options for how to change these defaults.

=cut

{   # closure
    my $round_robin = -1;
sub batch_by_round_robin {
    ++$round_robin;
    if ($round_robin * 3 >= @Outfiles) {
        $round_robin = 0;
    }
    return $round_robin;
}}  # end closure

{   # closure
    my $reads = 0;
    my $nbatch = 0;
sub batch_by_number_of_reads {
    ++$reads;
    if ($reads >= $Opt{number}) {
        $reads = 0;
        ++$nbatch;
        for my $r (0..2) {
            push @Outfiles, sprintf $OutfileFormat, $nbatch, $r;
        }
    }
    return $nbatch;
}}  # end closure

#------------
# Begin MAIN
#------------

process_commandline();
my $fh = Open("samtools view $Opt{filter} $ARGV[0] | ");
my $q20 = 0;
my $curr_pos = 0;
my %current;

my $unpaired = 0;
my $last_chr = q{};
while (<$fh>) {
    chomp;
    my @s = split "\t";
    if ($Opt{flush} && $s[SAM_CHROM] ne $last_chr) {
        my $unp = flush_buffer();
        $unpaired += $unp;
        if ($last_chr) {
            warn "Switch to $s[SAM_CHROM]: drop $unp unpaired mates\n";
        }
        %Mate = ();
        $last_chr = $s[SAM_CHROM];
    }
    $s[SAM_QNAME] =~ s{/[12]$}{};   # strip read number (will use flag field)
    if ($Opt{minQ20}) {
        # count the number of bases with qual between 20 and 93
        my $qual_count = ( $s[SAM_QUAL] =~ tr/\065-\176// );
        if ($qual_count < $Opt{minQ20}) {
            ++$q20;
            if (!$Opt{single} && ($s[SAM_FLAG] & FLAG_PAIRED)) {
                $Mate{$s[SAM_QNAME]} ||= -1;
            }
            next;
        }
    }
    if ($s[SAM_FLAG] & FLAG_REVERSE) {
        $s[SAM_SEQ]  = reverse $s[SAM_SEQ];
        $s[SAM_SEQ]  =~ tr/ACGTacgt/TGCAtgca/;
        $s[SAM_QUAL] = reverse $s[SAM_QUAL];
    }
    my $read = 0;
    if ($s[SAM_FLAG] & FLAG_PAIRED) {
        $read = ($s[SAM_FLAG] & FLAG_SECOND) ? 2 : 1;
    }
    if ($Opt{illumina}) {
        if (($s[SAM_QUAL] =~ tr/\041-\176/\100-\235/) != length $s[SAM_QUAL]) {
            warn "Qualities for read $s[SAM_QNAME]/$read out of "
                . "expected (phred+33) range\n";
        }
    }
    if ($Opt{fixquals}) {
        if (($s[SAM_QUAL] =~ tr/\100-\235/\041-\176/) != length $s[SAM_QUAL]) {
            warn "Qualities for read $s[SAM_QNAME]/$read out of "
                . "expected (but wrong, phred+64) range\n";
        }
    }
    if ($Opt{unmapped}) {
        if ($s[SAM_FLAG] & FLAG_UNMAPPED) {
            # check for mapped mate (at current position)
            if ($read && $current{$s[SAM_QNAME]}) {
                output_read($current{$s[SAM_QNAME]});
            }
        }
        else {
            # retain reads at current pos, for pairing with unmapped mates
            if ($curr_pos ne "$s[SAM_CHROM]:$s[SAM_POS]") {
                # before throwing away, check to see whether mate is waiting
                while (my ($name, $ra_s) = each %current) {
                    if ($Mate{$name}) {
                        output_read($ra_s);
                    }
                }
                %current = ();
                $curr_pos = "$s[SAM_CHROM]:$s[SAM_POS]";
            }
            $current{$s[SAM_QNAME]} = \@s;
            next;
        }
    }
    output_read(\@s);
}
$unpaired += flush_buffer();
# finish and report read counts
my @type = ("single-end", "paired-end read 1", "paired-end read 2");
if ($unpaired) {
    $type[0] = "unpaired";
}
my @total;
for (my $i = 0; $i < @Outfiles; ++$i) {
    if ($Ofh[$i]) {
        if ( !$Opt{quiet} ) {
            print "Wrote $Count[$i] $type[$i % 3] reads to $Outfiles[$i]\n";
        }
        $total[$i % 3] += $Count[$i];
        close $Ofh[$i] or warn "Failed to close $Outfiles[$i], $!\n";
    }
}
if ( $Opt{counts} ) {
    if ($total[1] != $total[2]) {
        warn "WARNING: Inconsistent numbers of Read1 ($total[1]) and "
            . "Read2 ($total[2]) written!\n";
    }
    if (uc($Opt{counts}) eq "STDERR") {
        print STDERR "$total[1] read-pairs, $total[0] unpaired reads written\n";
    }
    else {
        my $cfh = Open($Opt{counts}, 'w');
        print $cfh "$total[1] read-pairs, $total[0] unpaired reads written\n";
    }
}
if ( !$Opt{quiet} ) {
    if ($Opt{batches} || $Opt{number}) {
        print "\n";
        for (my $i = 0; $i < 3; ++$i) {
            print "Total $total[$i] $type[$i] reads\n" if $total[$i];
        }
    }
    if ($Opt{minQ20}) {
        printf "Dropped %d reads with fewer than %d Q20 bases\n",
            $q20 || 0, $Opt{minQ20};
    }
    if ($unpaired) {
        print "There were $unpaired paired reads where mate was never seen\n";
    }
}

# END of MAIN

sub output_read {
    my ($ra) = @_;
    my $read = $ra->[SAM_FLAG] & FLAG_PAIRED
                    ? ($ra->[SAM_FLAG] & FLAG_SECOND ? 2 : 1)
                    : 0;
    my $mate_read = $read ? 3-$read : 0;
    my $match = "$ra->[SAM_QNAME]/$mate_read";
    my $data = "\@$ra->[SAM_QNAME]/$read\n$ra->[SAM_SEQ]\n+\n$ra->[SAM_QUAL]\n";
    my $fileno = $Opt{single} ? 0 : $read;
    my $batch;
    my $mate;
    # do we expect a mate for this read?
    if ($fileno) {
        $mate = $Mate{$match};
        if (defined $mate) {
            if ($mate =~ s/^(\d+)//) {
                $batch = $1;
            }
            elsif ($mate =~ /^-1/) {
                $fileno = 0;
                $batch  = 0;
            }
            elsif ($mate =~ /^-2/) {
                return;
            }
            else {
                die "FATAL: Should never get here: mate = $mate ";
            }
            if ($Opt{unique}) {
                $Mate{$match} = -2;
                $Mate{"$ra->[SAM_QNAME]/$read"} = -2;
            }
            else {
                delete $Mate{$match};
            }
        }
    }
    elsif ($Opt{unique}) {
        if ($Mate{$match}) {
            return;
        }
        $Mate{$match} = -2;
    }
    if (!defined $batch) {
        $batch = $Opt{number} ? batch_by_number_of_reads()
                              : batch_by_round_robin();
    }
    # if found mate
    if ($mate) {
        # always print read 1 first, then read 2 (for named pipes)
        my $ofh = get_output_handle(1 + $batch * 3);
        print { $ofh } ($fileno == 1 ? $data : $mate);
        $ofh = get_output_handle(2 + $batch * 3);
        print { $ofh } ($fileno == 1 ? $mate : $data);
    }
    # if don't expect mate
    elsif (!$fileno) {
        my $ofh = get_output_handle($batch * 3);
        print { $ofh } $data;
    }
    # save until find mate
    else {
        $Mate{"$ra->[SAM_QNAME]/$read"} = $batch . $data;
    }
}

sub get_output_handle {
    my ($fileno) = @_;
    if (!$Ofh[$fileno]) {
        $Ofh[$fileno] = Open($Outfiles[$fileno], 'w');
    }
    ++$Count[$fileno];
    return $Ofh[$fileno];
}

sub flush_buffer {
    my $unp = 0;
    for my $data (values %Mate) {
        next if $data =~ /^-/;
        ++$unp;
        $data =~ s/^\d+//;
        my $ofh = get_output_handle(0);
        print { $ofh } $data;
    }
    return $unp;
}

sub my_basename {
    my ($infile) = @_;
    $infile =~ s{\.(?:bam|f(?:ast)?[aq](?:\.gz)?)?$}{};
    $infile =~ s{^.+/}{};
    return $infile;
}

sub process_commandline {
    %Opt = (filter  => '-F 0x700',
            suffix  => '.fq.gz',
            );
    GetOptions(\%Opt, qw(
          batches=i
          counts=s
          filter|f=s
          fixquals
          flush
          illumina
          minQ20=i
          number=s
          prefix=s
          quiet
          r0=s r1=s r2=s
          single
          suffix=s
          unique
          unmapped
          yes
          help+
          manual
          version
          )
    ) || pod2usage(0);
    if ( $Opt{manual} ) { pod2usage( verbose => 2 ); }
    if ( $Opt{help} )   { pod2usage( verbose => $Opt{help} - 1 ); }
    if ( $Opt{version} ) { die "$0, ", q$Revision: 3496 $, "\n"; }
    if ( @ARGV != 1 ) {
        pod2usage("Must provide a single BAM file");
    }
    if ( !$Opt{prefix} ) {
        $Opt{prefix} = my_basename($ARGV[0]);
    }
    elsif (-d $Opt{prefix} ) {
        $Opt{prefix} =~ s{/?$}{/};  # ensure ends with slash
        $Opt{prefix} .= my_basename($ARGV[0]);
    }
    if ( $Opt{suffix} ) {
        $Opt{suffix} =~ s{^\.?}{.}; # ensure begins with dot
    }
    if ( $Opt{batches} && $Opt{number} ) {
        pod2usage("Can't specify both --batches and --number at the same time");
    }
    elsif ( $Opt{batches} ) {
        if ($Opt{batches} < 1 || $Opt{batches} > 26) {
            die "Sorry, only -batches values from 1 to 26 are allowed.\n"
                . "If you need to split into more sets of files, see "
                . "the -number option.\n";
        }
        $OutfileFormat = $Opt{single} ? "$Opt{prefix}.%c$Opt{suffix}"
                                      : "$Opt{prefix}.%c.%d$Opt{suffix}";
        for my $b (1..$Opt{batches}) {
            for my $r (0..2) {
                push @Outfiles, sprintf $OutfileFormat, $b+96, $r;
            }
        }
    }
    elsif ( $Opt{number} ) {
        $Opt{number} = as_number($Opt{number});
        $OutfileFormat = $Opt{single} ? "$Opt{prefix}.%03d.$Opt{suffix}"
                                      : "$Opt{prefix}.%03d.%d$Opt{suffix}";
        for my $r (0..2) {
            push @Outfiles, sprintf $OutfileFormat, 0, $r;
        }
    }
    elsif ( $Opt{r1} && $Opt{r2} ) {
        @Outfiles = ($Opt{r0} || '/dev/null', $Opt{r1}, $Opt{r2});
    }
    elsif ( $Opt{r1} || $Opt{r2} ) {
        pod2usage(<< 'END_MSG');
Must specify both -r1 and -r2 if either is used.  If a single output file
is desired, use -single.  If only read 1 or read 2 is desired, supply 
appropriate flags to -filter, e.g. -filter '-f 0x40 -F 0x700' (read 1,
no secondary alignments, no QC failures, no duplicates)
END_MSG
    }
    else {
        if ( $Opt{single} ) {
            @Outfiles = ( "$Opt{prefix}$Opt{suffix}" );
        }
        else {
            @Outfiles = map { "$Opt{prefix}.$_$Opt{suffix}" } (0..2);
        }
    }
    if (!$Opt{r1}) {
        my @exist = grep { -f $_ } @Outfiles;
        if ( @exist ) {
            if ( !$Opt{yes} ) {
                print STDERR "File(s) @exist exist.  Overwrite? (y/N) ";
                my $ok = <STDIN>;
                if ($ok !~ /^Y/i) {
                    die "Aborted; files unchanged\n";
                }
            }
            # Not sure why I am deleting files; we'll overwrite them anyway
            # might be to ensure that they are writable? Interferes with
            # named pipes, so let's get rid of it (4/10/13)
            #unlink @exist or die "Could not remove existing file, $!\n";
        }
    }
}

=head1 OPTIONS

=over

=item B<--batches> N

Number of sets of FASTQ files to write.  When this option is given, the
filename is composed of the prefix (see below), a dot followed by a batch
letter (a-z), a dot followed by the read number (if applicable), and the
".fq.gz" suffix, e.g. "my_prefix.a.1.fq.gz".  The program doesn't handle
values greater than 26; this can be considered a bug, if you like.
Mutually exclusive with C<--number> option.

=item B<--counts> counts.txt

Optional file to receive simplified count of read pairs and unpaired reads
written.  Does not affect the more detailed accounting of reads per output
file written to STDOUT; use C<--quiet> to turn the latter off.

=item B<--filter> "<samtools view options>"

Command line options that will be used to filter BAM reads.  Defaults to
C<-F 0x700>, which excludes secondary alignments, duplicate reads and reads
that failed the chastity filter.  Use any combination of flags that samtools
view understands.

=item B<--fixquals>

If your BAM file is broken, and has Phred+64 qualities, use this option to
write your FASTQ with typical (Sanger) Phred+33 qualties.

=item B<--flush>

BAM files are now so large that even tracking the unmated reads becomes a
memory burden.  This option assumes that the BAM is sorted by chrom, and will
completely discard the mate information at the end of each chromosome.  Don't
use this option if your input BAM is not coordinate sorted, or if you expect
a significant number of read pairs to be split across chromosomes.

=item B<--illumina>

Write Illumina Phred+64, rather than typical Phred+33 quality scores.

=item B<--minQ20> N

Sets the minimum number of Q20 bases in a read.  Default is to allow any number
(including zero).

=item B<--number> N

The number of reads to write to each batch of output files.  May use K, M or
B suffix to specify large values (powers of 10).  The output files will be
named using the specified prefix and suffix, along with a zero-padded numeric
batch number: e.g "my_prefix.000.1.fq.gz".  This option is mutually exclusive
with C<--batches> option.

=item B<--prefix> output/prefix

Output prefix; may contain a path; all directories in the path must already
exist.  If --prefix is not specified, bam2fastq.pl will place output files in
the current directory, with a prefix based on the BAM file name.

=item B<--quiet>

Do not write count of reads to STDOUT at program end.

=item B<--r0> /dest/for/unpaired

=item B<--r1> /dest/for/read1

=item B<--r2> /dest/for/read2

When either C<-r1> or C<-r2> is specified, both must be; C<-r0> need not be
used, and can only be used when both C<-r1> and C<-r2> are specified.  All
first reads passing the filters will be written to the C<-r1> file, and
corresponding second reads will be written to the C<-r2> file.  Unpaired or
orphan reads will be written to the C<-r0> file, if one is specified.

These options are intended mainly for writing to named pipes, and as such no
testing for the existence of the files or addition of suffixes is performed.
If you just want to write read1 or just read2, or even just the unpaired
reads, use C<--single> with the appropriate flags to C<--filter>,

=item B<--single>

Output a single FASTQ file, called PREFIX.fq.gz.  If this option is not
specified, unpaired reads are output to PREFIX.0.fq.gz, and paired reads to
PREFIX.1.fq.gz and PREFIX.2.fq.gz.

=item B<--suffix> .fq.gz

Suffix for created files.  Use .fq or .fastq to avoid gzip compression.

=item B<--unique>

Output each pair only once.  Normally, a read pair appears in a BAM only
once, or if a read pair is aligned to more than one location, all but one
of the alignments is marked secondary (bam flag 0x100).  If this is the case,
then the standard filter '-F 0x700' will already remove them.  If you specify
your own C<--filter> option, just be sure to include this flag, if there are
duplicates.  Only use the C<--unique> option if the BAM file is not
well-formed, and has duplicate entries that are not flagged as secondary.
PLEASE NOTE: using this option will require a great deal of memory.

=item B<--unmapped>

Extract only unmapped reads (and their possibly-mapped mates, for those
stored in BAM using the position of their mate).  To extract only the
unmapped reads, and not the mates, use option C<--filter "-f 4"> instead.

=item B<--yes>

Force "yes" response to question to overwrite (clobber) output files.  When
this option is not specified, the program will pause to ask, if it detects
that files exist.

=item B<--help|--manual>

Display documentation: brief synopsis, or complete manual.

=back

=head1 AUTHORS

 Peter Chines <pchines@mail.nih.gov>
 
=head1 LEGAL

This software is "United States Government Work" under the terms of
the United States Copyright Act.  It was written as part of the authors'
official duties for the United States Government and thus cannot be
copyrighted.  This software is freely available to the public for
use without a copyright notice.  Restrictions cannot be placed on its present
or future use. 

Although all reasonable efforts have been taken to ensure the accuracy and
reliability of the software and data, the National Human Genome Research
Institute (NHGRI) and the U.S. Government does not and cannot warrant the
performance or results that may be obtained by using this software or data.
NHGRI and the U.S. Government disclaims all warranties as to performance,
merchantability or fitness for any particular purpose. 

In any work or product derived from this material, proper attribution of the
authors as the source of the software or data should be made, using "NHGRI
Genome Technology Branch" as the citation. 

=cut
