#!/usr/bin/perl -w
# $Id: run_star.pl 6612 2013-11-23 14:14:33Z pchines $

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Temp qw(tempfile tempdir);
use File::Path qw(mkpath);
use GTB::File qw(Open);
use GTB::File::Iter;
use GTB::Run::SGE;

our %Opt;
our $EMPTY = q{};
our $SPACE = q{ };
our $DPS_READS = '/cluster/ifs/projects/solexa/reads';
# use environment vars and defaults to find STAR and processed genomes
our $STAR_BASE = $ENV{STAR_BASE} || '/cluster/ifs/projects/star';
our $STAR_BIN = "$STAR_BASE/bin/STAR";
if (!-l $STAR_BIN) {
    $STAR_BIN  = `which STAR`;
    chomp $STAR_BIN;
}
our $STAR_GENOMES = $ENV{STAR_GENOMES} || "$STAR_BASE/genomes";
our @STAR_GFILES = qw(Genome SA SAindex chrLength.txt chrName.txt
        chrStart.txt);
our $Gdir;  # actual directory with genome index
our $Batch; # batch number, for STAR jobs

=head1 NAME

run_star.pl - run STAR RNAseq aligner

=head1 SYNOPSIS

To align one or many lanes of single- or paired-end data to hg19,
for one or several samples (see FILES in complete documentation
for format of samples.txt):

  run_star.pl -genome hg19 samples.txt

For complete documentation, run C<run_star.pl -man>

=head1 DESCRIPTION

Runs STAR aligner on an SGE cluster node, taking advantage of shared memory
by using multiple CPU slots.  Input files will be filtered to remove reads
that fail QC (e.g. Illumina chastity filter).  In order to filter FASTQ
files, deflines must be formatted as they are in Illumina Casava 1.8
pipeline.  Automatically determines amount of memory to request, based on
size of genome and number of CPU slots.

Will automatically merge new results with existing results, if target BAM
file exists.  BUG: instead of checking flowcell/lane, currently makes the
assumption that thee new results are part of new readgroups, and will add
them to output file.

=cut

#------------
# Begin MAIN 
#------------

process_commandline();
find_genome_dir();
my $sge = GTB::Run::SGE->new(dryrun => $Opt{dryrun});
my $ra_files = read_fq_dirs();
my %samp_dirs;
# TODO: use one or two fewer max_cpus to account for filter jobs?
while (@$ra_files) {
    my $n = $Opt{max_cpus};
    if (@$ra_files <= $Opt{max_cpus}) {
        $n = @$ra_files;
    }
    elsif (@$ra_files < $Opt{max_cpus} * 2) {
        $n = int(@$ra_files/2);
    }
    my @batch = splice @$ra_files, 0, $n;
    create_star_script(\@batch, \%samp_dirs);
    sleep 1;    # let one script start before sending next
}
for my $samp (keys %samp_dirs) {
    create_merge_script($samp, $samp_dirs{$samp});
}

#------------
# End MAIN
#------------

sub get_star_version {
    my $version;
    my $bin = $Opt{star};
    if (-l $Opt{star}) {
        $bin = readlink $Opt{star};
    }
    if ($bin =~ /_v?([\d.]+\w?)(?:_r\d+)?$/) {
        $version = $1;
    }
    else {
        unlink 'Log.out';
        system $bin;
        my $fh = Open('Log.out');
        while (<$fh>) {
            if (/STAR svn revision compiled=STAR_([\d\.]+\w?)/) {
                $version = $1;
                last;
            }
        }
        unlink 'Log.out', 'Log.timing.out', 'Log.progress.out';
    }
    if (!$version) {
        warn "Wasn't able to determine a version number for $bin\n";
    }
    return $version || $EMPTY;
}

sub find_genome_dir {
    my $version = get_star_version();
    if (-d $Opt{genome}) {
        $Gdir = $Opt{genome};
    }
    elsif (-d $STAR_GENOMES) {
        opendir D, $STAR_GENOMES
            or die "Can't read directory '$STAR_GENOMES', $!\n";
        my @dirs = readdir(D);
        my @vdir = grep { /\Q$Opt{genome}\E.+\Q$version\E/ } @dirs;
        if (!@vdir) {
            @vdir = grep { /\Q$Opt{genome}\E/ } @dirs;
            if (!@vdir) {
                my $known = join("\n", map { "    $_" } @dirs);
                die << "END_MSG";
No genome matching '$Opt{genome}' found in '$STAR_GENOMES'.
Known genomes include:
$known
Please select one of the above, or use -genome option to point to a
STAR-formatted genome of your choice.
END_MSG
            }
        }
        if (@vdir == 1) {
            $Gdir = "$STAR_GENOMES/$vdir[0]";
        }
        else {
            @vdir = reverse sort @vdir;
            $Gdir = "$STAR_GENOMES/$vdir[0]";
            my $matching = join("\n", map { "    $_" } @vdir);
            warn << "END_MSG";
Found more than one STAR genome dir matching '$Opt{genome}':
$matching
Selected '$Gdir'; if this is not what you wanted, specify full path to
genome directory, e.g. -genome $STAR_GENOMES/$Gdir
END_MSG
        }
    }
    else {
        die << 'END_MSG';
Can't find STAR genome directory.  Please set STAR_BASE or STAR_GENOMES
environment variable, or specify the complete path to the genome directory
using the -genome option.
END_MSG
    }
    my @missing = grep { !-f "$Gdir/$_" } @STAR_GFILES;
    if (@missing) {
        die << "END_MSG";
STAR genome directory '$Gdir' does not contain all of the necessary
files.  Missing: @missing.
Please specify the complete path to the genome directory using the
-genome option.
END_MSG
    }
    $Gdir = File::Spec->rel2abs($Gdir);
}

sub create_star_script {
    my ($ra_dir_info, $rh_samp_dirs) = @_;
    ++$Batch;
    my ($fh, $tmp_script) = tempfile("star.$Batch.sh.XXXX",
            DIR => $Opt{tmpdir}, SUFFIX => '.sh');
    my $slots = @$ra_dir_info;
    my $tot_mem = calc_memory_needed();
    my $mem_per_sort = 1024 * 1024 * 1024;
    $tot_mem += 2 * $mem_per_sort * $slots; # samtools sort uses 2x mem
    my $mem_per_slot = int( $tot_mem / $slots );
    my $first = 1;
    print $fh "#!/bin/bash\n";
    print $fh "# run_star.pl-generated script to run $slots parallel STAR jobs\n";
    printf $fh "cd $Opt{tmpdir}\n";
    for my $rh (@$ra_dir_info) {
        $rh->{dir} = File::Spec->rel2abs(
                tempdir( "star.$Batch.$rh->{rg_name}.$rh->{sample}.XXXX",
                    DIR => $Opt{tmpdir} )
                );
        # copy config file to dir (TODO: is this still the right name?)
        if ($Opt{config}) {
            system "cp $Opt{config} $rh->{dir}/Parameters.1.In.txt";
        }
        print $fh "cd $rh->{dir}\n";
        print $fh "mkfifo Read1 Read2\n";
        $rh->{fastq_counts} = "Read1.counts";
        my $fparam = $rh->{filter_params} || $Opt{filter_params} || $EMPTY;
        if ($rh->{in_bam}) {
            print $fh "bam2fastq.pl -r1 Read1 -r2 Read2 $fparam "
                . "-c $rh->{fastq_counts} $rh->{in_bam} &\n";
        }
        elsif ($rh->{r1} && $rh->{r2}) {
            print $fh "filter_one_fastq.pl -counts $rh->{fastq_counts} "
                . "$fparam $rh->{r1} > Read1 &\n";
            print $fh "filter_one_fastq.pl -counts Read2.counts "
                . "$fparam $rh->{r2} > Read2 &\n";
        }
        elsif ($rh->{r1}) {
            print $fh "filter_one_fastq.pl -counts $rh->{fastq_counts} "
                . "$fparam $rh->{r1} > Read1 &\n";
            $rh->{star_params} .= " --readFilesIn Read1";
        }
        else {
            die "FATAL: No 'in_bam', nor 'r1' and 'r2' for $rh->{sample} "
                . "$rh->{rg_name}";
        }
        my $sparam = $rh->{star_params} || $Opt{star_params} || $EMPTY;
        print $fh "$Opt{star} --genomeDir $Gdir --genomeLoad LoadAndKeep "
            . "$sparam --outStd SAM -outReadsUnmapped Fastx "
            . qq{ | awk '{print} /^\@HD/{print "$rh->{sam_header}"}' } #'for vim
            . " | samtools view -uS - "
            . " | samtools sort -m $mem_per_sort - aligned &\n";
        if ($first) {
            print $fh "sleep 10\n";  # allow first STAR to begin loading genome
            $first = 0;
        }
        print $fh "cd $Opt{tmpdir}\n";
    }
    print $fh "wait\n"; # for STAR aligners to finish
    print $fh "$Opt{star} --genomeDir $Gdir --genomeLoad Remove\n";
    my %samples;
    my $sleep = 1;
    for my $rh (@$ra_dir_info) {
        print $fh "check_star_run.pl -dir $rh->{dir} "
            . " -counts $rh->{dir}/$rh->{fastq_counts} "
            . ($Opt{report} ? " -out $Opt{report}.align.txt " : $EMPTY)
            . ($rh->{in_bam} ? " -f bam=$rh->{in_bam} " : $EMPTY)
            . " -f r1='$rh->{r1}' "
            . " -f sample=$rh->{sample} -f readgroup=$rh->{rg_name} "
            . " || mv aligned.bam incomplete.bam\n";
        $rh_samp_dirs->{$rh->{sample}} ||= {};
        my $rh_s = $rh_samp_dirs->{$rh->{sample}};
        $rh_s->{out_file} ||= $rh->{out_file};
        if (exists $rh->{md_params}) {
            $rh_s->{md_params} ||= $rh->{md_params};
            if ($rh_s->{md_params} ne $rh->{md_params}) {
                warn "Warning: two different 'md_params' settings for "
                    . "sample '$rh->{sample}': '$rh_s->{md_params}' and "
                    . "'$rh->{md_params}'\n";
            }
        }
        push @{ $rh_s->{files} }, $rh;
        $samples{ $rh->{sample} } = 1;
        ++$sleep;
    }
    #print $fh "wait\n"; # for SAM -> BAM conversions to finish
    # TODO: would be nice to test that everything worked as planned
    #print $fh "rm $tmp_script\n";   # cleanup this script
    close $fh or die "Error closing temp script file '$tmp_script', $!\n";
    my $jid = $sge->submit_job(
            cmd      => $tmp_script,
            sge_opts => { '-l'  => "mf=$mem_per_slot,h_vmem=$tot_mem",
                          '-pe' => "make-dedicated $slots",
                          '-V'  => $EMPTY,
                          '-cwd'=> $EMPTY,
                        },
            );
    for my $s (keys %samples) {
        push @{ $rh_samp_dirs->{$s}{jobs} }, $jid;
    }
}

sub calc_memory_needed {
    my $suffix_array = "$Gdir/SA";
    my $ndx_size = -s $suffix_array || die "Can't find $suffix_array";
    return int(1.25 * $ndx_size);
}

sub read_fq_dirs {
    my @files;
    my @jids;
    my $error;
    my %seen;
    for my $p (@ARGV) {
        if (-f $p) {
            if ($p =~ /\.bam$/) {
                # must first create FASTQs
                die "BAM file conversion not yet implemented";
            }
            else {
                my $iter = GTB::File::Iter->new( file => $p );
                while (my $rh = $iter->next_hash()) {
                    $rh->{sample} ||= $rh->{SAMPLE};
                    $rh->{in_bam} ||= $rh->{BAM};
                    if (!$rh->{sample}) {
                        $error = warn "Input file $p does not include "
                            . "a 'sample' column.\n";
                    }
                    if ($rh->{in_dir}) {
                        if (-d $rh->{in_dir}) {
                            my $key = join('.', $rh->{in_dir}, $rh->{lane}||0);
                            if ($seen{$key}) {
                                $error = warn "Directory $rh->{in_dir} appears more than once, without a specific lane indicator.\n";
                                next;
                            }
                            $seen{$key} = 1;
                            my ($ra_g, $err) =
                                get_read_pair_files($rh->{in_dir}, $rh->{lane});
                            $error ||= $err;
                            # batch files within each group
                            for my $rh_g (@$ra_g) {
                                my $i = 0;
                                my $r1 = $EMPTY;
                                my $r2 = $EMPTY;
                                for my $raf (@{ $rh_g->{files} }) {
                                    $r1 .= $raf->[0] . $SPACE;
                                    $r2 .= $raf->[1] . $SPACE;
                                    ++$i;
                                    if ($i >= $Opt{batch}) {
                                        push @files, { %$rh,
                                            lane => $rh_g->{lane} || 0,
                                            r1   => $r1,
                                            r2   => $r2,
                                        };
                                        $i = 0;
                                        $r1 = $r2 = $EMPTY;
                                    }
                                }
                                if ($r1) {
                                    push @files, { %$rh,
                                         lane => $rh_g->{lane} || 0,
                                         r1   => $r1,
                                         r2   => $r2,
                                        };
                                }
                            }
                        }
                        else {
                            $error = warn "Dir '$rh->{in_dir}' not present "
                                . "for sample '$rh->{sample}'\n";
                        }
                    }
                    elsif ($rh->{in_fq1} && $rh->{in_fq2}) {
                        if (-f $rh->{in_fq1} && -f $rh->{in_fq2}) {
                            $rh->{lane} ||= 0;
                            push @files, { %$rh,
                                r1 => File::Spec->rel2abs($rh->{in_fq1}),
                                r2 => File::Spec->rel2abs($rh->{in_fq2}),
                                };
                        }
                        else {
                            $error = warn "FASTQ files '$rh->{in_fq1}' and "
                                . "'$rh->{in_fq2}' are not both present\n";
                        }
                    }
                    elsif ($rh->{in_fq1}) { # single-end reads
                        if (-f $rh->{in_fq1}) {
                            $rh->{lane} ||= 0;
                            push @files, { %$rh,
                                r1 => File::Spec->rel2abs($rh->{in_fq1}),
                                };
                        }
                        else {
                            $error = warn "FASTQ file '$rh->{in_fq1}' "
                                . "is not present\n";
                        }
                    }
                    elsif ($rh->{in_bam}) {
                        if ($rh->{in_bam} =~
                                /(\d{6})_(\w+)_(\w+)\.(\d)\.(\d{6,8})\.bam$/) {
                            $rh->{run_date} ||= $1;
                            $rh->{machine}  ||= $2;
                            $rh->{flowcell} ||= $3;
                            $rh->{lane}     ||= $4;
                            $rh->{library}  ||= $5;
                            $rh->{run_name} = sprintf("%s_%s_%s",
                                map { $rh->{$_} } qw(run_date machine flowcell)
                                );
                        }
                        else {
                            # TODO: gather missing info from BAM header
                        }
                        $rh->{r1} = $EMPTY;
                        if (-f $rh->{in_bam}) {
                            push @files, { %$rh };
                        }
                        elsif (-f "$DPS_READS/$rh->{run_name}/$rh->{in_bam}") {
                            $rh->{in_bam} =
                                "$DPS_READS/$rh->{run_name}/$rh->{in_bam}";
                            push @files, { %$rh };
                        }
                        else {
                            $error = warn "BAM file '$rh->{in_bam}' not "
                                . "present for sample '$rh->{sample}'\n";
                        }
                    }
                    elsif (!grep { !$rh->{$_} } qw(run_date machine flowcell
                                lane library)) {
                        $rh->{r1} = $EMPTY;
                        $rh->{run_name} = sprintf("%s_%s_%s",
                                map { $rh->{$_} } qw(run_date machine flowcell)
                                );
                        $rh->{in_bam} = sprintf("$DPS_READS/%s/%s.%d.%d.bam",
                                    map { $rh->{$_} }
                                    qw(run_name run_name lane library)
                                    );
                        if (-f $rh->{in_bam}) {
                            push @files, { %$rh };
                        }
                        else {
                            $error = warn "Constructed BAM path "
                                . "'$rh->{in_bam}' does not exist for "
                                . "sample '$rh->{sample}'\n";
                        }
                    }
                    else {
                        $error = warn "No 'in_dir' and "
                            . "no 'in_fq1' and 'in_fq2' "
                            . "specified for sample '$rh->{sample}'\n";
                    }
                }
            }
        }
        elsif (-d $p) {
            my $s = $Opt{sample};
            if (!$Opt{sample}) {
                $s = $p;
                $s =~ s{.*/}{};
            }
            my ($ra_g, $err) = get_read_pair_files($p);
            $error ||= $err;
            for my $rh_g (@$ra_g) {
                push @files, map {
                    {in_dir => $p,
                     lane   => $rh_g->{lane},
                     sample => $s,
                     r1     => $_->[0],
                     r2     => $_->[1],
                    } } @{ $rh_g->{files} };
            }
        }
        else {
            $error = warn "'$p' is neither a file nor a directory\n";
        }
    }
    if ($error) {
        die "Aborting.\n";
    }
    my %samp_rgid;
    my %indl_rgid;
    for my $rh (@files) {
        my $samp = $rh->{sample};
        my $name = $Opt{prefix} ? "$Opt{prefix}.$samp" : $samp;
        my $path = File::Spec->join($Opt{outdir}, $name);
        $rh->{out_file} ||= $path;
        $rh->{out_file} =~ s/\.bam$//;
        $rh->{library}  ||= $samp;
        if (!defined $samp_rgid{ $samp }) {
            $samp_rgid{ $samp } = get_max_rgid("$rh->{out_file}.bam");
        }
        $rh->{rg_name} = $rh->{run_name};
        if (!$rh->{rg_name} && $rh->{flowcell}) {
            $rh->{rg_name} = join '_', map { $rh->{$_}||0 }
                            qw(run_date machine flowcell);
        }
        if (!$rh->{rg_name}) {
            $rh->{rg_name} = $rh->{in_dir} || $rh->{in_fq1};
            $rh->{rg_name} =~ s/(?:\.f(?:ast)?q)?(?:\.gz)?$//;
        }
        $rh->{rg_name} .= ".$rh->{lane}";
        if (!defined $indl_rgid{ $rh->{rg_name} }) {
            $indl_rgid{ $rh->{rg_name} } = ++$samp_rgid{ $samp };
        }
        $rh->{rgid} ||= $indl_rgid{ $rh->{rg_name} };
        $rh->{sam_header} ||= join "\t", '@RG', "ID:$rh->{rgid}",
            "SM:$samp", "LB:$rh->{library}", "PU:$rh->{rg_name}",
            "PL:Illumina", "CN:NISC";
    }
    # If long-running data preparation jobs were submitted, wait for them
    if (@jids) {
        $sge->wait(jobs => \@jids);
    }
    return \@files;
}

sub get_read_pair_files {
    my ($dir, $lane) = @_;
    $dir = File::Spec->rel2abs($dir);
    my @groups;
    my $error;
    my %pair;
    opendir(D, $dir) or die "Can't read directory '$dir', $!\n";
    while (my $f = readdir(D)) {
        # ignore invisibles, dirs, sample sheet, non-FASTQ
        if ($f =~ /^\./ || !-f "$dir/$f" || $f !~ /\.f(?:ast)?q(?:\.gz)?$/) {
            next;
        }
        if ($f =~ /^(\w+?(?:L00(\d))?)_R([12])_(\d+)\.f(?:ast)?q(?:\.gz)?$/) {
            if ($2 && $lane && $2 ne $lane) {
                next;
            }
            $pair{$2||0}{ "$1.$4" }[ $3 ] = "$dir/$f";
        }
        elsif ($f =~ /^(\S+?)\.([12])\.f(?:ast)?q(?:\.gz)?$/) {
            my ($fc, $rd) = ($1, $2);
            my $ln = 0;
            if ($fc =~ /\.(\d)\./) {
                $ln = $1;
            }
            if ($ln && $lane && $ln ne $lane) {
                next;
            }
            $pair{ $ln }{ $fc }[ $rd ] = "$dir/$f";
        }
        else {
            $error = warn "Can't tell whether '$f' is R1 or R2, "
                . "skipping file.\n";
        }
    }
    for my $lane (keys %pair) {
        my @pairs;
        push @groups, { lane => $lane, files => \@pairs };
        for my $ra (values %{ $pair{$lane} }) {
            if (!$ra->[1] || !$ra->[2]) {
                $error = warn "Incomplete pair: "
                    . ($ra->[1]||$ra->[2]) . "\n";
                next;
            }
            push @pairs, [ $ra->[1], $ra->[2] ];
        }
    }
    return (\@groups, $error);
}

sub get_max_rgid {
    my ($bam) = @_;
    my $rgid = 0;
    if (-f $bam) {
        my $fh = Open("samtools view -H $bam|");
        while (<$fh>) {
            if (/^\@RG.*ID:(\d+)/ && $1 > $rgid) {
                $rgid = $1;
            }
        }
    }
    return $rgid;
}

sub create_merge_script {
    my ($samp, $rh_samp) = @_;
    my $job = $samp;
    $job =~ s{[^\w.-]+}{_}g;
    my ($fh, $script) = tempfile("merge.$samp.XXXX",
            DIR => $Opt{tmpdir}, SUFFIX => '.sh' );
    print $fh "#!/bin/bash\n";
    print $fh "# run_star.pl-generated script to merge and cleanup $samp files\n";
    print $fh "cd $Opt{tmpdir}\n";
    my @bams;
    # archive existing file
    my $exist = $EMPTY;
    if ($rh_samp->{out_file} && -f "$rh_samp->{out_file}.bam") {
        $exist = $rh_samp->{out_file};
        if ($Opt{archivedir}) {    # use archive dir, same file name
            $exist =~ s#^.*/#$Opt{archivedir}/#;
        }
        elsif ($exist !~ s/aligned/archived/) {  # follow convention
            $exist .= ".archive";               # same dir, .archive suffix
        }
        print $fh "mv $rh_samp->{out_file}.bam $exist.bam\n";
        print $fh "mv $rh_samp->{out_file}\{,.bam}.bai $exist.bai\n";
        push @bams, "$exist.bam";
    }
    # merge BAM files
    my @dirs = map { $_->{dir} } @{ $rh_samp->{files} };
    push @bams, map { "$_/aligned.bam" } @dirs;
    my $inputs = join($SPACE, map { "I=$_" } @bams);
    # ..use MergeSamFiles to merge RG headers
    print $fh "java -Xmx4g -Djava.io.tmpdir=$Opt{tmpdir} "
        . "-jar $Opt{picard}/MergeSamFiles.jar "
        . "ASSUME_SORTED=true SORT_ORDER=coordinate CREATE_INDEX=true "
        . "OUTPUT=$rh_samp->{out_file}.tmp.bam $inputs\n";
    # flag duplicates, remove optical dups
    my $mdparams = $rh_samp->{md_params} || $Opt{md_params}
        || "READ_NAME_REGEX=':[0-9]:([0-9]+):([0-9]+):([0-9]+).*\$' REMOVE_OPTICAL=true";
    print $fh "java -Xmx4g -Djava.io.tmpdir=$Opt{tmpdir} "
        . "-jar $Opt{picard}/MarkDuplicates.jar "
        . "ASSUME_SORTED=true CREATE_INDEX=true $mdparams "
        . "METRICS_FILE=$rh_samp->{out_file}.metrics "
        . "OUTPUT=$rh_samp->{out_file}.bam "
        . "INPUT=$rh_samp->{out_file}.tmp.bam "
        . " && rm $rh_samp->{out_file}.tmp{,.bam}.ba?\n";
    # TODO: add QC step, to finish using temp dirs before deleting them?
    close $fh or die "Error closing $samp merge script '$script', $!\n";
    $sge->submit_job(
            cmd      => $script,
            hold_jid => $rh_samp->{jobs},
            memory   => '5g',
            sge_opts => { '-V'   => $EMPTY,
                          '-cwd' => $EMPTY,
                          '-N'   => "mrg.$job",
                        },
            );
}

sub process_commandline {
    # Set defaults here
    %Opt = (
            batch       => 99,
            debug       => 0,
            filter      => 1,
            genome      => 'hg19',
            max_cpus    => 12,
            picard      => '/cluster/ifs/projects/collins/gatk/picard',
            star        => $STAR_BIN,
            );
    GetOptions(\%Opt, qw(archivedir=s batch=i config=s debug+ dryrun
                filter! filter_params=s genome=s
                iknow library=s max_cpus=i md_params=s outfile=s picard=s
                report=s rgid=s sample=s star=s star_params=s
                tmpdir|tempdir=s
                manual help+ version)) || pod2usage(0);
    if ($Opt{manual})  { pod2usage(verbose => 2); }
    if ($Opt{help})    { pod2usage(verbose => $Opt{help}-1); }
    if ($Opt{version}) { die "run_star.pl, ", q$Revision: 6612 $, "\n"; }
    if (!$Opt{iknow}) {
        warn "This script is under development; things may change.  Use at your own risk.\n";
    }
    # If non-option arguments are required, uncomment next line
    if (@ARGV == 0) {
        pod2usage("Expected FASTQ directory or text file with multiple dirs");
    }
    # Add criteria to test option arguments below
    if ($Opt{max_cpus} < 1) {
        pod2usage("--max_cpus must be an integer greater than or equal to one");
    }
    if (!$Opt{star}) {
        pod2usage("STAR was not found; please specify binary location with "
                . "-star option,\nor set STAR_BASE environment variable");
    }
    $Opt{outfile} ||= './';
    if ($Opt{outfile} =~ s/\.bam$//) {
        (undef, $Opt{outdir}, $Opt{prefix}) =
            File::Spec->splitpath($Opt{outfile});
    }
    else {
        $Opt{outdir} = $Opt{outfile};
    }
    $Opt{tmpdir} ||= ".";
    for my $i (qw(archivedir outdir tmpdir report star picard)) {
        if ($Opt{$i}) {
            $Opt{$i} = File::Spec->rel2abs($Opt{$i});
        }
    }
    for my $i (qw(archivedir outdir tmpdir)) {
        if ($Opt{$i} && !-d $Opt{$i}) {
            mkpath($Opt{$i});
        }
    }
    if ($Opt{debug}) {
        $GTB::Run::VERBOSE = 1;
    }
    $Opt{library} ||= $Opt{sample};
    if ($Opt{report}) {
        $Opt{report} =~ s/\.txt(?:\.gz)?$//;
    }
}

__END__

=head1 OPTIONS

=over 4

=item B<--archivedir> /path

Optional directory where previous version (if any) of the sample BAM file is
moved to before creating new sample BAM.  If not provided, will create
.archive.bam in output directory.

=item B<--batch> N

Maximum number of FASTQ files in an input directory to bundle together as a
single batch.  Defaults to 99.

=item B<--config> /path/to/Parameters.1.In.txt

Optional STAR parameters file to be copied to working directory.

=item B<--debug>

Report commands that are being executed.

=item B<--dryrun>

Don't submit SGE commands, just generate batch scripts.

=item B<--genome> hg18

=item B<--genome> /path/to/star/genome/

Specify genome build to align against.  STAR genome built with same version
of program must be present.

=item B<--library> library_name

For SAM/BAM @RG header LB tag.  Defaults to same as sample name.

=item B<--max_cpus> 12

Number of STAR processes to run on a given node; defaults to 12.

=item B<--md_params> STRING

Parameters for MarkDuplicates.  Defaults to
"READ_NAME_REGEX=':[0-9]:([0-9]+):([0-9]+):([0-9]+).*\$' REMOVE_OPTICAL=true".

=item B<--outfile> /path/to/output.bam

=item B<--outfile> /path/to/output/dir/

=item B<--outfile> output.bam

Specify output location (directory and/or file prefix).  The directory must
be accessible to all nodes on the cluster via the same path (no local disks).
Temporary files will be written to the same directory.  The default is to
write to the current directory, using a filename composed of the prefix to
the first FASTQ file name and the prefix to the basename of the genome
reference aligned to.

=item B<--picard> /path/to/picard/jars

Use a different version of picard, if you like, but then some parameters
might not work (see C<md_params>).

=item B<--report> report_prefix

Where consolidated reports should be written.

=item B<--rgid> ID

For SAM/BAM @RG header ID tag.  Defaults to 1.

=item B<--sample> sample_name

For SAM/BAM @RG header SM tag.  Defaults to 'unknown'.

=item B<--star> /path/to/bin/STAR

=item B<--star> "/path/to/bin/STAR --opts"

Specify particular version of STAR to use; otherwise follows usual rules to
run the first one in your PATH.  May include commandline options.  Another 
way to include options is by supplying a C<-config> file.  A final way is to
provide "star_params", either within the tab-delimited input file, or on the
command line.

=item B<--tmpdir> /path/to/large/temp/dir

Path to use to write temporary files, including scripts, working dirs, and
picard temporary files.  There should be at least three times the size of
your raw data available at this location, and it must be writable by you, of
course.  Defaults to current working directory.

=item B<--help|--manual>

Display documentation.  One C<--help> gives a brief synopsis, C<-h -h> shows
all options, C<--manual> provides complete documentation.

=back

=head1 FILES

=head2 Samples File

A tab-delimited file that identifies the sample, the location of input files,
and various optional information about the library, sequencing run, and
parameters to use for STAR and Picard programs.  The file must begin with a
header line to identify the columns.  One of the columns must be named
'sample'.

Input files may be specified either by:

=over 4

=item 1 in_dir

For a directory full of pairs of FASTQ files, named as if directly from the
Illumina Casava pipeline.

=item 2 in_fq1 and in_fq2

For a single pair of FASTQ files.  For single-end reads, specify only in_fq1.

=item 3 in_bam

For a lane BAM file.  If no path is given, the program will look for the BAM
file in the (NISC) standard location:
/cluster/ifs/projects/solexa/reads/FLOWCELL

BAM file conversion was written for paired-end reads, but it may work with
single-end reads if you include an additional column "star_params" and
include the following option: "--readFilesIn Read1".  This is untested.

=back

Other fields that are interpreted and use include: run_name (or run_date,
machine, and flowcell), lane, library, star_params, md_params.

A single sample may appear on multiple lines of this file; indeed this is
expected, since RNAseq data is rarely run in just a single lane.

=head1 BUGS

Only runs under Sun Grid Engine.  Has minimal error checking.

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
