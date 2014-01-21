run-star
========

Scripts for submitting STAR aligner jobs to a cluster

So far, these scripts are only known to work on the trek SGE cluster at
NHGRI.  We hope to modify them to make them work in other environments,
including LSF clusters.

### Prerequisites:
  - STAR aligner (v.2.3.0 or later recommended)
        https://code.google.com/p/rna-star/
  - STAR indexed genome (created with STAR --runMode genomeGenerate and,
    ideally, an appropriate GTF file to define known splice junctions)
  - Sun Grid Engine (or Son of Grid Engine) cluster
  - perl (v.5.8 or later)
  - java (v.1.6 or later)

### Environment Variables:
  - STAR_BASE - directory including bin/STAR (and possibly genomes dir)
  - STAR_GENOMES - directory with a subdirectory for each STAR-indexed
    genome; defaults to $STAR_BASE/genomes

### Files:

The program expects a tab-delimited input file with a header line.  The file
format is described in detail in the manual page, but in its simplest form,
it looks like this:

```
sample      in_bam
SampleName1 /path/to/s1/readgroup1.bam
SampleName1 /path/to/s1/readgroup2.bam
SampleName2 /path/to/s2/readgroup1.bam
```

### Running script:

Writes several batch scripts, as well as temporary files to the current
working directory, if no location is specified.  So I typically create a
run-specific temporary directory to work in, and send the final aligned
output to a separate, higher level directory:

```
$ mkdir tmp_TODAY && cd tmp_TODAY
$ run_star.pl -genome hg19 -out ../aligned/ -picard /path/to/jars samples.txt
```

*Note the trailing slash on the output directory; this is required.*

To create batch scripts, but not submit them, add the --dryrun option.
This will lose the dependencies between the scripts, but in general you
can always first run all of the star.*.sh scripts, then run all of the
merge.*.sh scripts.  The star scripts require a dedicated machine with lots
of RAM and at least 12 CPUs (set with --max_cpus); the merge scripts
typically require only 4-6 Gb RAM, and even this can be further reduced.

To see all of the options, use run_star.pl --man.


Ideas
-----

Proposed (not yet implemented) Environment Variables:
  - PICARD - directory in which to find PICARD jar files, including
    MergeSamFiles, and MarkDuplicates (a special updated version is
    included here, in javai subdirectory).  Currently defaults to
    /cluster/ifs/projects/Collins/gatk/picard, but can be set with
    --picard command line option.
  - IN_BAM_DIR - directory in which to look for BAM files without a full
    path.  Hard coded; only alternative for now is to use full paths.
