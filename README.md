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

### Running script:

Writes several batch scripts, as well as temporary files if not otherwise
specified, to the current working directory.

  run_star.pl -genome hg19 -picard /path/to/jars samples.txt

To create batch scripts, but not submit them, use the --dryrun option.
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
