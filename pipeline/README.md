# SaaRclust pipeline

> **Ongoing work**

Clustering PacBio reads into their original chromosomes and directionalities using Strand-seq data - summarized in a [Snakemake](https://bitbucket.org/snakemake/snakemake) pipeline.


### How to use it

  1. **Download the pipeline**  
    Download the pipeline directory including Snakefile, the perl script processMinimaptab.pl, and a subdirectory "utils" containing the following Rscripts:
    - "install_SaaRclust.R"  : Rscript for installing SaaRclust R package
    - "SaaRclust_hardclust_pipeline.R" : Rscript for running hard clustering
    - "SaaRclust_softclust_pipeline.R" : Rscript for running soft clustering


  2. **Install required software:**

    * Install [minimap](https://github.com/lh3/minimap) for aligning Strand-seq reads to PacBio reads.
    
    
  3. **Prepare the input data**
    
    * Move to the directory containg the Snakefile
    * Create a directory with name "raw_reads"
    * Copy or (soft-link) the fastq file containing the raw Strand-seq reads in "raw_reads" directory
    * Create a directory "chunks" in "raw_reads" directory
    * Copy or (soft-link) the fasta formatted files containing PacBio reads in "raw_reads/chunks" directory. The names of files should have the format "{common_prefix}chunk{chunk_id}". Note that the files do not have to have the ".fasta" or ".fa" extension
    

  4. **Set up the parameters of the snakemake pipeline**

    * Open `Snakefile` and set the parameters of minimap and SaaRclust at the beginning of the file
      (such as Mosaicatcher) and to the R scripts.
    

  5. **Run Snakemake**

    * Run `snakemake` to compute all tasks locally
    
### The output files

Snakemake creates the following output files and directories:

  1. **Minimap alignment results in directory "aligns"**
    * The alignment files with names {PacBio_read_file_name}.maf.gz
    * Minimap output summaries with names {PacBio_read_file_name}.maf.gz.log
    
  2. **SaaRclust results**
    * SaaRclust summaries in directory "log"
    * SaaRclust hard and soft clustering (in different chunks) results saved as RData objects in directory "SaaRclust_results_{sample_name}/Clusters"
    * Some supplementary data like blacklisted alignments and some data quality measures

