# shellcheck shell=bash
# Machine-specific settings for the preprocessing scripts.
# Copy to paths.sh (ignored by git) and edit; every script sources it via common.sh.

# Reference files
REF_DIR=/path/to/ref
BOWTIE2_INDEX="$REF_DIR/hg38"                       # bowtie2 index prefix
CHROM_SIZES="$REF_DIR/hg38.chr.sizes"               # two columns: chrom, length
ADAPTERS="$REF_DIR/Sequencing_adaptors.fasta"       # Trimmomatic ILLUMINACLIP file
STAR_INDEX="$REF_DIR/hg38/STAR_genome"              # only for align_star.sh
HOMER_GENOME=hg38

THREADS=8

# Sample naming: ChIP BAMs are <SAMPLE_PREFIX>_<condition>_<replicate>.dedup.bam,
# controls are <CONTROL_PREFIX[condition]>_<replicate>.dedup.bam.
SAMPLE_PREFIX=RING1
CONDITIONS=(BF HS)          # bdgdiff compares CONDITIONS[0] (t1) with CONDITIONS[1] (t2)
REPLICATES=(rep1 rep2)
declare -A CONTROL_PREFIX=(
  [BF]=RING1_input_BF
  [HS]=RING1_control_HS
)

# Suffixes stripped from raw FASTQ names to get the sample prefix
RAW_CHIP_SUFFIX=_ChipSeq_biological_replicate_1_and_2_Homo_sapiens_ChIP-Seq
RAW_RNA_SUFFIX=_RNASeq_biological_replicate_1_and_2_Homo_sapiens_RNA-Seq
