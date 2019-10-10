# GenomeAnnotation

### Purpose

This repository contains various scripts that I have created that are useful for genome annotation (repeats & proteins). They may exist in some form already elsewhere, but they are here in case others need them. In particular, they are used in a [MAKER tutorial](https://github.com/darencard/darencard.github.io/blob/master/_blog/2017-05-16-maker-genome-annotation.md) I created, which I have heard others reference a lot.

### Brief Descriptions

I try to provide more details and usage information at the beginning of each script. Here is a quick description of what each script does for reference. User beware: these scripts have not been extensively tested and may not perform as desired, and as such, they are provided as-is, with no support and no guarantee of proper or desirable functioning.

`genestats`: calculates several count/length statistics for gene models in a GFF3 file

`ncbi2uniprot`: renames the header lines of an NCBI protein FASTA file so that they are very similar to the standard UniProt FASTA headers

`ncbi_CDS_creator.sh`: creates a CDS file from a NCBI genome using the feature table

`orthorbb`: infers homology between a query and a reference set of protein sequences using reciprocal and one-way best BLAST searches

`rmOutToGFF3custom`: converts the .out file from RepeatMasker to a "better" GFF3 file than default script in RepeatMasker
