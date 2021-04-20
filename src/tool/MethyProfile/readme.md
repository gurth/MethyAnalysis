# MethyProfile: 

---

## What is MethyProfile:

MethyProfile is focus on generating methylation profile based on methylation information (BED file) and genome information (GFF3 file), which is upstream of the methylation analysis process.

---

## Command Line:

usage:

```bash
methyprofile [options] input.bed input.gff3 (output.methyprofile.txt)
```
options:
```
   -P, --promoter=n    			- Analysing promoters methylation information with n bp. The default is 2000bp.
   -l, --single-list=genes.txt	- Gathering single gene information for a gene list behind.
   -h, --help          			- Show this message.
```

---

## Description:

MethyProfile import gene information from GFF3 file and build profile data structure. Then it find and traverse entry in BED file which corresponds to a particular gene to calculate mathylation ratio.

The user may specify the -P option to get promoter mathylation ratio with an option argument n. If there is no argument,  the default length of promoter is 2000bp before a gene sequence. If a gene is near the telomere with the length from the first nucleic acid on this chromosome to this gene sequence's beginnng smaller than 2000bp, which may be impossible, MethyProfile will regard the beginning this promoter as the first nucleic acid on this chromosome.

The -l option is specified for single gene analysis. It load genes in an required argument `genes.txt` with genes' ID line by line.  Mathylation ratio of every  subgene structure of given genes will be gathered and saved in a file. Besides,  entry in BED file which corresponds to the given genes will be saved in `single` directory separately by genes' ID which are the files' names.

---

## Advance:

`config.h`:

You can change macros in `config.h` and rebuild the project to meet the actual needs. For more information, please see code comment in `config.h`.

`plig-in`:

If ALLOW_PLUG_IN_SAVE is defined in `config.h`, this program allows you to save profile manually. The data structure of profile entry is defined in `profile_node.h`. See example of `plug-in/libsave.so/save.c`. 