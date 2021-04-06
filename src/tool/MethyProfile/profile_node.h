//
// Created by gurth on 4/6/21.
//

#ifndef METHYPROFILE_PROFILE_NODE_H
#define METHYPROFILE_PROFILE_NODE_H

struct ProfileNode
    /* Profile list item structure. Show the information of genes (or other
     * items in gff3 file such as exon). */
{
    char ID[0x20] = {0};        /* ID of the gene. */
    int chr = 0;                /* The chromosome where the gene is located. */
    bool chain = -1;            /* Positive chain (true) or negative chain (false). */
    unsigned long Start= 0;
    unsigned long End= 0;       /* Location on chromosome. */
#ifdef CG_NUMBER
    unsigned long NumCG = 0;    /* CG methylation number. */
    unsigned long NumCG_promoter = 0;   /* CG methylation number in promoter. */
#endif// !CG_NUMBER
    double methy_ratio = 0.0f;  /* Methylation ratio of this gene. */
    double methy_ratio_promoter = 0.0f;
    /* Methylation ratio of the promoter of this gene. */
};

#endif //METHYPROFILE_PROFILE_NODE_H
