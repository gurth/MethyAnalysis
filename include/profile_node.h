//
// Created by gurth on 4/6/21.
//

#ifndef METHYPROFILE_PROFILE_NODE_H
#define METHYPROFILE_PROFILE_NODE_H

enum class ProfileNodeType {none, gene, external};

struct BasicEntry
        /* Basic information of an entry.*/
{
    char str_type[ENTRY_TYPE_BUFFER_LENGTH] = {0};  /* Type in an entry. */
    int chr = 0;                /* The chromosome where the gene is located. */
    bool chain = -1;            /* Positive chain (true) or negative chain (false). */
    unsigned long long Start= 0;
    unsigned long long End= 0;       /* Location on chromosome. */
    double methy_ratio = 0.0f;  /* Methylation ratio of this gene. */
};

struct ExternNode : BasicEntry
        /* Extern information node. Show external information of a gene such as
         * exon and so on.*/
{
    ProfileNodeType type = ProfileNodeType::external;   /* Type of this profile node. */
    struct ExternNode* next = nullptr;                  /* Next node on the link list.*/
#ifdef _DEBUG_PROFILE_NODE
    unsigned long long depth = 0;
    unsigned long long mCdep = 0;
#endif // !_DEBUG_PROFILE_NODE
};

struct ProfileNode : BasicEntry
    /* Profile list item structure. Show the information of genes (or other
     * items in gff3 file such as exon). */
{
    ProfileNodeType type = ProfileNodeType::gene;   /* Type of this profile node. */
    char ID[GENE_ID_BUFFER_LENGTH] = {0};           /* ID of the gene. */
    bool single_tag = false;                        /* TRUE when this gene will be analysed singly. */
    unsigned long long Start_promoter= 0;
    unsigned long long End_promoter= 0;
#ifdef TYPE_NUMBER
    unsigned long long NumCG = 0;            /* CG methylation number. */
    unsigned long long NumCG_promoter = 0;   /* CG methylation number in promoter. */
    unsigned long long NumCHG = 0;            /* CG methylation number. */
    unsigned long long NumCHG_promoter = 0;   /* CG methylation number in promoter. */
    unsigned long long NumCHH = 0;            /* CG methylation number. */
    unsigned long long NumCHH_promoter = 0;   /* CG methylation number in promoter. */
#endif// !TYPE_NUMBER
#ifdef _DEBUG_PROFILE_NODE
    unsigned long long depth = 0;
    unsigned long long mCdep = 0;
    unsigned long long depth_promoter = 0;
    unsigned long long mCdep_promoter = 0;
#endif// !_DEBUG_PROFILE_NODE
    double methy_ratio_promoter = 0.0f;
        /* Methylation ratio of the promoter of this gene. */
    struct ExternNode* next = nullptr;  /* Next node on the link list.*/
};

#endif //METHYPROFILE_PROFILE_NODE_H
