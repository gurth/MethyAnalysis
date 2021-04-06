//
// Created by gurth on 3/11/21.
// This file define a class bed which realize the main function of the program.
// Please ensure only one bed progress can be run at a time, otherwise, static
// valuable pThis will be destroyed.
//

#ifndef METHYPROFILE_BED_H
#define METHYPROFILE_BED_H

#include <sys/stat.h>
#include <limits.h>
#include <indicators/block_progress_bar.hpp>
#include <indicators/cursor_control.hpp>

namespace bed
{
    #include "config.h"

    #define MUTEX_LOCK(__CODE, __MUTEX)  pthread_mutex_lock(&__MUTEX); \
        __CODE \
        pthread_mutex_unlock(&__MUTEX);                                                               \
    /* Lock and unlock __CODE in thread function. */

    enum class Method {raw, tag, profile};
    /*  Define the operation method. */

    struct BlockListNode
            /* Block list item structure. Show the block's information. */
    {
        unsigned long base=0;   /* Base offset of the block. */
        unsigned long length=0; /* Length of the block. */
    };

    #include "profile_node.h"

    class BED
    {
    private:
        char bedname[NAME_MAX];         /* Bed file name. */
        char *mapped= nullptr;          /* Bed file mapped location. */
        char *mappedIndex = nullptr;    /* GFF3 file mapped location. */
        unsigned long base_offset=0;    /* Base offset of every loop. */
        unsigned long sum =0;           /* Sum of bytes which have processed. */
        struct stat sb{};               /* Attributes of bed file. */
        struct stat sbIndex{};          /* Attributes of gff3 file. */
        int threadNum = 0;              /* Thread number per loop. */
        int nodeNum = 0;                /* First node number which is pending process. */
        int geneNum = 0;                /* Sum of gene number. */
        double progress = 0.0f;         /* Processing progress. */
        double progUnit=0.0f;           /* Progress unit. */
        indicators::BlockProgressBar bar     /* Progress bar. For details, please refer to
                                        * https://github.com/p-ranav/indicators.*/
        {
            indicators::option::BarWidth{50},
            indicators::option::Start{"["},
            //indicators::option::Fill{"■"},
            //indicators::option::Lead{"■"},
            //indicators::option::Remainder{" "},
            indicators::option::End{"]"},
            indicators::option::PrefixText{"Progress: "},
            indicators::option::ForegroundColor{indicators::Color::green},
            indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}
        };
        int fileHandle = -1;            /* File handle of bed file. */
        int indexHandle = -1;           /* File handle of gff3 file. */
        std::vector<BlockListNode>blockList;    /* Block list of bed file. */
        BlockListNode chrList[MAX_CHR];         /* Chromosome list in the form of blockList. */
        ProfileNode* profileList[MAX_GENE]={nullptr};
                                        /* Profile list, you can see ProfileNode definition above. */
        pthread_mutex_t mutex;          /* Global mutex of thread. */
    private:
        inline void processRaw();
            /* Process bed file with Method::raw, which means byte manipulation directly.
             * Always use to test the performance of this program. */

        inline void processTag();
            /* Process bed file and generate chromosome list.
             * Must run before generation of profile.
             * You can save chromosome list with savechrList().
             * Please notice that if you want to import chromosome list from file directly,
             * make sure the file is named as $(YOUR_BED_FILE_NAME).tag. */

        inline void processProfile(char*& gff3file);
            /* Process bed file and generate profile based on gff3 file.
             * Only support gff3 format.
             * The profile can be saved with saveProfile(). */

        static void dichotomySearchChr(char* m_beg, char* m_end);
            /* Search cut-off point between two block of chromosome data in bed file using dichotomy.
             * m_beg and m_end are the beginning an ending of a block.
             * Please make sure that m_beg and m_end point to the beginning of a line.
             * This function will be executed recursively. */

        static void dichotomySearchOffset(char* m_beg, char* m_end, char*& ppos, unsigned long pos, bool isBeg);
            /* Search interval of a gene in bed file using dichotomy.
             * This function will be executed recursively.
             * m_beg and m_end are the beginning an ending of a block.
             * ppos: The closest entry of bed file which is lower than pos if ifbeg true, or larger if ifbeg false.
             *       Function will make sure ppos points to the beginning of an entry.
             * pos: The point pending search.
             * isBeg: Indicates the beginning of an interval or the ending. */

        static void *pthFuncRaw(void *args);
            /* Thread function.
             * Implementation of processing bed file with Method::raw.
             * This function will process the block which is divided form bed file. */

        static void *pthFuncTag(void *args);
            /* Thread function.
             * Implementation of processing bed file to generate chromosome list.
             * Make sure block list has been build before.
             * The block which processing refers to the block list.
             * */

        static void *pthFuncBlockList(void *args);
            /* Thread function.
             * Implementation of processing bed file to generate block list.
             * This function will make sure a block is begin with the beginning of
             * an entry, and end with the beginning with another entry. */

        static void *pthFuncProfileList(void *args);
            /* Thread function.
             * Implementation of generation of profile list from gff3 file.
             * This function can gather basic information of entry in profile.
             * */

        static void *pthFuncProfile(void *args);
            /* Thread function.
             * Implementation of fulfill the profile list.
             * This function can calculate methylation ratio of every entry based on bed file.
             * */

        static inline double getMethyRatio(char* m_beg, char* m_end, size_t p_start, size_t p_end, bool chain
        #ifdef CG_NUMBER
                    , unsigned long& cg_numb
        #endif // CG_NUMBER
        );
            /* Get methy ratio.
             * m_beg and m_end are pointers that show the block edge.
             * p_start and p_end are the sequence in chromosome.
             * chain shows the positive chain (true) or negative chain (false).
             * This function will return methyratio from m_beg and m_end form p_start to p_end on chromosome.
             * */

        static inline char* goFrontItem(char* p, int n);
            /* Reach forward item in a line with item separated by '\t'.
             * p: The item now;
             * n: Go to the n_th item forward.
             * ex:
             * char* c= "abcd\tefgh\tzzzz\n";
             * p=c;
             * p= goFrontItem(p, 2);   // p = "zzzz\n"  */

        static inline void setValuePfNode(char*&pgoback, char*& pID, ProfileNode*& pPfN);
            /* Fulfill an profile entry. */

        static inline void methyMining(ProfileNode*& pGene);
            /* Calculate methylation ratio of a particular entry of profile list. */

        static int atoiChr(const char *nptr);
            /* Extension of atoi() in stdlib which can process chromosome number like X, Y, MT.*/

        static bool isspace(int x);
            /* Determine whether it is a space. Only used in atoiChr(). */

        static bool isdigit(int x);
            /* Determine whether it is a digit. Only used in atoiChr(). */

#ifdef ALLOW_PLUG_IN_SAVE
        static inline void LoadSavePlugInAndJmp(const char* foutput);
            /* Load plug in save add execute it.*/
#endif //!ALLOW_PLUG_IN_SAVE

    public:
        static BED* pThis;
            /* static pointer points to this class to ensure the thread function can
             * access data in this class. */
        bool have_promoter = false;
            /* Analysing promoters methylation information if true.*/
        size_t promoterLen = PROMOTER_LENGTH;
            /* Search range of promoters. */
    public:
        BED();                  /* Constructor without opening and mapping file.*/
        BED(char *bedfile);     /* Constructor with opening and mapping file.*/
        void bedfileOpen(const char* bedfile);
            /* Open bed file an map it to memory. */

        void bedfileClose();
            /* Close bed file unmap it from memory. */

        void process(char* outputfile, Method m);
            /* Main process of this class.
             * outputfile: Profile will save with outputfile if m is Method::profile.
             *             Default or in other method this value is nullptr.
             * m: Process Method.*/

        void process(char* gff3file, char* outputfile, Method m);
            /* Main process of this class.
             * gff3file: GFF3 file name.
             * outputfile: Profile will save with outputfile if m is Method::profile.
             *             Default or in other method this value is nullptr.
             * m: Process Method.*/

        void savechrList();
            /* Save chromosome list in file as $(YOUR_BED_FILE_NAME).tag.*/

        void saveProfile(const char* nameProfile);
            /* Save profile list in file as nameProfile.
             * If nameProfile == nullptr, profile list will be saved as
             * $(YOUR_BED_FILE_NAME).methyprofile.txt*/

        ~BED(); /* Call bedfileClose() to deconstruct class. */
    };

}
#endif //!METHYPROFILE_BED_H
