//
// Created by gurth on 3/11/21.
// This file define a class bed which realize the main function of the program.
// Please ensure only one bed progress can be run at a time, otherwise, static
// valuable pThis will be destroyed.
//

#ifndef METHYPROFILE_BED_H
#define METHYPROFILE_BED_H

#include <sys/stat.h>
#include <climits>
#include <thread>
#include <mutex>
#include <vector>

#ifdef _UNIX_PLATFORM_

#define BED_MAX_PATH NAME_MAX

#elif defined(_WIN32_PLATFORM_)

# include <Windows.h>
#define BED_MAX_PATH MAX_PATH

#endif //!_UNIX_PLATFORM_

#define MACRO_STR(x) #x
#define XMACRO_STR(s) MACRO_STR(s)

#include <config.h>

#include XMACRO_STR(MY_CONFIG_PATH)

#ifdef INDICATOR_PROGRESS_BAR

#include <indicators/block_progress_bar.hpp>

#endif //! INDICATOR_PROGRESS_BAR

#ifdef ENABLE_LOG

#include <zlog.h>

#endif//!ENABLE_LOG

namespace bed
{
    #define MUTEX_LOCK(__CODE, __MUTEX)  (__MUTEX).lock(); \
        __CODE \
        (__MUTEX).unlock();                                                               \
    /* Lock and unlock __CODE in thread function. */

    enum class Method {raw, tag, profile};
    /*  Define the operation method. */

    typedef void (*p_setProgress)(double);
    typedef void (*p_initProgress)(double, const char*);

    typedef void (*p_errorExit)(int);

    struct BlockListNode
            /* Block list item structure. Show the block's information. */
    {
        unsigned long long base=0;   /* Base offset of the block. */
        unsigned long long length=0; /* Length of the block. */
    };

    struct CHRList : BlockListNode
    {
        unsigned long long edge_base=0;
    };

    #include <profile_node.h>

    class BED
    {
    private:
        char bedname[BED_MAX_PATH] = {0};         /* Bed file name. */
        char *mapped= nullptr;          /* Bed file mapped location. */
        char *mappedIndex = nullptr;    /* GFF3 file mapped location. */
        unsigned long long base_offset=0;    /* Base offset of every loop. */
        unsigned long long sum =0;           /* Sum of bytes which have processed. */
        size_t size_file = 0;           /* Size of bed file. */
        size_t size_fileIndex = 0;      /* Size of gff3 file. */
#ifdef _UNIX_PLATFORM_
        int fileHandle = -1;            /* File handle of bed file. */
        int indexHandle = -1;           /* File handle of gff3 file. */
#elif defined(_WIN32_PLATFORM_)
        HANDLE fileHandle = INVALID_HANDLE_VALUE;       /* File handle of bed file. */
        HANDLE indexHandle = INVALID_HANDLE_VALUE;      /* File handle of gff3 file. */
        HANDLE bedMapHandle = INVALID_HANDLE_VALUE;     /* BED mapped handle. */
        HANDLE gff3MapHandle = INVALID_HANDLE_VALUE;    /* GFF3 mapped handle. */
#endif //!_UNIX_PLATFORM_
        int threadNum = 0;              /* Thread number per loop. */
        int nodeNum = 0;                /* First node number which is pending process. */
        int geneNum = 0;                /* Sum of gene number. */
        double progress = 0.0f;         /* Processing progress. */
        double progUnit=0.0f;           /* Progress unit. */
#ifdef INDICATOR_PROGRESS_BAR
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
#endif // !INDICATOR_PROGRESS_BAR
#ifdef ENABLE_LOG
        zlog_category_t *zc = nullptr;  /* zlog category, for log file output*/
#endif //!ENABLE_LOG
        std::vector<BlockListNode>blockList;    /* Block list of bed file. */
        CHRList chrList[MAX_CHR];         /* Chromosome list in the form of blockList. */
        ProfileNode* profileList[MAX_GENE]={nullptr};
                                        /* Profile list, you can see ProfileNode definition above. */
        std::mutex mtx;          /* Global mutex of thread. */
    private:
        void init();
            /* General initialization operation. */

        void processInit();
            /* Inittation when process start.*/

        char* open_map(const char* filename, size_t& length,
#ifdef _UNIX_PLATFORM_
        int& m_handle
#elif defined(_WIN32_PLATFORM_)
        HANDLE& m_handle, HANDLE& m_handleMap
#endif //!_UNIX_PLATFORM_
        ) const;
            /* Open and map file in memory.  */

        void unmap_close(size_t& length, char* p_m,
#ifdef _UNIX_PLATFORM_
        int& m_handle
#elif defined(_WIN32_PLATFORM_)
        HANDLE& m_handle, HANDLE& m_handleMap
#endif //!_UNIX_PLATFORM_
                         ) const;
            /* Unmap and close file*/

        inline void processRaw();
            /* Process bed file with Method::raw, which means byte manipulation directly.
             * Always use to test the performance of this program. */

        inline void processTag();
            /* Process bed file and generate chromosome list.
             * Must run before generation of profile.
             * You can save chromosome list with savechrList().
             * Please notice that if you want to import chromosome list from file directly,
             * make sure the file is named as $(YOUR_BED_FILE_NAME).tag. */

        inline void processProfile(const char*& gff3file);
            /* Process bed file and generate profile based on gff3 file.
             * Only support gff3 format.
             * The profile can be saved with saveProfile(). */

        static void dichotomySearchChr(char* m_beg, char* m_end);
            /* Search cut-off point between two block of chromosome data in bed file using dichotomy.
             * m_beg and m_end are the beginning an ending of a block.
             * Please make sure that m_beg and m_end point to the beginning of a line.
             * This function will be executed recursively. */

        static void dichotomySearchChain(int m_chr, char* m_beg, char* m_end);

        static void dichotomySearchOffset(char* m_beg, char* m_end, char*& ppos, unsigned long long pos, bool isBeg);
            /* Search interval of a gene in bed file using dichotomy.
             * This function will be executed recursively.
             * m_beg and m_end are the beginning an ending of a block.
             * ppos: The closest entry of bed file which is lower than pos if ifbeg true, or larger if ifbeg false.
             *       Function will make sure ppos points to the beginning of an entry.
             * pos: The point pending search.
             * isBeg: Indicates the beginning of an interval or the ending. */

        static void pthFuncRaw();
            /* Thread function.
             * Implementation of processing bed file with Method::raw.
             * This function will process the block which is divided form bed file. */

        static void pthFuncTag();
            /* Thread function.
             * Implementation of processing bed file to generate chromosome list.
             * Make sure block list has been build before.
             * The block which processing refers to the block list.
             * */

        static void pthFuncBlockList();
            /* Thread function.
             * Implementation of processing bed file to generate block list.
             * This function will make sure a block is begin with the beginning of
             * an entry, and end with the beginning with another entry. */

        static void pthFuncProfileList();
            /* Thread function.
             * Implementation of generation of profile list from gff3 file.
             * This function can gather basic information of entry in profile.
             * */

        static void pthFuncProfile();
            /* Thread function.
             * Implementation of fulfill the profile list.
             * This function can calculate methylation ratio of every entry based on bed file.
             * */

        static inline double getMethyRatio(char* m_beg, char* m_end, size_t p_start, size_t p_end, char* ID,
                                           bool single_tag, bool ispromoter, bool chain
        #ifdef TYPE_NUMBER
                    , unsigned long long& cg_numb
                    , unsigned long long& chg_numb
                    , unsigned long long& chh_numb
        #endif // TYPE_NUMBER
        #ifdef _DEBUG_PROFILE_NODE
                    ,unsigned long long& m_depth
                    ,unsigned long long& m_mCdep
        #endif // !_DEBUG_PROFILE_NODE
        );
            /* Get methy ratio.
             * m_beg and m_end are pointers that show the block edge.
             * p_start and p_end are the sequence in chromosome.
             * chain shows the positive chain (true) or negative chain (false).
             * ID is the ID of this entry.
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

        static inline void setValueBasic(char*& p, BasicEntry*& pEntry);
            /* Fulfill basic information in an entry. */

        static inline void setValuePfNode(char*& pgoback, char*& pID, ProfileNode*& pPfN);
            /* Fulfill an profile entry. */

        static inline void setValueExtern(char*& p, ExternNode*& pEntry);
            /* Fulfill basic information in an entry. */

        static inline void copyItem(char*& d, char*& s);
            /* Copy a item as string to destination. */

        static inline char* goNextEntry(char* p);
            /* Get to the beginning of next entry. */

        static inline void methyMining(ProfileNode*& pGene);
            /* Calculate methylation ratio of a particular entry of profile list. */

        static int atoiChr(const char *nptr);
            /* Extension of atoi() in stdlib which can process chromosome number like X, Y, MT.*/

        static bool isspace(int x);
            /* Determine whether it is a space. Only used in atoiChr(). */

        static bool isdigit(int x);
            /* Determine whether it is a digit. Only used in atoiChr(). */

        static void saveSingleData(char* m_beg, char* m_end, char* ID);
            /* Save single gene information as BED format. */

        static void saveSingleData(char* m_beg, char* m_end, char* ID, const char* suffix);
            /* Save single gene information as BED format. */

        static void profileNodeDump(const char* name_dump, int n, ProfileNode** pnlist);
            /* Save raw data of profile node. */

        static void m_setProgress(double m_progress);
            /* Setting progress of progress bar.*/

        static void m_initProgress(double m_progress, const char* info);
            /* Initiation of progress bar. */

        static void m_errorExit(int m_error);
            /* Exit program with error code m_error. */

#ifdef ALLOW_PLUG_IN_SAVE
        static inline void LoadSavePlugInAndJmp(const char* foutput);
            /* Load plug in save add execute it.*/
#endif //!ALLOW_PLUG_IN_SAVE

        bool LoadTag();
            /* Load chromosome list from file. */

        inline void TagChain();

    public:
        static BED* pThis;
            /* static pointer points to this class to ensure the thread function can
             * access data in this class. */
        bool have_promoter = false;
            /* Analysing promoters methylation information if true.*/
        bool do_single_analyse = false;
            /* Analysing single gene information. */    
        size_t promoterLen = PROMOTER_LENGTH;
            /* Search range of promoters. */
        void* sglist = nullptr;
            /* Single gene list which is pending analysis.*/  
        p_initProgress initProgress = nullptr;
            /* Pointer to re-call function for initiation of progress bar. */
        p_setProgress setProgress = nullptr; 
            /* Pointer to re-call function for setting progress of progress bar. */
        p_errorExit errorExit = nullptr;
            /* Pointer to re-call function for exiting with error code. */
        double latest_time_cost = 0;
            /* Time cost of latest operation. */
    public:
        BED();                  /* Constructor without opening and mapping file.*/
        explicit BED(char *bedfile);     /* Constructor with opening and mapping file.*/
        explicit BED(p_errorExit m_errorEx); /* Constructor with outside error parsing. */
        void bedfileOpen(const char* bedfile);
            /* Open bed file an map it to memory. */

        void bedfileClose();
            /* Close bed file unmap it from memory. */

        void process(const char* outputfile, Method m);
            /* Main process of this class.
             * outputfile: Profile will save with outputfile if m is Method::profile.
             *             Default or in other method this value is nullptr.
             * m: Process Method.*/

        void process(const char* gff3file, const char* outputfile, Method m);
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

        void saveExternProfile(const char* nameProfileEx);
            /* Save extern list of specific genes.
             * If nameProfileEx == nullptr, profile list will be saved as
             * $(YOUR_BED_FILE_NAME).gene.txt */

        void loadSingleList(const char* listfile);
            /* Load single list from file.
             * This function will load information from file and init relative data structure.
             * */
            
        ~BED(); /* Call bedfileClose() to deconstruct class. */
    };

}
#endif //!METHYPROFILE_BED_H
