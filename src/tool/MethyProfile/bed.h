//
// Created by gurth on 3/11/21.
//

#ifndef UNTITLED_BED_H
#define UNTITLED_BED_H

#include <sys/stat.h>
#include <limits.h>
#include <indicators/progress_bar.hpp>
#include <indicators/cursor_control.hpp>

namespace bed
{

    //#define ALLOW_LOWERCASE

    #define MAXTHREAD 16
    #define SEARCH_RANGE 512
    #define SEARCH_RANGE_INDEX 102400  // 100 K

    #define MAX_GENE 65536
    #define MAX_CHR 128

    #define PROC_GENE_SIZE 0x80
    #define PROC_GENE_READ 0x800

    //#define _FLAG_TEST

    #ifdef _FLAG_TEST

        #define TEST_CHAR 'G'
        #define BLOCK_SIZE 0x4000000 // 64 M
        #define BLOCK_READ  0x40000000 // 1 G

    #else

        #define BLOCK_SIZE 0x400000 // 4 M
        #define BLOCK_READ  0x4000000 // 64 M

    #endif //!_FLAG_TEST

    #define BLOCK_SIZE_INDEX 0x80000 // 512 K
    #define BLOCK_READ_INDEX 0x800000 // 4 M

    //#define SHOW_ALL_INFO
    #define SHOW_PROGRESSBAR

    enum class Method {raw, tag, profile};

    struct BlockListNode
    {
        unsigned long base=0;
        unsigned long length=0;
    };

    struct ProfileNode
    {
        char ID[0x20] = {0};
        int chr = 0;
        bool chain = -1;
        unsigned long Start= 0;
        unsigned long End= 0;
        double methy_ratio = 0.0f;
    };

    class BED
    {
    private:
        char bedname[NAME_MAX];
        char *mapped= nullptr;
        char *mappedIndex = nullptr;
        unsigned long base_offset=0;
        unsigned long sum =0;
        struct stat sb{};
        struct stat sbIndex{};
        int threadNum = 0;
        int nodeNum = 0;
        int geneNum = 0;
        double progress = 0.0f;
        double progUnit=0.0f;
        indicators::ProgressBar bar
        {
            indicators::option::BarWidth{50},
            indicators::option::Start{"["},
            indicators::option::Fill{"■"},
            indicators::option::Lead{"■"},
            indicators::option::Remainder{" "},
            indicators::option::End{"]"},
            indicators::option::PrefixText{"Progress: "},
            indicators::option::ForegroundColor{indicators::Color::green},
            indicators::option::FontStyles{std::vector<indicators::FontStyle>{indicators::FontStyle::bold}}
        };
        int fileHandle = -1;
        int indexHandle = -1;
        std::vector<BlockListNode>blockList;
        std::vector<BlockListNode>blockListIndex;
        BlockListNode chrList[MAX_CHR];
        ProfileNode* profileList[MAX_GENE]={nullptr};
        pthread_mutex_t mutex;
    private:
        inline void processRaw();
        inline void processTag();
        inline void processProfile(char*& gff3file);
        static void dichotomySearchChr(char* m_beg, char* m_end);
        static void dichotomySearchOffset(char* m_beg, char* m_end,char*& ppos, unsigned long pos, bool isBeg);
        static void *pthFuncRaw(void *args);
        static void *pthFuncTag(void *args);
        static void *pthFuncBlockList(void *args);
        static void *pthFuncProfileList(void *args);
        static void *pthFuncProfile(void *args);
        static inline char* goFrontItem(char* p, int n);
        static inline void setValuePfNode(char*&pgoback, char*& pID, ProfileNode*& pPfN);
        static inline void methyMining(ProfileNode*& pGene);
        static int atoiChr(const char *nptr);
        static bool isspace(int x);
        static bool isdigit(int x);
    public:
        static BED* pThis;
    public:
        BED();
        BED(char *bedfile);
        void bedfileOpen(char*& bedfile);
        void bedfileClose();
        void process(char* outputfile, Method m);
        void process(char* gff3file, char* outputfile, Method m);
        void savechrList();
        void saveProfile(const char* nameProfile);
        ~BED();
    };

}
#endif //!UNTITLED_BED_H
