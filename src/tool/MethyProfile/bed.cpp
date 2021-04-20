//
// Created by gurth on 3/11/21.
//

#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <thread>
#include <mutex>
#include <unistd.h>
#include <time.h>
#include <vector>
#include <string>
#include <set>
#include <fstream>
#include <algorithm>
#include "bed.h"

#ifdef ALLOW_PLUG_IN_SAVE
#include <dlfcn.h>
#endif //!ALLOW_PLUG_IN_SAVE

#ifdef _WIN32_PLATFORM_

#include <direct.h>
#include <io.h>

#elif defined(_UNIX_PLATFORM_)

#include <sys/stat.h>

#endif // !_WIN32_PLATFORM_

using namespace std;
using namespace bed;
using namespace indicators;

bool cmpChrList(BlockListNode* pa, BlockListNode* pb)
{
    return pa->base < pb->base;
}

bool cmpProfileList(ProfileNode* pa, ProfileNode* pb)
{
    if(pa->chr==pb->chr)
        return pa->Start < pb->Start;
    return pa->chr < pb->chr;
}

BED* BED::pThis = nullptr;

void BED::dichotomySearchChr(char* m_beg, char* m_end)
{
    int chrBeg=atoiChr(m_beg+3);
    int chrEnd=atoiChr(m_end+3);

    if(chrBeg!=chrEnd)
    {
        char* p=(char*)(size_t)(((size_t)m_beg + (size_t)m_end) >> 1);

        for(int j=0; j<SEARCH_RANGE;j++, p--)
            if(*(p)=='\n') break;
        p++;

        if(p>m_beg)  // if p == m_beg, dichotomySearch(p,m_end) will forever unchanged
        {
            dichotomySearchChr(m_beg, p);
            dichotomySearchChr(p,m_end);
        }
        else if(p==m_beg)
        {
            for(int j=0; j<SEARCH_RANGE;j++, p++)
                if(*(p)=='\n') break;
            p++;
            if(p < m_end)
            {
                dichotomySearchChr(m_beg, p);
                dichotomySearchChr(p,m_end);
            }
            else if(p == m_end)
                dichotomySearchChr(p,m_end);
        }
    }

    if(pThis->chrList[chrBeg].base > m_beg - pThis->mapped)
        pThis->chrList[chrBeg].base = m_beg - pThis->mapped;
}

void BED::dichotomySearchOffset(char *m_beg, char *m_end, char*& ppos, unsigned long pos, bool isBeg)
{
    int pos_beg=atoiChr(goFrontItem(m_beg, 2));
    int pos_end=atoiChr(goFrontItem(m_end, 1));
    char* p=(char*)(size_t)(((size_t)m_beg + (size_t)m_end) >> 1);
    for(int j=0; j<SEARCH_RANGE;j++, p--)
        if(*(p)=='\n') break;
    p++;
    int pos_mid_beg=atoiChr(goFrontItem(p, 1));
    int pos_mid_end=atoiChr(goFrontItem(p, 2));
    if(p>m_beg)
    {
        if(pos_beg <= pos && pos <= pos_mid_beg)
            dichotomySearchOffset(m_beg, p, ppos, pos, isBeg);
        if(pos_mid_end <= pos && pos <= pos_end)
            dichotomySearchOffset(p, m_end, ppos, pos, isBeg);
    }
    else if(p==m_beg)
    {
        for(int j=0; j<SEARCH_RANGE;j++, p++)
            if(*(p)=='\n') break;
        p++;
        if(p < m_end)
        {
            pos_mid_beg=atoiChr(goFrontItem(p, 1));
            pos_mid_end=atoiChr(goFrontItem(p, 2));
            if(pos_beg <= pos && pos <= pos_mid_beg)
                dichotomySearchOffset(m_beg, p, ppos, pos, isBeg);
            if(pos_mid_end <= pos && pos <= pos_end)
                dichotomySearchOffset(p, m_end, ppos, pos, isBeg);
        }
        else if(p == m_end)
        {
            if(isBeg)
                ppos = m_end;
            else
                ppos = m_beg;
        }
    }
}

void BED::methyMining(ProfileNode *& pGene)
{
    if(!(size_t)(pThis->chrList[pGene->chr].base+1)) return;
    // We can think about using a queue to shrink search range
    /* Prepare search range. */
    char* m_beg=pThis->mapped+pThis->chrList[pGene->chr].base;
    char* m_end=m_beg+pThis->chrList[pGene->chr].length;
    double dtemp=0.0f;

    /* Kill empty line before EOF. */
    m_end--;
    while(true){if((*m_end)!='\n') break; m_end--;}
    for(int j=0; j<SEARCH_RANGE;j++, m_end--)
        if(*(m_end)=='\n') break;
    m_end++;

    dtemp = getMethyRatio(m_beg, m_end, pGene->Start, pGene->End, pGene->ID, pGene->single_tag, pGene->chain
#ifdef CG_NUMBER
            ,pGene->NumCG
#endif //!CG_NUMBER
    );
    if(dtemp == -1)
    {
        pGene->methy_ratio = 0.0f;
        zlog_notice(pThis->zc,"getMethyRatio() return -1 when processing gene %s.", pGene->ID);
    }
    else pGene->methy_ratio=dtemp;

    if(pThis->have_promoter)
    {
        dtemp = getMethyRatio(m_beg, m_end
                              , (pGene->Start < pThis->promoterLen) ? 0 : pGene -> Start - pThis -> promoterLen
                              , pGene->Start, pGene->ID, false, pGene->chain
                    #ifdef CG_NUMBER
                            ,pGene->NumCG_promoter
                    #endif //!CG_NUMBER
                              );
        if(dtemp == -1)
        {
            pGene->methy_ratio_promoter = 0.0f;
            zlog_notice(pThis->zc,"getMethyRatio() return -1 when processing promoter of gene %s.", pGene->ID);
        }
        else pGene->methy_ratio_promoter=dtemp;
    }

    if(pThis->do_single_analyse && pGene->single_tag)
    {
        ExternNode* pEx = pGene->next;
        while (true)
        {
            if(pEx== nullptr) break;
#ifdef CG_NUMBER
            size_t CG_tmp = 0;
#endif //!CG_NUMBER
            dtemp = getMethyRatio(m_beg, m_end, pEx->Start, pEx->End, nullptr, false, pEx->chain
#ifdef CG_NUMBER
                    ,CG_tmp
#endif //!CG_NUMBER
            );
            if(dtemp == -1)
            {
                pEx->methy_ratio=0.0f;
                zlog_notice(pThis->zc,"getMethyRatio() return -1 when processing %s form %ld to %ld of gene %s.",
                            pEx->str_type, pEx->Start, pEx->End, pGene->ID);
            }
            else pEx->methy_ratio=dtemp;
            pEx=pEx->next;
        }
    }
}

double BED::getMethyRatio(char *m_beg, char *m_end, size_t p_start, size_t p_end, char* ID, bool single_tag, bool chain
#ifdef CG_NUMBER
, unsigned long& cg_numb
#endif // CG_NUMBER
)
{
    char* pGbeg= nullptr, *pGend= nullptr;
    char* p=nullptr;

    cg_numb=0;

    /* Search position. */
    dichotomySearchOffset(m_beg,m_end,pGbeg,p_start,true);
    dichotomySearchOffset(m_beg,m_end,pGend,p_end,false);

    if(!pGbeg && pGend)
        pGbeg=m_beg;
    else if(pGbeg && !pGend)
        pGend=m_end;
    else if(!pGbeg && !pGend)
        return -1;

    if(pThis->do_single_analyse && single_tag) saveSingleData(pGbeg, pGend, ID);

    /* Calculate methylation ratio. */
    size_t depth = 0, mCdep = 0;
    p=pGbeg;

    while(true)
    {
        if(p>pGend) break;
        p=goFrontItem(p,3);
        if((*p == '+') == (chain))
        {
#ifdef CG_NUMBER
            p=goFrontItem(p, 1);
            if((*p)=='C' && (*(p+1))=='G') cg_numb++;
            p=goFrontItem(p, 1);
#else
            p=goFrontItem(p, 2);
#endif //!CG_NUMBER

            depth+=atoi(p);
            p=goFrontItem(p,1);
            mCdep+=atoi(p);
        }
        for(int j=0;j<SEARCH_RANGE;j++, p++)
            if((*p)=='\n') break;
        p++;
    }
    if(depth!=0)
        return (double)((double)mCdep / (double)depth);
    if(mCdep) return -1;
    return 0;
}

void BED::pthFuncRaw()
{
    MUTEX_LOCK(
            size_t i = pThis->threadNum++;
            ,
            pThis->mtx
    )
    size_t begOffset = pThis->base_offset + i*BLOCK_SIZE;

#ifdef _FLAG_TEST
    if(pThis->sb.st_size - begOffset >= BLOCK_SIZE)
    {
        memset(pThis->mapped + begOffset, TEST_CHAR, BLOCK_SIZE);
        pThis->sum=begOffset+BLOCK_SIZE;
    #ifdef SHOW_ALL_INFO
        printf("Base: 0x%016lx Size: 0x%016lx\n",begOffset, BLOCK_SIZE);
    #endif //!SHOW_ALL_INFO
    }
    else
    {
        memset(pThis->mapped + begOffset, TEST_CHAR,pThis->sb.st_size - begOffset);
        pThis->sum=pThis->sb.st_size;
    #ifdef SHOW_ALL_INFO
        printf("Base: 0x%016lx (Last block)\n",begOffset);
    #endif //!SHOW_ALL_INFO
    }
#endif //!_FLAG_TEST

#ifdef SHOW_PROGRESSBAR
    char buff[0x30];
    sprintf(buff,"Processed %ld / %ld",pThis->sum,pThis->sb.st_size);
    MUTEX_LOCK(
        pThis->progress+=pThis->progUnit;
        pThis->bar.set_option(indicators::option::PostfixText{buff});
        pThis->bar.set_progress(pThis->progress);
        ,
        pThis->mtx
    )
#endif //!SHOW_PROGRESSBAR
}

void BED::pthFuncTag()
{
    MUTEX_LOCK(
        size_t k= pThis->nodeNum++;
        ,
        pThis->mtx
    )
    size_t begOffset = pThis->blockList[k].base;
    char* q=&(pThis->mapped[begOffset]);
    char* p=&(pThis->mapped[begOffset+pThis->blockList[k].length -1]);

    for(int j=0; j<SEARCH_RANGE;j++, p--)
        if(*(p)=='\n') break;
    p++;

    dichotomySearchChr(q, p);
#ifdef SHOW_PROGRESSBAR
    MUTEX_LOCK(
        pThis->progress+=pThis->progUnit;
        pThis->bar.set_progress(pThis->progress);
        ,
        pThis->mtx
    )
#endif //!SHOW_PROGRESSBAR
}

void BED::pthFuncBlockList()
{
    MUTEX_LOCK(
            size_t i = pThis->threadNum++;
            size_t k= pThis->nodeNum++;
            ,
            pThis->mtx
    )
    size_t begOffset = pThis->base_offset + i*BLOCK_SIZE;
    int j=0;

    if(k)
    {
        for (; j < SEARCH_RANGE; j++)
            if (pThis->mapped[begOffset + j] == '\n') break;
        begOffset = begOffset + j + 1;
    }
    pThis->blockList[k].base=begOffset;
#ifdef SHOW_PROGRESSBAR
    MUTEX_LOCK(
            pThis->progress+=pThis->progUnit;
            pThis->bar.set_progress(pThis->progress);
            ,
            pThis->mtx
        )
#endif //!SHOW_PROGRESSBAR
}

void BED::pthFuncProfile()
{
    MUTEX_LOCK(
            size_t i = pThis->threadNum++;
            ,
            pThis->mtx
    )
    size_t begOffset = pThis->base_offset + i*PROC_GENE_SIZE;
    size_t length = PROC_GENE_SIZE;
    if(pThis->geneNum - begOffset < PROC_GENE_SIZE)
        length = pThis->geneNum - begOffset;

    for(size_t j=0;j<length;j++)
    {
        methyMining(pThis->profileList[begOffset+j]);
#ifdef SHOW_PROGRESSBAR
        MUTEX_LOCK(
            pThis->progress+=pThis->progUnit;
            pThis->bar.set_progress(pThis->progress);
            ,
            pThis->mtx
        )
#endif //!SHOW_PROGRESSBAR
    }
}

void BED::pthFuncProfileList()
{
    MUTEX_LOCK(
        size_t i = pThis->threadNum++;
        ,
        pThis->mtx
    )
    size_t begOffset = pThis->base_offset + i*BLOCK_SIZE_INDEX;
    size_t length = BLOCK_SIZE_INDEX;
    if(pThis->sbIndex.st_size - begOffset < BLOCK_SIZE_INDEX)
        length=pThis->sbIndex.st_size - begOffset;

    char* p=pThis->mappedIndex+begOffset;
    char* pend=p+length;

    for(; p<pend; p++)
    {
        if((*p) == '#')   // ### tag a item
        {
            if(*(p+1) == '#' && *(p+2) == '#')
            {
                p=p+3;
                char* q=p+SEARCH_RANGE_INDEX;
                if(q>pend) q=pend;
                for(; p<q; p++)
                {
                    if((*p)=='I')          // Find "ID=gene:" which means a gene
                    {
                        if(memcmp(p, "ID=gene:", 8) == 0)
                        {
                            char* pgoback= p;
                            p+=8;
                            for(int j=0;j<SEARCH_RANGE;j++, pgoback--)
                                if((*pgoback)=='\n') break;
                            pgoback++;

                            auto* pfNtmp=(ProfileNode*) malloc(sizeof(ProfileNode));
                            setValuePfNode(pgoback,p,pfNtmp);

                            if(pThis->do_single_analyse)
                            {
                                auto* pSet = (set<string>*)pThis->sglist;
                                string id_now = pfNtmp->ID;

                                MUTEX_LOCK(
                                        auto iter = pSet->find(id_now);
                                        ,
                                        pThis->mtx
                                )

                                if(iter != pSet->end())
                                {
                                    ExternNode* pExRoot = nullptr;
                                    ExternNode* pExN = pExRoot;
                                    ExternNode* pExBuff = nullptr;
                                    pfNtmp->single_tag=true;
                                    p = goNextEntry(p);

                                    while(true)
                                    {
                                        if(*p=='#')
                                        {
                                            if(*(p+1) == '#' && *(p+2) == '#')
                                                break;
                                            else goto NEXT__;
                                        }
                                        if(p>pend) break;
                                        pExBuff=(ExternNode*)malloc(sizeof(ExternNode));
                                        setValueExtern(p, pExBuff);
                                        if(!pExRoot)
                                        {
                                            pExRoot=pExBuff;
                                            pExN=pExRoot;
                                        }
                                        else
                                        {
                                            pExN->next = pExBuff;
                                            pExN = pExN->next;
                                        }
                                    NEXT__:
                                        p = goNextEntry(p);
                                    } //! if(iter != pSet->end())

                                    p-=GO_BACK_LENGTH;
                                    pfNtmp->next=pExRoot;
                                }
                            } //! if(iter != pSet->end())

                            MUTEX_LOCK(
                                    pThis->profileList[pThis->geneNum ++]=pfNtmp;
                                    ,
                                    pThis->mtx
                                    )
                            break;
                        }
                    }//!if((*p)=='I')
                }//! for
            }//! if(*(p+1) == '#' && *(p+2) == '#')
        } //! if((*p) == '#')
    }//! for
#ifdef SHOW_PROGRESSBAR
    MUTEX_LOCK(
            pThis->progress+=pThis->progUnit;
            pThis->bar.set_progress(pThis->progress);
            ,
            pThis->mtx
            )
#endif //!SHOW_PROGRESSBAR
}

void BED::copyItem(char *&d, char *&s)
{
    char* p=s;
    for(int j=0;j<SEARCH_RANGE;j++, p++)
        if((*p)=='\t') break;
    memcpy(d, s, p-s);
    d[(p-s)]='\0';
}

char *BED::goNextEntry(char *p)
{
    for(int j=0;j<SEARCH_RANGE;j++, p++)
        if((*p)=='\n') break;
    p++;
    return p;
}

void BED::setValueBasic(char *&p, BasicEntry *&pEntry)
{
    pEntry->chr=atoiChr(p);
    p=goFrontItem(p, 2);
    if(!p) exit(SET_VALUE_ERROR + 1);
    char* ptmp=pEntry->str_type;
    copyItem(ptmp, p);
    p=goFrontItem(p, 1);
    if(!p) exit(SET_VALUE_ERROR + 2);
    pEntry->Start=atoll(p);
    p=goFrontItem(p, 1);
    if(!p) exit(SET_VALUE_ERROR + 3);
    pEntry->End=atoll(p);
    p=goFrontItem(p, 2);
    if(!p) exit(SET_VALUE_ERROR + 4);
    pEntry->chain=((*p)=='+');
}

void BED::setValuePfNode(char *& pgoback, char *& pID, ProfileNode *& pPfN)
{
    auto *ptmp = (BasicEntry *)pPfN;
    setValueBasic(pgoback, ptmp);

    pgoback=pID;
    for(int j=0;j<SEARCH_RANGE;j++, pID++)
        if((*pID)==';') break;
    memcpy(pPfN->ID, pgoback, (pID-pgoback));
    pPfN->ID[(pID-pgoback)]='\0';
    pPfN->next= nullptr;
    pPfN->single_tag=false;
}

void BED::setValueExtern(char *& p, ExternNode*& pEntry)
{
    auto *ptmp = (BasicEntry *)pEntry;
    setValueBasic(p, ptmp);
}

char *BED::goFrontItem(char *p, int n)
{
    int i=0;
    for(int j=0;j<SEARCH_RANGE;j++, p++)
    {
        if((*p) == 0x9 || (*p) == ' ')
        {
            for(int k=0;k<SEARCH_RANGE;k++, p++)
                if((*p) != 0x9 && (*p) != ' ') break;
            p--;
            i++;
            if(i>=n) return p+1;
        }
    }
    return nullptr;
}

void BED::init()
{
    pThis= this;
    for(auto & i : chrList)
        i={(unsigned long)-1,0};
    string zlogConf = XMACRO_STR(CMAKE_SOURCE_DIR);
    zlogConf+="/etc/zlog_methyprofile.conf";
    if(zlog_init(zlogConf.c_str()))
    {
        printf("\033[31m[Error]\033[0m:Zlog init failed.\n");
        exit(ZLOG_ERROR);
    }
    zc = zlog_get_category("MethyProfile");
}

BED::BED()
{
    init();
}

BED::BED(char* bedfile)
{
    init();
    if(bedfile== nullptr)
        strcpy(bedname,"test.bed");
    else
        strcpy(bedname,bedfile);
    bedfileOpen(bedname);
}

void BED::bedfileOpen(const char * bedfile)
{
    // Open file
    if((fileHandle = open(bedfile, O_RDWR)) < 0)
    {
        perror("open(): ") ;
        exit(FILE_OPEN_ERROR+1);
    }

    // Get file stat
    if((fstat(fileHandle, &sb)) == -1 )
    {
        perror("fstat(): ") ;
        exit(FILE_OPEN_ERROR+2);
    }

    strcpy(this->bedname, bedfile);

    // Map file in memory
    mapped = (char*)mmap(nullptr, sb.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, fileHandle, 0);
    if(mapped == (char*)-1)
    {
        perror("mmap(): ") ;
        exit(FILE_OPEN_ERROR+3);
    }
}

void BED::bedfileClose()
{
    munmap(mapped,sb.st_size);
    close(fileHandle);
}

BED::~BED()
{
    if(do_single_analyse)
        delete(set<string>*)sglist;
    bedfileClose();
    zlog_fini();
}

void BED::process(char *outputfile, Method m)
{
    process(nullptr, outputfile, m);
}

void BED::process(char *gff3file, char *outputfile, Method m)
{
    clock_t t;
    t=clock();
    string _outputfile;

    switch (m)
    {
        case Method::raw:
            processRaw();
            break;
        case Method::tag:
            processTag();
            break;
        case Method::profile:
            if(!outputfile)
                _outputfile=string(bedname)+string(".methyprofile.txt");
            else
                _outputfile=outputfile;
            processProfile(gff3file);
            saveProfile(_outputfile.c_str());
            _outputfile=string(bedname)+string(".gene.txt");
            saveExternProfile(_outputfile.c_str());
            break;
        default:
            break;
    }

    printf("\n");
    t = clock() - t;
    printf("Total cost: %lf ms\n", ((float)t) / CLOCKS_PER_SEC * 1000);
}

void BED::processRaw()
{
    size_t readNum=0;
    size_t restSize=0;
    size_t blockNum=0;
    thread* pth[MAXTHREAD];
    char buff[0x30];

    printf("\033[33m[Warning]\033[0m: Start reading ...\n");

    readNum=sb.st_size / BLOCK_READ;
    restSize=sb.st_size % BLOCK_READ;
#ifdef SHOW_PROGRESSBAR
    progUnit=(double)(100.0 / (sb.st_size / BLOCK_SIZE +1));
    sprintf(buff,"Processed 0 / %ld",sb.st_size);
    bar.set_option(indicators::option::PostfixText{buff});
    bar.set_progress(0);
#endif //!SHOW_PROGRESSBAR
    base_offset=0;

    for(int j=0; j<=readNum; j++, base_offset+=BLOCK_READ)
    {
        threadNum = 0;
        if(j==readNum)
            blockNum=(restSize / BLOCK_SIZE) +1;
        else
            blockNum=MAXTHREAD;

        // Create thread
        for (int i = 0; i < blockNum; i++)
            pth[i]=new thread(this->pthFuncRaw);

        // Join thread
        for (int i = 0; i < blockNum; i++)
        {
            pth[i]->join();
            delete pth[i];
        }
    }
}

void BED::processTag()
{
    size_t readNum=0;
    size_t restSize=0;
    size_t blockNum=0;
    thread* pth[MAXTHREAD];
    char buff[0x30];
    char* p= nullptr;

    printf("\033[33m[Warning]\033[0m: Start tagging ...\n");
#ifdef SHOW_PROGRESSBAR
    sprintf(buff,"Building block list");
    bar.set_option(indicators::option::PostfixText{buff});
    bar.set_progress(0);
#endif //!SHOW_PROGRESSBAR
    readNum=sb.st_size / BLOCK_READ;
    restSize=sb.st_size % BLOCK_READ;
    blockList.resize(sb.st_size / BLOCK_SIZE +1);

    base_offset=0;
    nodeNum = 0;
    progUnit=(double)(25.0 / blockList.size());

    for(int j=0; j<=readNum; j++, base_offset+=BLOCK_READ)
    {
        threadNum = 0;
        if(j==readNum)
            blockNum=(restSize / BLOCK_SIZE) +1;
        else
            blockNum=MAXTHREAD;

        // Create thread
        for (int i = 0; i < blockNum; i++)
            pth[i]=new thread(pthFuncBlockList);

        // Join thread
        for (int i = 0; i < blockNum; i++)
        {
            pth[i]->join();
            delete pth[i];
        }

    }

    p=mapped;
    while (true)
    {
        if(*(p+3) >= '0' && *(p+3) <= '9') break;
        for(int j=0; j<SEARCH_RANGE;j++, p++)
            if(*(p)=='\n') break;
        p++;
    }
    blockList[0].base=p-mapped;

    for(size_t j=0;j<blockList.size() - 1;j++)
        blockList[j].length=blockList[j+1].base-blockList[j].base;

    p=mapped + sb.st_size - 1;
    if(*(p)=='\n')
        for(int j=0; j<SEARCH_RANGE;j++, p--)
            if(*(p)!='\n') break;
    blockList[blockList.size()-1].length=p - mapped - blockList[blockList.size() - 1].base;

#ifdef SHOW_PROGRESSBAR
    sprintf(buff,"Tagging");
    bar.set_option(indicators::option::PostfixText{buff});
    pThis->bar.set_progress(25.0);
#endif //!SHOW_PROGRESSBAR

    nodeNum = 0;
    base_offset=0;
    progUnit=(double)((100.0 - progress) / blockList.size());

    for(int j=0; j<=readNum; j++, base_offset+=BLOCK_READ)
    {
        if(j==readNum)
            blockNum=(restSize / BLOCK_SIZE) +1;
        else
            blockNum=MAXTHREAD;

        // Create thread
        for (int i = 0; i < blockNum; i++)
            pth[i]=new thread(pthFuncTag);

        // Join thread
        for (int i = 0; i < blockNum; i++)
        {
            pth[i]->join();
            delete pth[i];
        }
    }

    vector<BlockListNode* >pTmpList;
    for(int i=1;i<128;i++)
        if(chrList[i].base + 1)
            pTmpList.push_back(&(chrList[i]));
    sort(pTmpList.begin(),pTmpList.end(),cmpChrList);
    for(int i=0;i<pTmpList.size()-1;i++)
        pTmpList[i]->length=pTmpList[i+1]->base-pTmpList[i]->base;
    pTmpList[pTmpList.size()-1]->length=sb.st_size-pTmpList[pTmpList.size()-1]->base;
}

void BED::savechrList()
{
    char buff[NAME_MAX];
    strcpy(buff,bedname);
    strcat(buff,".tag");
    ofstream out(buff,ios::out);
    if(!out.is_open())
    {
        perror("ofstream: ") ;
        exit(FILE_SAVE_ERROR);
    }
    out << "# chr base  length" << endl;
    for(int i=1;i<MAX_CHR-3;i++)
    {
        if(chrList[i].base + 1)
            out << i << " " << chrList[i].base << " " << chrList[i].length << endl;
#ifdef __CHR_CHL
        if(chrList[MAX_CHR-4].base + 1)
        out << "CH" << " " << chrList[MAX_CHR-3].base << " " << chrList[MAX_CHR-3].length << endl;
#endif //!__CHR_CHL
    }
    if(chrList[MAX_CHR-3].base + 1)
        out << "MT" << " " << chrList[MAX_CHR-3].base << " " << chrList[MAX_CHR-3].length << endl;
#ifdef __CHR_ZW
    if(chrList[MAX_CHR-2].base + 1)
        out << "Z" << " " << chrList[MAX_CHR-2].base << " " << chrList[MAX_CHR-2].length << endl;
    if(chrList[MAX_CHR-1].base + 1)
        out << "W" << " " << chrList[MAX_CHR-1].base << " " << chrList[MAX_CHR-1].length << endl;
#else
    if(chrList[MAX_CHR-2].base + 1)
        out << "X" << " " << chrList[MAX_CHR-2].base << " " << chrList[MAX_CHR-2].length << endl;
    if(chrList[MAX_CHR-1].base + 1)
        out << "Y" << " " << chrList[MAX_CHR-1].base << " " << chrList[MAX_CHR-1].length << endl;
#endif //__CHR_ZW
    out.close();
}

void BED::processProfile(char *&gff3file)
{
    if((indexHandle = open(gff3file, O_RDWR)) < 0)
    {
        perror("gff3 open(): ") ;
        exit(FILE_OPEN_ERROR + 0xA);
    }

    if((fstat(indexHandle, &sbIndex)) == -1 )
    {
        perror("gff3 fstat(): ") ;
        exit(FILE_OPEN_ERROR + 0xB);
    }

    mappedIndex = (char*)mmap(nullptr, sbIndex.st_size, PROT_READ | PROT_WRITE, MAP_SHARED, indexHandle, 0);
    if(mappedIndex == (char*)-1)
    {
        perror("gff3 mmap(): ") ;
        exit(FILE_OPEN_ERROR + 0xC);
    }
    // Build up profile list

    size_t readNum=0;
    size_t restSize=0;
    size_t blockNum=0;
    thread* pth[MAXTHREAD];
    char buff[0x30];

    printf("\033[33m[Warning]\033[0m: Start generating profile ...\n");

#ifdef SHOW_PROGRESSBAR
    sprintf(buff,"Indexing");
    bar.set_option(indicators::option::PostfixText{buff});
    bar.set_progress(0);
    progress=0;
    progUnit=(double)(25.0 / (sbIndex.st_size / BLOCK_SIZE_INDEX +1));
#endif //!SHOW_PROGRESSBAR

    readNum=sbIndex.st_size / BLOCK_READ_INDEX;
    restSize=sbIndex.st_size % BLOCK_READ_INDEX;

    base_offset=0;

    for(int j=0; j<=readNum; j++, base_offset+=BLOCK_READ_INDEX)
    {
        threadNum = 0;
        if(j==readNum)
            blockNum=(restSize / BLOCK_SIZE_INDEX) + 1;
        else
            blockNum=MAXTHREAD;

        // Create thread
        for (int i = 0; i < blockNum; i++)
            pth[i]=new thread(pthFuncProfileList);

        // Join thread
        for (int i = 0; i < blockNum; i++)
        {
            pth[i]->join();
            delete pth[i];
        }
    }
    // Profile generation

    sort(profileList, profileList+geneNum, cmpProfileList);
    base_offset=0;
    readNum=geneNum / PROC_GENE_READ;
    restSize=geneNum % PROC_GENE_READ;

#ifdef SHOW_PROGRESSBAR
    sprintf(buff,"Processing");
    bar.set_option(indicators::option::PostfixText{buff});
    progUnit=(double)((100.0 - progress) / geneNum);
#endif //!SHOW_PROGRESSBAR

    for(int j=0; j<=readNum; j++, base_offset+=PROC_GENE_READ)
    {
        threadNum = 0;
        if(j==readNum)
            blockNum=(restSize / PROC_GENE_SIZE) + 1;
        else
            blockNum=MAXTHREAD;

        // Create thread
        for (int i = 0; i < blockNum; i++)
            pth[i]=new thread(pthFuncProfile);

        // Join thread
        for (int i = 0; i < blockNum; i++)
        {
            pth[i]->join();
            delete pth[i];
        }
    }

    munmap(mappedIndex,sbIndex.st_size);
    close(indexHandle);
}

void BED::saveProfile(const char *nameProfile)
{
    FILE* fout=fopen(nameProfile,"w");
    if(fout== nullptr)
    {
        perror("fopen(): ");
        exit(FILE_SAVE_ERROR + 0xA);
    }
    fprintf(fout,"chr\tID\tStart\tEnd\tStrand\tMethy_ratio");
#ifdef CG_NUMBER
    fprintf(fout,"\tCG");
#endif //!CG_number
    if(have_promoter)
    {
        fprintf(fout, "\tPromoter_methy_ratio");
#ifdef CG_NUMBER
        fprintf(fout,"\tCG_promoter");
#endif //!CG_number
    }
    fprintf(fout,"\n");
    for(int i=0;i<geneNum;i++)
    {
        if (fabs(profileList[i]->methy_ratio) <= 1e-15) continue;
        string str_chr;
        switch (profileList[i]->chr)
        {
#ifdef __CHR_CHL
            case MAX_CHR-4:
                str_chr="CH";
                break;
#endif //!__CHR_CHL
            case MAX_CHR - 3:
                str_chr = "MT";
                break;
#ifdef __CHR_ZW
                case MAX_CHR-2:
                    str_chr="Z";
                    break;
                case MAX_CHR-1:
                    str_chr="W";
                    break;
#else
            case MAX_CHR - 2:
                str_chr = "X";
                break;
            case MAX_CHR - 1:
                str_chr = "Y";
                break;
#endif //!__CHR_ZW
            default:
                str_chr = to_string(profileList[i]->chr);
                break;
        }
        fprintf(fout, "%s\t%s\t%ld\t%ld\t%c\t%.15lf", str_chr.c_str(), profileList[i]->ID,
                profileList[i]->Start, profileList[i]->End, (profileList[i]->chain) ? '+' : '-',
                profileList[i]->methy_ratio);
#ifdef CG_NUMBER
        fprintf(fout, "\t%ld", profileList[i]->NumCG);
#endif //!CG_number
        if (have_promoter)
        {
            fprintf(fout, "\t%.15lf", profileList[i]->methy_ratio_promoter);
#ifdef CG_NUMBER
            fprintf(fout, "\t%ld", profileList[i]->NumCG_promoter);
#endif //!CG_number
        }
        fprintf(fout, "\n");
    }
#ifdef ALLOW_PLUG_IN_SAVE
    LoadSavePlugInAndJmp(nameProfile);
#endif //!ALLOW_PLUG_IN_SAVE
    fclose(fout);
}

void BED::saveExternProfile(const char *nameProfileEx)
{
    FILE* fout=fopen(nameProfileEx,"w");
    if(fout== nullptr)
    {
        perror("fopen(): ");
        exit(FILE_SAVE_ERROR + 0xC);
    }
    fprintf(fout,"chr\tGeneID\tType\tStart\tEnd\tStrand\tMethy_ratio");
    fprintf(fout,"\n");
    for(int i=0;i<geneNum;i++)
    {
        if (fabs(profileList[i]->methy_ratio) <= 1e-15 || !profileList[i]->single_tag) continue;
        string str_chr;
        switch (profileList[i]->chr)
        {
#ifdef __CHR_CHL
            case MAX_CHR-4:
                str_chr="CH";
                break;
#endif //!__CHR_CHL
            case MAX_CHR - 3:
                str_chr = "MT";
                break;
#ifdef __CHR_ZW
                case MAX_CHR-2:
                    str_chr="Z";
                    break;
                case MAX_CHR-1:
                    str_chr="W";
                    break;
#else
            case MAX_CHR - 2:
                str_chr = "X";
                break;
            case MAX_CHR - 1:
                str_chr = "Y";
                break;
#endif //!__CHR_ZW
            default:
                str_chr = to_string(profileList[i]->chr);
                break;
        }
        fprintf(fout, "%s\t%s\t%s\t%ld\t%ld\t%c\t%.15lf", str_chr.c_str(), profileList[i]->ID, profileList[i]->str_type,
                profileList[i]->Start, profileList[i]->End, (profileList[i]->chain) ? '+' : '-',
                profileList[i]->methy_ratio);
        fprintf(fout, "\n");
        ExternNode* p=profileList[i]->next;
        while (true)
        {
            if(p== nullptr) break;
            fprintf(fout, "%s\t%s\t%s\t%ld\t%ld\t%c\t%.15lf", str_chr.c_str(), profileList[i]->ID, p->str_type,
                    p->Start, p->End, (p->chain) ? '+' : '-',
                    p->methy_ratio);
            fprintf(fout, "\n");
            p=p->next;
        }
    }
    fclose(fout);
}

void BED::saveSingleData(char *m_beg, char *m_end, char *ID)
{
    if(!ID) return;
    char name_buff[NAME_MAX]={0};
    sprintf(name_buff, "./single/%s", ID);

    FILE* sgf= fopen(name_buff, "wb");
    if(sgf == nullptr)
    {
        perror("fopen(): ");
        exit(FILE_SAVE_ERROR + 0xd);
    }

    fprintf(sgf, "#chr\tstart\tend\tstrand\tmCtype\tdepth\tmCdep\tlevel\n");
    fwrite(m_beg, m_end-m_beg, 1, sgf);
    fclose(sgf);
}

int BED::atoiChr(const char *nptr)
{
#ifdef __CHR_CHL
    if((*nptr)=='C'
    #ifdef ALLOW_LOWERCASE
        || (*nptr)=='c'
    #endif //!ALLOW_LOWERCASE
            )
        return MAX_CHR-4;
#endif //!__CHR_CHL

    if((*nptr)=='M'
    #ifdef ALLOW_LOWERCASE
        || (*nptr)=='m'
    #endif //!ALLOW_LOWERCASE
    )
        return MAX_CHR-3;
#ifdef __CHR_ZW
    if((*nptr)=='Z'
    #ifdef ALLOW_LOWERCASE
        || (*nptr)=='z'
    #endif //!ALLOW_LOWERCASE
            )
        return MAX_CHR-2;
    if((*nptr)=='W'
    #ifdef ALLOW_LOWERCASE
        || (*nptr)=='w'
    #endif //!ALLOW_LOWERCASE
            )
        return MAX_CHR-1;
#else
    if((*nptr)=='X'
    #ifdef ALLOW_LOWERCASE
        || (*nptr)=='x'
    #endif //!ALLOW_LOWERCASE
    )
        return MAX_CHR-2;
    if((*nptr)=='Y'
    #ifdef ALLOW_LOWERCASE
        || (*nptr)=='y'
    #endif //!ALLOW_LOWERCASE
    )
        return MAX_CHR-1;
#endif //!__CHR_ZW

    int c;              /* current char */
    int total;         /* current total */
    int sign;           /* if '-', then negative, otherwise positive */

    /* skip whitespace */
    while ( isspace((int)(unsigned char)*nptr) )
        ++nptr;

    c = (int)(unsigned char)*nptr++;
    sign = c;           /* save sign indication */
    if (c == '-' || c == '+')
        c = (int)(unsigned char)*nptr++;    /* skip sign */

    total = 0;

    while (isdigit(c))
    {
        total = 10 * total + (c - '0');     /* accumulate digit */
        c = (int)(unsigned char)*nptr++;    /* get next char */
    }

    if (sign == '-')
        return -total;
    else
        return total;   /* return result, negated if necessary */
}

bool BED::isspace(int x)
{
    if(x==' '||x=='\t'||x=='\n'||x=='\f'||x=='\b'||x=='\r')
        return true;
    else
        return false;
}

bool BED::isdigit(int x)
{
    if(x<='9'&&x>='0')
        return true;
    else
        return false;
}

void BED::LoadSavePlugInAndJmp(const char* foutput)
{
    string plug_in_name = XMACRO_STR(CMAKE_SOURCE_DIR);
    plug_in_name+="/plug-in/libsave.so";
    void *handle = dlopen(plug_in_name.c_str(), RTLD_LAZY);

    if(!handle)
    {
        printf("\n\033[31m[Error]\033[0m: %s", dlerror());
        exit(FILE_OPEN_ERROR+0x20);
    }
    typedef void (*SaveAs)(const char*, ProfileNode**, int);
    SaveAs saveas=(SaveAs)dlsym(handle,"ProfileSave");
    if(!saveas)
    {
        printf("%s", dlerror());
        dlclose(handle);
        exit(MUTEX_ERROR+0x10);
    }
    (*saveas)(foutput, pThis->profileList, pThis->geneNum);
    dlclose(handle);
}

void BED::loadSingleList(char *listfile)
{
    ifstream list_in(listfile, ios::in);

    auto* pset=new set<string>;
    string str_buff;

    if(!list_in.is_open())
    {
        perror("ifstream: ") ;
        exit(FILE_OPEN_ERROR+0x11);
    }

    while(getline(list_in, str_buff))
    {
        if(str_buff.empty()) continue;
        pset->insert(str_buff);
    }

    if(!pset->empty())
    {
        do_single_analyse=true;
        sglist=pset;
        /* To do. */
        int acc=
#ifdef _WIN32_PLATFORM_
                _access("single", 0);
#elif defined(_UNIX_PLATFORM_)
                access("single", 0);
        if(acc==-1)
        {
#endif // !_WIN32_PLATFORM_
            int mk =
#ifdef _WIN32_PLATFORM_
                    _mkdir("single");
#elif defined(_UNIX_PLATFORM_)
                    mkdir("single", 0755);
#endif // !_WIN32_PLATFORM_
            if (mk == -1) {
                printf("\033[31m[Error]\033[0m: Cannot create folder for single genes.\n");
                exit(DIRECTORY_CREATE_ERROR);
            }
        }
    }

    list_in.close();
}
