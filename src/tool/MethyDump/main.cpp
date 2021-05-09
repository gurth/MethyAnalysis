#include <iostream>
#include <config.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include "my_config.h"
#include <profile_node.h>


using namespace std;

void ProfileSave(const char* nameProfile, ProfileNode* profileList, int n)
{
    char name_buff[PATH_MAX];
    sprintf(name_buff, "%s.pro.txt", nameProfile);
    FILE* fout=fopen(name_buff,"w");
    if(fout== nullptr)
    {
        perror("fopen(): ");
        exit(FILE_SAVE_ERROR + 0xA);
    }
    fprintf(fout,"chr\tID\tStart\tEnd\tStrand\tPromoter_methy_ratio");
#ifdef CG_NUMBER
    fprintf(fout, "\tCG_promoter");
#endif //!CG_number
    fprintf(fout, "\n");

    for(int i=0;i<n;i++)
    {
        if (fabs(profileList[i].methy_ratio) <= 1e-15) continue;
        char str_chr[0x10];
        switch (profileList[i].chr)
        {
#ifdef __CHR_CHL
            case MAX_CHR-4:
                    strcpy(str_chr, "CH");
                    break;
#endif //!__CHR_CHL
            case MAX_CHR - 3:
                strcpy(str_chr, "MT");
                break;
#ifdef __CHR_ZW
                case MAX_CHR-2:
                        strcpy(str_chr, "Z");
                        break;
                    case MAX_CHR-1:
                        strcpy(str_chr, "W");
                        break;
#else
            case MAX_CHR - 2:
                strcpy(str_chr, "X");
                break;
            case MAX_CHR - 1:
                strcpy(str_chr, "Y");
                break;
#endif //!__CHR_ZW
            default:
                sprintf(str_chr, "%d",profileList[i].chr);
                break;
        }
        char buff[0x30];
        sprintf(buff, "%s_promoter", profileList[i].ID);
        fprintf(fout, "%s\t%s\t%lld\t%lld\t%c\t%.15lf", str_chr, buff,
                profileList[i].Start, profileList[i].End, (profileList[i].chain) ? '+' : '-',
                profileList[i].methy_ratio_promoter);
#ifdef CG_NUMBER
        fprintf(fout, "\t%lld", profileList[i].NumCG_promoter);
#endif //!CG_number
        fprintf(fout, "\n");
    }
    fclose(fout);
}

int main(int argc, char*argv[])
{
    if(argc<=1)
    {
        printf("\033[31m[Error]\033[0m:No input file.\n");
        exit(1);
    }

    int len=0;
    FILE * pf= fopen(argv[1], "rb");
    ProfileNode* profileList = nullptr;

    if(pf== nullptr)
    {
        printf("\033[31m[Error]\033[0m:Open file error.\n");
        exit(1);
    }
    fread(&len, sizeof(int), 1, pf);
    profileList=(ProfileNode*) malloc(sizeof(ProfileNode)*len);
    fread(profileList, sizeof(ProfileNode)*len, 1, pf);
    fclose(pf);
    ProfileSave("a.txt",profileList, len);
    return 0;
}
