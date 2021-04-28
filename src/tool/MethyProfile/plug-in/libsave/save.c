#include <config.h>
#include "../../profile_node.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <limits.h>
#include <math.h>

extern "C"
{

    void __stdcall ProfileSave(const char* nameProfile, ProfileNode** profileList, int n)
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
        fprintf(fout, "\tCG");
    #endif //!CG_number
        fprintf(fout, "\n");

        for(int i=0;i<n;i++)
        {
            if (fabs(profileList[i]->methy_ratio) <= 1e-15) continue;
            char str_chr[0x10];
            switch (profileList[i]->chr)
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
                    sprintf(str_chr, "%d",profileList[i]->chr);
                    break;
            }
            char buff[0x30];
            sprintf(buff, "%s_promoter", profileList[i]->ID);
            fprintf(fout, "%s\t%s\t%ld\t%ld\t%c\t%.15lf", str_chr, buff,
                    profileList[i]->Start, profileList[i]->End, (profileList[i]->chain) ? '+' : '-',
                    profileList[i]->methy_ratio_promoter);
    #ifdef CG_NUMBER
            fprintf(fout, "\t%ld", profileList[i]->NumCG);
    #endif //!CG_number
            fprintf(fout, "\n");
        }
        fclose(fout);
    }

}