#include "../config.h"
#include "../profile_node.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

extern "C"
{

    void ProfileSave(const char* nameProfile, ProfileNode** profileList, int n)
    {
        FILE* fout=fopen(nameProfile,"w");
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
            if (abs(profileList[i]->methy_ratio) <= 1e-15) continue;
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
            fprintf(fout, "%s\t%s\t%ld\t%ld\t%c\t%.15lf", str_chr, profileList[i]->ID,
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