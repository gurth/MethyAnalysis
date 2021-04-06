#include <iostream>
#include "bed.h"
#include <getopt.h>
#include <fstream>

using namespace std;
using namespace bed;

void usage()
/* Explain usage. */
{
    printf(
            "MethyProfile v1 - Generating methylation profile.\n"
            "usage:\n"
            "   bigBedToBed input.bed input.gff3 (output.methyprofile.txt)\n"
            "options:\n"
            "   -P, --promoter=n    - Analysing promoters methylation information with n bp. The default is 2000bp.\n"
            "   -l, --single-list   - Gathering single gene information for a gene list behind.\n"
    );
}

int main(int argc,char **argv)
{
    if(argc<1)
    {
        printf("\033[31m[Error]\033[0m:No input file.\n");
        usage();
        exit(ARGUMENT_ERROR);
    }

    BED bedfile;
    static option long_options[]={
            {"promoter", optional_argument, nullptr, 'P'},
            {"single-list", required_argument, nullptr, 'l'},
            {0, 0, 0, 0}
    };
    static const char simple_options[]="P::l:";
    int longindex = -1;
    int opt;
    ifstream flist;

    while (true)
    {
        opt = getopt_long(argc, argv, simple_options, long_options, &longindex);
        if(opt == -1) break;
        switch (opt)
        {
            case 'P':
                bedfile.have_promoter = true;
                if(optarg)
                    bedfile.promoterLen = atoi(optarg+1);
                break;
            case 'l':
                flist.open(optarg+1, ios::in);
                if(!flist.is_open())
                {
                    perror("fstream::open(): ");
                    exit(FILE_OPEN_ERROR+0x10);
                }
                flist.close();
                break;
            default:
                printf("\033[31m[Error]\033[0m:Wrong argument.\n");
                usage();
                exit(ARGUMENT_ERROR);
                break;
        }
    }

    bedfile.bedfileOpen(argv[optind]);
    bedfile.process(nullptr, Method::tag);
    bedfile.savechrList();
    bedfile.process(argv[optind+1], argv[optind+2], Method::profile);
    return 0;
}
