#include <iostream>
#include "bed.h"
#include <getopt.h>
#include <fstream>
#include <ctype.h>
#include <string.h>

#include <string>
#include <set>

using namespace std;
using namespace bed;

void usage()
/* Explain usage. */
{
    printf(
            "MethyProfile v1.0 - Generating methylation profile.\n"
            "usage:\n"
            "   methyprofile [options] input.bed input.gff3 (output.methyprofile.txt)\n"
            "options:\n"
            "   -P, --promoter=n    - Analysing promoters methylation information with n bp. The default is 2000bp.\n"
            "   -l, --single-list   - Gathering single gene information for a gene list behind.\n"
            "   -h, --help          - Show this message.\n"
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
            {"help", no_argument, nullptr, 'h'},
            {0, 0, 0, 0}
    };
    static const char simple_options[]="P::l:h";
    int longindex = -1;
    int opt;

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
                bedfile.loadSingleList(optarg+1);
                break;
            case 'h':
                usage();
                return 0;
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
