#include <iostream>
#include "bed.h"

using namespace std;
using namespace bed;

int main(int argc,char **argv)
{
    BED bedfile(argv[1]);
    bedfile.process(nullptr, Method::tag);
    bedfile.savechrList();
    bedfile.process(argv[2], nullptr, Method::profile);
    return 0;
}
