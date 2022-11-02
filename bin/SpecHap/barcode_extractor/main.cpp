#include <iostream>
#include "BAMReader.h"
#include "worker.h"



int main(int argc, char *argv[])
{
    if (argc >= 3)
    {
        worker *worker1 = new worker(argv[1], argv[2]);
        if (argc == 3)
            worker1->extract();
        else if (argc == 4)
            worker1->extract_one_target(argv[3]);
        delete worker1;
    }
    else
        std::cerr << "Usage: BarcodeExtract in_bam out <chr>" << std::endl;
    return 0;
}