#include <iostream>
#include <fstream>
#include <getopt.h>
#include "sdsl/bit_vectors.hpp"
#include "stream.hpp"
#include "cpp/vg.pb.h"
#include "succinct.hpp"

using namespace std;
using namespace sdsl;
using namespace vg;
using namespace scg;

void help_main(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <graph.vg> >[graph.scg]" << endl
         << "Compressed graph, emits on stdout." << endl
         << endl
         << "options:" << endl
         << "    -h, --help     This text" << endl;
}

int main(int argc, char** argv) {

    if (argc == 1) {
        help_main(argv);
        return 1;
    }

    int c;
    optind = 1; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "h",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'h':
        case '?':
            help_main(argv);
            exit(1);
            break;

        default:
            abort ();
        }
    }

    SuccinctGraph* graph;
    string file_name = argv[optind];
    if (file_name == "-") {
        graph = new SuccinctGraph(std::cin);
    } else {
        ifstream in;
        in.open(file_name.c_str());
        graph = new SuccinctGraph(in);
    }

    //graph->serialize_to_ostream(std::cout);

    // fail...
    //delete graph;

    return 0;
}
