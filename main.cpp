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
         << "    -v, --vg FILE   compress graph in vg file" << endl
         << "    -o, --out FILE  serialize graph to file" << endl
         << "    -h, --help      This text" << endl;
}

int main(int argc, char** argv) {

    if (argc == 1) {
        help_main(argv);
        return 1;
    }

    string vg_name;
    string out_name;

    int c;
    optind = 1; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"vg", required_argument, 0, 'v'},
                {"out", required_argument, 0, 'o'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hv:o:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'v':
            vg_name = optarg;
            break;

        case 'o':
            out_name = optarg;
            break;

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
    //string file_name = argv[optind];
    assert(!vg_name.empty());
    if (vg_name == "-") {
        graph = new SuccinctGraph;
        graph->from_vg(std::cin);
    } else {
        ifstream in;
        in.open(vg_name.c_str());
        graph = new SuccinctGraph;
        graph->from_vg(in);
    }


    if (out_name.size()) {
        if (out_name == "-") {
            graph->serialize(std::cout);
            std::cout.flush();
        } else {
            ofstream out;
            out.open(out_name.c_str());
            graph->serialize(out);
            out.flush();
        }
    }

    // fail...
    //delete graph;

    return 0;
}
