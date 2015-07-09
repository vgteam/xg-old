#include <iostream>
#include <fstream>
#include <getopt.h>
#include "sdsl/bit_vectors.hpp"
#include "stream.hpp"
#include "cpp/vg.pb.h"
#include "xg.hpp"

using namespace std;
using namespace sdsl;
using namespace vg;
using namespace xg;

void help_main(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <graph.vg> >[graph.scg]" << endl
         << "Compressed graph, emits on stdout." << endl
         << endl
         << "options:" << endl
         << "    -v, --vg FILE        compress graph in vg FILE" << endl
         << "    -o, --out FILE       serialize graph to FILE" << endl
         << "    -i, --in FILE        use index in FILE" << endl
         << "    -n, --node ID        graph neighborhood around node with ID" << endl
         << "    -c, --context N      steps of context to extract when building neighborhood" << endl
         << "    -s, --node-seq ID    provide node sequence for ID" << endl
         << "    -f, --edges-from ID  list edges from node with ID" << endl
         << "    -t, --edges-to ID    list edges to node with ID" << endl
         << "    -p, --path TARGET    gets the region of the graph @ TARGET (chr:start-end)" << endl
         << "    -h, --help           this text" << endl;
}

int main(int argc, char** argv) {

    if (argc == 1) {
        help_main(argv);
        return 1;
    }

    string vg_name;
    string out_name;
    string in_name;
    int64_t node_id;
    bool edges_from = false;
    bool edges_to = false;
    bool node_sequence = false;
    int context_steps = 0;
    bool node_context = false;
    string target;
    
    int c;
    optind = 1; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
            {
                {"help", no_argument, 0, 'h'},
                {"vg", required_argument, 0, 'v'},
                {"out", required_argument, 0, 'o'},
                {"in", required_argument, 0, 'i'},
                {"node", required_argument, 0, 'n'},
                //{"range", required_argument, 0, 'r'},
                {"context", required_argument, 0, 'c'},
                {"edges-from", required_argument, 0, 'f'},
                {"edges-to", required_argument, 0, 't'},
                {"node-seq", required_argument, 0, 's'},
                {"path", required_argument, 0, 'p'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hv:o:i:f:t:s:c:n:p:",
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

        case 'i':
            in_name = optarg;
            break;

        case 'n':
            node_id = atol(optarg);
            node_context = true;
            break;

        case 'c':
            context_steps = atoi(optarg);
            break;

        case 'f':
            node_id = atol(optarg);
            edges_from = true;
            break;
            
        case 't':
            node_id = atol(optarg);
            edges_to = true;
            break;

        case 's':
            node_id = atol(optarg);
            node_sequence = true;
            break;

        case 'p':
            target = optarg;
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

    XG* graph;
    //string file_name = argv[optind];
    if (in_name.empty()) assert(!vg_name.empty());
    if (vg_name == "-") {
        graph = new XG;
        graph->from_vg(std::cin);
    } else if (vg_name.size()) {
        ifstream in;
        in.open(vg_name.c_str());
        graph = new XG;
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

    if (in_name.size()) {
        graph = new XG;
        if (in_name == "-") {
            graph->load(std::cin);
        } else {
            ifstream in;
            in.open(in_name.c_str());
            graph->load(in);
        }
    }

    // queries
    if (node_sequence) {
        cout << node_id << ": " << graph->node_sequence(node_id) << endl;
    }
    if (edges_from) {
        vector<Edge> edges = graph->edges_from(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << " -> " << edge.to() << endl;
        }
    }
    if (edges_to) {
        vector<Edge> edges = graph->edges_to(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << " -> " << edge.to() << endl;
        }
    }

    if (node_context) {
        Graph g;
        graph->neighborhood(node_id, context_steps, g);
        vector<Graph> gb = { g };
        stream::write_buffered(cout, gb, 0);
    }

    if (!target.empty()) {
        string name;
        int64_t start, end;
        Graph g;
        parse_region(target, name, start, end);
        graph->get_path_range(name, start, end, g);
        if (context_steps > 0) {
            graph->expand_context(g, context_steps);
        }
        vector<Graph> gb = { g };
        stream::write_buffered(cout, gb, 0);
    }

    // fail...
    //delete graph;

    return 0;
}
