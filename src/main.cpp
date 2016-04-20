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
    cerr << "usage: " << argv[0] << " [options]" << endl
         << "Succinct representation of a queryable sequence graph" << endl
         << endl
         << "options:" << endl
         << "    -v, --vg FILE        compress graph in vg FILE" << endl
         << "    -V, --validate       validate compression" << endl
         << "    -o, --out FILE       serialize graph to FILE" << endl
         << "    -i, --in FILE        use index in FILE" << endl
         << "    -n, --node ID        graph neighborhood around node with ID" << endl
         << "    -c, --context N      steps of context to extract when building neighborhood" << endl
         << "    -s, --node-seq ID    provide node sequence for ID" << endl
         << "    -P, --char POS       give the character at a given position in the graph" << endl
         << "    -F, --substr POS:LEN extract the substr of LEN on the node at the position" << endl
         << "    -f, --edges-from ID  list edges from node with ID" << endl
         << "    -t, --edges-to ID    list edges to node with ID" << endl
         << "    -O, --edges-of ID    list all edges related to node with ID" << endl
         << "    -S, --edges-on-start ID    list all edges on start of node with ID" << endl
         << "    -E, --edges-on-end ID      list all edges on start of node with ID" << endl
         << "    -p, --path TARGET    gets the region of the graph @ TARGET (chr:start-end)" << endl
         << "    -x, --extract-threads      extract succinct threads as paths" << endl
         << "    -r, --store-threads  store perfect match paths as succinct threads" << endl
         << "    -R, --report FILE    save an HTML space usage report to FILE when serializing" << endl
         << "    -D, --debug          show debugging output" << endl
         << "    -T, --text-output    write text instead of vg protobuf" << endl
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
    bool edges_of = false;
    bool edges_on_start = false;
    bool edges_on_end = false;
    bool node_sequence = false;
    string pos_for_char;
    string pos_for_substr;
    int context_steps = 0;
    bool node_context = false;
    string target;
    bool print_graph = false;
    bool text_output = false;
    bool validate_graph = false;
    bool extract_threads = false;
    bool store_threads = false;
    string report_name;
    
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
                {"char", required_argument, 0, 'P'},
                {"substr", required_argument, 0, 'F'},
                //{"range", required_argument, 0, 'r'},
                {"context", required_argument, 0, 'c'},
                {"edges-from", required_argument, 0, 'f'},
                {"edges-to", required_argument, 0, 't'},
                {"edges-of", required_argument, 0, 'O'},
                {"edges-on-start", required_argument, 0, 'S'},
                {"edges-on-end", required_argument, 0, 'E'},
                {"node-seq", required_argument, 0, 's'},
                {"path", required_argument, 0, 'p'},
                {"extract-threads", no_argument, 0, 'x'},
                {"store-threads", no_argument, 0, 'r'},
                {"report", required_argument, 0, 'R'},
                {"debug", no_argument, 0, 'D'},
                {"text-output", no_argument, 0, 'T'},
                {"validate", no_argument, 0, 'V'},
                {0, 0, 0, 0}
            };

        int option_index = 0;
        c = getopt_long (argc, argv, "hv:o:i:f:t:s:c:n:p:DxrTO:S:E:VR:P:F:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {

        case 'v':
            vg_name = optarg;
            break;

        case 'V':
            validate_graph = true;
            break;

        case 'o':
            out_name = optarg;
            break;

        case 'D':
            print_graph = true;
            break;

        case 'T':
            text_output = true;
            break;
            
        case 'x':
            extract_threads = true;
            break;
            
        case 'r':
            store_threads = true;
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

        case 'O':
            node_id = atol(optarg);
            edges_of = true;
            break;

        case 'S':
            node_id = atol(optarg);
            edges_on_start = true;
            break;

        case 'E':
            node_id = atol(optarg);
            edges_on_end = true;
            break;

        case 's':
            node_id = atol(optarg);
            node_sequence = true;
            break;

        case 'p':
            target = optarg;
            break;

        case 'P':
            pos_for_char = optarg;
            break;
            
        case 'F':
            pos_for_substr = optarg;
            break;
            
        case 'R':
            report_name = optarg;
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

    XG* graph = nullptr;
    //string file_name = argv[optind];
    if (in_name.empty()) assert(!vg_name.empty());
    if (vg_name == "-") {
        graph = new XG;
        graph->from_stream(std::cin, validate_graph, print_graph, store_threads);
    } else if (vg_name.size()) {
        ifstream in;
        in.open(vg_name.c_str());
        graph = new XG;
        graph->from_stream(in, validate_graph, print_graph, store_threads);
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

    // Prepare structure tree for serialization
    unique_ptr<sdsl::structure_tree_node> structure;
    
    if(!report_name.empty()) {
        // We need to make a report, so we need the structure. Make a real tree
        // node. The unique_ptr handles deleting.
        structure = unique_ptr<sdsl::structure_tree_node>(new sdsl::structure_tree_node("name", "type"));
    }

    if (out_name.size()) {
        if (out_name == "-") {
            graph->serialize(std::cout, structure.get(), "xg");
            std::cout.flush();
        } else {
            ofstream out;
            out.open(out_name.c_str());
            graph->serialize(out, structure.get(), "xg");
            out.flush();
        }
    }

    if(!report_name.empty()) {
        // Save the report
        ofstream out;
        out.open(report_name.c_str());
        sdsl::write_structure_tree<HTML_FORMAT>(structure.get(), out, 0);
    }

    // queries
    if (node_sequence) {
        cout << node_id << ": " << graph->node_sequence(node_id) << endl;
    }
    if (!pos_for_char.empty()) {
        // extract the position from the string
        int64_t id;
        bool is_rev;
        size_t off;
        extract_pos(pos_for_char, id, is_rev, off);
        // then pick it up from the graph
        cout << graph->pos_char(id, is_rev, off) << endl;
    }
    if (!pos_for_substr.empty()) {
        int64_t id;
        bool is_rev;
        size_t off;
        size_t len;
        extract_pos_substr(pos_for_substr, id, is_rev, off, len);
        cout << graph->pos_substr(id, is_rev, off, len) << endl;
    }
    
    if (edges_from) {
        vector<Edge> edges = graph->edges_from(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_to) {
        vector<Edge> edges = graph->edges_to(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_of) {
        vector<Edge> edges = graph->edges_of(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_on_start) {
        vector<Edge> edges = graph->edges_on_start(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }
    if (edges_on_end) {
        vector<Edge> edges = graph->edges_on_end(node_id);
        for (auto& edge : edges) {
            cout << edge.from() << (edge.from_start()?"-":"+")
                 << " -> " << edge.to() << (edge.to_end()?"-":"+") << endl;
        }
    }

    if (node_context) {
        Graph g;
        graph->neighborhood(node_id, context_steps, g);
        if (text_output) {
            to_text(cout, g);
        } else {
            vector<Graph> gb = { g };
            stream::write_buffered(cout, gb, 0);
        }
    }

    if (!target.empty()) {
        string name;
        int64_t start, end;
        Graph g;
        parse_region(target, name, start, end);
        graph->get_path_range(name, start, end, g);
        graph->expand_context(g, context_steps);
        if (text_output) {
            to_text(cout, g);
        } else {
            vector<Graph> gb = { g };
            stream::write_buffered(cout, gb, 0);
        }
    }
    
    if (extract_threads) {
        
    
        size_t thread_number = 0;
        for(Path& path : graph->extract_threads()) {
            // Give each thread a name
            path.set_name("_thread_" + to_string(thread_number++));
            
            // We need a Graph for serialization purposes. We do one chunk per
            // thread in case the threads are long.
            Graph g;
            
            *(g.add_path()) = path;
            
            // Dump the graph with its mappings. TODO: can we restrict these to
            // mappings to nodes we have already pulled out? Or pull out the
            // whole compressed graph?
            if (text_output) {
                to_text(cout, g);
            } else {
                vector<Graph> gb = { g };
                stream::write_buffered(cout, gb, 0);
            }
            
        }
    }

    // clean up
    if (graph) delete graph;

    return 0;
}
