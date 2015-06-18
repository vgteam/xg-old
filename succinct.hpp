#ifndef SUCCINCT_GRAPH_SG_HPP
#define SUCCINCT_GRAPH_SG_HPP

#include <iostream>
#include <fstream>
#include <map>
#include "cpp/vg.pb.h"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/enc_vector.hpp"
#include "sdsl/vlc_vector.hpp"

namespace scg {

using namespace std;
using namespace sdsl;
using namespace vg;

class SuccinctGraph {
public:
    SuccinctGraph(void) { }
    ~SuccinctGraph(void) { }
    SuccinctGraph(istream& file);
    // build up interface here
    void add_node(Node& node);
    void add_edge(Edge& edge);
    Node node(int64_t rank); // id?
    Edge edge(int64_t rank);
    Path path(string& name);
    Graph neighborhood(int64_t rank, int32_t steps);
    Graph range(int64_t rank1, int64_t rank2);
    Graph region(string& path_name, int64_t start, int64_t stop);
private:
    // for construction
    map<int64_t, string> node_label;
    map<int64_t, set<int64_t> > from_to;
    map<int64_t, set<int64_t> > to_from;
    map<string, set<int64_t> > path_nodes;
    size_t seq_length;
    size_t node_count;
    size_t edge_count;
    size_t path_entry_count;
    
    // sequence/integer vector
    int_vector<> s_iv;
    // maintain old ids from input
    //int_vector<> i_iv;
    // node starts in sequence, provides id schema
    // rank_1(i) = id
    // select_1(id) = i
    bit_vector s_bv; // node positions in siv
    int_vector<> f_iv;
    bit_vector f_bv;
    int_vector<> t_iv;
    bit_vector t_bv;
    map<string, bit_vector> p_bv;

};

}

#endif
