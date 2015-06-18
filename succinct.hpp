#ifndef SUCCINCT_GRAPH_SG_HPP
#define SUCCINCT_GRAPH_SG_HPP

#include <iostream>
#include <fstream>
#include <map>
#include "cpp/vg.pb.h"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/enc_vector.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/vlc_vector.hpp"
#include "sdsl/wavelet_trees.hpp"

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
    Node node(int64_t id); // gets node sequence
    string& node_sequence(int64_t id);
    Edge edges_to(int64_t id);
    Edge edges_from(int64_t id);
    Path path(string& name);
    Graph neighborhood(int64_t rank, int32_t steps);
    Graph range(int64_t rank1, int64_t rank2);
    Graph region(string& path_name, int64_t start, int64_t stop);
private:
    // sequence/integer vector
    int_vector<> s_iv;
    // node starts in sequence, provides id schema
    // rank_1(i) = id
    // select_1(id) = i
    bit_vector s_bv; // node positions in siv
    rank_support_v<1> s_bv_rank;
    bit_vector::select_1_type s_bv_select;
    rrr_vector<> s_cbv;
    rrr_vector<>::rank_1_type s_cbv_rank;
    rrr_vector<>::select_1_type s_cbv_select;
    // maintain old ids from input, ranked as in s_iv and s_bv
    int_vector<> i_iv;
    // maintain forward links
    int_vector<> f_iv;
    bit_vector f_bv;
    rank_support_v<1> f_bv_rank;
    bit_vector::select_1_type f_bv_select;
    // and the same data in the reverse direction
    int_vector<> t_iv;
    bit_vector t_bv;
    rank_support_v<1> t_bv_rank;
    bit_vector::select_1_type t_bv_select;
    map<string, bit_vector> p_bv;

};

}

#endif
