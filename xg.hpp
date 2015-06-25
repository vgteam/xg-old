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
//#include "sdsl/csa_bitcompressed.hpp"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"


namespace xg {

using namespace std;
using namespace sdsl;
using namespace vg;

class XG {
public:
    XG(void) : path_name_marker('#'),
                          seq_length(0),
                          node_count(0),
                          edge_count(0),
                          path_count(0) { }
    ~XG(void) { }
    XG(istream& in);
    void from_vg(istream& in);
    void load(istream& in);
    size_t serialize(std::ostream& out, sdsl::structure_tree_node* v = NULL, std::string name = "");
    size_t seq_length;
    size_t node_count;
    size_t edge_count;
    size_t path_count;
    size_t id_to_rank(int64_t id);
    int64_t rank_to_id(size_t rank);
    size_t max_node_rank(void);
    size_t max_path_rank(void);
    Node node(int64_t id); // gets node sequence
    string node_sequence(int64_t id);
    vector<Edge> edges_to(int64_t id);
    vector<Edge> edges_from(int64_t id);
    size_t node_rank_as_entity(int64_t id);
    size_t edge_rank_as_entity(int64_t id1, int64_t id2);
    Path path(const string& name);
    size_t path_rank(const string& name);
    string path_name(size_t rank);
    vector<size_t> paths_of_entity(size_t rank);
    vector<size_t> paths_of_node(int64_t id);
    vector<size_t> paths_of_edge(int64_t id1, int64_t id2);
    map<string, Mapping> node_mappings(int64_t id);
    bool path_contains_node(const string& name, int64_t id);
    bool path_contains_edge(const string& name, int64_t id1, int64_t id2);
    bool path_contains_entity(const string& name, size_t rank);
    bool has_edge(int64_t id1, int64_t id2);
    void neighborhood(int64_t id, size_t steps, Graph& g);
    void range(int64_t rank1, int64_t rank2, Graph& g);
    void region(string& path_name, int64_t start, int64_t stop, Graph& g);
    char path_name_marker;
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
    // edge table, allows o(1) determination of edge existence
    int_vector<> e_iv;
    //csa_wt<> e_csa;
    csa_sada<> e_csa;
    // allows lookups of id->rank mapping
    wt_int<> i_wt;
    // paths: serialized as bitvectors over nodes and edges
    int_vector<> pn_iv; // path names
    csa_wt<> pn_csa; // path name compressed suffix array
    bit_vector pn_bv;  // path name starts in uncompressed version of csa
    rank_support_v<1> pn_bv_rank;
    bit_vector::select_1_type pn_bv_select;
    int_vector<> pi_iv; // path ids by rank in the path names
    vector<bit_vector> pe_v; // path entity membership (ordered by rank in pn_iv
    vector<int_vector<>> pp_v; // path relative positions to each entity
    // entity/path mapping
    // but is over the entire entity space
    int_vector<> pe_iv; // list of path ids for each entity
    bit_vector pe_bv; // entity starts in ep vector
    wt_int<> pe_wt; // allows quick lookup of all entities in a particular path
};

Mapping new_mapping(const string& name, int64_t id);

}

#endif
