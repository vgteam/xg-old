#ifndef SUCCINCT_GRAPH_SG_HPP
#define SUCCINCT_GRAPH_SG_HPP

#include <iostream>
#include <fstream>
#include <map>
#include <omp.h>
#include "cpp/vg.pb.h"
#include "sdsl/bit_vectors.hpp"
#include "sdsl/enc_vector.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/vlc_vector.hpp"
#include "sdsl/wavelet_trees.hpp"
//#include "sdsl/csa_bitcompressed.hpp"
#include "sdsl/csa_wt.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "dynamic.hpp"

namespace xg {

using namespace std;
using namespace sdsl;
using namespace vg;

class Traversal {
public:
    int64_t id;
    bool rev;
    Traversal(int64_t i, bool r) : id(i), rev(r) { }
};

class XGPath;
typedef pair<int64_t, bool> Side;

class XG {
public:
    
    XG(void) : start_marker('#'),
               end_marker('$'),
               seq_length(0),
               node_count(0),
               edge_count(0),
               path_count(0) { }
    ~XG(void) { }
    XG(istream& in);
    XG(Graph& graph);
    void from_stream(istream& in, bool validate_graph = false, bool print_graph = false);
    void from_graph(Graph& graph, bool validate_graph = false, bool print_graph = false);
    void build(map<int64_t, string>& node_label,
               map<Side, set<Side> >& from_to,
               map<Side, set<Side> >& to_from,
               map<string, map<int, Mapping>>& path_nodes,
               bool validate_graph,
               bool print_graph);
    void load(istream& in);
    size_t serialize(std::ostream& out,
                     sdsl::structure_tree_node* v = NULL,
                     std::string name = "");
    size_t seq_length;
    size_t node_count;
    size_t edge_count;
    size_t path_count;

    size_t id_to_rank(int64_t id) const;
    int64_t rank_to_id(size_t rank) const;
    size_t max_node_rank(void) const;
    Node node(int64_t id) const; // gets node sequence
    string node_sequence(int64_t id) const;
    vector<Edge> edges_of(int64_t id) const;
    vector<Edge> edges_to(int64_t id) const;
    vector<Edge> edges_from(int64_t id) const;
    vector<Edge> edges_on_start(int64_t id) const;
    vector<Edge> edges_on_end(int64_t id) const;
    size_t node_rank_as_entity(int64_t id) const;
    size_t edge_rank_as_entity(int64_t id1, bool from_start, int64_t id2, bool to_end) const;
    bool entity_is_node(size_t rank) const;
    size_t entity_rank_as_node_rank(size_t rank) const;
    bool has_edge(int64_t id1, bool is_start, int64_t id2, bool is_end) const;

    Path path(const string& name) const;
    size_t path_rank(const string& name) const;
    size_t max_path_rank(void) const;
    string path_name(size_t rank) const;
    vector<size_t> paths_of_entity(size_t rank) const;
    vector<size_t> paths_of_node(int64_t id) const;
    vector<size_t> paths_of_edge(int64_t id1, bool from_start, int64_t id2, bool to_end) const;
    map<string, vector<Mapping>> node_mappings(int64_t id) const;
    bool path_contains_node(const string& name, int64_t id) const;
    bool path_contains_edge(const string& name,
                            int64_t id1, bool from_start,
                            int64_t id2, bool to_end) const;
    bool path_contains_entity(const string& name, size_t rank) const;
    void add_paths_to_graph(map<int64_t, Node*>& nodes, Graph& g) const;
    size_t node_occs_in_path(int64_t id, const string& name) const;
    vector<size_t> node_ranks_in_path(int64_t id, const string& name) const;
    vector<size_t> node_positions_in_path(int64_t id, const string& name) const;
    int64_t node_at_path_position(const string& name, size_t pos) const;
    Mapping mapping_at_path_position(const string& name, size_t pos) const;
    size_t path_length(const string& name) const;

    void neighborhood(int64_t id, size_t steps, Graph& g) const;
    //void for_path_range(string& name, int64_t start, int64_t stop, function<void(Node)> lambda);
    void get_path_range(string& name, int64_t start, int64_t stop, Graph& g) const;
    // basic method to query regions of the graph
    // add_paths flag allows turning off the (potentially costly, and thread-locking) addition of paths
    // when these are not necessary
    void expand_context(Graph& g, size_t steps, bool add_paths = true) const;
    void get_connected_nodes(Graph& g) const;
    void get_id_range(int64_t id1, int64_t id2, Graph& g) const;

    // gPBWT interface
    // Insert a thread. Path name must be unique or empty.
    void insert_thread(const Path& t);
    // Read all the threads embedded in the graph.
    list<Path> extract_threads() const;
    // Extract a particular thread by name. Name may not be empty.
    Path extract_thread(const string& name) const;
    // Count matches to a subthread among embedded threads
    size_t count_matches(const Path& t) const;
    
    /**
     * Represents the search state for the graph PBWT, so that you can continue
     * a search with more of a thread, or backtrack.
     *
     * By default, represents an un-started search (with no first visited side)
     * that can be extended to the whole collection of visits to a side.
     */
    struct ThreadSearchState {
        // What side have we just arrived at in the search?
        int64_t currentSide = 0;
        // What is the first visit at that side that is selected?
        int64_t rangeStart = 0;
        // And what is the past-the-last visit that is selected?
        int64_t rangeEnd = numeric_limits<int64_t>::max();
        
        // How many visits are selected?
        inline int64_t count() {
            return rangeEnd - rangeStart;
        }
        
        // Return true if the range has nothing selected.
        inline bool isEmpty() {
            return rangeEnd <= rangeStart;
        }
    };
    
    // Extend a search with the given section of a thread.
    void extend(ThreadSearchState& state, const Path& t) const;

    
    char start_marker;
    char end_marker;
    
private:

    // sequence/integer vector
    int_vector<> s_iv;
    // node starts in sequence, provides id schema
    // rank_1(i) = id
    // select_1(id) = i
    bit_vector s_bv; // node positions in siv
    rank_support_v<1> s_bv_rank;
    bit_vector::select_1_type s_bv_select;
    // compressed version, unused...
    rrr_vector<> s_cbv;
    rrr_vector<>::rank_1_type s_cbv_rank;
    rrr_vector<>::select_1_type s_cbv_select;

    // maintain old ids from input, ranked as in s_iv and s_bv
    int_vector<> i_iv;
    int64_t min_id; // id ranges don't have to start at 0
    int64_t max_id;
    int_vector<> r_iv; // ids-id_min is the rank

    // maintain forward links
    int_vector<> f_iv;
    bit_vector f_bv;
    rank_support_v<1> f_bv_rank;
    bit_vector::select_1_type f_bv_select;
    bit_vector f_from_start_bv;
    bit_vector f_to_end_bv;
    sd_vector<> f_from_start_cbv;
    sd_vector<> f_to_end_cbv;

    // and the same data in the reverse direction
    int_vector<> t_iv;
    bit_vector t_bv;
    rank_support_v<1> t_bv_rank;
    bit_vector::select_1_type t_bv_select;
    // these bit vectors are only used during construction
    // perhaps they should be moved?
    bit_vector t_from_start_bv;
    bit_vector t_to_end_bv;
    // used at runtime
    sd_vector<> t_from_start_cbv;
    sd_vector<> t_to_end_cbv;

    // edge table, allows o(1) determination of edge existence
    int_vector<> e_iv;

    //csa_wt<> e_csa;
    //csa_sada<> e_csa;

    // allows lookups of id->rank mapping
    //wt_int<> i_wt;

    // paths: serialized as bitvectors over nodes and edges
    int_vector<> pn_iv; // path names
    csa_wt<> pn_csa; // path name compressed suffix array
    bit_vector pn_bv;  // path name starts in uncompressed version of csa
    rank_support_v<1> pn_bv_rank;
    bit_vector::select_1_type pn_bv_select;
    int_vector<> pi_iv; // path ids by rank in the path names

    // probably these should get compressed, for when we have whole genomes with many chromosomes
    // the growth in required memory is quadratic but the stored matrix is sparse
    vector<XGPath*> paths; // path entity membership

    // entity->path membership
    int_vector<> ep_iv;
    bit_vector ep_bv; // entity delimiters in ep_iv
    rank_support_v<1> ep_bv_rank;
    bit_vector::select_1_type ep_bv_select;
    
    // Succinct thread storage
    
    // Threads are haplotype paths in the graph with no edits allowed, starting
    // and stopping at node boundaries.
    
    // We keep our strings in this dynamic succinct rank-select string from DYNAMIC.
    using dynamic_int_vector = dyn::rle_str;
    
    // TODO: Explain the whole graph PBWT extension here
    
    // Basically we keep usage counts for every element in the graph, and and
    // array of next-node-start sides for each side in the graph. We number
    // sides as 2 * xg internal node ID, +1 if it's a right side. This leaves us
    // 0 and 1 free for representing the null destination side and to use as a
    // per-side array run separator, respectively.
    
    // This holds, for each node and edge (with indexes as in the entity vector
    // f_iv), the usage count (i.e. the number of times it is visited by encoded
    // threads). This doesn't have to be dynamic since the length will never
    // change.
    int_vector<> h_iv;
    
    // This (as an extension to the algorithm described in the paper) holds the
    // number of threads beginning at each node. This isn't any extra
    // information relative to what's in the usage count array, but it's cheaper
    // (probably) to maintain this rather than to scan through all the edges on
    // a side every time.
    // ts stands for "thread start"
    int_vector<> ts_iv;
    
    // This holds the concatenated Benedict arrays (holding the next node's
    // entry side visited after the current entry side for the threads sorted in
    // reverse prefix lexicographic order). They are separated with 1s, with 0s
    // noting the null side (i.e. the thread ends at this node).
    dynamic_int_vector bs_iv;
    
    // We need the w function, which we call the "where_to" function. It tells
    // you, from a given visit at a given side, what visit offset if you go to
    // another side.
    int64_t where_to(int64_t current_side, int64_t visit_offset, int64_t new_side);
};

class XGPath {
public:
    XGPath(void) : member_count(0) { }
    ~XGPath(void) { }
    XGPath(const string& path_name,
           const vector<Mapping>& path,
           size_t entity_count,
           XG& graph,
           const map<int64_t, string>& node_label);
    string name;
    size_t member_count;
    sd_vector<> members;
    wt_int<> ids;
    sd_vector<> directions; // forward or backward through nodes
    int_vector<> positions;
    int_vector<> ranks;
    bit_vector offsets;
    rank_support_v<1> offsets_rank;
    bit_vector::select_1_type offsets_select;
    void load(istream& in);
    size_t serialize(std::ostream& out,
                     sdsl::structure_tree_node* v = NULL,
                     std::string name = "");
    Mapping mapping(size_t offset); // 0-based
};


Mapping new_mapping(const string& name, int64_t id, size_t rank, bool is_reverse);
void parse_region(const string& target, string& name, int64_t& start, int64_t& end);
void to_text(ostream& out, Graph& graph);

}

#endif
