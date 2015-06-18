#include "succinct.hpp"
#include "stream.hpp"

namespace scg {

int dna3bit(char c) {
    switch (c) {
    case 'A':
        return 0;
    case 'T':
        return 1;
    case 'C':
        return 2;
    case 'G':
        return 3;
    default:
        return 4;
    }
}

char revdna3bit(int i) {
    switch (i) {
    case 0:
        return 'A';
    case 1:
        return 'T';
    case 2:
        return 'C';
    case 3:
        return 'G';
    default:
        return 'N';
    }
}

SuccinctGraph::SuccinctGraph(istream& in) {
    //, size_t seq_length, size_t node_count, size_t edge_count) {
    // allocate space for construction
    seq_length=0; node_count=0; edge_count=0; path_entry_count=0;
    // we can always resize smaller, but not easily extend
    function<void(Graph&)> lambda = [this](Graph& graph) {
        for (int i = 0; i < graph.node_size(); ++i) {
            const Node& n = graph.node(i);
            if (node_label.find(n.id()) == node_label.end()) {
                ++node_count;
                seq_length += n.sequence().size();
                node_label[n.id()] = n.sequence();
            }
        }
        for (int i = 0; i < graph.edge_size(); ++i) {
            const Edge& e = graph.edge(i);
            if (from_to.find(e.from()) == from_to.end() || from_to[e.from()].count(e.to()) ==0) {
                ++edge_count;
                from_to[e.from()].insert(e.to());
                to_from[e.to()].insert(e.from());
            }
        }
        for (int i = 0; i < graph.path_size(); ++i) {
            const Path& p = graph.path(i);
            const string& name = p.name();
            for (int j = 0; j < p.mapping_size(); ++j) {
                const Mapping& m = p.mapping(j);
                path_nodes[name].insert(m.position().node_id());
                ++path_entry_count;
            }
        }
    };
    stream::for_each(in, lambda);

    size_t entity_count = node_count + edge_count;
    cerr << "graph has " << seq_length << "bp in sequence, "
         << node_count << " nodes, "
         << edge_count << " edges, and "
         << path_entry_count << " nodes in paths" << endl;

    // set up our compressed representation
    util::assign(s_iv, int_vector<>(seq_length, 0, 3));
    util::assign(s_bv, bit_vector(seq_length));
    util::assign(i_iv, int_vector<>(node_count));
    util::assign(f_iv, int_vector<>(entity_count));
    util::assign(f_bv, bit_vector(entity_count));
    util::assign(t_iv, int_vector<>(entity_count));
    util::assign(t_bv, bit_vector(entity_count));

    // for each node in the sequence
    // concatenate the labels into the s_iv
    size_t i = 0; // insertion point
    size_t r = 1;
    for (auto& p : node_label) {
        int64_t id = p.first;
        const string& l = p.second;
        s_bv[i] = 1; // record node start
        i_iv[r-1] = id;
        ++r;
        for (auto c : l) {
            s_iv[i++] = dna3bit(c); // store sequence
        }
    }
    //node_label.clear();

    // we have to process all the nodes before we do the edges
    // because we need to ensure full coverage of node space

    util::bit_compress(i_iv);
    enc_vector<> i_civ(i_iv);
    wt_int<> i_wt;
    construct_im(i_wt, i_iv);

    size_t f_itr = 0;
    for (size_t k = 0; k < node_count; ++k) {
        //cerr << k << endl;
        int64_t f_id = i_iv[k];
        size_t f_rank = k+1;
        f_iv[f_itr] = f_rank;
        f_bv[f_itr] = 1;
        ++f_itr;
        if (from_to.find(f_id) != from_to.end()) {
            for (auto& t_id : from_to[f_id]) {
                size_t t_rank = i_wt.select(1, t_id)+1;
                f_iv[f_itr] = t_rank;
                f_bv[f_itr] = 0;
                ++f_itr;
            }
        }
    }

    size_t t_itr = 0;
    for (size_t k = 0; k < i_iv.size(); ++k) {
        int64_t t_id = i_iv[k];
        size_t t_rank = k+1;
        t_iv[t_itr] = t_rank;
        t_bv[t_itr] = 1;
        ++t_itr;
        if (to_from.find(t_id) != to_from.end()) {
            for (auto& f_id : to_from[t_id]) {
                size_t f_rank = i_wt.select(1, f_id)+1;
                t_iv[t_itr] = f_rank;
                t_bv[t_itr] = 0;
                ++t_itr;
            }
        }
    }

    // to label the paths we'll need to compress and index our vectors
    util::bit_compress(s_iv);
    util::bit_compress(f_iv);
    util::bit_compress(t_iv);

    rank_support_v<1> s_bv_rank(&s_bv);
    bit_vector::select_1_type s_bv_select(&s_bv);
    rank_support_v<1> f_bv_rank(&f_bv);
    bit_vector::select_1_type f_bv_select(&f_bv);
    rank_support_v<1> t_bv_rank(&t_bv);
    bit_vector::select_1_type t_bv_select(&t_bv);

    
    // compressed vectors of the above
    //vlc_vector<> s_civ(s_iv);
    rrr_vector<> s_cbv(s_bv);
    rrr_vector<>::rank_1_type s_cbv_rank(&s_cbv);
    rrr_vector<>::select_1_type s_cbv_select(&s_cbv);

    /*
    vlc_vector<> f_civ(f_iv);
    rrr_vector<> f_cbv(f_bv);
    rrr_vector<>::rank_1_type f_cbv_rank(&f_cbv);
    rrr_vector<>::select_1_type f_cbv_select(&f_cbv);
    
    vlc_vector<> t_civ(t_iv);
    rrr_vector<> t_cbv(t_bv);
    rrr_vector<>::rank_1_type t_cbv_rank(&t_cbv);
    rrr_vector<>::select_1_type t_cbv_select(&t_cbv);
    */
    //map<string, sd_vector<> > p_cbv;
    
    
    cerr << "|s_iv| = " << size_in_mega_bytes(s_iv) << endl;
    cerr << "|i_iv| = " << size_in_mega_bytes(i_iv) << endl;
    cerr << "|f_iv| = " << size_in_mega_bytes(f_iv) << endl;
    cerr << "|t_iv| = " << size_in_mega_bytes(t_iv) << endl;

    cerr << "|s_bv| = " << size_in_mega_bytes(s_bv) << endl;
    cerr << "|f_bv| = " << size_in_mega_bytes(f_bv) << endl;
    cerr << "|t_bv| = " << size_in_mega_bytes(t_bv) << endl;

    //cerr << "|s_civ| = " << size_in_mega_bytes(s_civ) << endl;
    cerr << "|i_civ| = " << size_in_mega_bytes(i_civ) << endl;
    cerr << "|i_wt| = " << size_in_mega_bytes(i_wt) << endl;
    //cerr << "|f_civ| = " << size_in_mega_bytes(f_civ) << endl;
    //cerr << "|t_civ| = " << size_in_mega_bytes(t_civ) << endl;

    cerr << "|s_cbv| = " << size_in_mega_bytes(s_cbv) << endl;
    //cerr << "|f_cbv| = " << size_in_mega_bytes(f_cbv) << endl;
    //cerr << "|t_cbv| = " << size_in_mega_bytes(t_cbv) << endl;

    //cerr << s_civ << endl;
    /*
    for (int i = 0; i < s_civ.size(); ++i) {
        cerr << (char) s_civ[i];
    } cerr << endl;
    cerr << s_cbv << endl;
    cerr << f_civ << endl;
    cerr << f_cbv << endl;
    cerr << t_civ << endl;
    cerr << t_cbv << endl;
    cerr << i_civ << endl;
    */

    int max_id = s_cbv_rank(s_cbv.size());
    for (auto& p : node_label) {
        int64_t id = p.first;
        const string& l = p.second;
        //size_t rank = node_rank[id];
        size_t rank = i_wt.select(1, id)+1;
        //cerr << rank << endl;
        // find the node in the array
        //cerr << "id = " << id << " rank = " << s_cbv_select(rank) << endl;
        // this should be true given how we constructed things
        assert(rank == s_cbv_rank(s_cbv_select(rank)+1));
        // get the sequence from the s_iv
        size_t start = s_cbv_select(rank);
        size_t end = rank == max_id ? s_cbv.size() : s_cbv_select(rank+1);
        string s; s.resize(end-start);
        for (size_t i = start; i < s_cbv.size() && i < end; ++i) {
            s[i-start] = revdna3bit(s_iv[i]);
        }
        //cerr << id << " " << l << " =? " << s << endl;
        assert(l == s);
    }

    // todo, why -1?
    for (size_t j = 0; j < f_iv.size()-1; ++j) {
        //cerr << j << endl;
        if (f_bv[j] == 1) continue;
        // from id == rank
        size_t fid = i_civ[f_bv_rank(j)-1];
        // to id == f_cbv[j]
        size_t tid = i_civ[f_iv[j]-1];
        //cerr << fid << " " << tid << endl;
        assert(from_to[fid].count(tid));
    }
    
    for (size_t j = 0; j < t_iv.size()-1; ++j) {
        //cerr << j << endl;
        if (t_bv[j] == 1) continue;
        // from id == rank
        size_t tid = i_civ[t_bv_rank(j)-1];
        // to id == f_cbv[j]
        size_t fid = i_civ[t_iv[j]-1];
        //cerr << tid << " " << fid << endl;
        assert(to_from[tid].count(fid));
    }

    cerr << "graph ok" << endl;


}

}
