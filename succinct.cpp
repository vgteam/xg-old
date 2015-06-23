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
    load(in);
}

void SuccinctGraph::load(istream& in) {

    size_t sequence_length;
    size_t node_count;
    size_t edge_count;
    size_t path_count;
    
    sdsl::read_member(sequence_length, in);
    sdsl::read_member(node_count, in);
    sdsl::read_member(edge_count, in);
    size_t entity_count = node_count + edge_count;
    //cerr << sequence_length << ", " << node_count << ", " << edge_count << endl;

    i_iv.load(in);
    i_wt.load(in);

    s_iv.load(in);
    s_cbv.load(in);
    s_cbv_rank.load(in, &s_cbv);
    s_cbv_select.load(in, &s_cbv);

    f_iv.load(in);
    f_bv.load(in);
    f_bv_rank.load(in, &f_bv);
    f_bv_select.load(in, &f_bv);

    t_iv.load(in);
    t_bv.load(in);
    t_bv_rank.load(in, &t_bv);
    t_bv_select.load(in, &t_bv);

    pn_iv.load(in);
    pn_csa.load(in);
    pn_bv.load(in);
    pn_bv_rank.load(in, &pn_bv);
    pn_bv_select.load(in, &pn_bv);
    pi_iv.load(in);
    sdsl::read_member(path_count, in);
    for (size_t i = 0; i < path_count; ++i) {
        bit_vector bv;
        bv.load(in);
        pe_v.push_back(bv);
        int_vector<> iv;
        iv.load(in);
        pp_v.push_back(iv);
    }
}

size_t SuccinctGraph::serialize(ostream& out, sdsl::structure_tree_node* s, std::string name) {

    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
    size_t written = 0;

    written += sdsl::write_member(s_iv.size(), out, child, "sequence_length");
    written += sdsl::write_member(i_iv.size(), out, child, "node_count");
    written += sdsl::write_member(f_iv.size()-i_iv.size(), out, child, "edge_count");
    //written += sdsl::write_member(path_count, out, child, "path_count");

    written += i_iv.serialize(out, child, "id_vector");
    written += i_wt.serialize(out, child, "id_wavelet_tree");

    written += s_iv.serialize(out, child, "seq_vector");
    written += s_cbv.serialize(out, child, "seq_node_starts");
    written += s_cbv_rank.serialize(out, child, "seq_node_starts_rank");
    written += s_cbv_select.serialize(out, child, "seq_node_starts_select");

    written += f_iv.serialize(out, child, "from_vector");
    written += f_bv.serialize(out, child, "from_node_starts");
    written += f_bv_rank.serialize(out, child, "from_node_starts_rank");
    written += f_bv_select.serialize(out, child, "from_node_starts_select");
    
    written += t_iv.serialize(out, child, "to_vector");
    written += t_bv.serialize(out, child, "to_node_starts");
    written += t_bv_rank.serialize(out, child, "to_node_starts_rank");
    written += t_bv_select.serialize(out, child, "to_node_starts_select");

    written += pn_iv.serialize(out, child, "path_names");
    written += pn_csa.serialize(out, child, "path_names_csa");
    written += pn_bv.serialize(out, child, "path_names_starts");
    written += pn_bv_rank.serialize(out, child, "path_names_starts_rank");
    written += pn_bv_select.serialize(out, child, "path_names_starts_select");
    written += pi_iv.serialize(out, child, "path_ids");
    written += sdsl::write_member(pe_v.size(), out, child, "path_count");    
    for (size_t i = 0; i < pe_v.size(); ++i) {
        written += pe_v[i].serialize(out, child, "path_membership_" + path_name(i+1));
        written += pp_v[i].serialize(out, child, "path_positions_" + path_name(i+1));
    }

    sdsl::structure_tree::add_size(child, written);
    return written;
    
}

void SuccinctGraph::from_vg(istream& in) {

    // temporaries for construction
    map<int64_t, string> node_label;
    map<int64_t, set<int64_t> > from_to;
    map<int64_t, set<int64_t> > to_from;
    map<string, vector<int64_t> > path_nodes;
    size_t seq_length = 0;
    size_t node_count = 0;
    size_t edge_count = 0;
    size_t path_entry_count = 0;
    
    // we can always resize smaller, but not easily extend
    function<void(Graph&)> lambda = [this,
                                     &node_label,
                                     &from_to,
                                     &to_from,
                                     &path_nodes,
                                     &seq_length,
                                     &node_count,
                                     &edge_count,
                                     &path_entry_count](Graph& graph) {
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
                path_nodes[name].push_back(m.position().node_id());
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
    //util::assign(e_iv, int_vector<>(edge_count*3));

    // for each node in the sequence
    // concatenate the labels into the s_iv
    cerr << "storing node labels" << endl;
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
    construct_im(i_wt, i_iv);

    cerr << "storing forward edges and adjacency table" << endl;
    size_t f_itr = 0;
    size_t j_itr = 0; // edge adjacency pointer
    for (size_t k = 0; k < node_count; ++k) {
        //cerr << k << endl;
        int64_t f_id = i_iv[k];
        size_t f_rank = k+1;
        f_iv[f_itr] = f_rank;
        f_bv[f_itr] = 1;
        ++f_itr;
        if (from_to.find(f_id) != from_to.end()) {
            //e_iv[j_itr++] = 1; // indicates adjacency record start
            for (auto& t_id : from_to[f_id]) {
                size_t t_rank = i_wt.select(1, t_id)+1;
                f_iv[f_itr] = t_rank;
                f_bv[f_itr] = 0;
                ++f_itr;
                // store edge in adjacency index
                //e_iv[j_itr++] = f_rank+1;
                //e_iv[j_itr++] = t_rank+1;
            }
        }
    }
    
    //assert(e_iv.size() == edge_count*3);

    cerr << "storing reverse edges" << endl;
    size_t t_itr = 0;
    for (size_t k = 0; k < node_count; ++k) {
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

    /*
    csa_wt<wt_int<rrr_vector<63>>> csa;
    int_vector<> x = {1,8,15,23,1,8,23,11,8};
    construct_im(csa, x, 8);
    cout << " i SA ISA PSI LF BWT   T[SA[i]..SA[i]-1]" << endl;
    csXprintf(cout, "%2I %2S %3s %3P %2p %3B   %:3T", csa);
    */
    /*
    cerr << "building csa of edges" << endl;
    //string edges_file = "@edges.iv";
    //store_to_file(e_iv, edges_file);
    construct_im(e_csa, e_iv, 1);
    */

    // to label the paths we'll need to compress and index our vectors
    util::bit_compress(s_iv);
    util::bit_compress(f_iv);
    util::bit_compress(t_iv);
    //util::bit_compress(e_iv);

    //construct_im(e_csa, e_iv, 8);

    util::assign(s_bv_rank, rank_support_v<1>(&s_bv));
    util::assign(s_bv_select, bit_vector::select_1_type(&s_bv));
    util::assign(f_bv_rank, rank_support_v<1>(&f_bv));
    util::assign(f_bv_select, bit_vector::select_1_type(&f_bv));
    util::assign(t_bv_rank, rank_support_v<1>(&t_bv));
    util::assign(t_bv_select, bit_vector::select_1_type(&t_bv));
    
    // compressed vectors of the above
    //vlc_vector<> s_civ(s_iv);
    util::assign(s_cbv, rrr_vector<>(s_bv));
    util::assign(s_cbv_rank, rrr_vector<>::rank_1_type(&s_cbv));
    util::assign(s_cbv_select, rrr_vector<>::select_1_type(&s_cbv));


    cerr << "storing paths" << endl;
    // paths
    //path_nodes[name].push_back(m.position().node_id());
    string path_names;
    for (auto& pathpair : path_nodes) {
        // add path name
        const string& path_name = pathpair.first;
        //cerr << path_name << endl;
        const vector<int64_t>& path = pathpair.second;
        path_names += path_name_marker + path_name;
        pe_v.emplace_back();
        // path members (of nodes and edges ordered as per f_bv)
        bit_vector& pe_bv = pe_v.back();
        util::assign(pe_bv, bit_vector(entity_count));
        // TODO resize
        pp_v.emplace_back();
        // node positions in path
        int_vector<>& pp_iv = pp_v.back();
        util::assign(pp_iv, int_vector<>(path.size()));
        // TODO resize ... to fit
        size_t path_off = 0;
        size_t pe_off = 0;
        size_t pp_off = 0;
        for (size_t i = 0; i < path.size(); ++i) {
            //cerr << i << endl;
            auto& node_id = path[i];
            //cerr << node_id << endl;
            // record node
            pe_bv[node_rank_as_entity(node_id)] = 1;
            // and record node offset in path
            pp_iv[pp_off++] = path_off;
            path_off += node_label[node_id].size();
            // find the next edge in the path, and record it
            // (if there is a next node)
            if (i+1 < path.size()) {
                auto& next_node_id = path[i+1];
                if (has_edge(node_id, next_node_id)) {
                    pe_bv[edge_rank_as_entity(node_id, next_node_id)] = 1;
                }
            }
        }
        util::bit_compress(pp_iv);
    }

    // handle path names
    util::assign(pn_iv, int_vector<>(path_names.size()));
    util::assign(pn_bv, bit_vector(path_names.size()));
    // now record path name starts
    for (size_t i = 0; i < path_names.size(); ++i) {
        pn_iv[i] = path_names[i];
        if (path_names[i] == path_name_marker) {
            pn_bv[i] = 1; // register name start
        }
    }
    util::assign(pn_bv_rank, rank_support_v<1>(&pn_bv));
    util::assign(pn_bv_select, bit_vector::select_1_type(&pn_bv));
    
    //util::bit_compress(pn_iv);
    string path_name_file = "@pathnames.iv";
    store_to_file((const char*)path_names.c_str(), path_name_file);
    construct(pn_csa, path_name_file, 1);
    
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
    //cerr << "|i_iv| = " << size_in_mega_bytes(i_iv) << endl;
    cerr << "|f_iv| = " << size_in_mega_bytes(f_iv) << endl;
    cerr << "|t_iv| = " << size_in_mega_bytes(t_iv) << endl;

    //cerr << "|s_bv| = " << size_in_mega_bytes(s_bv) << endl;
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

    long double paths_mb_size = 0;
    cerr << "|pn_iv| = " << size_in_mega_bytes(pn_iv) << endl;
    paths_mb_size += size_in_mega_bytes(pn_iv);
    cerr << "|pn_csa| = " << size_in_mega_bytes(pn_csa) << endl;
    paths_mb_size += size_in_mega_bytes(pn_csa);
    cerr << "|pn_bv| = " << size_in_mega_bytes(pn_bv) << endl;
    paths_mb_size += size_in_mega_bytes(pn_bv);
    paths_mb_size += size_in_mega_bytes(pn_bv_rank);
    paths_mb_size += size_in_mega_bytes(pn_bv_select);
    paths_mb_size += size_in_mega_bytes(pi_iv);
    for (auto& bv : pe_v) {
        paths_mb_size += size_in_mega_bytes(bv);
    }
    for (auto& iv : pp_v) {
        paths_mb_size += size_in_mega_bytes(iv);
    }
    cerr << "total paths size " << paths_mb_size << endl;
    
    cerr << "total size [MB] = " << (
        size_in_mega_bytes(s_iv)
        + size_in_mega_bytes(f_iv)
        + size_in_mega_bytes(t_iv)
        //+ size_in_mega_bytes(s_bv)
        + size_in_mega_bytes(f_bv)
        + size_in_mega_bytes(t_bv)
        + size_in_mega_bytes(i_civ)
        + size_in_mega_bytes(i_wt)
        + size_in_mega_bytes(s_cbv)
        + paths_mb_size
        ) << endl;
    
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

    cerr << "validating graph sequence" << endl;
    int max_id = s_cbv_rank(s_cbv.size());
    for (auto& p : node_label) {
        int64_t id = p.first;
        const string& l = p.second;
        //size_t rank = node_rank[id];
        size_t rank = id_to_rank(id);
        //cerr << rank << endl;
        // find the node in the array
        //cerr << "id = " << id << " rank = " << s_cbv_select(rank) << endl;
        // this should be true given how we constructed things
        if (rank != s_cbv_rank(s_cbv_select(rank)+1)) {
            cerr << rank << " != " << s_cbv_rank(s_cbv_select(rank)+1) << " for node " << id << endl;
            assert(false);
        }
        // get the sequence from the s_iv
        string s = node_sequence(id);

        string ltmp, stmp;
        if (l.size() != s.size()) {
            cerr << l << " != " << endl << s << endl << " for node " << id << endl;
            assert(false);
        } else {
            int j = 0;
            for (auto c : l) {
                if (dna3bit(c) != dna3bit(s[j++])) {
                    cerr << l << " != " << endl << s << endl << " for node " << id << endl;
                    assert(false);
                }
            }
        }
    }
    node_label.clear();

    cerr << "validating forward edge table" << endl;
    // todo, why -1?
    for (size_t j = 0; j < f_iv.size()-1; ++j) {
        //cerr << j << endl;
        if (f_bv[j] == 1) continue;
        // from id == rank
        size_t fid = i_civ[f_bv_rank(j)-1];
        // to id == f_cbv[j]
        size_t tid = i_civ[f_iv[j]-1];
        //cerr << fid << " " << tid << endl;
        if (from_to[fid].count(tid) == 0) {
            cerr << "could not find edge (f) " << fid << " -> " << tid << endl;
            assert(false);
        }
    }

    cerr << "validating reverse edge table" << endl;
    for (size_t j = 0; j < t_iv.size()-1; ++j) {
        //cerr << j << endl;
        if (t_bv[j] == 1) continue;
        // from id == rank
        size_t tid = i_civ[t_bv_rank(j)-1];
        // to id == f_cbv[j]
        size_t fid = i_civ[t_iv[j]-1];
        //cerr << tid << " " << fid << endl;
        if (to_from[tid].count(fid) == 0) {
            cerr << "could not find edge (t) " << fid << " -> " << tid << endl;
            assert(false);
        }
    }

    cerr << "validating paths" << endl;
    for (auto& pathpair : path_nodes) {
        const string& name = pathpair.first;
        const vector<int64_t>& path = pathpair.second;
        size_t prank = path_rank(name);
        assert(path_name(prank) == name);
        bit_vector& pe_bv = pe_v[prank-1];
        int_vector<>& pp_iv = pp_v[prank-1];
        // check each entity in the nodes is present
        for (auto& id : path) {
            assert(pe_bv[node_rank_as_entity(id)]);
        }
        //cerr << path_name << " rank = " << prank << endl;
        // check membership now for each entity in the path
    }

    cerr << "graph ok" << endl;

}

Node SuccinctGraph::node(int64_t id) {
    Node n;
    n.set_id(id);
    n.set_sequence(node_sequence(id));
    return n;
}

string SuccinctGraph::node_sequence(int64_t id) {
    size_t rank = id_to_rank(id);
    size_t start = s_cbv_select(rank);
    size_t end = rank == max_rank() ? s_cbv.size() : s_cbv_select(rank+1);
    string s; s.resize(end-start);
    for (size_t i = start; i < s_cbv.size() && i < end; ++i) {
        s[i-start] = revdna3bit(s_iv[i]);
    }
    return s;
}

size_t SuccinctGraph::id_to_rank(int64_t id) {
    return i_wt.select(1, id)+1;
}

int64_t SuccinctGraph::rank_to_id(size_t rank) {
    return i_iv[rank];
}

vector<Edge> SuccinctGraph::edges_to(int64_t id) {
    size_t rank = id_to_rank(id);
    size_t t_start = t_bv_select(rank);
    size_t t_end = t_bv_select(rank+1);
    vector<Edge> edges;
    for (size_t i = t_start; i < t_end; ++i) {
        Edge edge;
        edge.set_from(rank_to_id(t_iv[i]));
        edge.set_to(id);
        edges.push_back(edge);
    }
    return edges;
}

vector<Edge> SuccinctGraph::edges_from(int64_t id) {
    size_t rank = id_to_rank(id);
    size_t f_start = f_bv_select(rank);
    size_t f_end = f_bv_select(rank+1);
    vector<Edge> edges;
    for (size_t i = f_start; i < f_end; ++i) {
        Edge edge;
        edge.set_from(id);
        edge.set_to(rank_to_id(f_iv[i]));
        edges.push_back(edge);
    }
    return edges;
}

int64_t SuccinctGraph::max_rank(void) {
    return s_cbv_rank(s_cbv.size());
}

size_t SuccinctGraph::node_rank_as_entity(int64_t id) {
    //cerr << id_to_rank(id) << endl;
    return f_bv_select(id_to_rank(id));
}

// snoop through the forward table to check if the edge exists
bool SuccinctGraph::has_edge(int64_t id1, int64_t id2) {
    size_t rank1 = id_to_rank(id1);
    size_t rank2 = id_to_rank(id2);
    size_t f_start = f_bv_select(rank1);
    size_t f_end = f_bv_select(rank1+1);
    for (size_t i = f_start; i < f_end; ++i) {
        if (rank2 == f_iv[i]) {
            return true;
        }
    }
    return false;
}

size_t SuccinctGraph::edge_rank_as_entity(int64_t id1, int64_t id2) {
    size_t rank1 = id_to_rank(id1);
    size_t rank2 = id_to_rank(id2);
    //cerr << "Finding rank for " << id1 << "(" << rank1 << ") " << " -> " << id2 << "(" << rank2 << ")"<< endl;
    size_t f_start = f_bv_select(rank1);
    size_t f_end = f_bv_select(rank1+1);
    //cerr << f_start << " to " << f_end << endl;
    for (size_t i = f_start; i < f_end; ++i) {
        //cerr << f_iv[i] << endl;
        if (rank2 == f_iv[i]) {
            return i;
        }
    }
    cerr << "edge does not exist: " << id1 << " -> " << id2 << endl;
    assert(false);
}

size_t SuccinctGraph::path_rank(const string& name) {
    // find the name in the csa
    auto occs = locate(pn_csa, name);
    if (occs.size() > 1) {
        cerr << "multiple hits for " << name << endl;
        assert(false);
    }
    return pn_bv_rank(occs[0]);
}

string SuccinctGraph::path_name(size_t rank) {
    size_t start = pn_bv_select(rank)+1; // step past '#'
    size_t end = (rank == pn_bv_rank(pn_bv.size()-1)) ? pn_bv.size() : pn_bv_select(rank+1);
    string name; name.resize(end-start);
    for (size_t i = start; i < end; ++i) {
        name[i-start] = pn_iv[i];
    }
    return name;
}

/*
Path SuccinctGraph::path(string& name) {
}
Graph SuccinctGraph::neighborhood(int64_t rank, int32_t steps) {
}
Graph SuccinctGraph::range(int64_t rank1, int64_t rank2) {
}
Graph SuccinctGraph::region(string& path_name, int64_t start, int64_t stop) {
}
*/

}
