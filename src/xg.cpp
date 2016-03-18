#include "xg.hpp"
#include "stream.hpp"

namespace xg {

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

XG::XG(istream& in)
    : start_marker('#'),
      end_marker('$'),
      seq_length(0),
      node_count(0),
      edge_count(0),
      path_count(0) {
    load(in);
}

XG::XG(Graph& graph)
    : start_marker('#'),
      end_marker('$'),
      seq_length(0),
      node_count(0),
      edge_count(0),
      path_count(0) {
    from_graph(graph);
}

XG::XG(function<void(function<void(Graph&)>)> get_chunks)
    : start_marker('#'),
      end_marker('$'),
      seq_length(0),
      node_count(0),
      edge_count(0),
      path_count(0) {
    from_callback(get_chunks);
}

XG::~XG(void) {
     if (bs_iv != nullptr) delete bs_iv;
}

void XG::load(istream& in) {

    if (!in.good()) {
        cerr << "[xg] error: index does not exist!" << endl;
        exit(1);
    }

    sdsl::read_member(seq_length, in);
    sdsl::read_member(node_count, in);
    sdsl::read_member(edge_count, in);
    sdsl::read_member(path_count, in);
    size_t entity_count = node_count + edge_count;
    //cerr << sequence_length << ", " << node_count << ", " << edge_count << endl;
    sdsl::read_member(min_id, in);
    sdsl::read_member(max_id, in);

    i_iv.load(in);
    r_iv.load(in);

    s_iv.load(in);
    s_cbv.load(in);
    s_cbv_rank.load(in, &s_cbv);
    s_cbv_select.load(in, &s_cbv);

    f_iv.load(in);
    f_bv.load(in);
    f_bv_rank.load(in, &f_bv);
    f_bv_select.load(in, &f_bv);
    f_from_start_cbv.load(in);
    f_to_end_cbv.load(in);

    t_iv.load(in);
    t_bv.load(in);
    t_bv_rank.load(in, &t_bv);
    t_bv_select.load(in, &t_bv);
    t_to_end_cbv.load(in);
    t_from_start_cbv.load(in);

    pn_iv.load(in);
    pn_csa.load(in);
    pn_bv.load(in);
    pn_bv_rank.load(in, &pn_bv);
    pn_bv_select.load(in, &pn_bv);
    pi_iv.load(in);
    sdsl::read_member(path_count, in);
    for (size_t i = 0; i < path_count; ++i) {
        auto path = new XGPath;
        path->load(in);
        paths.push_back(path);
    }
    ep_iv.load(in);
    ep_bv.load(in);
    ep_bv_rank.load(in, &ep_bv);
    ep_bv_select.load(in, &ep_bv);
    
    h_iv.load(in);
    ts_iv.load(in);
    bs_iv = deserialize(in);
}

void XGPath::load(istream& in) {
    members.load(in);
    ids.load(in);
    directions.load(in);
    ranks.load(in);
    positions.load(in);
    offsets.load(in);
    offsets_rank.load(in, &offsets);
    offsets_select.load(in, &offsets);
}

size_t XGPath::serialize(std::ostream& out,
                         sdsl::structure_tree_node* v,
                         std::string name) {
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
    size_t written = 0;
    written += members.serialize(out, child, "path_membership_" + name);
    written += ids.serialize(out, child, "path_node_ids_" + name);
    written += directions.serialize(out, child, "path_node_directions_" + name);
    written += ranks.serialize(out, child, "path_mapping_ranks_" + name);
    written += positions.serialize(out, child, "path_node_offsets_" + name);
    written += offsets.serialize(out, child, "path_node_starts_" + name);
    written += offsets_rank.serialize(out, child, "path_node_starts_rank_" + name);
    written += offsets_select.serialize(out, child, "path_node_starts_select_" + name);
    return written;
}

XGPath::XGPath(const string& path_name,
               const vector<Mapping>& path,
               size_t entity_count,
               XG& graph,
               const map<int64_t, string>& node_label) {

    name = path_name;
    member_count = 0;
    
    // path members (of nodes and edges ordered as per f_bv)
    bit_vector members_bv;
    util::assign(members_bv, bit_vector(entity_count));
    // node ids, the literal path
    int_vector<> ids_iv;
    util::assign(ids_iv, int_vector<>(path.size()));
    // directions of traversal (typically forward, but we allow backwards
    bit_vector directions_bv;
    util::assign(directions_bv, bit_vector(path.size()));
    // node positions in path
    util::assign(positions, int_vector<>(path.size()));
    // mapping ranks in path
    util::assign(ranks, int_vector<>(path.size()));

    size_t path_off = 0;
    size_t members_off = 0;
    size_t positions_off = 0;
    size_t path_length = 0;

    // determine total length
    for (size_t i = 0; i < path.size(); ++i) {
        auto node_id = path[i].position().node_id();
        auto label_itr = node_label.find(node_id);
        if (label_itr != node_label.end()) {
            path_length += label_itr->second.size();
            ids_iv[i] = node_id;
        } else {
            cerr << "[xg] error: when making paths could not find node label for " << node_id << endl;
            assert(false);
        }
    }

    // make the bitvector for path offsets
    util::assign(offsets, bit_vector(path_length));
    set<int64_t> uniq_nodes;
    set<pair<pair<int64_t, bool>, pair<int64_t, bool>>> uniq_edges;
    //cerr << "path " << path_name << " has " << path.size() << endl;
    for (size_t i = 0; i < path.size(); ++i) {
        //cerr << i << endl;
        auto& mapping = path[i];
        auto node_id = mapping.position().node_id();
        bool is_reverse = mapping.position().is_reverse();
        //cerr << node_id << endl;
        // record node
        members_bv[graph.node_rank_as_entity(node_id)-1] = 1;
        // record direction of passage through node
        directions_bv[i] = mapping.position().is_reverse();
        // and the external rank of the mapping
        ranks[i] = mapping.rank();
        // we've seen another entity
        uniq_nodes.insert(node_id);
        // and record node offset in path
        positions[positions_off++] = path_off;
        // record position of node
        offsets[path_off] = 1;
        // and update the offset counter
        auto label_itr = node_label.find(node_id);
        if (label_itr != node_label.end()) {
            path_off += label_itr->second.size();
        } else {
            cerr << "[xg] error: when recording offsets could not find node label for " << node_id << endl;
            assert(false);
        }

        // find the next edge in the path, and record it
        if (i+1 < path.size()) { // but only if there is a next node
            auto next_node_id = path[i+1].position().node_id();
            bool next_is_reverse = path[i+1].position().is_reverse();
            //cerr << "checking if we have the edge" << endl;
            int64_t id1, id2;
            bool rev1, rev2;
            if (is_reverse && next_is_reverse) {
                id1 = next_node_id; id2 = node_id;
                rev1 = false; rev2 = false;
            } else {
                id1 = node_id; id2 = next_node_id;
                rev1 = is_reverse; rev2 = next_is_reverse;
            }
            if (graph.has_edge(id1, rev1, id2, rev2)) {
                members_bv[graph.edge_rank_as_entity(id1, rev1, id2, rev2)-1] = 1;
                uniq_edges.insert(
                    make_pair(
                        make_pair(id1, rev1),
                        make_pair(id2, rev2)));
            } else if (graph.has_edge(id2, !rev2, id1, !rev1)) {
                members_bv[graph.edge_rank_as_entity(id2, !rev2, id1, !rev1)-1] = 1;
                uniq_edges.insert(
                    make_pair(
                        make_pair(id2, !rev2),
                        make_pair(id1, !rev1)
                        ));
            } else {
                cerr << "graph does not have edge from "
                     << node_id << (path[i].position().is_reverse()?"+":"-")
                     << " to "
                     << next_node_id << (path[i+1].position().is_reverse()?"-":"+")
                     << endl;
            }
        }
    }
    // set member count as the unique entities that are in the path
    //cerr << uniq_nodes.size() << " vs " << path.size() << endl;
    member_count = uniq_nodes.size() + uniq_edges.size();
    // compress path membership vectors
    util::assign(members, sd_vector<>(members_bv));
    // and traversal information
    util::assign(directions, sd_vector<>(directions_bv));
    // handle entity lookup structure (wavelet tree)
    util::bit_compress(ids_iv);
    construct_im(ids, ids_iv);
    // bit compress the positional offset info
    util::bit_compress(positions);
    // bit compress mapping ranks
    util::bit_compress(ranks);

    util::assign(offsets_rank, rank_support_v<1>(&offsets));
    util::assign(offsets_select, bit_vector::select_1_type(&offsets));
}

Mapping XGPath::mapping(size_t offset) {
    // TODO actually store the "real" mapping
    Mapping m;
    // store the starting position and series of edits
    m.mutable_position()->set_node_id(ids[offset]);
    m.mutable_position()->set_is_reverse(directions[offset]);
    m.set_rank(ranks[offset]);
    return m;
}

size_t XG::serialize(ostream& out, sdsl::structure_tree_node* s, std::string name) {

    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(s, name, sdsl::util::class_name(*this));
    size_t written = 0;

    written += sdsl::write_member(s_iv.size(), out, child, "sequence_length");
    written += sdsl::write_member(i_iv.size(), out, child, "node_count");
    written += sdsl::write_member(f_iv.size()-i_iv.size(), out, child, "edge_count");
    written += sdsl::write_member(path_count, out, child, "path_count");
    written += sdsl::write_member(min_id, out, child, "min_id");
    written += sdsl::write_member(max_id, out, child, "max_id");

    written += i_iv.serialize(out, child, "id_rank_vector");
    written += r_iv.serialize(out, child, "rank_id_vector");

    written += s_iv.serialize(out, child, "seq_vector");
    written += s_cbv.serialize(out, child, "seq_node_starts");
    written += s_cbv_rank.serialize(out, child, "seq_node_starts_rank");
    written += s_cbv_select.serialize(out, child, "seq_node_starts_select");

    written += f_iv.serialize(out, child, "from_vector");
    written += f_bv.serialize(out, child, "from_node");
    written += f_bv_rank.serialize(out, child, "from_node_rank");
    written += f_bv_select.serialize(out, child, "from_node_select");
    written += f_from_start_cbv.serialize(out, child, "from_is_from_start");
    written += f_to_end_cbv.serialize(out, child, "from_is_to_end");
    
    written += t_iv.serialize(out, child, "to_vector");
    written += t_bv.serialize(out, child, "to_node");
    written += t_bv_rank.serialize(out, child, "to_node_rank");
    written += t_bv_select.serialize(out, child, "to_node_select");
    written += t_to_end_cbv.serialize(out, child, "to_is_to_end");
    written += t_from_start_cbv.serialize(out, child, "to_is_from_start");

    written += pn_iv.serialize(out, child, "path_names");
    written += pn_csa.serialize(out, child, "path_names_csa");
    written += pn_bv.serialize(out, child, "path_names_starts");
    written += pn_bv_rank.serialize(out, child, "path_names_starts_rank");
    written += pn_bv_select.serialize(out, child, "path_names_starts_select");
    written += pi_iv.serialize(out, child, "path_ids");
    written += sdsl::write_member(paths.size(), out, child, "path_count");    
    for (auto path : paths) {
        written += path->serialize(out, child, path->name);
    }
    
    written += ep_iv.serialize(out, child, "entity_path_mapping");
    written += ep_bv.serialize(out, child, "entity_path_mapping_starts");
    written += ep_bv_rank.serialize(out, child, "entity_path_mapping_starts_rank");
    written += ep_bv_select.serialize(out, child, "entity_path_mapping_starts_select");

    written += h_iv.serialize(out, child, "thread_usage_count");
    written += ts_iv.serialize(out, child, "thread_start_count");
    written += xg::serialize(bs_iv, out, child, "benedict_arrays");

    sdsl::structure_tree::add_size(child, written);
    return written;
    
}

void XG::from_stream(istream& in, bool validate_graph, bool print_graph) {

    from_callback([&](function<void(Graph&)> handle_chunk) {
        // TODO: should I be bandying about function references instead of
        // function objects here?
        stream::for_each(in, handle_chunk);
    }, validate_graph, print_graph);
}

void XG::from_graph(Graph& graph, bool validate_graph, bool print_graph) {

    from_callback([&](function<void(Graph&)> handle_chunk) {
        // There's only one chunk in this case.
        handle_chunk(graph);
    }, validate_graph, print_graph);

}

void XG::from_callback(function<void(function<void(Graph&)>)> get_chunks, 
    bool validate_graph, bool print_graph) {

    // temporaries for construction
    map<int64_t, string> node_label;
    // need to store node sides
    map<Side, set<Side> > from_to;
    map<Side, set<Side> > to_from;
    map<string, map<int, Mapping> > path_nodes;

    // This takes in graph chunks and adds them into our temporary storage.
    function<void(Graph&)> lambda = [this,
                                     &node_label,
                                     &from_to,
                                     &to_from,
                                     &path_nodes](Graph& graph) {
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
            if (from_to.find(Side(e.from(), e.from_start())) == from_to.end()
                || from_to[Side(e.from(), e.from_start())].count(Side(e.to(), e.to_end())) == 0) {
                ++edge_count;
                from_to[Side(e.from(), e.from_start())].insert(Side(e.to(), e.to_end()));
                to_from[Side(e.to(), e.to_end())].insert(Side(e.from(), e.from_start()));
            }
        }

        // Print out all the paths in the graph we are loading
        for (int i = 0; i < graph.path_size(); ++i) {
            const Path& p = graph.path(i);
            const string& name = p.name();
#ifdef VERBOSE_DEBUG
            cerr << "Path " << name << ": ";
#endif
            for (int j = 0; j < p.mapping_size(); ++j) {
                const Mapping& m = p.mapping(j);
                path_nodes[name][m.rank()] = m;
#ifdef VERBOSE_DEBUG
                cerr << m.position().node_id() * 2 + m.position().is_reverse() << "; ";
#endif
            }
#ifdef VERBOSE_DEBUG
            cerr << endl;
#endif
        }
    };
    
    // Get all the chunks via the callback, and have them called back to us.
    // The other end handles figuring out how much to loop.
    get_chunks(lambda);
    
    path_count = path_nodes.size();

    build(node_label, from_to, to_from, path_nodes, validate_graph, print_graph);
    
}

void XG::build(map<int64_t, string>& node_label,
               map<Side, set<Side> >& from_to,
               map<Side, set<Side> >& to_from,
               map<string, map<int, Mapping>>& path_nodes,
               bool validate_graph,
               bool print_graph) {

    size_t entity_count = node_count + edge_count;
#ifdef VERBOSE_DEBUG
    cerr << "graph has " << seq_length << "bp in sequence, "
         << node_count << " nodes, "
         << edge_count << " edges, and "
         << path_count << " paths "
         << "for a total of " << entity_count << " entities" << endl;
#endif

    // for mapping of ids to ranks using a vector rather than wavelet tree
    min_id = node_label.begin()->first;
    max_id = node_label.rbegin()->first;
    
    // set up our compressed representation
    util::assign(s_iv, int_vector<>(seq_length, 0, 3));
    util::assign(s_bv, bit_vector(seq_length));
    util::assign(i_iv, int_vector<>(node_count));
    util::assign(r_iv, int_vector<>(max_id-min_id+1)); // note possibly discontiguous
    util::assign(f_iv, int_vector<>(entity_count));
    util::assign(f_bv, bit_vector(entity_count));
    util::assign(f_from_start_bv, bit_vector(entity_count));
    util::assign(f_to_end_bv, bit_vector(entity_count));
    util::assign(t_iv, int_vector<>(entity_count));
    util::assign(t_bv, bit_vector(entity_count));
    util::assign(t_to_end_bv, bit_vector(entity_count));
    util::assign(t_from_start_bv, bit_vector(entity_count));
    
    // Prepare empty vectors for path indexing
#ifdef VERBOSE_DEBUG
    cerr << "creating empty succinct thread store" << endl;
#endif
    util::assign(h_iv, int_vector<>(entity_count * 2, 0));
    util::assign(ts_iv, int_vector<>((node_count + 1) * 2, 0));
    bs_iv = new dynamic_int_vector;
    for(int64_t i = 0; i < node_count * 2; i++) {
        // Add in a separator marking the start of the B_s[] array for each
        // side. TODO: can we make the compressed representation expose a batch
        // append API?
        bs_iv->push_back(BS_SEPARATOR);
    }

    // for each node in the sequence
    // concatenate the labels into the s_iv
#ifdef VERBOSE_DEBUG
    cerr << "storing node labels" << endl;
#endif
    size_t i = 0; // insertion point
    size_t r = 1;
    for (auto& p : node_label) {
        int64_t id = p.first;
        const string& l = p.second;
        s_bv[i] = 1; // record node start
        i_iv[r-1] = id;
        // store ids to rank mapping
        r_iv[id-min_id] = r;
        ++r;
        for (auto c : l) {
            s_iv[i++] = dna3bit(c); // store sequence
        }
    }
    //node_label.clear();

    // we have to process all the nodes before we do the edges
    // because we need to ensure full coverage of node space

    util::bit_compress(i_iv);
    util::bit_compress(r_iv);

#ifdef VERBOSE_DEBUG    
    cerr << "storing forward edges and adjacency table" << endl;
#endif
    size_t f_itr = 0;
    size_t j_itr = 0; // edge adjacency pointer
    for (size_t k = 0; k < node_count; ++k) {
        int64_t f_id = i_iv[k];
        size_t f_rank = k+1;
        f_iv[f_itr] = f_rank;
        f_bv[f_itr] = 1;
        ++f_itr;
        for (auto end : { false, true }) {
            if (from_to.find(Side(f_id, end)) != from_to.end()) {
                auto t_side_itr = from_to.find(Side(f_id, end));
                if (t_side_itr != from_to.end()) {
                    for (auto& t_side : t_side_itr->second) {
                        size_t t_rank = id_to_rank(t_side.first);
                        // store link
                        f_iv[f_itr] = t_rank;
                        f_bv[f_itr] = 0;
                        // store side for start of edge
                        f_from_start_bv[f_itr] = end;
                        f_to_end_bv[f_itr] = t_side.second;
                        ++f_itr;
                    }
                }
            }
        }
    }

    // compress the forward direction side information
    util::assign(f_from_start_cbv, sd_vector<>(f_from_start_bv));
    util::assign(f_to_end_cbv, sd_vector<>(f_to_end_bv));
    
    //assert(e_iv.size() == edge_count*3);
#ifdef VERBOSE_DEBUG
    cerr << "storing reverse edges" << endl;
#endif

    size_t t_itr = 0;
    for (size_t k = 0; k < node_count; ++k) {
        //cerr << k << endl;
        int64_t t_id = i_iv[k];
        size_t t_rank = k+1;
        t_iv[t_itr] = t_rank;
        t_bv[t_itr] = 1;
        ++t_itr;
        for (auto end : { false, true }) {
            if (to_from.find(Side(t_id, end)) != to_from.end()) {
                auto f_side_itr = to_from.find(Side(t_id, end));
                if (f_side_itr != to_from.end()) {
                    for (auto& f_side : f_side_itr->second) {
                        size_t f_rank = id_to_rank(f_side.first);
                        // store link
                        t_iv[t_itr] = f_rank;
                        t_bv[t_itr] = 0;
                        // store side for end of edge
                        t_to_end_bv[t_itr] = end;
                        t_from_start_bv[t_itr] = f_side.second;
                        ++t_itr;
                    }
                }
            }
        }
    }

    // compress the reverse direction side information
    util::assign(t_to_end_cbv, sd_vector<>(t_to_end_bv));
    util::assign(t_from_start_cbv, sd_vector<>(t_from_start_bv));


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

#ifdef VERBOSE_DEBUG
    cerr << "storing paths" << endl;
#endif
    // paths
    //path_nodes[name].push_back(m.position().node_id());
    string path_names;
    size_t path_entities = 0; // count of nodes and edges
    for (auto& pathpair : path_nodes) {
        // add path name
        const string& path_name = pathpair.first;
        //cerr << path_name << endl;
        vector<Mapping> walk;
        for (auto& m : pathpair.second) {
            walk.push_back(m.second);
        }
        path_names += start_marker + path_name + end_marker;
        XGPath* path = new XGPath(path_name, walk, entity_count, *this, node_label);
        paths.push_back(path);
        path_entities += path->member_count;
    }

    // handle path names
    util::assign(pn_iv, int_vector<>(path_names.size()));
    util::assign(pn_bv, bit_vector(path_names.size()));
    // now record path name starts
    for (size_t i = 0; i < path_names.size(); ++i) {
        pn_iv[i] = path_names[i];
        if (path_names[i] == start_marker) {
            pn_bv[i] = 1; // register name start
        }
    }
    util::assign(pn_bv_rank, rank_support_v<1>(&pn_bv));
    util::assign(pn_bv_select, bit_vector::select_1_type(&pn_bv));
    
    //util::bit_compress(pn_iv);
    string path_name_file = "@pathnames.iv";
    store_to_file((const char*)path_names.c_str(), path_name_file);
    construct(pn_csa, path_name_file, 1);

    // entity -> paths
    util::assign(ep_iv, int_vector<>(path_entities+entity_count));
    util::assign(ep_bv, bit_vector(path_entities+entity_count));
    size_t ep_off = 0;
    for (size_t i = 0; i < entity_count; ++i) {
        ep_bv[ep_off] = 1;
        ep_iv[ep_off] = 0; // null so we can detect entities with no path membership
        ++ep_off;
        for (size_t j = 0; j < paths.size(); ++j) {
            if (paths[j]->members[i] == 1) {
                ep_iv[ep_off++] = j+1;
            }
        }
    }

    util::bit_compress(ep_iv);
    //cerr << ep_off << " " << path_entities << " " << entity_count << endl;
    assert(ep_off <= path_entities+entity_count);
    util::assign(ep_bv_rank, rank_support_v<1>(&ep_bv));
    util::assign(ep_bv_select, bit_vector::select_1_type(&ep_bv));

#ifdef VERBOSE_DEBUG
    cerr << "storing threads" << endl;
#endif
    
    // Just store all the paths that are all perfect mappings as threads.
    // We end up converting *back* into Path objects.
    for (auto& pathpair : path_nodes) {
        Path reconstructed;
        
        // Grab the name
        reconstructed.set_name(pathpair.first);
        
        bool all_perfect = true;
        
        // Grab the Mappings, which are now sorted by rank
        for (auto& m : pathpair.second) {
            *reconstructed.add_mapping() = m.second;
            
            // Make sure the mapping is perfect
            // TODO: handle offsets
            if(m.second.edit_size() > 1) {
                // Not a perfect mapping
                all_perfect = false;
                break;
            } else if(m.second.edit_size() == 0) {
                // Is a perfect mapping
                continue;
            } else {
                // We have exactly one edit. Is it perfect?
                auto edit = m.second.edit(0);
                
                if(edit.from_length() != edit.to_length() || edit.sequence() != "") {
                    // The edit calls for actual editing
                    all_perfect = false;
                    break;
                }
            }
        }
        
        if(all_perfect) {
            // We actually want to insert this path as a thread
            insert_thread(reconstructed);
        }
        
    }
    

#ifdef DEBUG_CONSTRUCTION
    cerr << "|s_iv| = " << size_in_mega_bytes(s_iv) << endl;
    cerr << "|f_iv| = " << size_in_mega_bytes(f_iv) << endl;
    cerr << "|t_iv| = " << size_in_mega_bytes(t_iv) << endl;

    cerr << "|f_from_start_cbv| = " << size_in_mega_bytes(f_from_start_cbv) << endl;
    cerr << "|t_to_end_cbv| = " << size_in_mega_bytes(t_to_end_cbv) << endl;

    cerr << "|f_bv| = " << size_in_mega_bytes(f_bv) << endl;
    cerr << "|t_bv| = " << size_in_mega_bytes(t_bv) << endl;

    cerr << "|i_iv| = " << size_in_mega_bytes(i_iv) << endl;
    //cerr << "|i_wt| = " << size_in_mega_bytes(i_wt) << endl;

    cerr << "|s_cbv| = " << size_in_mega_bytes(s_cbv) << endl;
    
    cerr << "|h_iv| = " << size_in_mega_bytes(h_iv) << endl;
    cerr << "|ts_iv| = " << size_in_mega_bytes(ts_iv) << endl;

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
    cerr << "|ep_iv| = " << size_in_mega_bytes(ep_iv) << endl;
    paths_mb_size += size_in_mega_bytes(ep_iv);
    cerr << "|ep_bv| = " << size_in_mega_bytes(ep_bv) << endl;
    paths_mb_size += size_in_mega_bytes(ep_bv);
    paths_mb_size += size_in_mega_bytes(ep_bv_rank);
    paths_mb_size += size_in_mega_bytes(ep_bv_select);
    cerr << "total paths size " << paths_mb_size << endl;
    // TODO you are missing the rest of the paths size in xg::paths
    // but this fragment should be factored into a function anyway
    
    cerr << "total size [MB] = " << (
        size_in_mega_bytes(s_iv)
        + size_in_mega_bytes(f_iv)
        + size_in_mega_bytes(t_iv)
        //+ size_in_mega_bytes(s_bv)
        + size_in_mega_bytes(f_bv)
        + size_in_mega_bytes(t_bv)
        + size_in_mega_bytes(i_iv)
        //+ size_in_mega_bytes(i_wt)
        + size_in_mega_bytes(s_cbv)
        + size_in_mega_bytes(h_iv)
        + size_in_mega_bytes(ts_iv)
        // TODO: add in size of the dynamic bs_iv
        + paths_mb_size
        ) << endl;

#endif

    if (print_graph) {
        cerr << "printing graph" << endl;
        cerr << s_iv << endl;
        for (int i = 0; i < s_iv.size(); ++i) {
            cerr << revdna3bit(s_iv[i]);
        } cerr << endl;
        cerr << s_bv << endl;
        cerr << i_iv << endl;
        cerr << f_iv << endl;
        cerr << f_bv << endl;
        cerr << t_iv << endl;
        cerr << t_bv << endl;
        cerr << "paths" << endl;
        for (auto& path : paths) {
            cerr << path->name << endl;
            cerr << path->members << endl;
            cerr << path->ids << endl;
            cerr << path->ranks << endl;
            cerr << path->directions << endl;
            cerr << path->positions << endl;
            cerr << path->offsets << endl;
        }
        cerr << ep_bv << endl;
        cerr << ep_iv << endl;
    }

    if (validate_graph) {
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

        // -1 here seems weird
        // what?
        cerr << "validating forward edge table" << endl;
        for (size_t j = 0; j < f_iv.size()-1; ++j) {
            if (f_bv[j] == 1) continue;
            // from id == rank
            size_t fid = i_iv[f_bv_rank(j)-1];
            // to id == f_cbv[j]
            size_t tid = i_iv[f_iv[j]-1];
            bool from_start = f_from_start_bv[j];
            // get the to_end
            bool to_end = false;
            for (auto& side : from_to[Side(fid, from_start)]) {
                if (side.first == tid) {
                    to_end = side.second;
                }
            }
            if (from_to[Side(fid, from_start)].count(Side(tid, to_end)) == 0) {
                cerr << "could not find edge (f) "
                     << fid << (from_start ? "+" : "-")
                     << " -> "
                     << tid << (to_end ? "+" : "-")
                     << endl;
                assert(false);
            }
        }

        cerr << "validating reverse edge table" << endl;
        for (size_t j = 0; j < t_iv.size()-1; ++j) {
            //cerr << j << endl;
            if (t_bv[j] == 1) continue;
            // from id == rank
            size_t tid = i_iv[t_bv_rank(j)-1];
            // to id == f_cbv[j]
            size_t fid = i_iv[t_iv[j]-1];
            //cerr << tid << " " << fid << endl;

            bool to_end = t_to_end_bv[j];
            // get the to_end
            bool from_start = false;
            for (auto& side : to_from[Side(tid, to_end)]) {
                if (side.first == fid) {
                    from_start = side.second;
                }
            }
            if (to_from[Side(tid, to_end)].count(Side(fid, from_start)) == 0) {
                cerr << "could not find edge (t) "
                     << fid << (from_start ? "+" : "-")
                     << " -> "
                     << tid << (to_end ? "+" : "-")
                     << endl;
                assert(false);
            }
        }
    
        cerr << "validating paths" << endl;
        for (auto& pathpair : path_nodes) {
            const string& name = pathpair.first;
            auto& path = pathpair.second;
            size_t prank = path_rank(name);
            //cerr << path_name(prank) << endl;
            assert(path_name(prank) == name);
            sd_vector<>& pe_bv = paths[prank-1]->members;
            int_vector<>& pp_iv = paths[prank-1]->positions;
            sd_vector<>& dir_bv = paths[prank-1]->directions;
            // check each entity in the nodes is present
            // and check node reported at the positions in it
            size_t pos = 0;
            size_t in_path = 0;
            for (auto& t : path) {
                auto& m = t.second;
                int64_t id = m.position().node_id();
                bool rev = m.position().is_reverse();
                // todo rank
                assert(pe_bv[node_rank_as_entity(id)-1]);
                assert(dir_bv[in_path] == rev);
                Node n = node(id);
                //cerr << id << " in " << name << endl;
                auto p = node_positions_in_path(id, name);
                assert(std::find(p.begin(), p.end(), pos) != p.end());
                for (size_t k = 0; k < n.sequence().size(); ++k) {
                    //cerr << "id " << id << " ==? " << node_at_path_position(name, pos+k) << endl;
                    assert(id == node_at_path_position(name, pos+k));
                    assert(id == mapping_at_path_position(name, pos+k).position().node_id());
                }
                pos += n.sequence().size();
                ++in_path;
            }
            //cerr << path_name << " rank = " << prank << endl;
            // check membership now for each entity in the path
        }
        
        cerr << "validating threads" << endl;
        
        // How many thread orientations are in the index?
        size_t threads_found = 0;
        // And how many shoukd we have inserted?
        size_t threads_expected = 0;
        
        for(auto path : extract_threads()) {
#ifdef VERBOSE_DEBUG
            cerr << "Path: ";
            for(size_t i = 0; i < path.mapping_size(); i++) {
                Mapping mapping = path.mapping(i);
                cerr << mapping.position().node_id() * 2 + mapping.position().is_reverse() << "; ";
            }
            cerr << endl;
#endif
            // Make sure we can search all the threads we find present in the index
            assert(count_matches(path) > 0);
            
            threads_found++;
        }
        
        for (auto& pathpair : path_nodes) {
            Path reconstructed;
            
            // Grab the name
            reconstructed.set_name(pathpair.first);
            
            bool all_perfect = true;
            
            // Grab the Mappings, which are now sorted by rank
            for (auto& m : pathpair.second) {
                *reconstructed.add_mapping() = m.second;
                
                // Make sure the mapping is perfect
                // TODO: handle offsets
                if(m.second.edit_size() > 1) {
                    // Not a perfect mapping
                    all_perfect = false;
                    break;
                } else if(m.second.edit_size() == 0) {
                    // Is a perfect mapping
                    continue;
                } else {
                    // We have exactly one edit. Is it perfect?
                    auto edit = m.second.edit(0);
                    
                    if(edit.from_length() != edit.to_length() || edit.sequence() != "") {
                        // The edit calls for actual editing
                        all_perfect = false;
                        break;
                    }
                }
            }
            
            if(all_perfect) {
                // This path should have been inserted. Look for it.
                assert(count_matches(reconstructed) > 0);
                
                threads_expected += 2;
            }
            
        }
        
        // Make sure we have the right number of threads.
        assert(threads_found == threads_expected);

        cerr << "graph ok" << endl;
    }
}

Node XG::node(int64_t id) const {
    Node n;
    n.set_id(id);
    //cerr << omp_get_thread_num() << " looks for " << id << endl;
    n.set_sequence(node_sequence(id));
    return n;
}

string XG::node_sequence(int64_t id) const {
    size_t rank = id_to_rank(id);
    size_t start = s_cbv_select(rank);
    size_t end = rank == node_count ? s_cbv.size() : s_cbv_select(rank+1);
    string s; s.resize(end-start);
    for (size_t i = start; i < s_cbv.size() && i < end; ++i) {
        s[i-start] = revdna3bit(s_iv[i]);
    }
    return s;
}

size_t XG::node_length(int64_t id) const {
    size_t rank = id_to_rank(id);
    size_t start = s_cbv_select(rank);
    size_t end = rank == node_count ? s_cbv.size() : s_cbv_select(rank+1);
    return end-start;
}

char XG::pos_char(int64_t id, bool is_rev, size_t off) const {
    assert(off < node_length(id));
    if (!is_rev) {
        size_t rank = id_to_rank(id);
        size_t pos = s_cbv_select(rank) + off;
        assert(pos < s_iv.size());
        char c = revdna3bit(s_iv[pos]);
        return c;
    } else {
        size_t rank = id_to_rank(id);
        size_t pos = s_cbv_select(rank+1) - (off+1);
        assert(pos < s_iv.size());
        char c = revdna3bit(s_iv[pos]);
        return reverse_complement(c);
    }
}

string XG::pos_substr(int64_t id, bool is_rev, size_t off, size_t len) const {
    if (!is_rev) {
        size_t rank = id_to_rank(id);
        size_t start = s_cbv_select(rank) + off;
        assert(start < s_iv.size());
        // get until the end position, or the end of the node, which ever is first
        size_t end;
        if (!len) {
            end = s_cbv_select(rank+1);
        } else {
            end = min(start + len, (size_t)s_cbv_select(rank+1));
        }
        assert(end < s_iv.size());
        string s; s.resize(end-start);
        for (size_t i = start; i < s_cbv.size() && i < end; ++i) {
            s[i-start] = revdna3bit(s_iv[i]);
        }
        return s;
    } else {
        size_t rank = id_to_rank(id);
        size_t end = s_cbv_select(rank+1) - off;
        assert(end < s_iv.size());
        // get until the end position, or the end of the node, which ever is first
        size_t start;
        if (len > end || !len) {
            start = s_cbv_select(rank);
        } else {
            start = max(end - len, (size_t)s_cbv_select(rank));
        }
        assert(end < s_iv.size());
        string s; s.resize(end-start);
        for (size_t i = start; i < s_cbv.size() && i < end; ++i) {
            s[i-start] = revdna3bit(s_iv[i]);
        }
        return reverse_complement(s);
    }
}

size_t XG::id_to_rank(int64_t id) const {
    return r_iv[id-min_id];
}

int64_t XG::rank_to_id(size_t rank) const {
    //assert(rank > 0);
    return i_iv[rank-1];
}

vector<Edge> XG::edges_of(int64_t id) const {
    auto e1 = edges_to(id);
    auto e2 = edges_from(id);
    e1.reserve(e1.size() + distance(e2.begin(), e2.end()));
    e1.insert(e1.end(), e2.begin(), e2.end());
    // now get rid of duplicates
    vector<Edge> e3;
    set<string> seen;
    for (auto& edge : e1) {
        string s; edge.SerializeToString(&s);
        if (!seen.count(s)) {
            e3.push_back(edge);
            seen.insert(s);
        }
    }
    return e3;
}

vector<Edge> XG::edges_to(int64_t id) const {
    vector<Edge> edges;
    size_t rank = id_to_rank(id);
    size_t t_start = t_bv_select(rank)+1;
    size_t t_end = rank == node_count ? t_bv.size() : t_bv_select(rank+1);
    for (size_t i = t_start; i < t_end; ++i) {
        Edge edge;
        edge.set_to(id);
        edge.set_from(rank_to_id(t_iv[i]));
        edge.set_from_start(t_from_start_cbv[i]);
        edge.set_to_end(t_to_end_cbv[i]);
        edges.push_back(edge);
    }
    return edges;
}

vector<Edge> XG::edges_from(int64_t id) const {
    vector<Edge> edges;
    size_t rank = id_to_rank(id);
    size_t f_start = f_bv_select(rank)+1;
    size_t f_end = rank == node_count ? f_bv.size() : f_bv_select(rank+1);
    for (size_t i = f_start; i < f_end; ++i) {
        Edge edge;
        edge.set_from(id);
        edge.set_to(rank_to_id(f_iv[i]));
        edge.set_from_start(f_from_start_cbv[i]);
        edge.set_to_end(f_to_end_cbv[i]);
        edges.push_back(edge);
    }
    return edges;
}

vector<Edge> XG::edges_on_start(int64_t id) const {
    vector<Edge> edges;
    for (auto& edge : edges_of(id)) {
        if((edge.to() == id && !edge.to_end()) || (edge.from() == id && edge.from_start())) {
            edges.push_back(edge);
        }
    }
    return edges;
}

vector<Edge> XG::edges_on_end(int64_t id) const {
    vector<Edge> edges;
    for (auto& edge : edges_of(id)) {
        if((edge.to() == id && edge.to_end()) || (edge.from() == id && !edge.from_start())) {
            edges.push_back(edge);
        }
    }
    return edges;
}

size_t XG::max_node_rank(void) const {
    return s_cbv_rank(s_cbv.size());
}

size_t XG::max_path_rank(void) const {
    //cerr << pn_bv << endl;
    //cerr << "..." << pn_bv_rank(pn_bv.size()) << endl;
    return pn_bv_rank(pn_bv.size());
}

size_t XG::node_rank_as_entity(int64_t id) const {
    //cerr << id_to_rank(id) << endl;
    return f_bv_select(id_to_rank(id))+1;
}

bool XG::entity_is_node(size_t rank) const {
    return 1 == f_bv[rank-1];
}

size_t XG::entity_rank_as_node_rank(size_t rank) const {
    return !entity_is_node(rank) ? 0 : f_iv[rank-1];
}

// snoop through the forward table to check if the edge exists
bool XG::has_edge(int64_t id1, bool from_start, int64_t id2, bool to_end) const {
    // invert the edge if we are doubly-reversed
    // this has the same meaning
    // ...confused
    /*
    if (from_start && to_end) {
        int64_t tmp = id1;
        id1 = id2; id2 = tmp;
        from_start = false;
        to_end = false;
    }
    */
    size_t rank1 = id_to_rank(id1);
    size_t rank2 = id_to_rank(id2);
    size_t f_start = f_bv_select(rank1);
    size_t f_end = rank1 == node_count ? f_bv.size() : f_bv_select(rank1+1);
    for (size_t i = f_start; i < f_end; ++i) {
        if (rank2 == f_iv[i]
            && f_from_start_cbv[i] == from_start
            && f_to_end_cbv[i] == to_end) {
            return true;
        }
    }
    return false;
}

size_t XG::edge_rank_as_entity(int64_t id1, bool from_start, int64_t id2, bool to_end) const {
    size_t rank1 = id_to_rank(id1);
    size_t rank2 = id_to_rank(id2);
#ifdef VERBOSE_DEBUG
    cerr << "Finding rank for "
         << id1 << (from_start?"+":"-") << " (" << rank1 << ") " << " -> "
         << id2 << (to_end?"-":"+") << " (" << rank2 << ")"<< endl;
#endif
    size_t f_start = f_bv_select(rank1) + 1;
    size_t f_end = rank1 == node_count ? f_bv.size() : f_bv_select(rank1+1);
#ifdef VERBOSE_DEBUG
    cerr << f_start << " to " << f_end << endl;
#endif
    for (size_t i = f_start; i < f_end; ++i) {
#ifdef VERBOSE_DEBUG
        cerr << f_iv[i] << endl;
#endif
        if (rank2 == f_iv[i]
            && f_from_start_cbv[i] == from_start
            && f_to_end_cbv[i] == to_end) {
#ifdef VERBOSE_DEBUG
            cerr << i << endl;
#endif
            return i+1;
        }
    }
    //cerr << "edge does not exist: " << id1 << " -> " << id2 << endl;
    assert(false);
}

size_t XG::edge_rank_as_entity(const Edge& edge) const {
    if(has_edge(edge.from(), edge.from_start(), edge.to(), edge.to_end())) {
        int64_t rank = edge_rank_as_entity(edge.from(), edge.from_start(), edge.to(), edge.to_end());
#ifdef VERBOSE_DEBUG
        cerr << "Found rank " << rank << endl;
#endif
        assert(!entity_is_node(rank));
        return rank;
    } else if(has_edge(edge.to(), !edge.to_end(), edge.from(), !edge.from_start())) {
        // Handle the case where the edge is spelled backwards; get the rank of the forwards version.
        int64_t rank = edge_rank_as_entity(edge.to(), !edge.to_end(), edge.from(), !edge.from_start());
#ifdef VERBOSE_DEBUG
        cerr << "Found rank " << rank << endl;
#endif
        assert(!entity_is_node(rank));
        return rank;
    } else {
        // Someone gave us an edge that doesn't exist.
        assert(false);
    }
}

Path XG::path(const string& name) const {
    // Extract a whole path by name
    
    // First find the XGPath we're using to store it.
    XGPath& xgpath = *(paths[path_rank(name)-1]);
    
    // Make a new path to fill in
    Path to_return;
    // Fill in the name
    to_return.set_name(name);
    
    for(size_t i = 0; i < xgpath.member_count; i++) {
        // For everything on the XGPath, put a Mapping on the real path.
        *(to_return.add_mapping()) = xgpath.mapping(i);
    }
    
    return to_return;
    
}

size_t XG::path_rank(const string& name) const {
    // find the name in the csa
    string query = start_marker + name + end_marker;
    auto occs = locate(pn_csa, query);
    if (occs.size() > 1) {
        cerr << "multiple hits for " << query << endl;
        assert(false);
    }
    if(occs.size() == 0) {
        // This path does not exist. Give back 0, which can never be a real path
        // rank.
        return 0;
    }
    //cerr << "path named " << name << " is at " << occs[0] << endl;
    return pn_bv_rank(occs[0])+1; // step past '#'
}

string XG::path_name(size_t rank) const {
    //cerr << "path rank " << rank << endl;
    size_t start = pn_bv_select(rank)+1; // step past '#'
    size_t end = rank == path_count ? pn_iv.size() : pn_bv_select(rank+1);
    end -= 1;  // step before '$'
    string name; name.resize(end-start);
    for (size_t i = start; i < end; ++i) {
        name[i-start] = pn_iv[i];
    }
    return name;
}

bool XG::path_contains_entity(const string& name, size_t rank) const {
    return 1 == paths[path_rank(name)-1]->members[rank-1];
}

bool XG::path_contains_node(const string& name, int64_t id) const {
    return path_contains_entity(name, node_rank_as_entity(id));
}

bool XG::path_contains_edge(const string& name, int64_t id1, bool from_start, int64_t id2, bool to_end) const {
    return path_contains_entity(name, edge_rank_as_entity(id1, from_start, id2, to_end));
}

vector<size_t> XG::paths_of_entity(size_t rank) const {
    size_t off = ep_bv_select(rank);
    assert(ep_bv[off++]);
    vector<size_t> path_ranks;
    while (off < ep_bv.size() && ep_bv[off] == 0) {
        path_ranks.push_back(ep_iv[off++]);
    }
    return path_ranks;
}

vector<size_t> XG::paths_of_node(int64_t id) const {
    return paths_of_entity(node_rank_as_entity(id));
}

vector<size_t> XG::paths_of_edge(int64_t id1, bool from_start, int64_t id2, bool to_end) const {
    return paths_of_entity(edge_rank_as_entity(id1, from_start, id2, to_end));
}

map<string, vector<Mapping>> XG::node_mappings(int64_t id) const {
    map<string, vector<Mapping>> mappings;
    // for each time the node crosses the path
    for (auto i : paths_of_entity(node_rank_as_entity(id))) {
        // get the path name
        string name = path_name(i);
        // get reference to the offset of the mapping in the path
        // to get the direction and (stored) rank
        for (auto j : node_ranks_in_path(id, name)) {
            // nb: path rank is 1-based, path index is 0-based
            mappings[name].push_back(paths[i-1]->mapping(j));
        }
    }
    return mappings;
}

void XG::neighborhood(int64_t id, size_t steps, Graph& g) const {
    *g.add_node() = node(id);
    expand_context(g, steps);
}

void XG::expand_context(Graph& g, size_t steps, bool add_paths) const {
    map<int64_t, Node*> nodes;
    map<pair<Side, Side>, Edge*> edges;
    set<int64_t> to_visit;
    // start with the nodes in the graph
    for (size_t i = 0; i < g.node_size(); ++i) {
        to_visit.insert(g.node(i).id());
        // handles the single-node case: we should still get the paths
        Node* np = g.mutable_node(i);
        nodes[np->id()] = np;
    }
    for (size_t i = 0; i < g.edge_size(); ++i) {
        auto& edge = g.edge(i);
        to_visit.insert(edge.from());
        to_visit.insert(edge.to());
        edges[make_pair(Side(edge.from(), edge.from_start()),
                        Side(edge.to(), edge.to_end()))] = g.mutable_edge(i);
    }
    // and expand
    for (size_t i = 0; i < steps; ++i) {
        set<int64_t> to_visit_next;
        for (auto id : to_visit) {
            // build out the graph
            // if we have nodes we haven't seeen
            if (nodes.find(id) == nodes.end()) {
                Node* np = g.add_node();
                nodes[id] = np;
                *np = node(id);
            }
            for (auto& edge : edges_of(id)) {
                auto sides = make_pair(Side(edge.from(), edge.from_start()),
                                       Side(edge.to(), edge.to_end()));
                if (edges.find(sides) == edges.end()) {
                    Edge* ep = g.add_edge(); *ep = edge;
                    edges[sides] = ep;
                }
                if (edge.from() == id) {
                    to_visit_next.insert(edge.to());
                } else {
                    to_visit_next.insert(edge.from());
                }
            }
        }
        to_visit = to_visit_next;
    }
    // then add connected nodes
    for (auto& e : edges) {
        auto& edge = e.second;
        // get missing nodes
        int64_t f = edge->from();
        if (nodes.find(f) == nodes.end()) {
            Node* np = g.add_node();
            nodes[f] = np;
            *np = node(f);
        }
        int64_t t = edge->to();
        if (nodes.find(t) == nodes.end()) {
            Node* np = g.add_node();
            nodes[t] = np;
            *np = node(t);
        }
    }
    if (add_paths) {
        add_paths_to_graph(nodes, g);
    }
}

    
// if the graph ids partially ordered, this works no prob
// otherwise... owch
// the paths become disordered due to traversal of the node ids in order
void XG::add_paths_to_graph(map<int64_t, Node*>& nodes, Graph& g) const {
    // map from path name to (map from mapping rank to mapping)
    map<string, map<size_t, Mapping>> paths;
    // mappings without 
    map<string, vector<Mapping>> unplaced;
    // use:
    //size_t node_position_in_path(int64_t id, const string& name) const;

    // pick up current paths in the graph
    for (size_t i = 0; i < g.path_size(); ++i) {
        auto& path = g.path(i);
        for (size_t j = 0; j < path.mapping_size(); ++j) {
            auto& m = path.mapping(j);
            if (m.rank()) {
                paths[path.name()][m.rank()] = m;
            } else {
                unplaced[path.name()].push_back(m);
            }                
        }
    }
    // do the same for the mappings in the list of nodes
    for (auto& n : nodes) {
        auto& id = n.first;
        auto& node = n.second;
        for (auto& n : node_mappings(id)) {
            auto& name = n.first;
            for (auto& m : n.second) {
                if (m.rank()) {
                    paths[name][m.rank()] = m;
                } else {
                    unplaced[name].push_back(m);
                }
            }
        }
    }
    // rebuild graph's paths
    // NB: mapping ranks allow us to remove this bit
    // only adding what we haven't seen before
    g.clear_path();
    for (auto& p : paths) {
        auto& name = p.first;
        auto& mappings = p.second;
        Path* path = g.add_path();
        path->set_name(name);
        for (auto& n : mappings) {
            *path->add_mapping() = n.second;
        }
        if (unplaced.find(name) != unplaced.end()) {
            auto& unp = unplaced[name];
            for (auto& m : unp) {
                *path->add_mapping() = m;
            }
        }
    }
}

void XG::get_id_range(int64_t id1, int64_t id2, Graph& g) const {
    id1 = max(min_id, id1);
    id2 = min(max_id, id2);
    for (auto i = id1; i <= id2; ++i) {
        *g.add_node() = node(i);
    }
}

/*
void XG::get_connected_nodes(Graph& g) {
}
*/

size_t XG::path_length(const string& name) const {
    return paths[path_rank(name)-1]->offsets.size();
}

// TODO, include paths
void XG::get_path_range(string& name, int64_t start, int64_t stop, Graph& g) const {
    // what is the node at the start, and at the end
    auto& path = *paths[path_rank(name)-1];
    size_t plen = path.offsets.size();
    if (start > plen) return; // no overlap with path
    size_t pr1 = path.offsets_rank(start+1)-1;
    // careful not to exceed the path length
    if (stop >= plen) stop = plen-1;
    size_t pr2 = path.offsets_rank(stop+1)-1;
    set<int64_t> nodes;
    set<pair<Side, Side> > edges;
    auto& pi_wt = path.ids;
    for (size_t i = pr1; i <= pr2; ++i) {
        int64_t id = rank_to_id(pi_wt[i]);
        nodes.insert(id);
        for (auto& e : edges_from(id)) {
            edges.insert(make_pair(Side(e.from(), e.from_start()), Side(e.to(), e.to_end())));
        }
        for (auto& e : edges_to(id)) {
            edges.insert(make_pair(Side(e.from(), e.from_start()), Side(e.to(), e.to_end())));
        }
    }
    for (auto& n : nodes) {
        *g.add_node() = node(n);
    }
    map<string, Path*> local_paths;
    for (auto& n : nodes) {
        for (auto& m : node_mappings(n)) {
            if (local_paths.find(m.first) == local_paths.end()) {
                Path* p = g.add_path();
                local_paths[m.first] = p;
                p->set_name(m.first);
            }
            Path* new_path = local_paths[m.first];
            // TODO output mapping direction
            //if () { m.second.is_reverse(true); }
            for (auto& n : m.second) {
                *new_path->add_mapping() = n;
            }
        }
    }
    for (auto& e : edges) {
        Edge edge;
        edge.set_from(e.first.first);
        edge.set_from_start(e.first.second);
        edge.set_to(e.second.first);
        edge.set_to_end(e.second.second);
        *g.add_edge() = edge;
    }
}

size_t XG::node_occs_in_path(int64_t id, const string& name) const {
    size_t p = path_rank(name)-1;
    auto& pi_wt = paths[p]->ids;
    return pi_wt.rank(pi_wt.size(), id);
}

vector<size_t> XG::node_ranks_in_path(int64_t id, const string& name) const {
    vector<size_t> ranks;
    size_t p = path_rank(name)-1;
    for (size_t i = 1; i <= node_occs_in_path(id, name); ++i) {
        ranks.push_back(paths[p]->ids.select(i, id));
        auto m = paths[p]->mapping(ranks.back());
    }
    return ranks;
}

vector<size_t> XG::node_positions_in_path(int64_t id, const string& name) const {
    auto& path = *paths[path_rank(name)-1];
    auto rank = id_to_rank(id);
    vector<size_t> pos_in_path;
    for (auto i : node_ranks_in_path(id, name)) {
        pos_in_path.push_back(path.positions[i]);
    }
    return pos_in_path;
}

int64_t XG::node_at_path_position(const string& name, size_t pos) const {
    size_t p = path_rank(name)-1;
    return paths[p]->ids[paths[p]->offsets_rank(pos+1)-1];
}

Mapping XG::mapping_at_path_position(const string& name, size_t pos) const {
    size_t p = path_rank(name)-1;
    return paths[p]->mapping(paths[p]->offsets_rank(pos+1)-1);
}

Mapping new_mapping(const string& name, int64_t id, size_t rank, bool is_reverse) {
    Mapping m;
    m.mutable_position()->set_node_id(id);
    m.mutable_position()->set_is_reverse(is_reverse);
    m.set_rank(rank);
    return m;
}

void parse_region(const string& target, string& name, int64_t& start, int64_t& end) {
    start = -1;
    end = -1;
    size_t foundFirstColon = target.find(":");
    // we only have a single string, use the whole sequence as the target
    if (foundFirstColon == string::npos) {
        name = target;
    } else {
        name = target.substr(0, foundFirstColon);
	    size_t foundRangeDash = target.find("-", foundFirstColon);
        if (foundRangeDash == string::npos) {
            start = atoi(target.substr(foundFirstColon + 1).c_str());
            end = start;
        } else {
            start = atoi(target.substr(foundFirstColon + 1, foundRangeDash - foundRangeDash - 1).c_str());
            end = atoi(target.substr(foundRangeDash + 1).c_str());
        }
    }
}

void to_text(ostream& out, Graph& graph) {
    out << "H" << "\t" << "HVN:Z:1.0" << endl;
    for (size_t i = 0; i < graph.node_size(); ++i) {
        auto& node = graph.node(i);
        out << "S" << "\t" << node.id() << "\t" << node.sequence() << endl;
    }
    for (size_t i = 0; i < graph.path_size(); ++i) {
        auto& path = graph.path(i);
        for (size_t j = 0; j < path.mapping_size(); ++j) {
            auto& mapping = path.mapping(i);
            string orientation = mapping.position().is_reverse() ? "-" : "+";
            out << "P" << "\t" << mapping.position().node_id() << "\t" << path.name() << "\t"
                << mapping.rank() << "\t" << orientation << "\n";
        }
    }
    for (int i = 0; i < graph.edge_size(); ++i) {
        auto& edge = graph.edge(i);
        out << "L" << "\t" << edge.from() << "\t"
            << (edge.from_start() ? "-" : "+") << "\t"
            << edge.to() << "\t"
            << (edge.to_end() ? "-" : "+") << endl;
    }
}

int64_t XG::where_to(int64_t current_side, int64_t visit_offset, int64_t new_side) const {
    // Given that we were at visit_offset on the current side, where will we be
    // on the new side? 
    
    // What will the new visit offset be?
    int64_t new_visit_offset = 0;
    
    // Work out where we're going as a node and orientation
    int64_t new_node_id = rank_to_id(new_side / 2);
    bool new_node_is_reverse = new_side % 2;
    
    // Work out what edges are going into the place we're going into.
    vector<Edge> edges = new_node_is_reverse ? edges_on_end(new_node_id) : edges_on_start(new_node_id);
    
    // Work out what node and orientation we came from
    int64_t old_node_id = rank_to_id(current_side / 2);
    bool old_node_is_reverse = current_side % 2;
    
    // What edge are we following
    Edge edge_taken = make_edge(old_node_id, old_node_is_reverse, new_node_id, new_node_is_reverse);
    
    // Make sure we find it
    bool edge_found = false;
    
    for(auto& edge : edges) {
        // Look at every edge in order.
        
        if(edges_equivalent(edge, edge_taken)) {
            // If we found the edge we're taking, break.
            edge_found = true;
            break;
        }
        
        // Otherwise add in the threads on this edge to the offset

        // Get the orientation number for this edge (which is *not* the edge we
        // are following and so doesn't involve old_node_anything). We only
        // treat it as reverse if the reverse direction is the only direction we
        // can take to get here.
        int64_t edge_orientation_number = ((edge_rank_as_entity(edge) - 1) * 2) + arrive_by_reverse(edge, new_node_id, new_node_is_reverse);
        
        int64_t contribution = h_iv[edge_orientation_number];
#ifdef VERBOSE_DEBUG
        cerr << contribution << " (from prev edge " << edge.from() << (edge.from_start() ? "L" : "R") << "-" << edge.to() << (edge.to_end() ? "R" : "L") << " at " << edge_orientation_number << ") + ";
#endif
        new_visit_offset += contribution;
    }
    
    assert(edge_found);
    
    // What edge out of all the edges we can take are we taking?
    int64_t edge_taken_index = -1;
    
    // Look at the edges we could have taken next
    vector<Edge> edges_out = old_node_is_reverse ? edges_on_start(old_node_id) : edges_on_end(old_node_id);
    
    for(int64_t i = 0; i < edges_out.size(); i++) {
        if(edges_equivalent(edges_out[i], edge_taken)) {
            // i is the index of the edge we took, of the edges available to us.
            edge_taken_index = i;
            break;
        }
    }
    
    assert(edge_taken_index != -1);
    
    // We need to un-const this int vector because DYNAMIC has no consts anywhere.
    //auto nonconst_bs_iv = const_cast<dynamic_int_vector*>(&bs_iv);
    dynamic_int_vector& nonconst_bs_iv = *bs_iv;
    
    // Where does the B_s[] range for the side we're leaving start?
    int64_t bs_start = nonconst_bs_iv.select(current_side - 2, BS_SEPARATOR) + 1;
    
    // Get the rank in B_s[] for our current side of our visit offset among
    // B_s[] entries pointing to the new node and add that in. Make sure to +2
    // to account for the nulls and separators.
    int64_t contribution = nonconst_bs_iv.rank(bs_start + visit_offset, edge_taken_index + 2) - nonconst_bs_iv.rank(bs_start, edge_taken_index + 2);
#ifdef VERBOSE_DEBUG
    cerr << contribution << " (via this edge) + ";
#endif
    new_visit_offset += contribution;
    
    // Get the number of threads starting at the new side and add that in.
    contribution = ts_iv[new_side];
#ifdef VERBOSE_DEBUG
    cerr << contribution << " (starting here) = ";
#endif
    new_visit_offset += contribution;
#ifdef VERBOSE_DEBUG
    cerr << new_visit_offset << endl;
#endif
    
    // Now we know where it actually ends up: after all the threads that start,
    // all the threads that come in via earlier edges, and all the previous
    // threads going there that come via this edge.
    return new_visit_offset;
}

void XG::insert_thread(const Path& t) {
    // We're going to insert this thread
    
    auto insert_thread_forward = [&](const Path& thread) {
    
#ifdef VERBOSE_DEBUG
        cerr << "Inserting thread " << thread.name() << " with " << thread.mapping_size() << " mappings" << endl;
#endif
        // Where does the current visit fall on its node? On the first node we
        // arbitrarily decide to be first of all the threads starting there.
        // TODO: Make sure that we actually end up ordering starts based on path
        // name or something later.
        int64_t visit_offset = 0;
        for(size_t i = 0; i < thread.mapping_size(); i++) {
            // For each visit to a node...
        
            // What side are we visiting?
            int64_t node_id = thread.mapping(i).position().node_id();
            bool node_is_reverse = thread.mapping(i).position().is_reverse();
            int64_t node_side = id_to_rank(node_id) * 2 + node_is_reverse;
            
#ifdef VERBOSE_DEBUG
            cerr << "Visit side " << node_side << endl;
#endif

            // Where are we going next?
            
            if(i == thread.mapping_size() - 1) {
                // This is the last visit. Send us off to null
                
#ifdef VERBOSE_DEBUG
                cerr << "End the thread." << endl;
#endif
                // Stick a new entry in the B array at the place where it belongs.
                bs_iv->insert(bs_iv->select(node_side - 2, BS_SEPARATOR) + 1 + visit_offset, BS_NULL);
            } else {
                // This is not the last visit. Send us off to the next place, and update the count on the edge.
                
                // Work out where we're actually going next
                int64_t next_id = thread.mapping(i + 1).position().node_id();
                bool next_is_reverse = thread.mapping(i + 1).position().is_reverse();
                int64_t next_side = id_to_rank(next_id) * 2 + next_is_reverse;

                // What edge do we take to get there? We're going to search for
                // the actual edge articulated in the official orientation.
                Edge edge_wanted = make_edge(node_id, node_is_reverse, next_id, next_is_reverse);
                
                // And what is its index in the edges we can take?
                int64_t edge_taken_index = -1;
                
                // Look at the edges we could have taken
                vector<Edge> edges_out = node_is_reverse ? edges_on_start(node_id) : edges_on_end(node_id);
                for(int64_t i = 0; i < edges_out.size(); i++) {
                    if(edge_taken_index == -1 && edges_equivalent(edges_out[i], edge_wanted)) {
                        // i is the index of the edge we took, of the edges available to us.
                        edge_taken_index = i;
                    }
#ifdef VERBOSE_DEBUG
                    cerr << "Edge #" << i << ": "  << edges_out[i].from() << (edges_out[i].from_start() ? "L" : "R") << "-" << edges_out[i].to() << (edges_out[i].to_end() ? "R" : "L") << endl;
                    
                    // Work out what its number is for the orientation it goes
                    // out in. We know we have this edge, so we just see if we
                    // have to depart in reveres and if so take the edge in
                    // reverse.
                    int64_t edge_orientation_number = (edge_rank_as_entity(edges_out[i]) - 1) * 2 + depart_by_reverse(edges_out[i], node_id, node_is_reverse);
                    
                    cerr << "\tOrientation rank: " << edge_orientation_number << endl;
#endif
                }
                
                assert(edge_taken_index != -1);

#ifdef VERBOSE_DEBUG
                cerr << "Proceed to " << next_side << " via edge #" << edge_taken_index << "/" << edges_out.size() << endl;
#endif
            
                // Make a nice reference to the edge we're taking in its real orientation.
                auto& edge_taken = edges_out[edge_taken_index];
            
                // Where do we insert the B value?
                int64_t target_entry = bs_iv->select(node_side - 2, BS_SEPARATOR) + 1 + visit_offset;
                
                // Where is the first spot not in our range at the moment?
                int64_t next_range_start = (node_side - 2 == bs_iv->rank(bs_iv->size(), BS_SEPARATOR) - 1 ? bs_iv->size() : bs_iv->select(node_side - 2 + 1, BS_SEPARATOR));
    
#ifdef VERBOSE_DEBUG
                cerr << "Will go in entry " << target_entry << " of " << bs_iv->size() << endl;
#endif

                // Make sure our new entry won't end up owned by the next side over.
                assert(target_entry <= next_range_start);

                // Stick a new entry in the B array at the place where it belongs.
                // Make sure to +2 to leave room in the number space for the separators and null destinations.
                bs_iv->insert(bs_iv->select(node_side - 2, BS_SEPARATOR) + 1 + visit_offset, edge_taken_index + 2);
                
                // Update the usage count for the edge going form here to the next node
                // Make sure that edge storage direction is correct.
                // Which orientation of an edge are we crossing?
                int64_t edge_orientation_number = (edge_rank_as_entity(edge_taken) - 1) * 2 + depart_by_reverse(edge_taken, node_id, node_is_reverse);
                
#ifdef VERBOSE_DEBUG
                cerr << "We need to add 1 to the traversals of oriented edge " << edge_orientation_number << endl;
#endif
                
                // Increment the count for the edge
                h_iv[edge_orientation_number]++;
                
#ifdef VERBOSE_DEBUG
                cerr << "h = [";
                for(int64_t k = 0; k < h_iv.size(); k++) {
                    cerr << h_iv[k] << ", ";
                }
                cerr << "]" << endl;
#endif
                
                // Now where do we go to on the next visit?
                visit_offset = where_to(node_side, visit_offset, next_side);
                
#ifdef VERBOSE_DEBUG

                cerr << "Offset " << visit_offset << " at side " << next_side << endl;
#endif
                
            }
            
#ifdef VERBOSE_DEBUG
            
            for(size_t j = 0; j < bs_iv->size(); j++) {
                cerr << bs_iv->at(j) << "; ";
            }
            cerr << endl;

            cerr << "Node " << node_id << " orientation " << node_is_reverse <<
                " has rank " <<
                ((node_rank_as_entity(node_id) - 1) * 2 + node_is_reverse) <<
                " of " << h_iv.size() << endl;
#endif
            
            // Increment the usage count of the node in this orientation
            h_iv[(node_rank_as_entity(node_id) - 1) * 2 + node_is_reverse]++;
            
            if(i == 0) {
                // The thread starts here
                ts_iv[node_side]++;
            }
        }
        
    };
    
    // We need a simple reverse that works only for perfect match paths
    auto simple_reverse = [&](const Path& thread) {
        // Clone the thread (TODO: this is going to copy all the thread data anyway)
        Path reversed = thread;
        
        // TODO: give it a reversed name or something
        
        for(size_t i = 0; i < thread.mapping_size(); i++) { 
            // Copy the mappings from back to front, flipping the is_reverse on their positions.
            Mapping reversing = thread.mapping(thread.mapping_size() - 1 - i);
            reversing.mutable_position()->set_is_reverse(!reversing.position().is_reverse());
            *reversed.mutable_mapping(i) = reversing;
        }
        
        return reversed;        
    };
    
    // Insert forward
    insert_thread_forward(t);
    
    // Insert reverse
    insert_thread_forward(simple_reverse(t));
    
    // TODO: name annotation
    
}

list<Path> XG::extract_threads() const {

    // Fill in a lsut of paths found
    list<Path> found;

#ifdef VERBOSE_DEBUG
    cerr << "Extracting threads" << endl;
#endif

    for(int64_t i = 1; i < ts_iv.size(); i++) {
        // For each real side
    
#ifdef VERBOSE_DEBUG
        cerr << ts_iv[i] << " threads start at side " << i << endl;
#endif
    
        // For each side
        if(ts_iv[i] == 0) {
            // Skip it if no threads start at it
            continue;
        }
        
        for(int64_t j = 0; j < ts_iv[i]; j++) {
            // For every thread starting there
      
#ifdef debug    
            cerr << "Extracting thread " << j << endl;
#endif
            
            // make a new path
            Path path;
            
            // Start the side at i and the offset at j
            int64_t side = i;
            int64_t offset = j;
            
            // We need to un-const this int vector because DYNAMIC has no consts anywhere.
            dynamic_int_vector& nonconst_bs_iv = *bs_iv;
            
            while(true) {
                // Unpack the side into a node traversal
                Mapping m;
                m.mutable_position()->set_node_id(rank_to_id(side / 2));
                m.mutable_position()->set_is_reverse(side % 2);
                
                // Add the mapping to the path
                *path.add_mapping() = m;
                
#ifdef VERBOSE_DEBUG
                cerr << "At side " << side << endl;
                
                cerr << "Want B group " << side - 2 << " of " << nonconst_bs_iv.rank(nonconst_bs_iv.size(), BS_SEPARATOR) << " at " << nonconst_bs_iv.size() << endl;
#endif
                // Work out where we go
                
                // What edge of the available edges do we take?
                int64_t edge_index = nonconst_bs_iv.at(nonconst_bs_iv.select(side - 2, BS_SEPARATOR) + 1 + offset);
                
#ifdef VERBOSE_DEBUG
                cerr << "Group starts at " << nonconst_bs_iv.select(side - 2, BS_SEPARATOR) << endl;
                int64_t outbound_count = ((side - 2 == nonconst_bs_iv.rank(nonconst_bs_iv.size(), BS_SEPARATOR) - 1) ? nonconst_bs_iv.size() : nonconst_bs_iv.select(side - 2 + 1, BS_SEPARATOR)) - nonconst_bs_iv.select(side - 2, BS_SEPARATOR);
                cerr << "Local B entry " << offset << " of " << outbound_count << endl;
                assert(offset < outbound_count);
#endif
                
                
                // If we find a separator, we're very broken.
                assert(edge_index != BS_SEPARATOR);
                
                if(edge_index == BS_NULL) {
                    // Path ends here.
                    break;
                } else {
                    // Convert to an actual edge index
                    edge_index -= 2;
                }
                
#ifdef VERBOSE_DEBUG
                cerr << "Taking edge #" << edge_index << " from " << side << endl;
#endif

                // We also should not have negative edges.
                assert(edge_index >= 0);
                
                // Look at the edges we could have taken next
                vector<Edge> edges_out = side % 2 ? edges_on_start(rank_to_id(side / 2)) : edges_on_end(rank_to_id(side / 2));
                
                assert(edge_index < edges_out.size());
                
                Edge& taken = edges_out[edge_index];
                
#ifdef VERBOSE_DEBUG
                cerr << edges_out.size() << " edges possible." << endl;
#endif
                // Follow the edge
                int64_t other_node = taken.from() == rank_to_id(side / 2) ? taken.to() : taken.from();
                bool other_orientation = (side % 2) != taken.from_start() != taken.to_end();
                
                // Get the side 
                int64_t other_side = id_to_rank(other_node) * 2 + other_orientation;
                
#ifdef VERBOSE_DEBUG
                cerr << "Go to side " << other_side << endl;
#endif
                
                // Go there with where_to
                offset = where_to(side, offset, other_side);
                side = other_side;
    
            }
            
            found.push_back(path);
            
        }
    }
    
    return found;
}

size_t XG::count_matches(const Path& t) const {
    // This is just a really simple wrapper that does a single extend
    ThreadSearchState state;
    extend_search(state, t);
    return state.count();
}

void XG::extend_search(ThreadSearchState& state, const Path& t) const {
    
#ifdef VERBOSE_DEBUG
    cerr << "Looking for path: ";
    for(int64_t i = 0; i < t.mapping_size(); i++) {
        // For each item in the path
        const Mapping& mapping = t.mapping(i);
        int64_t next_id = mapping.position().node_id();
        bool next_is_reverse = mapping.position().is_reverse();
        int64_t next_side = id_to_rank(next_id) * 2 + next_is_reverse;
        cerr << next_side << "; ";
    }
    cerr << endl;
#endif
    
    
    for(int64_t i = 0; i < t.mapping_size(); i++) {
        // For each item in the path
        const Mapping& mapping = t.mapping(i);
        
        if(state.is_empty()) {
            // Don't bother trying to extend empty things.
            
#ifdef VERBOSE_DEBUG
            cerr << "Nothing selected!" << endl;
#endif
            
            break;
        }
        
        // TODO: make this mapping to side thing a function
        int64_t next_id = mapping.position().node_id();
        bool next_is_reverse = mapping.position().is_reverse();
        int64_t next_side = id_to_rank(next_id) * 2 + next_is_reverse;
        
#ifdef VERBOSE_DEBUG
        cerr << "Extend mapping to " << state.current_side << " range " << state.range_start << " to " << state.range_end << " with " << next_side << endl;
#endif
        
        if(state.current_side == 0) {
            // If the state is a start state, just select the whole node using
            // the node usage count in this orientation. TODO: orientation not
            // really important unless we're going to search during a path
            // addition.
            state.range_start = 0;
            state.range_end = h_iv[(node_rank_as_entity(next_id) - 1) * 2 + next_is_reverse];
        } else {
            // Else, look at where the path goes to and apply the where_to function to shrink the range down.
            state.range_start = where_to(state.current_side, state.range_start, next_side);
            state.range_end = where_to(state.current_side, state.range_end, next_side);
        }
        
        // Update the side that the state is on
        state.current_side = next_side;
    }
}

size_t serialize(XG::dynamic_int_vector* to_serialize, ostream& out, sdsl::structure_tree_node* child, const std::string name) {
    size_t written = 0;
    
    // Convert the dynamic int vector to an SDSL int vector
    int_vector<> converted(to_serialize->size());
    
    for(size_t i = 0; i < to_serialize->size(); i++) {
        converted[i] = to_serialize->at(i);
    }
    
    written += converted.serialize(out, child, name);
    
    return written;
}

XG::dynamic_int_vector* deserialize(istream& in) {

    int_vector<> to_convert;

    to_convert.load(in);
    
    XG::dynamic_int_vector* converted = new XG::dynamic_int_vector;
    
    for(size_t i = 0; i < to_convert.size(); i++) {
        converted->push_back(to_convert[i]);
    }
    
    return converted;
}

bool edges_equivalent(const Edge& e1, const Edge& e2) {
    return ((e1.from() == e2.from() && e1.to() == e2.to() && e1.from_start() == e2.from_start() && e1.to_end() == e2.to_end()) ||
        (e1.from() == e2.to() && e1.to() == e2.from() && e1.from_start() == !e2.to_end() && e1.to_end() == !e2.from_start()));
}

bool relative_orientation(const Edge& e1, const Edge& e2) {
    assert(edges_equivalent(e1, e2));
    
    // Just use the reverse-equivalence check from edges_equivalent.
    return e1.from() == e2.to() && e1.to() == e2.from() && e1.from_start() == !e2.to_end() && e1.to_end() == !e2.from_start();
}

bool arrive_by_reverse(const Edge& e, int64_t node_id, bool node_is_reverse) {
    if(e.to() == node_id && (node_is_reverse == e.to_end())) {
        // We can follow the edge forwards and arrive at the correct side of the node
        return false;
    } else if(e.to() == e.from() && e.from_start() != e.to_end()) {
        // Reversing self loop
        return false;
    }
    // Otherwise, since we know the edge goes to the node, we have to take it backward.
    return true;
}

bool depart_by_reverse(const Edge& e, int64_t node_id, bool node_is_reverse) {
    if(e.from() == node_id && (node_is_reverse == e.from_start())) {
        // We can follow the edge forwards and arrive at the correct side of the node
        return false;
    } else if(e.to() == e.from() && e.from_start() != e.to_end()) {
        // Reversing self loop
        return false;
    }
    // Otherwise, since we know the edge goes to the node, we have to take it backward.
    return true;
}

Edge make_edge(int64_t from, bool from_start, int64_t to, bool to_end) {
    Edge e;
    e.set_from(from);
    e.set_to(to);
    e.set_from_start(from_start);
    e.set_to_end(to_end);
    
    return e;
}

char reverse_complement(const char& c) {
    switch (c) {
        case 'A': return 'T'; break;
        case 'T': return 'A'; break;
        case 'G': return 'C'; break;
        case 'C': return 'G'; break;
        case 'N': return 'N'; break;
        // Handle the GCSA2 start/stop characters.
        case '#': return '$'; break;
        case '$': return '#'; break;
        default: return 'N';
    }
}

string reverse_complement(const string& seq) {
    string rc;
    rc.assign(seq.rbegin(), seq.rend());
    for (auto& c : rc) {
        switch (c) {
        case 'A': c = 'T'; break;
        case 'T': c = 'A'; break;
        case 'G': c = 'C'; break;
        case 'C': c = 'G'; break;
        case 'N': c = 'N'; break;
        // Handle the GCSA2 start/stop characters.
        case '#': c = '$'; break;
        case '$': c = '#'; break;
        default: break;
        }
    }
    return rc;
}

void extract_pos(const string& pos_str, int64_t& id, bool& is_rev, size_t& off) {
    // format is id:off for forward, and id:-off for reverse
    // find the colon
    auto s = pos_str.find(":");
    assert(s != string::npos);
    id = stol(pos_str.substr(0, s));
    auto r = pos_str.find("-", s);
    if (r == string::npos) {
        is_rev = false;
        off = stoi(pos_str.substr(s+1, pos_str.size()));
    } else {
        is_rev = true;
        off = stoi(pos_str.substr(r+1, pos_str.size()));
    }
}

void extract_pos_substr(const string& pos_str, int64_t& id, bool& is_rev, size_t& off, size_t& len) {
    // format is id:off:len on forward and id:-off:len on reverse
    auto s = pos_str.find(":");
    assert(s != string::npos);
    id = stol(pos_str.substr(0, s));
    auto r = pos_str.find("-", s);
    if (r == string::npos) {
        is_rev = false;
        // find second colon
        auto t = pos_str.find(":", s+1);
        assert(t != string::npos);
        off = stoi(pos_str.substr(s+1, t-s));
        len = stoi(pos_str.substr(t+1, pos_str.size()));
    } else {
        is_rev = true;
        auto t = pos_str.find(":", r+1);
        assert(t != string::npos);
        off = stoi(pos_str.substr(r+1, t-r+1));
        len = stoi(pos_str.substr(t+1, pos_str.size()));
    }
}

}
