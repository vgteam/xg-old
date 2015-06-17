#include "succinct.hpp"
#include "stream.hpp"

namespace scg {

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

    util::assign(s_iv, int_vector<>(seq_length));
    util::assign(s_bv, bit_vector(seq_length));
    util::assign(f_iv, int_vector<>(entity_count));
    util::assign(f_bv, bit_vector(entity_count));
    util::assign(t_iv, int_vector<>(entity_count));
    util::assign(t_bv, bit_vector(entity_count));
    
}

}
