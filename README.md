# succinct graphs

## challenge

`vg`'s current graph memory model is weak and extremely bloated. It relies on fixed-width 64-bit integer ids and large hash tables mapping these to other entities. This makes it difficult to store in memory, and a general-purpose key-value store (rocksdb) is used to allow low-memory access to the entire graph. Although this design has some advantages, querying the graph requires costly IO operations, and thus use must be managed carefully when developing high-performance applications.

Fully-indexed graphs should be cheap to store and hold in memory, but it doesn't seem there is a standard approach that can be used just for high-performance access to the sequence and identifier space of the graph. Most work has gone into improving performance for querying the text of such a graph [GCSA](https://github.com/jltsiren/gcsa2) or generating one out of sequencing reads (assemblers such as [SGA](https://github.com/jts/sga) or [fermi2](https://github.com/lh3/fermi2)).

The basic requirement is a system that a minimal amount of memory to store the sequence of the graph, its edges, and paths in the graph, but still allows constant-time access to the essential features of the graph. The system should support accessing:

* the node's label (a DNA sequence, for instance, or URL)
* the node's neighbors (inbound and outbound edges)
* the node's region in the graph (ranges of node id space that are within some distance of the node)
* node locations relative to stored paths in the graph
* node and edge path membership

## sketch

A mutable system might be constructable with [wavelet tries](http://arxiv.org/abs/1204.3581), but research in this area is very new, and I have not found readily-available code for working with these systems. It should be possible to construct mutable wavelet tries using sdsl-lite as a basis, but at present this may be too complex an objective. An immutable system seems like a straightforward thing to do.

The basic idea is to store the sequences of all elements in the graph in a [compressed integer vector](https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/enc_vector.hpp#L48-L58). A second [compressed bitvector](https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/rrr_vector.hpp) of the same length flags node starts, and implicitly provides identifiers to the system--- the rank of this vector for a given position in the sequence vector is the node id. A given node's position in the sequence vector can be found by using select on this bit vector.

To store edges we keep compressed integer vectors of node ids. Edges are directed and we'll doubly index them, so we keep an edge-from and edge-to index. The edge-from and edge-to integer vectors have the same length as the number of nodes and edges in the graph. For each node we first write the node id, then for each edge going from or coming to the node, we write the id of the node that would be reached by traversing the edge. We mark the node "self-links" using a bitvector with a 1 for each node start and a 0 for each edge entry. The order of the node in these edge vectors must correspond to its order in the sequence vector. If a node has no to or from edges, it will only be represented by a "self" link in the edge vector.

We can represent the path space of the graph using a bitvector marking which entities in the edge-from integer vector lie in a path. For each traversed node or edge, we mark a 1 in a new bitvector, leaving the un-traversed nodes and edges 0. Each path thus maps a label to a list of nodes and edges.
