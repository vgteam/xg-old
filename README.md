# succinct graphs

## challenge

`vg`'s current graph memory model is weak and extremely bloated. It relies on fixed-width 64-bit integer ids and large hash tables mapping these to other entities. This makes it difficult to store in memory, and a general-purpose key-value store (rocksdb) is used to allow low-memory access to the entire graph. Although this design has some advantages, querying the graph requires costly IO operations, and thus use must be managed carefully when developing high-performance applications.

Fully-indexed graphs should be cheap to store and hold in memory, but it doesn't seem there is a standard approach that can be used just for high-performance access to the sequence and identifier space of the graph. Most work has gone into improving performance for querying the text of such a graph (GCSA) or generating one out of sequencing reads (assemblers such as SGA, fermi2, and minia).

The basic requirement is a system that a minimal amount of memory to store the sequence of the graph, its edges, and paths in the graph, but still allows constant-time access to the essential features of the graph. The system should support accessing:

* the node's label (a DNA sequence, for instance, or URL)
* the node's neighbors (inbound and outbound edges)
* the node's region in the graph (ranges of node id space that are within some distance of the node)
* node locations relative to stored paths in the graph
* node and edge path membership

## sketch

A mutable system might be constructable with [wavelet tries](http://arxiv.org/abs/1204.3581), but research in this area is very new, and I have not found readily-available code for working with these systems. It should be possible to construct mutable wavelet tries using sdsl-lite as a basis, but at present this may be too complex an objective.

An immutable system seems like a straightforward thing to do. The basic idea is to store the sequences of all elements in the graph in a [compressed integer vector](https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/enc_vector.hpp#L48-L58). A second [compressed bitvector](https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/rrr_vector.hpp) of the same length flags node starts, and implicitly provides identifiers to the system--- the rank of this vector for a given position in the sequence vector is the node id. A given node's position in the sequence vector can be found by using select(id) on this bit vector.

Next come edges. For each to and from of an edge we store an integer vector of node ids (edges are directed and we'll doubly index them). All the edges from a particular node must be contiguous in this array, and the order of this node in the edge vector must correspond to its order in the sequence vector. If a node has no to or from edges, we give it a blank entry in this vector. We now add a new bit vector of the same length as the number of edges. We register a 1 in this bit vector for each position that is the first outbound (or inbound) edge from a particular node. As for the sequence space, we can take positions in the edge table and determine which node they correspond to, and we can select all the edges of a particular node.

Finally, paths can be represented as a collection of compressed bit vectors over the node space (which allow constant-time determination of node path membership) or as lists of node ids (which may compress better). If nodes are single-based, then path rank operations can be used to provide relative positioning of nodes in the path.
