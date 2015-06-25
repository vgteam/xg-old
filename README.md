---
generator: pandoc
...

xg
==

*succinct labeled graphs with collections and paths*
----------------------------------------------------

[![Build
Status](https://travis-ci.org/ekg/xg.svg)](https://travis-ci.org/ekg/xg)

what?
-----

This library provides an interface to the construction and query of
succinct/compressed labeled graphs, labeled collections of nodes and
edges, and paths through the graph. The primary motivation is to index
[variation graphs](https://github.com/ekg/vg#vg) of the type produced by
[`vg`](https://github.com/ekg/vg), but in principle the library could be
used for any labeled, directed graph.

Graphs indexed with `xg` can be loaded and queried using variety of
high-performance interfaces backed by efficient, succinct data
structures from [sdsl-lite](https://github.com/simongog/sdsl-lite).

usage
-----

First build:

``` {.shell}
make
make test
```

Now you can index a graph constructed with vg, then obtain a subgraph
starting at a particular node and extending up to 5 steps away:

``` {.shell}
./xg -v test/data/z.vg -o z.xg
./xg -i z.xg -n 235 -c 5 | vg view -
```

`vg view` allows us to see the graph in GFA format.

See [`main.cpp`](https://github.com/ekg/xg/blob/master/main.cpp) for
example usage.

challenge
---------

`vg`'s current graph memory model is weak and extremely bloated. It
relies on fixed-width 64-bit integer ids and large hash tables mapping
these to other entities. This makes it difficult to store in memory, and
a general-purpose key-value store (rocksdb) is used to allow low-memory
access to the entire graph. Although this design has some advantages,
querying the graph requires costly IO operations, and thus use must be
managed carefully when developing high-performance applications.

Fully-indexed graphs should be cheap to store and hold in memory, but it
doesn't seem there is a standard approach that can be used just for
high-performance access to the sequence and identifier space of the
graph. Most work has gone into improving performance for querying the
text of such a graph ([GCSA](https://github.com/jltsiren/gcsa2)) or
generating one out of sequencing reads (assemblers such as
[SGA](https://github.com/jts/sga) or
[fermi2](https://github.com/lh3/fermi2)).

The basic requirement is a system that a minimal amount of memory to
store the sequence of the graph, its edges, and paths in the graph, but
still allows constant-time access to the essential features of the
graph. The system should support accessing:

-   the node's label (a DNA sequence, for instance, or URL)
-   the node's neighbors (inbound and outbound edges)
-   the node's region in the graph (ranges of node id space that are
    within some distance of the node)
-   node locations relative to stored paths in the graph
-   node and edge path membership

sketch
------

In theory we could construct a mutable system based on [wavelet
tries](http://arxiv.org/abs/1204.3581), but research in this area is
very new, and I have not found readily-available code for working with
these systems. It should be possible to construct mutable wavelet tries
using sdsl-lite as a basis, but at present this may be too complex an
objective. An immutable system seems like a straightforward thing to do.

First some definitions. We have a graph *G* = *N*, *E*, *P* with nodes
*N* = *n*~1~, …, *n*~∣*N*∣~, directed edges *E* = *e*~1~, …, *e*~∣*E*∣~,
and paths *P* = *p*~1~, …, *p*~∣*P*∣~. Nodes match labels *l*~*n*~*i*~~
to ranks *i* in the collection of node labels:
*n*~*i*~ = *l*~*n*~*i*~~, *i*. Edges go from one node to another
*e*~*j*~ = *n*~*x*~, *n*~*y*~. Paths match labels *l*~*p*~*k*~~ to sets
of nodes and edges
*p*~*k*~ = *l*~*p*~*k*~~, {*n*~1~, *e*~3~, *n*~4~, *e*~5~, …}.

We first store the concatenated sequences of all elements,
*S* = *l*~*n*~1~~*l*~*n*~2~~*l*~*n*~3~~…*l*~*n*~∣*N*∣~~, in the graph in
a [compressed integer
vector](https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/enc_vector.hpp#L48-L58),
*S*~*iv*~. A second [compressed
bitvector](https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/rrr_vector.hpp),
*S*~*bv*~ : ∣*S*~*iv*~∣ = ∣*S*~*bv*~∣, flags node starts, providing a
system of node identifiers. We can apply *rank*~1~(*S*~*bv*~, *x*) to
determine the node rank/id at a given position in *S*~*iv*~, and we can
use *select*~1~(*S*~*bv*~, *x*) to find the positions in *S*~*iv*~
corresponding to node with rank/id *x*, thus allowing basic navigation
of the nodes and their labels.

To store edges we keep compressed integer vectors of node ids for the
forward *F*~*iv*~ and reverse *T*~*iv*~ link directions, where
*F*~*iv*~ = *f*~1~, …, *f*~∣*N*∣~ and
*f*~*i*~ = *i*, *to*~*i*~1~~, …, *to*~*i*~∣*to*~*i*~∣~~. *T*~*iv*~
inverts this relationship, providing *T*~*iv*~ = *t*~1~, …, *t*~∣*N*∣~
and *t*~*i*~ = *i*, *from*~*i*~1~~, …, *from*~*i*~∣*from*~*i*~∣~~.
Recall that *i* is the rank of the node. Using another bitvector
*F*~*bv*~ : ∣*F*~*bv*~∣ = ∣*F*~*iv*~∣ and
*T*~*bv*~ : ∣*T*~*bv*~∣ = ∣*T*~*iv*~∣ for we record the first position
of each node's entries in *F*~*iv*~ and *T*~*iv*~. This first position
simply records the rank *i* in *S*~*iv*~. The rest of the positions in
the node's range record the ranks/ids of the nodes on the other end of
the edge--- on the "to" end in the *F*~*iv*~ and the "from" end in
*T*~*iv*~. If a node has no edges either coming from or going to it, it
will only be represented by reference to its own rank in the
correspending edge integer vector.

We can represent the path space of the graph using a bitvector marking
which entities in the edge-from integer vector *F*~*iv*~ lie in a path.
For each traversed node or edge, we mark a 1 in a new bitvector
*P*~*i*~**~*bv*~ : ∣*P*~*i*~*bv*~~∣ = ∣*F*~*iv*~∣. We mark contained
entries with 1 and set the un-traversed nodes and edges to 0. Each path
thus maps a label to a list of nodes and edges.
