# xg
## *succinct labeled graphs with collections and paths*

[![Build Status](https://travis-ci.org/ekg/xg.svg)](https://travis-ci.org/ekg/xg)

## what?

This library provides an interface to the construction and query of succinct/compressed labeled graphs, labeled collections of nodes and edges, and paths through the graph. The primary motivation is to index [variation graphs](https://github.com/ekg/vg#vg) of the type produced by [`vg`](https://github.com/ekg/vg), but in principle the library could be used for any labeled, directed graph.

Graphs indexed with `xg` can be loaded and queried using variety of high-performance interfaces backed by efficient, succinct data structures from [sdsl-lite](https://github.com/simongog/sdsl-lite).

## usage

First build:

```shell
make
make test
```

Now you can index a graph constructed with vg, then obtain a subgraph starting at a particular node and extending up to 5 steps away:

```shell
./xg -v test/data/z.vg -o z.xg
./xg -i z.xg -n 235 -c 5 | vg view -
```

`vg view` allows us to see the graph in GFA format.

See [`main.cpp`](https://github.com/ekg/xg/blob/master/main.cpp) for example usage.

## challenge

`vg`'s current graph memory model is weak and extremely bloated. It relies on fixed-width 64-bit integer ids and large hash tables mapping these to other entities. This makes it difficult to store in memory, and a general-purpose key-value store (rocksdb) is used to allow low-memory access to the entire graph. Although this design has some advantages, querying the graph requires costly IO operations, and thus use must be managed carefully when developing high-performance applications.

Fully-indexed graphs should be cheap to store and hold in memory, but it doesn't seem there is a standard approach that can be used just for high-performance access to the sequence and identifier space of the graph. Most work has gone into improving performance for querying the text of such a graph ([GCSA](https://github.com/jltsiren/gcsa2)) or generating one out of sequencing reads (assemblers such as [SGA](https://github.com/jts/sga) or [fermi2](https://github.com/lh3/fermi2)).

The basic requirement is a system that a minimal amount of memory to store the sequence of the graph, its edges, and paths in the graph, but still allows constant-time access to the essential features of the graph. The system should support accessing:

* the node's label (a DNA sequence, for instance, or URL)
* the node's neighbors (inbound and outbound edges)
* the node's region in the graph (ranges of node id space that are within some distance of the node)
* node locations relative to stored paths in the graph
* node and edge path membership

## sketch

In theory we could construct a mutable system based on <a href="http://arxiv.org/abs/1204.3581">wavelet tries</a>, but research in this area is very new, and I have not found readily-available code for working with these systems. It should be possible to construct mutable wavelet tries using sdsl-lite as a basis, but at present this may be too complex an objective. An immutable system seems like a straightforward thing to do.

First some definitions. We have a graph <span class="math"><em>G</em> = <em>N</em>, <em>E</em>, <em>P</em></span> with nodes <span class="math"><em>N</em> = <em>n</em><sub>1</sub>, …, <em>n</em><sub>∣<em>N</em>∣</sub></span>, directed edges <span class="math"><em>E</em> = <em>e</em><sub>1</sub>, …, <em>e</em><sub>∣<em>E</em>∣</sub></span>, and paths <span class="math"><em>P</em> = <em>p</em><sub>1</sub>, …, <em>p</em><sub>∣<em>P</em>∣</sub></span>. Nodes match labels <span class="math"><em>l</em><sub><em>n</em><sub><em>i</em></sub></sub></span> to ranks <span class="math"><em>i</em></span> in the collection of node labels: <span class="math"><em>n</em><sub><em>i</em></sub> = <em>l</em><sub><em>n</em><sub><em>i</em></sub></sub>, <em>i</em></span>. Edges go from one node to another <span class="math"><em>e</em><sub><em>j</em></sub> = <em>n</em><sub><em>x</em></sub>, <em>n</em><sub><em>y</em></sub></span>. Paths match labels <span class="math"><em>l</em><sub><em>p</em><sub><em>k</em></sub></sub></span> to sets of nodes and edges <span class="math"><em>p</em><sub><em>k</em></sub> = <em>l</em><sub><em>p</em><sub><em>k</em></sub></sub>, {<em>n</em><sub>1</sub>, <em>e</em><sub>3</sub>, <em>n</em><sub>4</sub>, <em>e</em><sub>5</sub>, …}</span>.

We first store the concatenated sequences of all elements, <span class="math"><em>S</em> = <em>l</em><sub><em>n</em><sub>1</sub></sub><em>l</em><sub><em>n</em><sub>2</sub></sub><em>l</em><sub><em>n</em><sub>3</sub></sub>…<em>l</em><sub><em>n</em><sub>∣<em>N</em>∣</sub></sub></span>, in the graph in a <a href="https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/enc_vector.hpp#L48-L58">compressed integer vector</a>, <span class="math"><em>S</em><sub><em>i</em><em>v</em></sub></span>. A second <a href="https://github.com/simongog/sdsl-lite/blob/master/include/sdsl/rrr_vector.hpp">compressed bitvector</a>, <span class="math"><em>S</em><sub><em>b</em><em>v</em></sub> : ∣<em>S</em><sub><em>i</em><em>v</em></sub>∣ = ∣<em>S</em><sub><em>b</em><em>v</em></sub>∣</span>, flags node starts, providing a system of node identifiers. We can apply <span class="math"><em>r</em><em>a</em><em>n</em><em>k</em><sub>1</sub>(<em>S</em><sub><em>b</em><em>v</em></sub>, <em>x</em>)</span> to determine the node rank/id at a given position in <span class="math"><em>S</em><sub><em>i</em><em>v</em></sub></span>, and we can use <span class="math"><em>s</em><em>e</em><em>l</em><em>e</em><em>c</em><em>t</em><sub>1</sub>(<em>S</em><sub><em>b</em><em>v</em></sub>, <em>x</em>)</span> to find the positions in <span class="math"><em>S</em><sub><em>i</em><em>v</em></sub></span> corresponding to node with rank/id <span class="math"><em>x</em></span>, thus allowing basic navigation of the nodes and their labels.

To store edges we keep compressed integer vectors of node ids for the forward <span class="math"><em>F</em><sub><em>i</em><em>v</em></sub></span> and reverse <span class="math"><em>T</em><sub><em>i</em><em>v</em></sub></span> link directions, where <span class="math"><em>F</em><sub><em>i</em><em>v</em></sub> = <em>f</em><sub>1</sub>, …, <em>f</em><sub>∣<em>N</em>∣</sub></span> and <span class="math"><em>f</em><sub><em>i</em></sub> = <em>i</em>, <em>t</em><em>o</em><sub><em>i</em><sub>1</sub></sub>, …, <em>t</em><em>o</em><sub><em>i</em><sub>∣<em>t</em><em>o</em><sub><em>i</em></sub>∣</sub></sub></span>. <span class="math"><em>T</em><sub><em>i</em><em>v</em></sub></span> inverts this relationship, providing <span class="math"><em>T</em><sub><em>i</em><em>v</em></sub> = <em>t</em><sub>1</sub>, …, <em>t</em><sub>∣<em>N</em>∣</sub></span> and <span class="math"><em>t</em><sub><em>i</em></sub> = <em>i</em>, <em>f</em><em>r</em><em>o</em><em>m</em><sub><em>i</em><sub>1</sub></sub>, …, <em>f</em><em>r</em><em>o</em><em>m</em><sub><em>i</em><sub>∣<em>f</em><em>r</em><em>o</em><em>m</em><sub><em>i</em></sub>∣</sub></sub></span>. Recall that <span class="math"><em>i</em></span> is the rank of the node. Using another bitvector <span class="math"><em>F</em><sub><em>b</em><em>v</em></sub> : ∣<em>F</em><sub><em>b</em><em>v</em></sub>∣ = ∣<em>F</em><sub><em>i</em><em>v</em></sub>∣</span> and <span class="math"><em>T</em><sub><em>b</em><em>v</em></sub> : ∣<em>T</em><sub><em>b</em><em>v</em></sub>∣ = ∣<em>T</em><sub><em>i</em><em>v</em></sub>∣</span> for we record the first position of each node's entries in <span class="math"><em>F</em><sub><em>i</em><em>v</em></sub></span> and <span class="math"><em>T</em><sub><em>i</em><em>v</em></sub></span>. This first position simply records the rank <span class="math"><em>i</em></span> in <span class="math"><em>S</em><sub><em>i</em><em>v</em></sub></span>. The rest of the positions in the node's range record the ranks/ids of the nodes on the other end of the edge--- on the &quot;to&quot; end in the <span class="math"><em>F</em><sub><em>i</em><em>v</em></sub></span> and the &quot;from&quot; end in <span class="math"><em>T</em><sub><em>i</em><em>v</em></sub></span>. If a node has no edges either coming from or going to it, it will only be represented by reference to its own rank in the correspending edge integer vector.

We can represent the path space <span class="math"><em>P</em><sub><em>i</em></sub>, …, <em>P</em><sub><em>n</em></sub></span> of the graph using a bitvector marking which entities in the edge-from integer vector <span class="math"><em>F</em><sub><em>i</em><em>v</em></sub></span> lie in a path. For each traversed node or edge, we mark a 1 in a new bitvector <span class="math"><em>P</em><em>e</em><sub><em>i</em><sub><em>b</em><em>v</em></sub></sub> : ∣<em>P</em><em>e</em><sub><em>i</em><sub><em>b</em><em>v</em></sub></sub>∣ = ∣<em>F</em><sub><em>i</em><em>v</em></sub>∣</span>, which is typically sparse and compressed. We mark contained entries with 1 and set the un-traversed nodes and edges to 0. Each path thus maps a label to a list of nodes and edges. To support cycles in which a single node may be traversed multiple times, we also store the path <span class="math"><em>P</em><sub><em>i</em></sub></span> as a vector of node ids <span class="math"><em>P</em><em>i</em><em>d</em><sub><em>i</em></sub></span>, and use a wavelet tree to provide rank and select operations on this structure to find node ranks in the path. These node ranks map into another vector <span class="math"><em>P</em><em>o</em><sub><em>i</em><sub><em>i</em><em>v</em></sub></sub></span> which lists the offset of each node relative to the path. Finally, a bit vector of the length of each path <span class="math"><em>P</em><em>p</em><sub><em>i</em><sub><em>b</em><em>v</em></sub></sub></span> in which we store a <span class="math">1</span> at each position in the path where a node starts allows us to find the node at a particular position in the path. Rank can be used to determine the node id at a given position in the path. In conjunction these structures allow us to store the paths and employ them as relativistic coordinate systems. Paths can overlap and serve as a form of annotation of features in the path space of the graph, in DNA for instance we might record a gene or exon as a path.
