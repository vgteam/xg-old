# xg
## *succinct labeled graphs with collections and paths*

[![Build Status](https://travis-ci.org/vgteam/xg.svg)](https://travis-ci.org/vgteam/xg)

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
