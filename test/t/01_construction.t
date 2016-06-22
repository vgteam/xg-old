#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for xg

plan tests 12

is $(xg -Vv data/l.vg 2>&1 | grep ok | wc -l) 1 "a small graph verifies"
is $(xg -Vv data/lg.vg 2>&1 | grep ok | wc -l) 1 "a small graph with two named paths verifies"
is $(xg -Vrv data/lg.vg 2>&1 | grep ok | wc -l) 1 "a small graph with threads verifies"
is $(xg -Vrdv data/lg.vg 2>&1 | grep ok | wc -l) 1 "a small graph with batch-inserted threads verifies"
is $(xg -Vrv data/l+.vg 2>&1 | grep ok | wc -l) 1 "node ids need not start at 1"
is $(xg -Vrdv data/z.vg 2>&1 | grep ok | wc -l) 1 "a 1mb graph verifies"
xg -Vv data/z.vg -o data/z.vg.idx 2>/dev/null
is $? 0 "serialization works"
rm -f data/z.vg.idx
xg -Vv data/with_m.vg 2>/dev/null
is $? 0 "graphs can be compressed even with M"

xg -Vv data/cyclic_all.vg -o c.idx 2>/dev/null
is $(xg -c 10 -n 1 -i c.idx -T | grep '-' | wc -l) 2 "graphs with cycles and edges from specific sides can be stored and queried"
rm c.idx

is $(xg -Vrv data/cyclic_path.vg 2>&1 | grep ok | wc -l) 1 "a graph with a path cycle validates"

is $(xg -Vrv data/self_loop_paths.vg 2>&1 | grep ok | wc -l) 1 "a small graph with all self loops validates"
is $(xg -Vrv data/b.vg 2>&1 | grep ok | wc -l) 1 "a large graph with doubly-reversing edges validates"
