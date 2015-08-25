#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for xg

plan tests 7

is $(xg -v data/l.vg 2>&1 | grep ok | wc -l) 1 "a small graph verifies"
is $(xg -v data/lg.vg 2>&1 | grep ok | wc -l) 1 "a small graph with two named paths verifies"
is $(xg -v data/l+.vg 2>&1 | grep ok | wc -l) 1 "node ids need not start at 1"
is $(xg -v data/z.vg 2>&1 | grep ok | wc -l) 1 "a 1mb graph verifies"
xg -v data/z.vg -o data/z.vg.idx 2>/dev/null
is $? 0 "serialization works"
rm -f data/z.vg.idx
xg -v data/with_m.vg 2>/dev/null
is $? 0 "graphs can be compressed even with M"

xg -v data/cyclic_all.vg -o c.idx 2>/dev/null
is $(xg -c 10 -n 1 -i c.idx -T | grep '-' | wc -l) 4 "graphs with cycles and edges from specific sides can be stored and queried"
rm c.idx
