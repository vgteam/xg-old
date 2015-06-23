#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for succinctg

plan tests 9

is $(succinctg -v data/l.vg 2>&1 | grep ok | wc -l) 1 "a small graph verifies"
is $(succinctg -v data/l+.vg 2>&1 | grep ok | wc -l) 1 "node ids need not start at 1"
is $(succinctg -v data/z.vg 2>&1 | grep ok | wc -l) 1 "a 1mb graph verifies"
succinctg -v data/z.vg -o data/z.vg.idx 2>/dev/null
is $? 0 "serialization works"
rm -f data/z.vg.idx
succinctg -v data/with_m.vg 2>/dev/null
is $? 0 "graphs can be compressed even with M"

succinctg -v data/z.vg -o z.idx 2>/dev/null
is $(succinctg -i z.idx -s 10331 | cut -f 2 -d\ ) "CAGCAGTGGAGCAGAAACAGAGGAGATGACACCATGGGGTAAGCACAGTC" "graph can be queried to obtain node labels"
is $(succinctg -i z.idx -f 10331 | md5sum | awk '{print $1}') "91e7cfdbeee58e6fb2dbd51d169f0277" "graph can be queried to get from nodes"
is $(succinctg -i z.idx -t 10331 | md5sum | awk '{print $1}') "6fafb73d8559760040253bc036ad0b41" "graph can be queried to get to nodes"

rm -f z.idx

cat data/ll.vg data/z.vg | succinctg -v - 2>/dev/null
is $? 0 "graphs can be constructed with multiple paths"
