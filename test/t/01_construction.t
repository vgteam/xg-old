#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for xg

plan tests 13

is $(xg -v data/l.vg 2>&1 | grep ok | wc -l) 1 "a small graph verifies"
is $(xg -v data/lg.vg 2>&1 | grep ok | wc -l) 1 "a small graph with two named paths verifies"
is $(xg -v data/l+.vg 2>&1 | grep ok | wc -l) 1 "node ids need not start at 1"
is $(xg -v data/z.vg 2>&1 | grep ok | wc -l) 1 "a 1mb graph verifies"
xg -v data/z.vg -o data/z.vg.idx 2>/dev/null
is $? 0 "serialization works"
rm -f data/z.vg.idx
xg -v data/with_m.vg 2>/dev/null
is $? 0 "graphs can be compressed even with M"

xg -v data/z.vg -o z.idx 2>/dev/null
is $(xg -i z.idx -s 10331 | cut -f 2 -d\ ) "CAGCAGTGGAGCAGAAACAGAGGAGATGACACCATGGGGTAAGCACAGTC" "graph can be queried to obtain node labels"
is $(xg -i z.idx -f 10331 | md5sum | awk '{print $1}') "26c04eefc1a3bdcb97eb3be3d65743be" "graph can be queried to get from nodes"
is $(xg -i z.idx -t 10331 | md5sum | awk '{print $1}') "1073816e90cdd7c1b5dacdcf0a980bb3" "graph can be queried to get to nodes"
is $(xg -i z.idx -n 10331 -c 10 | md5sum | awk '{print $1}') "ccf3a8b68b188ffe86d867eb58dadbed" "graph can be queried to get node context"
is $(xg -i z.idx -p z:500000-500500 | md5sum | awk '{print $1}') "9d1f5448af3a5342c37b8e5c9523a0de" "graph can be queried to get a region of a particular path"
rm -f z.idx

cat data/ll.vg data/z.vg | xg -v - 2>/dev/null
is $? 0 "graphs can be constructed with multiple paths"

xg -v data/ll.vg -o ll.idx 2>/dev/null
is $(xg -i ll.idx -n 1 -c 10 | md5sum | awk '{print $1}') "13294ce74910217ba041b65c7f775552" "a small graph can be exactly reconstructed from the index"
rm ll.idx
