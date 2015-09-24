#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for xg

plan tests 15

xg -v data/z.vg -o z.idx 2>/dev/null
is $(xg -i z.idx -s 10331 | cut -f 2 -d\ ) "CAGCAGTGGAGCAGAAACAGAGGAGATGACACCATGGGGTAAGCACAGTC" "graph can be queried to obtain node labels"
is $(xg -i z.idx -f 10331 | md5sum | awk '{print $1}') "b7a5dbb50a04c66c3f9e25afcfa987b6" "graph can be queried to get from nodes"
is $(xg -i z.idx -t 10331 | md5sum | awk '{print $1}') "f7d6410e597fd59eb9ccbc1d7bfe24d1" "graph can be queried to get to nodes"
is $(xg -i z.idx -n 10331 -c 10 | md5sum | awk '{print $1}') "9aceef363cd0419d1bdb154fc9418ef2" "graph can be queried to get node context"
is $(xg -i z.idx -p z:500000-500500 | md5sum | awk '{print $1}') "711247fc64aab5bc8bd745e5579e0d9e" "graph can be queried to get a region of a particular path"
rm -f z.idx

xg -v data/l.vg -o l.idx 2>/dev/null
xg -i l.idx -p z:0-100 >/dev/null
is $? 0 "path queries can exceed reference length without error"

is $(xg -i l.idx -p z:0-10 | md5sum | cut -f 1 -d\ ) "d5b83270620631112a138419e3395e76" "paths can be queried from the small graph"
is $(xg -i l.idx -p z:0-100 -c 2 | md5sum | cut -f 1 -d\ ) "ad8f8e0fb1c9fc628f6d312d5b56a9c9" "the entire graph can be extracted with a long query and context"
rm -f l.idx

xg -v data/xyz.vg -o xyz.idx 2>/dev/null
(xg -i xyz.idx -p x:10-20 && xg -i xyz.idx -p y:10-20 && xg -i xyz.idx -p z:10-20) >/dev/null
is $? 0 "a multi-path graph can be queried"
rm -f xyz.idx

xg -v data/mult.xyz.vg -o mult.xyz.idx 2>/dev/null
is $(xg -i mult.xyz.idx -p y:10-20 | md5sum | awk '{print $1}') $(xg -i mult.xyz.idx -p x:10-20 | md5sum | awk '{print $1}') "a multi-path graph with overlapping paths can be queried"
rm -f mult.xyz.idx

cat data/ll.vg data/z.vg | xg -v - 2>/dev/null
is $? 0 "graphs can be constructed with multiple paths"

xg -v data/ll.vg -o ll.idx 2>/dev/null
is $(xg -i ll.idx -n 1 -c 10 | md5sum | cut -f 1 -d\ ) $(md5sum data/ll.vg | cut -f 1 -d\ ) "a small graph can be exactly reconstructed from the index"
rm ll.idx

xg -v data/cyclic_all.vg -o c.idx 2>/dev/null
#is $(xg -c 10 -n 1 -i c.idx -T | grep '-' | wc -l) 4 "graphs with cycles and edges from specific sides can be stored and queried"
is $(xg -O 1 -i c.idx | grep '1- -> 1+' | wc -l) 1 "can obtain edges to/from self, which are unique"
is $(xg -S 1 -i c.idx | grep '4+ -> 1+' | wc -l) 1 "can obtain edges on start"
is $(xg -E 1 -i c.idx | grep '1+ -> 2+' | wc -l) 1 "can obtain edges on end"

rm c.idx
