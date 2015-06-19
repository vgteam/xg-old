#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for succinctg

plan tests 4

is $(succinctg -v data/l.vg 2>&1 | grep ok | wc -l) 1 "a small graph verifies"
is $(succinctg -v data/l+.vg 2>&1 | grep ok | wc -l) 1 "node ids need not start at 1"
is $(succinctg -v data/z.vg 2>&1 | grep ok | wc -l) 1 "a 1mb graph verifies"
succinctg -v data/z.vg -o data/z.vg.idx 2>/dev/null
is $? 0 "serialization works"
rm -f data/z.vg.idx
