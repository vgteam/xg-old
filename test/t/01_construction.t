#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=..:$PATH # for succinctg

plan tests 3

is $(succinctg data/l.vg 2>&1 | grep ok | wc -l) 1 "a small graph verifies"
is $(succinctg data/l+.vg 2>&1 | grep ok | wc -l) 1 "node ids need not start at 1"
is $(succinctg data/z.vg 2>&1 | grep ok | wc -l) 1 "a 1mb graph verifies"
