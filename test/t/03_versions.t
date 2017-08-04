#!/usr/bin/env bash

BASH_TAP_ROOT=../bash-tap
. ../bash-tap/bash-tap-bootstrap

PATH=../bin:$PATH # for xg

plan tests 4

# Make sure we can read various old versions of XG format.
xg -i data/versions/v00.xg -o /dev/null
is $? 0 "XG version 0 can be deserialized"

xg -i data/versions/v01.xg -o /dev/null
is $? 0 "XG version 1 can be deserialized"

is "$(xg -i data/versions/vLarge.xg -o /dev/null 2>&1 | grep 'too new' | wc -l)" "1" "Future XG versions are rejected"

xg -v data/l.vg -o serialized.xg
is "$(cat serialized.xg | head -c6 | tail -c4 | xxd | cut -d' ' -f2,3 | tr -d ' ')" "00000001" "New XG files are written in version 1 format"
rm -f serialized.xg


