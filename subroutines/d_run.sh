#!/usr/bin/bash

set -ex
trap "exit" INT TERM
trap "kill 0" EXIT

RUSTFLAGS="-Awarnings" cargo build --example dLigesis
BIN=../target/debug/examples/dLigesis
$BIN $i ./data/4 &

echo done