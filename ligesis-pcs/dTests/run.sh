#!/usr/bin/bash

set -ex
trap "exit" INT TERM
trap "kill 0" EXIT

cd "$(dirname "$0")/.."
RUSTFLAGS="-Awarnings" cargo build --example $1 --release
BIN=../target/release/examples/$1

PROCS=()
for i in 0 1 2 3
do
  $BIN $i ./dTests/data/4 &
  pid=$!
  PROCS+=("$pid")
done

for pid in "${PROCS[@]}"
do
  wait $pid
done

echo done