#!/usr/bin/env sh
echo_green(){
    echo -e "\e[1;32m$1\e[0m"
}

set -e

cur_dir=$(pwd)
for test in $(find . -name makefile); do
    echo --- Running test: $test test ---
    cd ${cur_dir}/$(dirname $test)
    make test
done

echo_green 'All tests pass!'
