#!/bin/bash

rm -rf ./build
rm -rf ./bin/*

cd ./3rdParty/zlog 
make clean
cd ../../
cd ./scripts/

echo Cleaning libzlog ...
sh ./clean_libzlog.sh

echo Done
