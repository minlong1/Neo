#!/bin/bash

git submodule update --init --recursive
cd ../contrib/sabcor/
make 
