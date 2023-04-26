#!/usr/bin/bash

tasks=$(ls -d final*)

for item in $tasks; do
    cp -a $item/${item}.orca.engrad genref
done
