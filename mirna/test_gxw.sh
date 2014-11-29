#!/bin/bash

regex=$(./gxw_to_seq.pl $1 1)

while read L; do
    tabless=$(echo $L | sed "s/ [0-9].*//")
    if echo $tabless | grep "$regex" >/dev/null; then
	echo $tabless | sed "s/ .*$//"
    fi
done < "$2"
