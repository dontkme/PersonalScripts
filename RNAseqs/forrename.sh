#!/bin/bash
while read old new; do
    if [[ -f "$old" ]]; then
        mv "$old" "$new"
        echo "Renamed: $old -> $new"
    else
        echo "File not found: $old"
    fi
done < renamelist.txt
