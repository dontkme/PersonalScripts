#!/bin/bash
grep -P -R -H -B 10 -A 4 '^\s+N50' ./*/ass.err >N50all.txt
