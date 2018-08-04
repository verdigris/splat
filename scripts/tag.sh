#!/bin/bash

set -e

last_tag=`git tag -l | grep build- | sort | tail -n 1`
n=`echo $last_tag | awk '{print substr($0,7)}'`
n="10#$n"
let 'n=n+1'
sed -i -e s/"^\(BUILD = \)\([0-9]*\)$"/"\1$n"/ splat/__init__.py
git add splat/__init__.py
m=`printf "%04d" $n`
new_tag="build-$m"
git commit -m "$new_tag"
echo "tag: $new_tag"
today=`date +%Y-%m-%d`
msg="Splat build #$m $today"
echo "msg: $msg"
git tag -a "$new_tag" -m "$msg"

exit 0
