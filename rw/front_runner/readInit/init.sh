#!/bin/sh

# count the line number of the water depth field
lines=`wc -l initH2D | awk '{print $1}'`

awk -v lines=$lines '
BEGIN {
  print lines " 0 0"
} {
  print $1 " " $2 " " $3;
}' < initH2D | delaunay > test.gts
