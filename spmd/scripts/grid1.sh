#!/bin/bash

./times.sh 1 ./run.sh ~/mpd1.hosts $2 $3 4 4 2>&1
./times.sh $1 ./run.sh ~/mpd1.hosts $2 $3 4 4 2>&1 | grep Min
echo

./times.sh 1 ./run.sh ~/mpd1.hosts $2 $3 4 3 2>&1
./times.sh $1 ./run.sh ~/mpd1.hosts $2 $3 4 3 2>&1 | grep Min
echo

./times.sh 1 ./run.sh ~/mpd1.hosts $2 $3 3 3 2>&1
./times.sh $1 ./run.sh ~/mpd1.hosts $2 $3 3 3 2>&1 | grep Min
echo

./times.sh 1 ./run.sh ~/mpd1.hosts $2 $3 3 2 2>&1
./times.sh $1 ./run.sh ~/mpd1.hosts $2 $3 3 2 2>&1 | grep Min
echo

./times.sh 1 ./run.sh ~/mpd1.hosts $2 $3 2 2 2>&1
./times.sh $1 ./run.sh ~/mpd1.hosts $2 $3 2 2 2>&1 | grep Min
echo

./times.sh 1 ./run.sh ~/mpd1.hosts $2 $3 2 1 2>&1
./times.sh $1 ./run.sh ~/mpd1.hosts $2 $3 2 1 2>&1 | grep Min
echo

./times.sh 1 ./run.sh ~/mpd1.hosts $2 $3 1 1 2>&1
./times.sh $1 ./run.sh ~/mpd1.hosts $2 $3 1 1 2>&1 | grep Min
