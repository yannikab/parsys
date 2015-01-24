#!/bin/bash

n=0;
m=$1;
shift 1;
while [ $n -lt $m ];
do
	$*;
	n=`expr $n + 1`;
done;
