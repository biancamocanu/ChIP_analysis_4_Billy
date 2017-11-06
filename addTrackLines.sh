#! /bin/sh

INFILE=$1
NAME=$2

PREFIX=$(echo $INFILE | cut -d "." -f 1)
echo Adding tracklines and this name: $NAME

awk -v NAME="$NAME" 'BEGIN { print "browser position chr11:5,289,521-5,291,937"
	print "track type=bedGraph name=\""NAME"\" description=\""NAME"\" visibility=full windowingFunction=maximum color=0,125,0"}
	{print $0}' $INFILE > ${PREFIX}_header.bedgraph

