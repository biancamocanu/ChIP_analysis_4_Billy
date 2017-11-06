#! /bin/sh

INFILE=$1
NAME=$2

PREFIX=$(echo $INFILE | cut -d "." -f 1)
echo Adding tracklines and this name: $NAME

awk -v NAME="$NAME" 'BEGIN { print "browser position chr11:5,289,521-5,291,937"
	 print "track type=bed name=\""NAME"\" description=\""NAME"\" visibility=squish autoScale=on colorByStrand=\"255,0,0 0,0,255\""}
{ print $0}' $INFILE > ${PREFIX}_header.bed
