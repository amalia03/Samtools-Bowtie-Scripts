#!/bin/bash

###Assuming that you have multiple cigar files, the script is written as a loop. 

sam_loc="/directory/to/sam/file"
out_loc="/output/directory/"


cd $star_loc

id="$(ls)"

for i in $id;
do
    mkdir $out_loc/$i
#Declares each element in the loop and removes the header with a simple grep command
    ls $i/*sam
    echo
    grep -v '^@' $i/*sam > $out_loc/$i/$i\_dehead.sam
    cd ..

done
