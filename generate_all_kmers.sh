#!/bin/bash

for i in {6..40}; 
do 
    # TODO make file structure make sense
    OUTDIR="counts_""$i""-mer"
    mkdir $OUTDIR
    python generate_AnnData.py $i $OUTDIR 
done