#!/bin/bash
# usage: sh auroc_golds_predictions.R
# To be run from the In_out directory (like In_out_degree_2_num_nodes_50)
# From there it exploits the Goldstands and the different subdirectories of Predictions
# via paths are relative (including ./ but not ending with /)

# GOLDS=$1
# PREDICTIONS=$2


GOLDS="./Goldstandards"
PREDICTIONS="./Predictions/dynGENIE3"

Rscript ~/Desktop/Benchmarking/In_out_degree_2_num_nodes_100/auroc/auroc_golds_predictions.R $GOLDS $PREDICTIONS

getwd()



