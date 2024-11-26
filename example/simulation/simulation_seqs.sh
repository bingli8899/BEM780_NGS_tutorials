#!/bin/bash
# Use it as: ./simulation_seqs.sh <executable_folder> <configuration_file> <output_directory>

# Assign arguments to variables
EXEC_FOLDER="$1"
SIMPHY_EXEC="$EXEC_FOLDER/simphy"
SEQGEN_EXEC="$EXEC_FOLDER/seq-gen"
CONFIG_DIR="$2"
CONFIG_FILE="$2/configuration.txt"
OUTPUT_DIR="$3"
SEQGEN_OUTPUT_DIR="$OUTPUT_DIR/seqgen-out"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$SEQGEN_OUTPUT_DIR"

echo "Running SimPhy..."
"$SIMPHY_EXEC" -i "$CONFIG_FILE" -o "$OUTPUT_DIR"
echo "Simphy running done! Let's run seq-gen"


for i in $(seq -w 1 10); do # 10 is hard-coded here! 
	
	for j in $(seq 1 2); do # Only two replicates, so this is hard-coded
	
	SIMPHY_TREE_DIR="$OUTPUT_DIR/$j" 
  	TREE_FILE="$SIMPHY_TREE_DIR/g_trees$i.trees"
  	SEQGEN_OUTPUT_DIR_REP="$SEQGEN_OUTPUT_DIR/$j"
    SEQGEN_OUTPUT="$SEQGEN_OUTPUT_DIR_REP/gene$i.nex"
  	
  	mkdir -p "$SEQGEN_OUTPUT_DIR_REP"

  	"$SEQGEN_EXEC" -l1000 -mHKY -a0.5 -t2.0 -f0.25,0.25,0.25,0.25 -z 1234 -on < "$TREE_FILE" > "$SEQGEN_OUTPUT"
	done
done

echo "All done."
