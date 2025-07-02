#!/bin/bash
# transform_files.sh - Prepares single-cell data files for Seurat compatibility

INPUT_DIR=$1
OUTPUT_DIR=$2

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Transforming files from $INPUT_DIR to $OUTPUT_DIR..."

# Process barcodes file
if [ -f "$INPUT_DIR/GSE165686_barcodes.tsv.gz" ]; then
  echo "Transforming barcodes file..."
  gunzip -c "$INPUT_DIR/GSE165686_barcodes.tsv.gz" | tail -n +2 > "$OUTPUT_DIR/barcodes.tsv"
  gzip -c "$OUTPUT_DIR/barcodes.tsv" > "$OUTPUT_DIR/barcodes.tsv.gz"
  rm "$OUTPUT_DIR/barcodes.tsv"
fi

# Process features file
if [ -f "$INPUT_DIR/GSE165686_features.tsv.gz" ]; then
  echo "Transforming features file..."
  gunzip -c "$INPUT_DIR/GSE165686_features.tsv.gz" | tail -n +2 > "$OUTPUT_DIR/features.tsv"
  gzip -c "$OUTPUT_DIR/features.tsv" > "$OUTPUT_DIR/features.tsv.gz"
  rm "$OUTPUT_DIR/features.tsv"
fi

# Copy and rename matrix file
if [ -f "$INPUT_DIR/GSE165686_matrix.mtx.gz" ]; then
  echo "Copying matrix file..."
  cp "$INPUT_DIR/GSE165686_matrix.mtx.gz" "$OUTPUT_DIR/matrix.mtx.gz"
fi

# Copy metadata file
if [ -f "$INPUT_DIR/GSE165686_metadata.tsv.gz" ]; then
  echo "Copying metadata file..."
  cp "$INPUT_DIR/GSE165686_metadata.tsv.gz" "$OUTPUT_DIR/metadata.tsv.gz"
fi

echo "File transformation complete."