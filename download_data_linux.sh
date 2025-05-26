#!/bin/bash

TARGET_DIR="./"
mkdir -p "$TARGET_DIR"

# Zenodo folder URLs (returning ZIP archives)
FOLDERS=(
  "https://zenodo.org/records/14967806/files/data.zip?download=1"
)

# Download folder ZIPs
for url in "${FOLDERS[@]}"; do
  echo "Downloading folder ZIP: $url"
  pushd "$TARGET_DIR" > /dev/null
  # Download as data.zip
  curl -L -o data.zip "$url"
  # Check if file is a valid zip
  if file data.zip | grep -q 'Zip archive data'; then
    echo "Successfully downloaded ZIP: data.zip"
  else
    echo "Error: Downloaded file is not a valid ZIP archive!"
    rm -f data.zip
    popd > /dev/null
    exit 1
  fi
  # Extract directly into current directory (overwrite if exists)
  echo "Extracting data.zip into current directory ..."
  unzip -o data.zip -d .
  # Remove the zip after extraction
  rm -f data.zip
  popd > /dev/null
done

