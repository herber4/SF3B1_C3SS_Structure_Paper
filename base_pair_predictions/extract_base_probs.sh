#!/bin/bash

# Loop through all *.ps files
for file in *.ps; do
    # Check if the file exists and is readable
    if [ -r "$file" ]; then
        # Extract lines between "%start of base pair probability data" and "showpage"
        awk '/%start of base pair probability data/{flag=1; next} /showpage/{flag=0} flag' "$file" > "${file%.ps}_extracted.txt"
        echo "Extracted lines from $file"
    else
        echo "Error: $file is not readable or does not exist"
    fi
done
