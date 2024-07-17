#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 2 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 input_file new_header [output_file] [-a|--append]"
    exit 1
fi

# Parse options
append=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -a|--append)
            append=true
            shift
            ;;
        *)
            break
            ;;
    esac
done

# Input file and new header
input_file="$1"
new_header="$2"

# Determine the output file
output_file="${3:-output.tsv}"

# Replace or append the header using awk
awk -v new_header="$new_header" -v append="$append" '
BEGIN {
    # Split the new header into an array
    split(new_header, new_cols);
}
NR == 1 {
    if (append) {
        # If appending, print the original header first
        print;
    }
    # Print the new header
    for (i = 1; i <= length(new_cols); i++) {
        printf "%s", new_cols[i];
        if (i < length(new_cols)) printf "\t";
    }
    if (!append) {
        printf "\n";
    }
    next;
}
{
    # Print the rest of the lines as is
    print;
}' "$input_file" > "$output_file"

echo "Header updated and saved to $output_file"
