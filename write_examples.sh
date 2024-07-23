#! /bin/bash

write_examples() {
    examples_path=$1
    output_path=$2

    echo >> $output_path

    find $examples_path -type f -iname '*.rs' -print0 | 
    while IFS= read -r -d '' line; do 
        example=$(basename ${line} .rs)
        echo "[[example]]" >> $output_path
        echo "name = \"$example\"" >> $output_path
        echo "path = \"examples/$example.rs\"" >> $output_path
        echo "doc-scrape-examples = true" >> $output_path
        echo >> $output_path
    done
}

write_examples ./crates/parry2d/examples ./crates/parry2d/Cargo.toml
write_examples ./crates/parry3d/examples ./crates/parry3d/Cargo.toml
