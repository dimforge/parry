#! /bin/bash

tmp=`mktemp -d`

echo $tmp

cp -r src $tmp/.
cp -r LICENSE README.md $tmp/.

### Publish the 2D version.
sed 's#\.\./\.\./src#src#g' crates/parry2d/Cargo.toml > $tmp/Cargo.toml
rm -rf $tmp/examples
cp -r crates/parry2d/examples $tmp/examples
currdir=`pwd`
cd $tmp && cargo publish
cd $currdir

### Publish the 2D f64 version.
sed 's#\.\./\.\./src#src#g' crates/parry2d-f64/Cargo.toml > $tmp/Cargo.toml
cd $tmp && cargo publish
cd $currdir


### Publish the 3D version.
sed 's#\.\./\.\./src#src#g' crates/parry3d/Cargo.toml > $tmp/Cargo.toml
rm -rf $tmp/examples
cp -r crates/parry3d/examples $tmp/examples
cp -r LICENSE README.md $tmp/.
cd $tmp && cargo publish
cd $currdir

### Publish the 3D f64 version.
sed 's#\.\./\.\./src#src#g' crates/parry3d-f64/Cargo.toml > $tmp/Cargo.toml
cp -r LICENSE README.md $tmp/.
cd $tmp && cargo publish

rm -rf $tmp

