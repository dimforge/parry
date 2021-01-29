#! /bin/bash

tmp=`mktemp -d`

echo $tmp

cp -r src $tmp/.
cp -r LICENSE README.md $tmp/.

### Publish the 2D version.
sed 's#\.\./\.\./src#src#g' build/parry2d/Cargo.toml > $tmp/Cargo.toml
currdir=`pwd`
cd $tmp && cargo publish
cd $currdir

### Publish the 2D f64 version.
sed 's#\.\./\.\./src#src#g' build/parry2d-f64/Cargo.toml > $tmp/Cargo.toml
cd $tmp && cargo publish
cd $currdir


### Publish the 3D version.
sed 's#\.\./\.\./src#src#g' build/parry3d/Cargo.toml > $tmp/Cargo.toml
cp -r LICENSE README.md $tmp/.
cd $tmp && cargo publish
cd $currdir

### Publish the 3D f64 version.
sed 's#\.\./\.\./src#src#g' build/parry3d-f64/Cargo.toml > $tmp/Cargo.toml
cp -r LICENSE README.md $tmp/.
cd $tmp && cargo publish

rm -rf $tmp

