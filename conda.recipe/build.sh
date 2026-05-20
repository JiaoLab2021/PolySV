#!/bin/bash
set -e -x # Exit on error, print commands

# Use $PREFIX which is the standard conda build environment variable
BIN_DIR=$PREFIX/bin

# 1. Install MUMmer4 from source
echo "[INFO] Installing MUMmer..."
wget -c https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar -zxvf mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure --prefix=$PREFIX # Install directly into the environment
make
make install
cd ..
rm -rf mummer-4.0.0rc1 mummer-4.0.0rc1.tar.gz
echo "[INFO] MUMmer4 installed!"

# 2. Install Assemblytics from source
echo "[INFO] Installing Assemblytics..."
git clone https://github.com/marianattestad/assemblytics.git
chmod a+x assemblytics/scripts/Assemblytics*
cp assemblytics/scripts/Assemblytics* $BIN_DIR/
rm -rf assemblytics
echo "[INFO] Assemblytics installed!"

# 3. Finally, install the hqsv package itself
echo "[INFO] Installing hqsv package..."
$PYTHON -m pip install . --no-deps --ignore-installed -vv
echo "[INFO] hqsv installed!"
