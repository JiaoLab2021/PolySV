#!/bin/bash
set -e -x

# Conda 的默认执行目录，放在这里的程序，系统随叫随到
BIN_DIR=$PREFIX/bin
mkdir -p $BIN_DIR

# 1. Install DeBreak
echo "[INFO] Cloning DeBreak..."
export GIT_SSL_NO_VERIFY=true
git clone --depth 1 https://github.com/Maggi-Chen/DeBreak.git
echo "[INFO] Installing DeBreak into Conda bin directory..."
chmod +x DeBreak/*.py
chmod +x DeBreak/debreak
cp DeBreak/debreak $BIN_DIR/
cp DeBreak/*.py $BIN_DIR/
rm -rf DeBreak
echo "[INFO] DeBreak installed successfully!"

# 2. polysv
echo "[INFO] Installing polysv package..."
$PYTHON -m pip install . --no-deps --ignore-installed -vv
echo "[INFO] polysv installed!"
