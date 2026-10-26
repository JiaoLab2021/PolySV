# PolySV

## Introduction
A high-quality SV detection, genotyping and optimization pipeline


## Requirements
Please note the following requirements before building and running the software:
- `linux` operating system
- Python3 (version: `3.9.20`)
- The following dependencies must also be installed: [samtools](https://github.com/samtools/samtools)(version: `1.18`), [Minimap2](https://github.com/lh3/minimap2) (version: `2.28-r1209`), [MuMmer](https://github.com/mummer4/mummer) (version: `4.0.0rc1`), [Sniffles2](https://github.com/fritzsedlazeck/Sniffles), (version: `2.2`), [SVIM](https://github.com/eldariont/svim) (version: `2.0.0`), [cuteSV2](https://github.com/tjiangHIT/cuteSV) (version: `2.1.0`), [debreak](https://github.com/Maggi-Chen/DeBreak) (version `1.0.2`), [svim-asm](https://github.com/eldariont/svim-asm) (version: `1.0.3`), [Assemblytics](https://github.com/marianattestad/assemblytics) (version: `1.2.1`), [Jasmine](https://github.com/mkirsche/Jasmine) (version `1.1.5`), [bcftools](https://github.com/samtools/bcftools) (version `1.15.1`)


## Installation
### Install via conda
The most recommended installation method is to install PolySV and all its dependencies at once through Anaconda
```
# Create a new Conda environment
conda create -n polysv_env
# Activate the environment
conda activate polysv_env
# Install PolySV
conda install -c hejiabao -c bioconda -c conda-forge polysv
polysv -h
```
### Build and Install from Source
Use the following script to build the software:
1. Make sure all required dependencies are installed before execution
```
# Install all dependencie
conda create -n polysv_env
conda activate polysv_env
conda install make=4.4.1 python=3.9 r-base=4.2.2 samtools=1.18 minimap2=2.28 mummer4=4.0.0rc1 \
cutesv sniffles svim svim-asm jasminesv assemblytics=1.2.1 -c bioconda -c conda-forge
git clone https://github.com/Maggi-Chen/DeBreak.git
cd DeBreak
echo 'export PATH="$PATH:'$(pwd)'"' >> ~/.bashrc
source ~/.bashrc
```
2. Get the PolySV source code
```
git clone git@github.com:JiaoLab2021/PolySV.git
cd PolySV && ln -s polysv.py polysv
echo 'export PATH="$PATH:'$(pwd)'"' >> ~/.bashrc
source ~/.bashrc
polysv -h
```

## Test
```
# Install
samtools
minimap2 -h
nucmer -h
delta-filter -h
show-coords -h
sniffles -h
cuteSV -h
svim -h
svim-asm -h
Assemblytics
jasmine -h
```

# Run
```
cd test
polysv -r refgenome.fa -tgs sample.fastq.gz -tgsdt ccs -genome sample.fa -n sample -o sample.vcf.gz --restart
```

## Citation
```
When using the following tools, please cite the corresponding articles:
- 

