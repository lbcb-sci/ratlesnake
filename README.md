# Ratlesnake

[![Build status for c++/clang++](https://travis-ci.org/lbcb-sci/ratlesnake.svg?branch=master)](https://travis-ci.org/lbcb-sci/ratlesnake)

Tool for assessing third generation sequencing data for genome assembly.

## Description

Ratlesnake was designed to aid development of new assembly algorithms. Given a reference genome, it calculates the most contiguous assembly possible for each chromosome separately. In addition, it classifies sequences into distinct classes and annotates related events, such as breaking points in chimeric sequences, inclusion intervals of contained sequences and repetitive genomic regions in sequences overlapping them.

Ratlesnake takes as input two files: third generation sequencing data in FASTA/FASTQ format and the corresponding reference genome, as well in FASTA/FASTQ format. All input files **can be compressed with gzip** (which will have impact on parsing time). Output is a set of files which are described bellow:

* ratlesnake.gfa
    * contains the assembly graph of each chromosome (multiple parts per chromosome possible)
* ratlesnake_solid.fasta
    * contains uncontained and non-chimeric sequences used for building the assembly graph
    * headers have the following custom tags:
        * XB:i:\<int\> - valid region begin
        * XE:i:\<int\> - valid region end
        * YB:i:\<int\> - chimeric region begin
        * YE:i:\<int\> - chimeric region end
        * ZB:i:\<int\> - repetitive region begin
        * ZE:i:\<int\> - repetitive region end
* ratlesnake_contained.fasta
    * contains sequences that are contained in others
    * headers have the following custom tags:
        * XB/XE - as above
        * XI:i:\<int\> - zero-based index of sequence in which it is contained
* ratlesnake_chimeric.fasta
    * contains sequences that are chimeric
    * headers have the following custom tags:
        * YB/YE - as above
* ratlesnake_repeats.fasta
    * contains sequences that are overlapping repetitive genomic regions
        * ZB/ZE - as above

## Dependencies

1. gcc 4.8+ or clang 3.4+
2. cmake 3.2+

## Installation

CmakeLists is provided in the project root folder. By running the following commands:

```bash
git clone --recursive https://github.com/lbcb-sci/ratlesnake ratlesnake
cd ratlesnake
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
a executable named `ratlesnake` will appear in the `build/lib` directory.

Optionally, you can run `sudo make install` to install Ratlesnake to your machine.

***Note***: if you omitted `--recursive` from `git clone`, run `git submodule init` and `git submodule update` before proceeding with compilation.

## Usage

Usage of ratlesnake is as following:

```bash
usage: ratlesnake [options ...] <sequences> <reference>

    <sequences>
        input file in FASTA/FASTQ format (can be compressed with gzip)
        containing sequences
    <reference>
        input file in FASTA/FASTQ format (can be compressed with gzip)
        containing a reference genome

    options:
        -t, --threads <int>
            default: 1
            number of threads
        --version
            prints the version number
        -h, --help
            prints the usage
```

## Contact information

For additional information, help and bug reports please send an email to one of the following: robert.vaser@fer.hr.

## Acknowledgement

This work has been supported in part by the European Regional Development Fund under the grant KK.01.1.1.01.0009 (DATACROSS) and in part by the Croatian Science Foundation under the project Single genome and metagenome assembly (IP-2018-01-5886).
