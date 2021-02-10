# Ratlesnake

[![Build status for c++/clang++](https://travis-ci.org/lbcb-sci/ratlesnake.svg?branch=master)](https://travis-ci.org/lbcb-sci/ratlesnake)

Tool for assessing third generation sequencing data for genome assembly.

## Usage

Ratlesnake was designed to aid development of new assembly algorithms. Given a reference genome, it calculates the most contiguous assembly possible for each chromosome separately. In addition, it classifies sequences into distinct classes and annotates related events, such as breaking points in chimeric sequences, inclusion intervals of contained sequences and repetitive genomic regions in sequences overlapping them.

Ratlesnake takes as input two files: third generation sequencing data in FASTA/FASTQ format and the corresponding reference genome, as well in FASTA/FASTQ format. Output is a set of files which are described bellow:

* `ratlesnake.gfa`
    * contains the assembly graph of each chromosome (multiple parts per chromosome possible)
* `ratlesnake_solid.fasta`
    * contains uncontained and non-chimeric sequences used for building the assembly graph
    * headers have the following custom tags:
        * XB:i:\<int\> - valid region begin
        * XE:i:\<int\> - valid region end
        * YB:i:\<int\> - chimeric region begin
        * YE:i:\<int\> - chimeric region end
        * ZB:i:\<int\> - repetitive region begin
        * ZE:i:\<int\> - repetitive region end
* `ratlesnake_contained.fasta`
    * contains sequences that are contained in others
    * headers have the following custom tags:
        * XB/XE - as above
        * XI:i:\<int\> - zero-based index of sequence in which it is contained
* `ratlesnake_chimeric.fasta`
    * contains sequences that are chimeric
    * headers have the following custom tags:
        * YB/YE - as above
* `ratlesnake_repeats.fasta`
    * contains sequences that are overlapping repetitive genomic regions
        * ZB/ZE - as above

To build ratlesnake run the following commands:

```bash
git clone https://github.com/lbcb-sci/ratlesnake && cd ratlesnake && mkdir build &&  cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
```

which will create the ratlesnake executable (installable with `make install`). Running it will display the following usage:

```bash
usage: ratlesnake [options ...] <sequences> <reference>

  <sequences>/<reference>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options:
    -t, --threads <int>
      default: 1
      number of threads
    --version
      prints the version number
    -h, --help
      prints the usage
```

#### Dependencies

- gcc 4.8+ | clang 3.5+
- cmake 3.11+
- zlib 1.2.8+

###### Hidden
- lbcb-sci/ram 2.0.0

## Acknowledgement

This work has been supported in part by the European Regional Development Fund under the grant KK.01.1.1.01.0009 (DATACROSS) and in part by the Croatian Science Foundation under the project Single genome and metagenome assembly (IP-2018-01-5886).
