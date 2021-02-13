# Ratlesnake

[![Build status for c++/clang++](https://travis-ci.org/lbcb-sci/ratlesnake.svg?branch=master)](https://travis-ci.org/lbcb-sci/ratlesnake)

Tool for assessing third generation sequencing data for genome assembly.

## Usage

Ratlesnake was designed to aid development of new assembly algorithms. Given a reference genome, it calculates the most contiguous assembly possible for each chromosome separately. In addition, it classifies sequences into distinct classes and annotates related events, such as breaking points in chimeric sequences, inclusion intervals of contained sequences and repetitive genomic regions in sequences overlapping them.

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
    -r, --result <int> (option can be used multiple times)
      default: 0
      result mode:
        0 - sequence statistics (histogram(s) on stdout)
        1 - sequence accuracy given reference (histogram on stdout)
        2 - sequence annotation given reference (FASTA files)
            - XI:i:<int> 0-based id
            - XC:i:<int> 0-based id of containing read
            - XB:i:<int> (XE:i:<int>) valid region begin (end)
            - YB:i:<int> (YE:i:<int>) chimeric region begin (end)
            - ZB:i:<int> (ZE:i:<int>) repetitive region begin (end)
        3 - reference reconstruction given sequences (GFA files)
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
