#include <getopt.h>

#include <cstdint>
#include <iostream>
#include <string>
#include <memory>
#include <unordered_map>

#include "bioparser/bioparser.hpp"

static const std::string version = "v0.0.0";

static struct option options[] = {
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

void help();

struct Sequence {
    Sequence(const char* name, std::uint32_t name_length,
        const char* data, std::uint32_t data_length)
            : name(name, name_length), data(data, data_length) {
    }

    Sequence(const char* name, std::uint32_t name_length,
        const char* data, std::uint32_t data_length,
        const char*, std::uint32_t)
            : Sequence(name, name_length, data, data_length) {
    }

    ~Sequence() {
    }

    std::string name;
    std::string data;
};

struct Annotation {
    Annotation() {
    }
    ~Annotation() {
    }

    std::uint64_t valid;
    std::vector<std::uint64_t> inclusion;
    std::vector<std::uint64_t> chimeric;
    std::vector<std::uint64_t> repeat;
};

struct Overlap {
    Overlap(
        const char* q_name, std::uint32_t q_name_length,
        std::uint32_t q_length,
        std::uint32_t q_begin,
        std::uint32_t q_end,
        char orientation,
        const char* t_name, std::uint32_t t_name_length,
        std::uint32_t t_length,
        std::uint32_t t_begin,
        std::uint32_t t_end,
        std::uint32_t,
        std::uint32_t,
        std::uint32_t)
            : q_id(),
            q_length(q_length),
            q_begin(q_begin),
            q_end(q_end),
            t_id(),
            t_length(t_length),
            t_begin(t_begin),
            t_end(t_end),
            strand(orientation == '-') {

        std::string q(q_name, q_name_length);
        if (name_to_id.find(q) != name_to_id.end()) {
            q_id = name_to_id[q];
        } else {
            q_id = name_to_id.size();
            name_to_id[q] = name_to_id.size();
        }

        std::string t(t_name, t_name_length);
        if (name_to_id.find(t) != name_to_id.end()) {
            t_id = name_to_id[t];
        } else {
            t_id = name_to_id.size();
            name_to_id[t] = name_to_id.size();
        }
    }

    static std::unordered_map<std::string, std::uint32_t> name_to_id;

    std::uint64_t q_id;
    std::uint32_t q_length;
    std::uint32_t q_begin;
    std::uint32_t q_end;
    std::uint64_t t_id;
    std::uint32_t t_length;
    std::uint32_t t_begin;
    std::uint32_t t_end;
    bool strand;
};

std::unordered_map<std::string, std::uint32_t> Overlap::name_to_id;

inline bool isSuffix(const std::string& src, const std::string& suffix) {
    return src.size() < suffix.size() ? false :
        src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
}

void annotate(std::vector<Annotation>& dst,
    const std::vector<std::unique_ptr<Overlap>>& overlaps) {
}

void reconstruct(std::vector<Annotation>& dst,
    std::vector<std::unique_ptr<Overlap>>& overlaps) {
}

int main(int argc, char** argv) {

    std::vector<std::string> input_paths;

    char argument;
    while ((argument = getopt_long(argc, argv, "h", options, nullptr)) != -1) {
        switch (argument) {
            case 'v': std::cout << version << std::endl; return 0;
            case 'h': help(); return 0;
            default: return 1;
        }
    }

    for (std::int32_t i = optind; i < argc; ++i) {
        input_paths.emplace_back(argv[i]);
    }

    if (input_paths.size() < 3) {
        std::cerr << "[ratlesnake::] error: missing input file(s)!" << std::endl;
        help();
        return 1;
    }

    std::unique_ptr<bioparser::Parser<Sequence>> sparser = nullptr;

    if (isSuffix(input_paths[0], ".fasta") || isSuffix(input_paths[0], ".fa") ||
        isSuffix(input_paths[0], ".fasta.gz") || isSuffix(input_paths[0], ".fa.gz")) {
        sparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            input_paths[0]);
    } else if (isSuffix(input_paths[0], ".fastq") || isSuffix(input_paths[0], ".fq") ||
               isSuffix(input_paths[0], ".fastq.gz") || isSuffix(input_paths[0], ".fq.gz")) {
        sparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            input_paths[0]);
    } else {
        std::cerr << "[ratlesnake::] error: file " << input_paths[0] <<
            " has unsupported format extension (valid extensions: .fasta, "
            ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!" <<
        std::endl;
        return 1;
    }

    std::unique_ptr<bioparser::Parser<Overlap>> oparser = nullptr;

    if (isSuffix(input_paths[1], ".paf") || isSuffix(input_paths[1], ".paf.gz")) {
        oparser = bioparser::createParser<bioparser::PafParser, Overlap>(
            input_paths[1]);
    } else {
        std::cerr << "[ratlesnake::] error: file " << input_paths[1] <<
            " has unsupported format extension (valid extensions: .fasta, "
            ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!" <<
        std::endl;
        return 1;
    }

    std::unique_ptr<bioparser::Parser<Sequence>> rparser = nullptr;

    if (isSuffix(input_paths[2], ".fasta") || isSuffix(input_paths[2], ".fa") ||
        isSuffix(input_paths[2], ".fasta.gz") || isSuffix(input_paths[2], ".fa.gz")) {
        rparser = bioparser::createParser<bioparser::FastaParser, Sequence>(
            input_paths[2]);
    } else if (isSuffix(input_paths[2], ".fastq") || isSuffix(input_paths[2], ".fq") ||
               isSuffix(input_paths[2], ".fastq.gz") || isSuffix(input_paths[2], ".fq.gz")) {
        rparser = bioparser::createParser<bioparser::FastqParser, Sequence>(
            input_paths[2]);
    } else {
        std::cerr << "[ratlesnake::] error: file " << input_paths[2] <<
            " has unsupported format extension (valid extensions: .fasta, "
            ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!" <<
        std::endl;
        return 1;
    }

    std::vector<std::unique_ptr<Overlap>> overlaps;
    oparser->parse(overlaps, -1);

    std::vector<Annotation> annotations(Overlap::name_to_id.size());

    annotate(annotations, overlaps);
    reconstruct(annotations, overlaps);

    return 0;
}

void help() {
    std::cout <<
        "usage: ratlesnake [options ...] <sequences> <overlaps> <reference>\n"
        "\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing third generation sequences\n"
        "    <overlaps>\n"
        "        input file in PAF format (can be compressed with gzip)\n"
        "        containing overlaps between sequences and reference\n"
        "    <reference>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing the reference genome\n"
        "\n"
        "    options:\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n";
}
