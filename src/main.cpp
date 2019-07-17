#include <getopt.h>

#include <cstdint>
#include <iostream>
#include <string>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <math.h>
#include <set>

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "ram/ram.hpp"

static const std::string version = "v0.0.1";

static struct option options[] = {
    {"threads", required_argument, nullptr, 't'},
    {"version", no_argument, nullptr, 'v'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, 0, nullptr, 0}
};

void help();

inline bool isSuffix(const std::string& src, const std::string& suffix) {
    return src.size() < suffix.size() ? false :
        src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::unique_ptr<bioparser::Parser<ram::Sequence>> createSequenceParser(
    const std::string& path) {

    if (isSuffix(path, ".fasta")    || isSuffix(path, ".fa") ||
        isSuffix(path, ".fasta.gz") || isSuffix(path, ".fa.gz")) {
        return bioparser::createParser<bioparser::FastaParser, ram::Sequence>(path);
    }
    if (isSuffix(path, ".fastq")    || isSuffix(path, ".fq") ||
        isSuffix(path, ".fastq.gz") || isSuffix(path, ".fq.gz")) {
        return bioparser::createParser<bioparser::FastqParser, ram::Sequence>(path);
    }

    std::cerr << "[ratlesnake::] error: file " << path
              << " has unsupported format extension (valid extensions: .fasta, "
              << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!"
              << std::endl;
    return nullptr;
}

// void explore(std::vector<std::unique_ptr<ram::Sequence>>& src);

// void transform(std::vector<std::unique_ptr<ram::Sequence>>& src);

struct Annotation {
    std::vector<std::uint64_t> inclusion_interval;
    std::vector<std::uint64_t> chimeric_regions;
    std::vector<std::uint64_t> repetitive_regions;
};

std::vector<Annotation> annotate(
    const std::vector<std::unique_ptr<ram::Sequence>>& src,
    const std::vector<std::unique_ptr<ram::Sequence>>& dst,
    std::uint32_t num_threads);

void reconstruct(const std::vector<Annotation>& annotations,
    const std::vector<std::unique_ptr<ram::Sequence>>& src,
    const std::vector<std::unique_ptr<ram::Sequence>>& dst,
    std::uint32_t num_threads);

int main(int argc, char** argv) {

    std::uint32_t num_threads = 1;

    std::vector<std::string> input_paths;

    char argument;
    while ((argument = getopt_long(argc, argv, "t:h", options, nullptr)) != -1) {
        switch (argument) {
            case 't': num_threads = atoi(optarg); break;
            case 'v': std::cout << version << std::endl; return 0;
            case 'h': help(); return 0;
            default: return 1;
        }
    }

    for (std::int32_t i = optind; i < argc; ++i) {
        input_paths.emplace_back(argv[i]);
    }

    if (input_paths.size() < 2) {
        std::cerr << "[ratlesnake::] error: missing input file(s)!" << std::endl;
        help();
        return 1;
    }

    std::unique_ptr<bioparser::Parser<ram::Sequence>> sparser = createSequenceParser(input_paths[0]);
    if (sparser == nullptr) {
        return 1;
    }

    std::unique_ptr<bioparser::Parser<ram::Sequence>> rparser = createSequenceParser(input_paths[1]);
    if (rparser == nullptr) {
        return 1;
    }

    std::vector<std::unique_ptr<ram::Sequence>> sequences;
    sparser->parse(sequences, -1);

    std::vector<std::unique_ptr<ram::Sequence>> references;
    rparser->parse(references, -1);

    auto annotations = annotate(sequences, references, num_threads);

    reconstruct(annotations, sequences, references, num_threads);

    return 0;
}

void help() {
    std::cout <<
        "usage: ratlesnake [options ...] <sequences> <reference>\n"
        "\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing sequences\n"
        "    <reference>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing a reference genome\n"
        "\n"
        "    options:\n"
        "        -t, --threads <int>\n"
        "            default: 1\n"
        "            number of threads\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n";
}

std::vector<Annotation> annotate(
    const std::vector<std::unique_ptr<ram::Sequence>>& src,
    const std::vector<std::unique_ptr<ram::Sequence>>& dst,
    std::uint32_t num_threads) {

    auto thread_pool = thread_pool::createThreadPool(num_threads);

    std::vector<Annotation> annotations(src.size() + dst.size());

    ram::MinimizerEngine minimizer_engine(15, 5, num_threads);
    minimizer_engine.minimize(dst.begin(), dst.end());
    minimizer_engine.filter(0.001);

    std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;
    for (std::uint32_t i = 0; i < src.size(); ++i) {
        thread_futures.emplace_back(thread_pool->submit(
            [&] (std::uint32_t i) -> std::vector<ram::Overlap> {
                auto overlaps = minimizer_engine.map(src[i], false, false);
                std::vector<ram::Overlap> repetitive_overlaps;
                if (overlaps.size() < 2) {
                    return repetitive_overlaps;
                }

                std::sort(overlaps.begin(), overlaps.end(),
                    [] (const ram::Overlap& lhs, const ram::Overlap& rhs) -> bool {
                        return lhs.q_begin < rhs.q_begin ||
                            (lhs.q_begin == rhs.q_begin && lhs.q_end > rhs.q_end);
                    });

                // annotate chimeric regions
                std::uint32_t k = 0;
                for (std::uint32_t j = 1; j < overlaps.size(); ++j) {
                    if (overlaps[k].q_end < overlaps[j].q_end) {
                        bool is_chimeric = true;
                        if (overlaps[k].t_id == overlaps[j].t_id) {
                            std::uint32_t q_gap = abs(overlaps[j].q_begin - overlaps[k].q_end);
                            std::uint32_t t_gap = abs(overlaps[j].t_begin - overlaps[k].t_end);
                            is_chimeric &= std::min(q_gap, t_gap) <
                                0.88 * std::max(q_gap, t_gap);
                        }
                        if (is_chimeric) {
                            annotations[src[i]->id].chimeric_regions.emplace_back(
                                overlaps[j].q_begin > overlaps[k].q_end ?
                                (static_cast<uint64_t>(overlaps[k].q_end) << 32 | overlaps[j].q_begin) :
                                (static_cast<uint64_t>(overlaps[j].q_begin) << 32 | overlaps[k].q_end));
                        }
                        k = j;
                    } else {
                        repetitive_overlaps.emplace_back(overlaps[j]);
                    }
                }

                if (repetitive_overlaps.empty()) {
                    return repetitive_overlaps;
                }

                // annotate repetitive regions
                k = 0;
                std::uint32_t q_end = repetitive_overlaps[k].q_end;
                for (std::uint32_t j = 1; j < repetitive_overlaps.size(); ++j) {
                    if (q_end < repetitive_overlaps[j].q_begin) {
                        annotations[src[i]->id].repetitive_regions.emplace_back(
                            static_cast<uint64_t>(repetitive_overlaps[k].q_begin) << 32 | q_end);
                        k = j;
                    }
                    q_end = std::max(q_end, repetitive_overlaps[j].q_end);
                }
                annotations[src[i]->id].repetitive_regions.emplace_back(
                    static_cast<uint64_t>(repetitive_overlaps[k].q_begin) << 32 | q_end);

                return repetitive_overlaps;
            }
        , i));
    }

    std::vector<ram::Overlap> repetitive_overlaps;
    for (std::uint32_t i = 0; i < thread_futures.size(); ++i) {
        thread_futures[i].wait();
        auto overlaps = thread_futures[i].get();
        repetitive_overlaps.insert(repetitive_overlaps.end(), overlaps.begin(),
            overlaps.end());
    }

    if (repetitive_overlaps.empty()) {
        return annotations;
    }

    // annotate repetitive regions on references
    std::sort(repetitive_overlaps.begin(), repetitive_overlaps.end(),
        [] (const ram::Overlap& lhs, const ram::Overlap& rhs) -> bool {
            return lhs.t_id < rhs.t_id ||
                (lhs.t_id == rhs.t_id && lhs.t_begin < rhs.t_begin) ||
                (lhs.t_id == rhs.t_id && lhs.t_begin == rhs.t_begin && lhs.t_end > rhs.t_end);
        });

    std::uint32_t j = 0;
    std::uint32_t t_end = repetitive_overlaps[j].q_end;
    for (std::uint32_t i = 1; i < repetitive_overlaps.size(); ++i) {
        if (repetitive_overlaps[i - 1].t_id != repetitive_overlaps[i].t_id ||
            t_end < repetitive_overlaps[i].t_begin) {

            annotations[repetitive_overlaps[i - 1].t_id].repetitive_regions.emplace_back(
                static_cast<uint64_t>(repetitive_overlaps[j].t_begin) << 32 | t_end);
            j = i;
            t_end = repetitive_overlaps[i].t_end;
        } else {
            t_end = std::max(t_end, repetitive_overlaps[i].t_end);
        }
    }
    annotations[repetitive_overlaps[j].t_id].repetitive_regions.emplace_back(
        static_cast<uint64_t>(repetitive_overlaps[j].t_begin) << 32 | t_end);

    return annotations;
}

void reconstruct(const std::vector<Annotation>& annotations,
    const std::vector<std::unique_ptr<ram::Sequence>>& src,
    const std::vector<std::unique_ptr<ram::Sequence>>& dst,
    std::uint32_t num_threads) {

    auto thread_pool = thread_pool::createThreadPool(num_threads);

    ram::MinimizerEngine minimizer_engine(15, 5, num_threads);
    minimizer_engine.minimize(dst.begin(), dst.end());
    minimizer_engine.filter(0.001);

    std::vector<std::future<ram::Overlap>> thread_futures;
    for (std::uint32_t i = 0; i < src.size(); ++i) {
        thread_futures.emplace_back(thread_pool->submit(
            [&] (std::uint32_t i) -> ram::Overlap {
                auto overlaps = minimizer_engine.map(src[i], false, false);
                if (overlaps.empty()) {
                    return ram::Overlap(-1, -1, -1, -1, -1, -1);
                }
                std::sort(overlaps.begin(), overlaps.end(),
                    [] (const ram::Overlap& lhs, const ram::Overlap& rhs) -> bool {
                        return lhs.q_end - lhs.q_begin > rhs.q_end - rhs.q_begin;
                    });
                return overlaps.front();
            }
        , i));
    }

    std::vector<std::uint64_t> sources;
    std::vector<ram::Overlap> overlaps;
    for (std::uint32_t i = 0; i < thread_futures.size(); ++i) {
        thread_futures[i].wait();
        auto overlap = thread_futures[i].get();
        if (overlap.t_id != -1ULL) {
            sources.emplace_back(src[i]->id);
            overlaps.emplace_back(overlap);
        }
    }

    // remove contained reads
    // TODO: sort sources!
    std::sort(overlaps.begin(), overlaps.end(),
        [] (const ram::Overlap& lhs, const ram::Overlap& rhs) -> bool {
            return lhs.t_id < rhs.t_id ||
                (lhs.t_id == rhs.t_id && lhs.t_begin < rhs.t_begin) ||
                (lhs.t_id == rhs.t_id && lhs.t_begin == rhs.t_begin && lhs.t_end > rhs.t_end);
        });

    for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
        for (std::uint32_t j = i + 1; j < overlaps.size(); ++j) {
            if (overlaps[i].t_id != overlaps[j].t_id ||
                overlaps[i].t_end < overlaps[j].t_end) {
                i = j - 1;
                break;
            }
            sources[j] = -1;
        }
    }

    {
        std::vector<uint64_t> s;
        std::vector<ram::Overlap> o;
        for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
            if (sources[i] != -1ULL) {
                s.emplace_back(sources[i]);
                o.emplace_back(overlaps[i]);
            }
        }
        s.swap(sources);
        o.swap(overlaps);
    }

    // create assembly graph
    using Node = std::vector<std::pair<std::uint64_t, std::uint32_t>>;
    std::vector<Node> nodes(src.size());

    for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
        std::cout << "S\t" << src[sources[i]]->name << "\t"
                  << src[sources[i]]->data << "\tLN:i:"
                  << src[sources[i]]->data.size() << "\tRC:i:1" << std::endl;
        for (std::uint32_t j = i + 1; j < overlaps.size(); ++j) {
            if (overlaps[i].t_id != overlaps[j].t_id ||
                overlaps[i].t_end < overlaps[j].t_begin) {
                break;
            }
            nodes[i].emplace_back(sources[j], overlaps[i].t_end - overlaps[j].t_begin);
            std::cout << "L\t" << src[sources[i]]->name << "\t"
                      << (overlaps[i].strand ? "+\t" : "-\t")
                      << src[sources[j]]->name << "\t"
                      << (overlaps[j].strand ? "+\t" : "-\t")
                      << overlaps[i].t_end - overlaps[j].t_begin << "M" << std::endl;
        }
    }

    for (const auto& it: nodes) {

    }
}
