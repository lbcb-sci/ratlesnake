#include <getopt.h>

#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <numeric>
#include <iomanip>
#include <algorithm>

#include "bioparser/bioparser.hpp"
#include "thread_pool/thread_pool.hpp"
#include "ram/ram.hpp"

static const std::string version = "v0.0.3";

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
    std::vector<std::uint64_t> chimeric_regions;
    std::vector<std::uint64_t> repetitive_regions;
    bool is_junk = false;
};

std::vector<Annotation> annotate(
    const std::vector<std::unique_ptr<ram::Sequence>>& src,
    const std::vector<std::unique_ptr<ram::Sequence>>& dst,
    std::uint32_t num_threads);

void reconstruct(std::vector<Annotation>& annotations,
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

    std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(num_threads);

    std::vector<Annotation> annotations(src.size() + dst.size());

    auto minimizer_engine = ram::createMinimizerEngine(15, 5, thread_pool);
    minimizer_engine->minimize(dst.begin(), dst.end());
    minimizer_engine->filter(0.001);

    std::vector<std::future<std::vector<ram::Overlap>>> thread_futures;
    for (std::uint32_t i = 0; i < src.size(); ++i) {
        thread_futures.emplace_back(thread_pool->submit(
            [&] (std::uint32_t i) -> std::vector<ram::Overlap> {
                auto overlaps = minimizer_engine->map(src[i], false, false);
                std::vector<ram::Overlap> repetitive_overlaps;
                if (overlaps.size() < 2) {
                    if (overlaps.size() == 1) {
                        if ((overlaps[0].q_end - overlaps[0].q_begin) < 0.5 * src[i]->data.size() ||
                            static_cast<double>(overlaps[0].matches) / (overlaps[0].q_end - overlaps[0].q_begin) < 0.27) {
                                annotations[src[i]->id].is_junk = true;
                            }
                    }
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
                            std::uint32_t q_gap = abs(static_cast<int32_t>(overlaps[j].q_begin - overlaps[k].q_end));
                            std::uint32_t t_gap = abs(static_cast<int32_t>(overlaps[j].t_begin - overlaps[k].t_end));
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

                if (repetitive_overlaps.empty() || !annotations[src[i]->id].chimeric_regions.empty()) {
                    return std::vector<ram::Overlap>();
                }

                // annotate repetitive regions
                auto is_valid_repeat = [] (std::uint32_t q_begin, std::uint32_t q_end, std::uint32_t q_length) -> bool {
                    return q_end - q_begin > 500 &&
                        (q_begin < 0.05 * q_length || q_end > 0.95 * q_length || q_end - q_begin > 1500);
                };
                k = 0;
                std::uint32_t q_length = src[i]->data.size();
                std::uint32_t q_begin;
                std::uint32_t q_end = repetitive_overlaps[k].q_end;
                for (std::uint32_t j = 1; j < repetitive_overlaps.size(); ++j) {
                    if (q_end < repetitive_overlaps[j].q_begin) {
                        q_begin = repetitive_overlaps[k].q_begin;
                        if (is_valid_repeat(q_begin, q_end, q_length)) {
                            annotations[src[i]->id].repetitive_regions.emplace_back(
                                static_cast<uint64_t>(q_begin) << 32 | q_end);
                        }
                        k = j;
                    }
                    q_end = std::max(q_end, repetitive_overlaps[j].q_end);
                }
                q_begin = repetitive_overlaps[k].q_begin;
                if (is_valid_repeat(q_begin, q_end, q_length)) {
                    annotations[src[i]->id].repetitive_regions.emplace_back(
                        static_cast<uint64_t>(q_begin) << 32 | q_end);
                }

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
    std::uint32_t t_end = repetitive_overlaps[j].t_end;
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

void reconstruct(std::vector<Annotation>& annotations,
    const std::vector<std::unique_ptr<ram::Sequence>>& src,
    const std::vector<std::unique_ptr<ram::Sequence>>& dst,
    std::uint32_t num_threads) {

    std::shared_ptr<thread_pool::ThreadPool> thread_pool = thread_pool::createThreadPool(num_threads);

    auto minimizer_engine = ram::createMinimizerEngine(15, 5, thread_pool);
    minimizer_engine->minimize(dst.begin(), dst.end());
    minimizer_engine->filter(0.001);

    std::vector<std::future<ram::Overlap>> thread_futures;
    for (std::uint32_t i = 0; i < src.size(); ++i) {
        thread_futures.emplace_back(thread_pool->submit(
            [&] (std::uint32_t i) -> ram::Overlap {
                auto overlaps = minimizer_engine->map(src[i], false, false);
                if (overlaps.empty()) {
                    return ram::Overlap(-1, -1, -1, -1, -1, -1, -1, -1);
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
        if (overlap.t_id != -1U) {
            sources.emplace_back(src[i]->id);
            overlaps.emplace_back(overlap);
        }
    }

    // remove contained reads
    std::vector<uint32_t> rank(overlaps.size());
    std::iota(rank.begin(), rank.end(), 0);
    std::sort(rank.begin(), rank.end(),
        [&overlaps] (std::uint32_t lhsr, std::uint32_t rhsr) -> bool {
            const auto& lhs = overlaps[lhsr];
            const auto& rhs = overlaps[rhsr];
            return lhs.t_id < rhs.t_id ||
                (lhs.t_id == rhs.t_id && lhs.t_begin < rhs.t_begin) ||
                (lhs.t_id == rhs.t_id && lhs.t_begin == rhs.t_begin && lhs.t_end > rhs.t_end);
        });

    for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
        for (std::uint32_t j = i + 1; j < overlaps.size(); ++j) {
            if (overlaps[rank[i]].t_id != overlaps[rank[j]].t_id ||
                overlaps[rank[i]].t_end < overlaps[rank[j]].t_end) {
                i = j - 1;
                break;
            }
            sources[rank[j]] = -1;
        }
    }

    // create assembly graph
    struct Node {
        std::uint64_t id;
        std::vector<std::uint64_t> edges;
    };

    std::uint32_t num_nodes = 0;
    for (const auto& it: sources) {
        if (it != -1ULL) {
            ++num_nodes;
        }
    }
    std::vector<Node> nodes(num_nodes);
    std::vector<std::uint64_t> rank_to_node(src.size());

    std::ofstream graph_s("ratlesnake.gfa");

    for (std::uint32_t i = 0, k = 0; i < overlaps.size(); ++i) {
        if (sources[rank[i]] == -1ULL) {
            continue;
        }

        graph_s << "S\t" << src[sources[rank[i]]]->name << "\t"
                << "*" << "\t"
                << "LN:i:" << src[sources[rank[i]]]->data.size() << "\t"
                << std::endl;

        nodes[k].id = rank[i];

        for (std::uint32_t j = i + 1; j < overlaps.size(); ++j) {
            if (sources[rank[j]] == -1ULL) {
                continue;
            }
            if (overlaps[rank[i]].t_id != overlaps[rank[j]].t_id ||
                overlaps[rank[i]].t_end < overlaps[rank[j]].t_begin) {
                break;
            }

            nodes[k].edges.emplace_back(rank[j]);

            graph_s << "L\t" << src[sources[rank[i]]]->name << "\t"
                    << (overlaps[rank[i]].strand ? "+\t" : "-\t")
                    << src[sources[rank[j]]]->name << "\t"
                    << (overlaps[rank[j]].strand ? "+\t" : "-\t")
                    << "*"
                    << std::endl;
        }
        rank_to_node[rank[i]] = k++;
    }

    graph_s.close();

    std::cerr << "[ratlesnake::] chromosome reconstruction ratios"
              << std::endl;
    for (std::uint32_t i = 0; i < nodes.size(); ++i) {
        std::uint32_t t_id = overlaps[nodes[i].id].t_id;
        if (i == 0 || t_id != overlaps[nodes[i - 1].id].t_id) {
            if (i != 0) {
                std::cerr << std::endl;
            }
            std::cerr << "[ratlesnake::] " << dst[t_id - src.size()]->name
                      << std::setprecision(3);
        }
        std::uint32_t t_begin = overlaps[nodes[i].id].t_begin;
        while (!nodes[i].edges.empty()) {
            i = rank_to_node[nodes[i].edges.front()];
        }
        std::uint32_t t_end = overlaps[nodes[i].id].t_end;
        std::cerr << " -> "
                  << static_cast<double>(t_end - t_begin) / dst[t_id - src.size()]->data.size();
    }
    std::cerr << std::endl;

    std::ofstream regular_s("ratlesnake_regular.fasta");
    std::ofstream chimeric_s("ratlesnake_chimeric.fasta");
    std::ofstream repetitive_s("ratlesnake_repeats.fasta");
    std::ofstream junk_s("ratlesnake_junk.fasta");

    std::uint32_t num_regular = 0, num_chimeric = 0, num_repetitive = 0, num_junk = 0;

    for (const auto& it: src) {
        if (annotations[it->id].is_junk) {
            num_junk++;
            junk_s << ">" << it->name
                   << " LN:i:" << it->data.size()
                   << " " << it->id;
            junk_s << std::endl
                   << it->data
                   << std::endl;
            continue;
        }
        if (!annotations[it->id].chimeric_regions.empty()) {
            ++num_chimeric;
            chimeric_s << ">" << it->name
                       << " LN:i:" << it->data.size();
            for (const auto& jt: annotations[it->id].chimeric_regions) {
                chimeric_s << " YB:i:" << (jt >> 32)
                           << " YE:i:" << (jt << 32 >> 32);
            }
            chimeric_s << " " << it->id;
            chimeric_s << std::endl
                       << it->data
                       << std::endl;
            continue;
        }
        if (!annotations[it->id].repetitive_regions.empty()) {
            ++num_repetitive;
            repetitive_s << ">" << it->name
                       << " LN:i:" << it->data.size();
            for (const auto& jt: annotations[it->id].repetitive_regions) {
                repetitive_s << " ZB:i:" << (jt >> 32)
                           << " ZE:i:" << (jt << 32 >> 32);
            }
            repetitive_s << " " << it->id;
            repetitive_s << std::endl
                       << it->data
                       << std::endl;
            continue;
        }
        ++num_regular;
        regular_s << ">" << it->name
                  << " LN:i:" << it->data.size()
                  << " " << it->id;
        regular_s << std::endl
                  << it->data
                  << std::endl;
    }

    repetitive_s.close();
    chimeric_s.close();
    regular_s.close();
    junk_s.close();

    std::cerr << "[ratlesnake::] sequence information" << std::endl
              << "[ratlesnake::] # -> " << src.size() << std::endl
              << "[ratlesnake::] # regular -> " << num_regular << std::endl
              << "[ratlesnake::] # chimeric -> " << num_chimeric << std::endl
              << "[ratlesnake::] # repetitive -> " << num_repetitive << std::endl
              << "[ratlesnake::] # junk -> " << num_junk << std::endl;
}
