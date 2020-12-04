// Copyright (c) 2020 Sara Bakic, Luka Pozega, Robert Vaser

#include <getopt.h>

#include <cstdint>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <numeric>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <utility>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "ram/minimizer_engine.hpp"

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace {

const char* ratlesnake_version = RATLESNAKE_VERSION;

static struct option options[] = {
  {"threads", required_argument, nullptr, 't'},
  {"version", no_argument, nullptr, 'v'},
  {"help", no_argument, nullptr, 'h'},
  {nullptr, 0, nullptr, 0}
};

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(
    const std::string& path) {
  auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta")    || is_suffix(path, ".fa") ||
      is_suffix(path, ".fasta.gz") || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq")    || is_suffix(path, ".fq") ||
      is_suffix(path, ".fastq.gz") || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[ratlesnake::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fa, .fa.gz, .fastq, .fastq.gz, .fq, .fq.gz)!"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout <<
      "usage: ratlesnake [options ...] <sequences> <reference>\n"
      "\n"
      "  <sequences>\n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "  <reference>\n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "\n"
      "  options:\n"
      "    -t, --threads <int>\n"
      "      default: 1\n"
      "      number of threads\n"
      "    --version\n"
      "      prints the version number\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

namespace ratlesnake {

struct Annotation {
  std::vector<std::pair<std::uint32_t, std::uint32_t>> inclusion_intervals;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> chimeric_regions;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> repetitive_regions;
  bool is_junk = false;
};

std::vector<Annotation> Annotate(
    const std::vector<std::unique_ptr<biosoup::Sequence>>& src,
    const std::vector<std::unique_ptr<biosoup::Sequence>>& dst,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool) {
  std::vector<Annotation> annotations(src.size() + dst.size());

  ram::MinimizerEngine minimizer_engine{15, 5, 500, 4, 100, 10000, thread_pool};
  minimizer_engine.Minimize(dst.begin(), dst.end());
  minimizer_engine.Filter(0.001);

  std::vector<std::future<std::vector<biosoup::Overlap>>> thread_futures;
  for (std::uint32_t i = 0; i < src.size(); ++i) {
    thread_futures.emplace_back(thread_pool->Submit(
        [&] (std::uint32_t i) -> std::vector<biosoup::Overlap> {
          auto overlaps = minimizer_engine.Map(src[i], false, false);

          if (overlaps.size() < 2) {
            if (overlaps.size() == 1) {
              std::size_t overlap_len =
                  overlaps.front().lhs_end - overlaps.front().lhs_begin;
              double overlap_score =
                overlaps.front().score / static_cast<double>(overlap_len);
              if (overlap_len < 0.4 * src[i]->data.size() ||
                  overlap_score < 0.1) {
                annotations[src[i]->id].is_junk = true;
              }
            }
            return std::vector<biosoup::Overlap>{};
          }

          std::sort(overlaps.begin(), overlaps.end(),
              [] (const biosoup::Overlap& lhs,
                  const biosoup::Overlap& rhs) -> bool {
                return lhs.lhs_begin < rhs.lhs_begin ||
                      (lhs.lhs_begin == rhs.lhs_begin && lhs.lhs_end > rhs.lhs_end);  // NOLINT
              });

          std::vector<biosoup::Overlap> repetitive_overlaps;

          // annotate chimeric regions
          for (std::uint32_t j = 1, k = 0; j < overlaps.size(); ++j) {
            if (overlaps[k].lhs_end < overlaps[j].lhs_end) {
              bool is_chimeric = true;
              if (overlaps[k].rhs_id == overlaps[j].rhs_id) {
                std::uint32_t lhs_gap = abs(static_cast<std::int32_t>(overlaps[j].lhs_begin) - overlaps[k].lhs_end);  // NOLINT
                std::uint32_t rhs_gap = abs(static_cast<std::int32_t>(overlaps[j].rhs_begin) - overlaps[k].rhs_end);  // NOLINT
                is_chimeric &= std::min(lhs_gap, rhs_gap) < 0.88 * std::max(lhs_gap, rhs_gap);  // NOLINT
              }
              if (is_chimeric) {
                if (overlaps[j].lhs_begin > overlaps[k].lhs_end) {
                  annotations[src[i]->id].chimeric_regions.emplace_back(
                      overlaps[k].lhs_end, overlaps[j].lhs_begin);
                } else {
                  annotations[src[i]->id].chimeric_regions.emplace_back(
                      overlaps[j].lhs_begin, overlaps[k].lhs_end);
                }
              }
              k = j;
            } else {
              repetitive_overlaps.emplace_back(overlaps[j]);
            }
          }

          if (repetitive_overlaps.empty() ||
              !annotations[src[i]->id].chimeric_regions.empty()) {
            return std::vector<biosoup::Overlap>();
          }

          // annotate repetitive regions
          auto is_valid_repeat = [] (
              std::uint32_t lhs_begin,
              std::uint32_t lhs_end,
              std::uint32_t lhs_length) -> bool {
            return lhs_end - lhs_begin > 500 &&
                (lhs_begin < 0.05 * lhs_length || lhs_end > 0.95 * lhs_length || lhs_end - lhs_begin > 2000);  // NOLINT
          };

          std::uint32_t lhs_length = src[i]->data.size();
          std::uint32_t lhs_begin;
          std::uint32_t lhs_end = repetitive_overlaps.front().lhs_end;
          // stop dummy
          repetitive_overlaps.emplace_back(-1, -1, -1, -1, -1, -1, -1, 0);
          for (std::uint32_t j = 1, k = 0; j < repetitive_overlaps.size(); ++j) {  // NOLINT
            if (lhs_end < repetitive_overlaps[j].lhs_begin) {
              lhs_begin = repetitive_overlaps[k].lhs_begin;
              if (is_valid_repeat(lhs_begin, lhs_end, lhs_length)) {
                annotations[src[i]->id].repetitive_regions.emplace_back(
                    lhs_begin, lhs_end);
              }
              k = j;
            }
            lhs_end = std::max(lhs_end, repetitive_overlaps[j].lhs_end);
          }

          return repetitive_overlaps;
        },
        i));
    }

    std::vector<biosoup::Overlap> repetitive_overlaps;
    for (auto& it : thread_futures) {
      auto overlaps = it.get();
      repetitive_overlaps.insert(
          repetitive_overlaps.end(),
          overlaps.begin(),
          overlaps.end());
    }

    if (repetitive_overlaps.empty()) {
      return annotations;
    }

    // annotate repetitive regions on references
    std::sort(repetitive_overlaps.begin(), repetitive_overlaps.end(),
        [] (const biosoup::Overlap& lhs, const biosoup::Overlap& rhs) -> bool {
          return lhs.rhs_id < rhs.rhs_id ||
              (lhs.rhs_id == rhs.rhs_id && lhs.rhs_begin < rhs.rhs_begin) ||
              (lhs.rhs_id == rhs.rhs_id && lhs.rhs_begin == rhs.rhs_begin && lhs.rhs_end > rhs.rhs_end);  // NOLINT
        });

    std::uint32_t rhs_end = repetitive_overlaps.front().rhs_end;
    for (std::uint32_t i = 1, j = 0; i < repetitive_overlaps.size(); ++i) {
      if (repetitive_overlaps[i - 1].rhs_id != repetitive_overlaps[i].rhs_id ||
          rhs_end < repetitive_overlaps[i].rhs_begin) {
        annotations[repetitive_overlaps[i - 1].rhs_id].repetitive_regions.emplace_back(  // NOLINT
            repetitive_overlaps[j].rhs_begin, rhs_end);
        j = i;
        rhs_end = repetitive_overlaps[i].rhs_end;
      } else {
        rhs_end = std::max(rhs_end, repetitive_overlaps[i].rhs_end);
      }
    }

    return annotations;
}

void Reconstruct(
    std::vector<Annotation>* annotations,
    const std::vector<std::unique_ptr<biosoup::Sequence>>& src,
    const std::vector<std::unique_ptr<biosoup::Sequence>>& dst,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool) {
  ram::MinimizerEngine minimizer_engine{15, 5, 500, 4, 100, 10000, thread_pool};
  minimizer_engine.Minimize(dst.begin(), dst.end());
  minimizer_engine.Filter(0.001);

  std::vector<std::future<biosoup::Overlap>> thread_futures;
  for (std::uint32_t i = 0; i < src.size(); ++i) {
    thread_futures.emplace_back(thread_pool->Submit(
        [&] (std::uint32_t i) -> biosoup::Overlap {
          auto overlaps = minimizer_engine.Map(src[i], false, false);
          if (overlaps.empty()) {
            return biosoup::Overlap(-1, -1, -1, -1, -1, -1, -1, 0);
          }
          std::sort(overlaps.begin(), overlaps.end(),
              [] (const biosoup::Overlap& lhs,
                  const biosoup::Overlap& rhs) -> bool {
                return lhs.lhs_end - lhs.lhs_begin > rhs.lhs_end - rhs.lhs_begin;  // NOLINT
              });
            return overlaps.front();
          },
          i));
  }

  std::vector<std::uint64_t> sources;
  std::vector<biosoup::Overlap> overlaps;
  for (std::uint32_t i = 0; i < thread_futures.size(); ++i) {
    auto overlap = thread_futures[i].get();
    if (overlap.rhs_id != -1U) {
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
        return lhs.rhs_id < rhs.rhs_id ||
            (lhs.rhs_id == rhs.rhs_id && lhs.rhs_begin < rhs.rhs_begin) ||
            (lhs.rhs_id == rhs.rhs_id && lhs.rhs_begin == rhs.rhs_begin && lhs.rhs_end > rhs.rhs_end);  // NOLINT
      });

  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    for (std::uint32_t j = i + 1; j < overlaps.size(); ++j) {
      if (overlaps[rank[i]].rhs_id != overlaps[rank[j]].rhs_id ||
          overlaps[rank[i]].rhs_end < overlaps[rank[j]].rhs_end) {
        i = j - 1;
        break;
      }
      (*annotations)[sources[rank[j]]].inclusion_intervals.emplace_back(
          0, sources[rank[i]]);
    }
  }

  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    if ((*annotations)[sources[rank[i]]].inclusion_intervals.empty() == false) {
      (*annotations)[sources[rank[i]]].inclusion_intervals.emplace_back(
          overlaps[rank[i]].lhs_begin, overlaps[rank[i]].lhs_end);
      sources[rank[i]] = -1;
    }
  }

  // create assembly graph
  struct Node {
    std::uint64_t id;
    std::vector<std::uint64_t> edges;
  };

  std::uint32_t num_nodes = 0;
  for (const auto& it : sources) {
    if (it != -1ULL) {
      ++num_nodes;
    }
  }
  std::vector<Node> nodes(num_nodes);
  std::vector<std::uint64_t> rank_to_node(src.size());

  std::ofstream graph_s("ratlesnake.gfa");
  std::ofstream solid_s("ratlesnake_solid.fasta");

  for (std::uint32_t i = 0, k = 0; i < overlaps.size(); ++i) {
    if (sources[rank[i]] == -1ULL) {
      continue;
    }

    graph_s << "S\t" << src[sources[rank[i]]]->name << "\t"
            << "*" << "\t"
            << "LN:i:" << src[sources[rank[i]]]->data.size() << "\t"
            << "UR:Z:ratlesnake_solid.fasta"
            << std::endl;

    solid_s << ">" << src[sources[rank[i]]]->name
            << " LN:i:" << src[sources[rank[i]]]->data.size()
            << " XB:i:" << overlaps[rank[i]].lhs_begin
            << " XE:i:" << overlaps[rank[i]].lhs_end;
    for (const auto& it : (*annotations)[sources[rank[i]]].chimeric_regions) {
      solid_s << " YB:i:" << it.first
              << " YE:i:" << it.second;
    }
    for (const auto& it : (*annotations)[sources[rank[i]]].repetitive_regions) {
      solid_s << " ZB:i:" << it.first
              << " ZE:i:" << it.second;
    }

    solid_s << " " << sources[rank[i]];
    solid_s << std::endl
            << src[sources[rank[i]]]->data
            << std::endl;

    nodes[k].id = rank[i];

    for (std::uint32_t j = i + 1; j < overlaps.size(); ++j) {
      if (sources[rank[j]] == -1ULL) {
        continue;
      }
      if (overlaps[rank[i]].rhs_id != overlaps[rank[j]].rhs_id ||
          overlaps[rank[i]].rhs_end < overlaps[rank[j]].rhs_begin) {
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
      std::uint32_t rhs_id = overlaps[nodes[i].id].rhs_id;
      if (i == 0 || rhs_id != overlaps[nodes[i - 1].id].rhs_id) {
        if (i != 0) {
          std::cerr << std::endl;
        }
        std::cerr << "[ratlesnake::] " << dst[rhs_id - src.size()]->name
                  << std::setprecision(3);
      }
      std::uint32_t rhs_begin = overlaps[nodes[i].id].rhs_begin;
      while (!nodes[i].edges.empty()) {
        i = rank_to_node[nodes[i].edges.front()];
      }
      std::uint32_t rhs_end = overlaps[nodes[i].id].rhs_end;
      std::cerr << " -> " << static_cast<double>(rhs_end - rhs_begin) / dst[rhs_id - src.size()]->data.size();  // NOLINT
    }
    std::cerr << std::endl;

    std::ofstream contained_s("ratlesnake_contained.fasta");
    std::ofstream chimeric_s("ratlesnake_chimeric.fasta");
    std::ofstream repetitive_s("ratlesnake_repetitive.fasta");
    std::ofstream junk_s("ratlesnake_junk.fasta");

    std::uint32_t num_contained = 0;
    std::uint32_t num_chimeric = 0;
    std::uint32_t num_repetitive = 0;
    std::uint32_t num_junk = 0;

    for (const auto& it : src) {
      if ((*annotations)[it->id].is_junk) {
        num_junk++;
        junk_s << ">" << it->name
               << " LN:i:" << it->data.size()
               << " " << it->id;
        junk_s << std::endl
               << it->data
               << std::endl;
        continue;
      }
      if (!(*annotations)[it->id].inclusion_intervals.empty()) {
        ++num_contained;
        const auto& intervals = (*annotations)[it->id].inclusion_intervals;
        contained_s << ">" << it->name
                    << " LN:i:" << it->data.size()
                    << " XB:i:" << intervals.back().first
                    << " XE:i:" << intervals.back().second;
        for (std::uint32_t i = 0; i < intervals.size() - 1; ++i) {
          contained_s << " XI:i:" << intervals[i].second;
        }
        contained_s << " " << it->id;
        contained_s << std::endl
                    << it->data
                    << std::endl;
      }
      if (!(*annotations)[it->id].chimeric_regions.empty()) {
        ++num_chimeric;
        chimeric_s << ">" << it->name
                   << " LN:i:" << it->data.size();
        for (const auto& jt : (*annotations)[it->id].chimeric_regions) {
          chimeric_s << " YB:i:" << jt.first
                     << " YE:i:" << jt.second;
        }
        chimeric_s << " " << it->id;
        chimeric_s << std::endl
                   << it->data
                   << std::endl;
      }
      if (!(*annotations)[it->id].repetitive_regions.empty()) {
        ++num_repetitive;
        repetitive_s << ">" << it->name
                     << " LN:i:" << it->data.size();
        for (const auto& jt : (*annotations)[it->id].repetitive_regions) {
          repetitive_s << " ZB:i:" << jt.first
                       << " ZE:i:" << jt.second;
        }
        repetitive_s << " " << it->id;
        repetitive_s << std::endl
                    << it->data
                    << std::endl;
      }
    }

    repetitive_s.close();
    chimeric_s.close();
    contained_s.close();
    junk_s.close();

    std::cerr << "[ratlesnake::] sequence information" << std::endl
              << "[ratlesnake::] # -> " << src.size() << std::endl
              << "[ratlesnake::] # contained -> " << num_contained << std::endl
              << "[ratlesnake::] # chimeric -> " << num_chimeric << std::endl
              << "[ratlesnake::] # repetitive -> " << num_repetitive << std::endl  // NOLINT
              << "[ratlesnake::] # junk -> " << num_junk << std::endl;
  }

}  // namespace ratlesnake

int main(int argc, char** argv) {
  std::uint32_t num_threads = 1;

  std::vector<std::string> input_paths;

  std::string optstr = "t:h";
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 't': num_threads = atoi(optarg); break;
      case 'v': std::cout << ratlesnake_version << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }

  for (std::int32_t i = optind; i < argc; ++i) {
    input_paths.emplace_back(argv[i]);
  }
  if (input_paths.size() < 2) {
    std::cerr << "[ratlesnake::] error: missing input file(s)!" << std::endl;
    Help();
    return 1;
  }

  auto sparser = CreateParser(input_paths[0]);
  if (sparser == nullptr) {
    return 1;
  }

  auto rparser = CreateParser(input_paths[1]);
  if (rparser == nullptr) {
    return 1;
  }

  auto sequences = sparser->Parse(-1);
  auto references = rparser->Parse(-1);

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

  auto annotations = ratlesnake::Annotate(sequences, references, thread_pool);

  ratlesnake::Reconstruct(&annotations, sequences, references, thread_pool);

  return 0;
}
