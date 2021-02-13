// Copyright (c) 2020 Sara Bakic, Luka Pozega, Robert Vaser

#include <getopt.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <map>
#include <utility>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/sequence.hpp"
#include "edlib.h"  // NOLINT
#include "ram/minimizer_engine.hpp"

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};
std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

namespace {

static struct option options[] = {
  {"threads", required_argument, nullptr, 't'},
  {"result", required_argument, nullptr, 'r'},
  {"version", no_argument, nullptr, 'v'},
  {"help", no_argument, nullptr, 'h'},
  {nullptr, 0, nullptr, 0}
};

template<typename T>
std::unique_ptr<bioparser::Parser<T>> CreateParser(const std::string& path) {
  auto is_suffix = [] (const std::string& s, const std::string& suff) {
    return s.size() < suff.size() ? false :
        s.compare(s.size() - suff.size(), suff.size(), suff) == 0;
  };

  if (is_suffix(path, ".fasta") || is_suffix(path, ".fasta.gz") ||
      is_suffix(path, ".fna")   || is_suffix(path, ".fna.gz") ||
      is_suffix(path, ".fa")    || is_suffix(path, ".fa.gz")) {
    try {
      return bioparser::Parser<T>::template Create<bioparser::FastaParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }
  if (is_suffix(path, ".fastq") || is_suffix(path, ".fastq.gz") ||
      is_suffix(path, ".fq")    || is_suffix(path, ".fq.gz")) {
    try {
      return bioparser::Parser<T>::template Create<bioparser::FastqParser>(path);  // NOLINT
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[ratlesnake::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fna, .fna.gz, .fa, .fa.gz, .fastq, .fastq.gz, "
            << ".fq, .fq.gz)!"
            << std::endl;
  return nullptr;
}

void Help() {
  std::cout <<
      "usage: ratlesnake [options ...] <sequences> [<reference>]\n"
      "\n"
      "  <sequences>/<reference>\n"
      "    input file in FASTA/FASTQ format (can be compressed with gzip)\n"
      "\n"
      "  options:\n"
      "    -t, --threads <int>\n"
      "      default: 1\n"
      "      number of threads\n"
      "    -r, --result <int> (option can be used multiple times)\n"
      "      default: 0\n"
      "      result mode:\n"
      "        0 - sequence statistics (histogram(s) on stdout)\n"
      "        1 - sequence accuracy given reference (histogram on stdout)\n"
      "        2 - sequence annotation given reference (FASTA files)\n"
      "            - XI:i:<int> 0-based id\n"
      "            - XC:i:<int> 0-based id of containing read\n"
      "            - XB:i:<int> (XE:i:<int>) valid region begin (end)\n"
      "            - YB:i:<int> (YE:i:<int>) chimeric region begin (end)\n"
      "            - ZB:i:<int> (ZE:i:<int>) repetitive region begin (end)\n"
      "        3 - reference reconstruction given sequences (GFA files)\n"
      "    --version\n"
      "      prints the version number\n"
      "    -h, --help\n"
      "      prints the usage\n";
}

}  // namespace

namespace ratlesnake {

std::uint32_t OverlapLength(const biosoup::Overlap& o) {
  return std::max(o.rhs_end - o.rhs_begin, o.lhs_end - o.lhs_begin);
}

void Analyse(const std::string& sequences_path) {
  auto sparser = CreateParser<biosoup::Sequence>(sequences_path);

  std::vector<std::size_t> lens, quals;
  while (true) {
    auto sequences = sparser->Parse(1ULL << 30);
    if (sequences.empty()) {
      break;
    }

    lens.reserve(lens.size() + sequences.size());
    for (const auto& it : sequences) {
      lens.emplace_back(it->data.size());
      if (!it->quality.empty()) {
        quals.emplace_back(
            std::accumulate(it->quality.begin(), it->quality.end(), 0ULL,
                [] (std::size_t s, char c) -> std::size_t {
                  return s + (c - '!');
                })
            / it->data.size());
      }
    }
  }
  if (lens.empty()) {
    return;
  }

  // length
  std::sort(lens.begin(), lens.end());
  std::size_t total = std::accumulate(lens.begin(), lens.end(), 0ULL);
  std::size_t td = log10(std::max(total, lens.size())) + 1;

  std::cout << "Num seq    = " << std::setw(td) << lens.size()
            << std::endl
            << "Total len  = " << std::setw(td) << total
            << std::endl
            << "Min len    = " << std::setw(td) << lens.front()
            << std::endl
            << "Median len = " << std::setw(td) << lens[lens.size() / 2]
            << std::endl
            << "Avg len    = " << std::setw(td) << total / lens.size()
            << std::endl
            << "Max len    = " << std::setw(td) << lens.back()
            << std::endl;

  std::map<std::size_t, std::size_t> hist = {
    {1000, 0}, {2000, 0}, {3000, 0}, {4000, 0}, {5000, 0}, {7500, 0},
    {10000, 0}, {12500, 0}, {15000, 0}, {17500, 0}, {20000, 0}, {25000, 0},
    {30000, 0}, {35000, 0}, {40000, 0}, {45000, 0}, {50000, 0}, {75000, 0},
    {100000, 0}, {125000, 0}, {150000, 0}, {175000, 0}, {200000, 0}
  };
  std::size_t whales = 0;
  for (const auto& it : lens) {
    if (it > 200000) {
      ++whales;
      continue;
    }
    for (auto& jt : hist) {
      if (it < jt.first) {
        ++jt.second;
        break;
      }
    }
  }
  std::size_t hist_max = 0;
  for (const auto& it : hist) {
    hist_max = std::max(hist_max, it.second);
  }
  hist_max = std::max(hist_max, whales);
  std::size_t nd = std::max(hist_max > 0 ? log10(hist_max) + 1 : 1, 3.);
  std::size_t ld = std::max(hist_max > 0 ? log(hist_max) : 0, 3.);

  std::cout << std::string(nd + ld + 14, '-') << std::endl;
  std::cout << "    Len" << " | "
            << std::string(nd - 3 , ' ') << "Num" << " | "
            << "Log" << std::endl;
  std::cout << std::string(nd + ld + 14, '-') << std::endl;
  for (const auto& it : hist) {
    std::cout << std::setw(7) << it.first
              << " | "
              << std::setw(nd) << it.second
              << " | "
              << std::string(it.second > 0 ? log(it.second) : 0, '/')
              << std::endl;
  }
  std::cout << "  _v_  "
            << " | "
            << std::setw(nd) << whales
            << " | "
            << std::string(whales > 0 ? log(whales) : 0, '/')
            << std::endl;
  std::cout << " (___\\/" << " | " << std::string(nd, ' ')<< " |" << std::endl;
  std::cout << std::string(nd + ld + 14, '-') << std::endl;

  if (quals.empty()) {
    return;
  }

  // Phred base quality
  std::sort(quals.begin(), quals.end());
  total = std::accumulate(quals.begin(), quals.end(), 0ULL);

  std::cout << "Min qual    = " << std::setw(2) << quals.front()
            << std::endl
            << "Median qual = " << std::setw(2) << quals[quals.size() / 2]
            << std::endl
            << "Avg qual    = " << std::setw(2) << total / quals.size()
            << std::endl
            << "Max qual    = " << std::setw(2) << quals.back()
            << std::endl;

  hist.clear();
  for (std::uint32_t i = quals.front(); i < quals.back(); ++i) {
    hist[i] = 0;
  }
  for (const auto& it : quals) {
    ++hist[it];
  }
  hist_max = 0;
  for (const auto& it : hist) {
    hist_max = std::max(hist_max, it.second);
  }
  nd = std::max(hist_max > 0 ? log10(hist_max) + 1 : 1, 3.);
  ld = std::max(hist_max > 0 ? log(hist_max) : 0, 3.);

  std::cout << std::string(nd + ld + 13, '-') << std::endl;
  std::cout << " Phred" << " | "
            << std::string(nd - 3 , ' ') << "Num" << " | "
            << "Log" << std::endl;
  std::cout << std::string(nd + ld + 13, '-') << std::endl;
  for (const auto& it : hist) {
    std::cout << std::setw(6) << it.first
              << " | "
              << std::setw(nd) << it.second
              << " | "
              << std::string(it.second > 0 ? log(it.second) : 0, '/')
              << std::endl;
  }
  std::cout << std::string(nd + ld + 13, '-') << std::endl;
}

void Analyse(
    const std::string& sequences_path,
    const std::string& reference_path,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool) {
  auto sparser = CreateParser<biosoup::NucleicAcid>(sequences_path);
  auto rparser = CreateParser<biosoup::NucleicAcid>(reference_path);

  ram::MinimizerEngine minimizer_engine{thread_pool};
  std::vector<biosoup::Overlap> overlaps;

  biosoup::NucleicAcid::num_objects = 0;
  while (true) {
    auto reference = rparser->Parse(1ULL << 30);
    if (reference.empty()) {
      break;
    }
    std::size_t num_references = reference.size();

    minimizer_engine.Minimize(reference.begin(), reference.end());
    minimizer_engine.Filter(0.001);

    biosoup::NucleicAcid::num_objects = 0;
    sparser->Reset();
    while (true) {
      auto sequences = sparser->Parse(1ULL << 30);
      if (sequences.empty()) {
        break;
      }
      if (sequences.back()->id >= overlaps.size()) {
        overlaps.resize(sequences.back()->id + 1);
      }

      std::vector<std::future<void>> futures;
      for (const auto& it : sequences) {
        futures.emplace_back(thread_pool->Submit(
            [&] (decltype(it) sequence) -> void {
              auto ovl = minimizer_engine.Map(sequence, false, false);
              if (ovl.empty()) {
                return;
              }
              std::sort(ovl.begin(), ovl.end(),
                  [&] (const biosoup::Overlap& lhs,
                       const biosoup::Overlap& rhs) -> bool {
                    return OverlapLength(lhs) > OverlapLength(rhs);
                  });
              if (OverlapLength(overlaps[sequence->id]) <
                  OverlapLength(ovl.front())) {
                overlaps[sequence->id] = ovl.front();
              }
            },
            std::ref(it)));
      }
      for (const auto& it : futures) {
        it.wait();
      }
    }

    biosoup::NucleicAcid::num_objects = num_references;
  }

  // accuracy
  std::vector<double> accs;

  biosoup::NucleicAcid::num_objects = 0;
  rparser->Reset();
  auto reference = rparser->Parse(-1);

  biosoup::NucleicAcid::num_objects = 0;
  sparser->Reset();
  while (true) {
    auto sequences = sparser->Parse(1ULL << 30);
    if (sequences.empty()) {
      break;
    }

    std::vector<std::future<double>> futures;
    for (const auto& it : sequences) {
      futures.emplace_back(thread_pool->Submit(
          [&] (decltype(it) sequence) -> double {
            auto& o = overlaps[sequence->id];
            if (OverlapLength(o) == 0) {
              return 0;
            }
            if (o.strand == 0) {
              sequence->ReverseAndComplement();
              auto lhs_begin = o.lhs_begin;
              o.lhs_begin = it->inflated_len - o.lhs_end;
              o.lhs_end = it->inflated_len - lhs_begin;
            }

            std::string s = sequence->InflateData(
                o.lhs_begin,
                o.lhs_end - o.lhs_begin);
            std::string r = reference[o.rhs_id]->InflateData(
                o.rhs_begin,
                o.rhs_end - o.rhs_begin);

            auto result = edlibAlign(
                s.c_str(), s.size(),
                r.c_str(), r.size(),
                edlibDefaultAlignConfig());
            double score = 0;
            if (result.status == EDLIB_STATUS_OK) {
              score = 1 - result.editDistance /
                  static_cast<double>(OverlapLength(o));
            }
            edlibFreeAlignResult(result);

            return score;
          },
          std::ref(it)));
    }
    for (auto& it : futures) {
      accs.emplace_back(it.get());
    }
  }
  std::sort(accs.begin(), accs.end());
  double total = std::accumulate(accs.begin(), accs.end(), 0.);
  std::size_t td = std::max(log10(accs.size()) + 1, 6.);

  std::cout << "Num seq    = " << std::setw(td) << accs.size()
            << std::endl
            << "Min acc    = " << std::fixed << std::setprecision(4)
            << std::setw(td) << accs.front()
            << std::endl
            << "Median acc = " << std::fixed << std::setprecision(4)
            << std::setw(td) << accs[accs.size() / 2]
            << std::endl
            << "Avg acc    = " << std::fixed << std::setprecision(4)
            << std::setw(td) << total / accs.size()
            << std::endl
            << "Max acc    = " << std::fixed << std::setprecision(4)
            << std::setw(td) << accs.back()
            << std::endl;

  std::map<double, std::size_t> hist;
  for (double i = 0.00000001; i < 1.001; i += 0.0125) {
    hist[i] = 0;
  }
  std::size_t num_junk = 0;
  for (const auto& it : accs) {
    if (it <= 0.0001) {
      ++num_junk;
      continue;
    }
    for (auto& jt : hist) {
      if (it <= jt.first) {
        ++jt.second;
        break;
      }
    }
  }
  std::size_t hist_max = num_junk;
  for (auto it = hist.cbegin(); it != hist.cend();) {
    if (it->second == 0) {
      it = hist.erase(it);
    } else {
      hist_max = std::max(hist_max, it->second);
      ++it;
    }
  }
  std::size_t nd = std::max(hist_max > 0 ? log10(hist_max) + 1 : 1, 3.);
  std::size_t ld = std::max(hist_max > 0 ? log(hist_max) : 0, 3.);

  std::cout << std::string(nd + ld + 14, '-') << std::endl;
  std::cout << "    Acc" << " | "
            << std::string(nd - 3 , ' ') << "Num" << " | "
            << "Log" << std::endl;
  std::cout << std::string(nd + ld + 14, '-') << std::endl;
  for (const auto& it : hist) {
    std::cout << std::fixed << std::setprecision(4) << std::setw(7) << it.first
              << " | "
              << std::setw(nd) << it.second
              << " | "
              << std::string(it.second > 0 ? log(it.second) : 0, '/')
              << std::endl;
  }
  std::cout << std::string(nd + ld + 14, '-') << std::endl;
  std::cout << "   Junk" << " | "
            << std::setw(nd) << num_junk
            << " | "
            << std::string(num_junk > 0 ? log(num_junk) : 0, '/')
            << std::endl;
  std::cout << std::string(nd + ld + 14, '-') << std::endl;
}

void Annotate(
    const std::string& sequences_path,
    const std::string& reference_path,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool) {
  auto sparser = CreateParser<biosoup::NucleicAcid>(sequences_path);
  auto rparser = CreateParser<biosoup::NucleicAcid>(reference_path);

  biosoup::NucleicAcid::num_objects = 0;
  auto reference = rparser->Parse(-1);

  ram::MinimizerEngine minimizer_engine{thread_pool};
  minimizer_engine.Minimize(reference.begin(), reference.end());
  minimizer_engine.Filter(0.001);

  std::vector<biosoup::Overlap> overlaps;
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> chimeric_regions;  // NOLINT
  std::vector<std::vector<std::pair<std::uint32_t, std::uint32_t>>> repetitive_regions;  // NOLINT

  biosoup::NucleicAcid::num_objects = 0;
  while (true) {
    auto sequences = sparser->Parse(1ULL << 30);
    if (sequences.empty()) {
      break;
    }
    if (sequences.back()->id >= chimeric_regions.size()) {
      chimeric_regions.resize(sequences.back()->id + 1);
      repetitive_regions.resize(chimeric_regions.size());
    }

    std::vector<std::future<biosoup::Overlap>> futures;
    for (const auto& it : sequences) {
      futures.emplace_back(thread_pool->Submit(
          [&] (decltype(it) sequence) -> biosoup::Overlap {
            auto overlaps = minimizer_engine.Map(sequence, false, false);

            if (overlaps.empty()) {
              return biosoup::Overlap{sequence->id, 0, 0, -1U, 0, 0, 0};
            } else if (overlaps.size() == 1) {
              std::size_t lhs_overlap_len =
                  overlaps.front().lhs_end -
                  overlaps.front().lhs_begin;
              double lhs_score =
                  overlaps.front().score /
                  static_cast<double>(lhs_overlap_len);
              if (lhs_overlap_len < 0.4 * sequence->inflated_len || lhs_score < 0.1) {  // NOLINT
                return biosoup::Overlap{sequence->id, 0, 0, -1U, 0, 0, 0};
              }
              return overlaps.front();
            }

            std::sort(overlaps.begin(), overlaps.end(),
                [] (const biosoup::Overlap& lhs,
                    const biosoup::Overlap& rhs) -> bool {
                  return (lhs.lhs_begin <  rhs.lhs_begin) ||
                         (lhs.lhs_begin == rhs.lhs_begin && lhs.lhs_end > rhs.lhs_end);  // NOLINT
                });

            std::vector<biosoup::Overlap> repeats;

            // annotate chimeric regions
            for (std::uint32_t i = 0, j = 1; j < overlaps.size(); ++j) {
              if (overlaps[i].lhs_end >= overlaps[j].lhs_end) {
                repeats.emplace_back(overlaps[j]);
              } else {
                bool is_chimeric = true;
                if (overlaps[i].rhs_id == overlaps[j].rhs_id &&
                    overlaps[i].strand == overlaps[j].strand) {
                  std::size_t ref_len =
                      reference[overlaps[i].rhs_id]->inflated_len;

                  if (( overlaps[i].strand && overlaps[j].rhs_begin < 256 && overlaps[i].rhs_end > ref_len - 256) ||  // NOLINT
                      (!overlaps[i].strand && overlaps[i].rhs_begin < 256 && overlaps[j].rhs_end > ref_len - 256)) {  // NOLINT
                    is_chimeric = false;
                  } else {
                    std::int32_t lhs_gap =
                        overlaps[j].lhs_begin - overlaps[i].lhs_end;
                    std::int32_t rhs_gap = overlaps[i].strand ?
                        overlaps[j].rhs_begin - overlaps[i].rhs_end :
                        overlaps[i].rhs_begin - overlaps[j].rhs_end;
                    is_chimeric = std::abs(lhs_gap - rhs_gap) > 1280;
                  }
                }
                if (is_chimeric) {
                  if (overlaps[i].lhs_end < overlaps[j].lhs_begin) {
                    chimeric_regions[sequence->id].emplace_back(
                        overlaps[i].lhs_end,
                        overlaps[j].lhs_begin);
                  } else {
                    chimeric_regions[sequence->id].emplace_back(
                        overlaps[j].lhs_begin,
                        overlaps[i].lhs_end);
                  }
                }
                i = j;
              }
            }

            // annotate repetitive regions
            repeats.emplace_back(-1, -1, -1, -1, -1, -1, -1, 0);
            for (std::uint32_t i = 0, j = 1; j < repeats.size(); ++j) {
              if (repeats[i].lhs_end < repeats[j].lhs_begin) {
                if ((repeats[i].lhs_end - repeats[i].lhs_begin > 500) &&
                    (repeats[i].lhs_begin < 0.05 * sequence->inflated_len ||
                     repeats[i].lhs_end   > 0.95 * sequence->inflated_len ||
                     repeats[i].lhs_end - repeats[i].lhs_begin > 2000)) {
                  repetitive_regions[sequence->id].emplace_back(
                      repeats[i].lhs_begin,
                      repeats[i].lhs_end);
                }
                i = j;
              }
              repeats[i].lhs_end = std::max(
                  repeats[i].lhs_end,
                  repeats[j].lhs_end);
            }

            // return longest overlap
            std::sort(overlaps.begin(), overlaps.end(),
                [] (const biosoup::Overlap& lhs,
                    const biosoup::Overlap& rhs) -> bool {
                  return OverlapLength(lhs) > OverlapLength(rhs);
                });

            return overlaps.front();
          },
          std::ref(it)));
    }
    for (auto& it : futures) {
      overlaps.emplace_back(it.get());
    }
  }

  // annotate inclusion regions
  std::vector<std::vector<std::uint32_t>> containment(overlaps.size());

  std::sort(overlaps.begin(), overlaps.end(),
      [] (const biosoup::Overlap& lhs, const biosoup::Overlap& rhs) -> bool {
        return
            (lhs.rhs_id <  rhs.rhs_id) ||
            (lhs.rhs_id == rhs.rhs_id && lhs.rhs_begin <  rhs.rhs_begin) ||
            (lhs.rhs_id == rhs.rhs_id && lhs.rhs_begin == rhs.rhs_begin && lhs.rhs_end > rhs.rhs_end);  // NOLINT
      });

  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    if (OverlapLength(overlaps[i]) == 0) {
      continue;
    }
    for (std::uint32_t j = i + 1; j < overlaps.size(); ++j) {
      if (overlaps[i].rhs_id != overlaps[j].rhs_id ||
          overlaps[i].rhs_end < overlaps[j].rhs_end) {
        i = j - 1;
        break;
      }
      if (containment[overlaps[j].lhs_id].empty()) {
        containment[overlaps[j].lhs_id].emplace_back(overlaps[i].lhs_id);
      }
    }
  }

  std::sort(overlaps.begin(), overlaps.end(),
      [] (const biosoup::Overlap& lhs, const biosoup::Overlap& rhs) -> bool {
        return lhs.lhs_id < rhs.lhs_id;
      });

  // store
  std::string prefix = sequences_path;
  {
    auto c = prefix.rfind('/');
    prefix = prefix.substr(c == std::string::npos ? 0 : c + 1);
    c = prefix.find('.');
    prefix = prefix.substr(0, c);
  }

  std::ofstream os_uncontained(prefix + "_uncontained.fasta");
  std::ofstream os_contained(prefix + "_contained.fasta");
  std::ofstream os_junk(prefix + "_junk.fasta");
  std::ofstream os_chimeric(prefix + "_chimeric.fasta");
  std::ofstream os_repetitive(prefix + "_repetitive.fasta");

  std::size_t num_uncontained = 0;
  std::size_t num_contained = 0;
  std::size_t num_junk = 0;
  std::size_t num_chimeric = 0;
  std::size_t num_repetitive = 0;

  biosoup::NucleicAcid::num_objects = 0;
  sparser->Reset();
  while (true) {
    auto sequences = sparser->Parse(1ULL << 30);
    if (sequences.empty()) {
      break;
    }

    for (const auto& it : sequences) {
      it->name += " XI:i:" + std::to_string(it->id);
      if (OverlapLength(overlaps[it->id]) != 0) {
        it->name += " XB:i:" + std::to_string(overlaps[it->id].lhs_begin);
        it->name += " XE:i:" + std::to_string(overlaps[it->id].lhs_end);
      }
      for (const auto& jt : containment[it->id]) {
        it->name += " XC:i:" + std::to_string(jt);
      }
      for (const auto& jt : chimeric_regions[it->id]) {
        it->name += " YB:i:" + std::to_string(jt.first);
        it->name += " YE:i:" + std::to_string(jt.second);
      }
      for (const auto& jt : repetitive_regions[it->id]) {
        it->name += " ZB:i:" + std::to_string(jt.first);
        it->name += " ZE:i:" + std::to_string(jt.second);
      }

      if (OverlapLength(overlaps[it->id]) == 0) {
        ++num_junk;
        os_junk << ">" << it->name << std::endl << it->InflateData() << std::endl;  // NOLINT
        continue;
      }
      if (!containment[it->id].empty()) {
        ++num_contained;
        os_contained << ">" << it->name << std::endl << it->InflateData() << std::endl;  // NOLINT
      } else {
        ++num_uncontained;
        os_uncontained << ">" << it->name << std::endl << it->InflateData() << std::endl;  // NOLINT
      }
      if (!chimeric_regions[it->id].empty()) {
        ++num_chimeric;
        os_chimeric << ">" << it->name << std::endl << it->InflateData() << std::endl;  // NOLINT
      }
      if (!repetitive_regions[it->id].empty()) {
        ++num_repetitive;
        os_repetitive << ">" << it->name << std::endl << it->InflateData() << std::endl;  // NOLINT
      }
    }
  }

  os_repetitive.close();
  os_chimeric.close();
  os_junk.close();
  os_contained.close();
  os_uncontained.close();

  std::size_t td = log10(overlaps.size()) + 1;
  std::cout << "Num seq    = " << std::setw(td) << overlaps.size()
            << std::endl
            << "Num uncont = " << std::setw(td) << num_uncontained
            << std::endl
            << "Num cont   = " << std::setw(td) << num_contained
            << std::endl
            << "Num junk   = " << std::setw(td) << num_junk
            << std::endl
            << "Num chim   = " << std::setw(td) << num_chimeric
            << std::endl
            << "Num repeat = " << std::setw(td) << num_repetitive
            << std::endl;
}

void Reconstruct(
    const std::string& sequences_path,
    const std::string& reference_path,
    std::shared_ptr<thread_pool::ThreadPool> thread_pool) {
  auto sparser = CreateParser<biosoup::NucleicAcid>(sequences_path);
  auto rparser = CreateParser<biosoup::NucleicAcid>(reference_path);

  biosoup::NucleicAcid::num_objects = 0;
  auto reference = rparser->Parse(-1);

  ram::MinimizerEngine minimizer_engine{thread_pool};
  minimizer_engine.Minimize(reference.begin(), reference.end());
  minimizer_engine.Filter(0.001);

  std::vector<biosoup::Overlap> overlaps;

  biosoup::NucleicAcid::num_objects = 0;
  while (true) {
    auto sequences = sparser->Parse(1ULL << 30);
    if (sequences.empty()) {
      break;
    }

    std::vector<std::future<biosoup::Overlap>> futures;
    for (const auto& it : sequences) {
      futures.emplace_back(thread_pool->Submit(
          [&] (decltype(it) sequence) -> biosoup::Overlap {
            auto overlaps = minimizer_engine.Map(sequence, false, false);

            if (overlaps.empty()) {
              return biosoup::Overlap{sequence->id, 0, 0, -1U, 0, 0, 0};
            } else if (overlaps.size() == 1) {
              std::size_t lhs_overlap_len =
                  overlaps.front().lhs_end -
                  overlaps.front().lhs_begin;
              double lhs_score =
                  overlaps.front().score /
                  static_cast<double>(lhs_overlap_len);
              if (lhs_overlap_len < 0.4 * sequence->inflated_len || lhs_score < 0.1) {  // NOLINT
                return biosoup::Overlap{sequence->id, 0, 0, -1U, 0, 0, 0};
              }
              return overlaps.front();
            }

            std::sort(overlaps.begin(), overlaps.end(),
                [] (const biosoup::Overlap& lhs,
                    const biosoup::Overlap& rhs) -> bool {
                  return OverlapLength(lhs) > OverlapLength(rhs);
                });
            return overlaps.front();
          },
          std::ref(it)));
    }
    for (auto& it : futures) {
      overlaps.emplace_back(it.get());
    }
  }

  // remove contained sequences
  std::vector<std::vector<std::uint32_t>> containment(overlaps.size());

  std::sort(overlaps.begin(), overlaps.end(),
      [] (const biosoup::Overlap& lhs, const biosoup::Overlap& rhs) -> bool {
        return
            (lhs.rhs_id <  rhs.rhs_id) ||
            (lhs.rhs_id == rhs.rhs_id && lhs.rhs_begin <  rhs.rhs_begin) ||
            (lhs.rhs_id == rhs.rhs_id && lhs.rhs_begin == rhs.rhs_begin && lhs.rhs_end > rhs.rhs_end);  // NOLINT
      });

  for (std::uint32_t i = 0; i < overlaps.size(); ++i) {
    if (OverlapLength(overlaps[i]) == 0) {
      continue;
    }
    for (std::uint32_t j = i + 1; j < overlaps.size(); ++j) {
      if (overlaps[i].rhs_id != overlaps[j].rhs_id ||
          overlaps[i].rhs_end < overlaps[j].rhs_end) {
        i = j - 1;
        break;
      }
      if (OverlapLength(overlaps[j]) != 0) {
        overlaps[j] = biosoup::Overlap{overlaps[j].lhs_id, 0, 0, -1U, 0, 0, 0};
      }
    }
  }

  std::sort(overlaps.begin(), overlaps.end(),
      [] (const biosoup::Overlap& lhs, const biosoup::Overlap& rhs) -> bool {
        return lhs.lhs_id < rhs.lhs_id;
      });

  biosoup::NucleicAcid::num_objects = 0;
  sparser->Reset();
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> sequences;
  while (true) {
    auto chunk = sparser->Parse(1ULL << 30);
    if (chunk.empty()) {
      break;
    }

    for (auto& it : chunk) {
      if (OverlapLength(overlaps[it->id]) == 0) {
        it.reset();
      }
    }
    sequences.insert(
        sequences.end(),
        std::make_move_iterator(chunk.begin()),
        std::make_move_iterator(chunk.end()));
  }

  // store
  std::sort(overlaps.begin(), overlaps.end(),
      [] (const biosoup::Overlap& lhs, const biosoup::Overlap& rhs) -> bool {
        return
            (lhs.rhs_id <  rhs.rhs_id) ||
            (lhs.rhs_id == rhs.rhs_id && lhs.rhs_begin <  rhs.rhs_begin) ||
            (lhs.rhs_id == rhs.rhs_id && lhs.rhs_begin == rhs.rhs_begin && lhs.rhs_end > rhs.rhs_end);  // NOLINT
      });

  std::size_t max_name_len = 0;
  std::uint32_t max_ref_len = 0;
  for (const auto& it : reference) {
    max_name_len = std::max(max_name_len, it->name.size());
    max_ref_len = std::max(max_ref_len, it->inflated_len);
  }
  std::size_t td = log10(max_ref_len) + 1;

  for (std::uint32_t i = 0, j = 0; i < reference.size(); ++i) {
    std::ofstream os(reference[i]->name + "_reconstruction.gfa");

    std::cout << std::string(max_name_len + 2 * td + 15, '-')
              << std::endl;

    std::cout << reference[i]->name
              << std::string(max_name_len - reference[i]->name.size(), ' ')
              << " =";
    bool is_first = true;

    std::uint32_t rhs_begin = -1, rhs_end;
    for (; j < overlaps.size(); ++j) {
      if (overlaps[j].rhs_id != i) {
        break;
      }
      if (rhs_begin == -1U) {
        rhs_begin = overlaps[j].rhs_begin;
      }
      rhs_end = overlaps[j].rhs_end;

      os << "S\t"
         << sequences[overlaps[j].lhs_id]->name << "\t"
         << sequences[overlaps[j].lhs_id]->InflateData(
              overlaps[j].lhs_begin,
              overlaps[j].lhs_end - overlaps[j].lhs_begin) << "\t"  // NOLINT
         << std::endl;

      for (std::uint32_t k = j + 1; k < overlaps.size(); ++k) {
        if (overlaps[k].rhs_id != i ||
            overlaps[k].rhs_begin > overlaps[j].rhs_end) {
          break;
        }

        os << "L\t"
           << sequences[overlaps[j].lhs_id]->name << "\t"
           << (overlaps[j].strand ? "+" : "-") << "\t"
           << sequences[overlaps[k].lhs_id]->name << "\t"
           << (overlaps[k].strand ? "+" : "-") << "\t"
           << overlaps[j].rhs_end - overlaps[k].rhs_begin << "M"
           << std::endl;
      }

      if (j == overlaps.size() - 1 ||
          overlaps[j + 1].rhs_id > i ||
          overlaps[j + 1].rhs_begin > overlaps[j].rhs_end) {
        if (!is_first) {
          std::cerr << std::string(max_name_len + 2, ' ');
        }
        is_first = false;
        std::cout << " [" << std::setw(td) << rhs_begin << ", "
                          << std::setw(td) << rhs_end << "]"
                  << std::fixed << std::setprecision(3)
                  << " (" << (rhs_end - rhs_begin) /
                      static_cast<double>(reference[i]->inflated_len) << ")"
                  << std::endl;
        rhs_begin = -1;
      }
    }
    os.close();
  }
  std::cout << std::string(max_name_len + 2 * td + 15, '-')
            << std::endl;
}

}  // namespace ratlesnake

int main(int argc, char** argv) {
  std::uint32_t num_threads = 1;
  std::vector<std::uint8_t> results = { 0 };

  std::vector<std::string> input_paths;

  std::string optstr = "t:r:h";
  int arg;
  while ((arg = getopt_long(argc, argv, optstr.c_str(), options, nullptr)) != -1) {  // NOLINT
    switch (arg) {
      case 't': num_threads = atoi(optarg); break;
      case 'r': results.emplace_back(atoi(optarg)); break;
      case 'v': std::cout << VERSION << std::endl; return 0;
      case 'h': Help(); return 0;
      default: return 1;
    }
  }
  if (results.size() > 1) {
    results.erase(results.begin());
  }

  for (std::int32_t i = optind; i < argc; ++i) {
    input_paths.emplace_back(argv[i]);
  }
  if ((input_paths.size() < 1) ||
      (input_paths.size() == 1 && std::accumulate(results.begin(), results.end(), 0) > 0)) {  // NOLINT
    std::cerr << "[ratlesnake::] error: missing input file(s)!" << std::endl;
    Help();
    return 1;
  }
  for (const auto& it : input_paths) {
    if (CreateParser<biosoup::NucleicAcid>(it) == nullptr) {
      return 1;
    }
  }

  auto thread_pool = std::make_shared<thread_pool::ThreadPool>(num_threads);

  for (const auto& it : results) {
    if (it == 0) {
      ratlesnake::Analyse(input_paths[0]);
    } else if (it == 1) {
      ratlesnake::Analyse(input_paths[0], input_paths[1], thread_pool);
    } else if (it == 2) {
      ratlesnake::Annotate(input_paths[0], input_paths[1], thread_pool);
    } else if (it == 3) {
      ratlesnake::Reconstruct(input_paths[0], input_paths[1], thread_pool);
    }
  }

  return 0;
}
