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

struct Vertex {
    Vertex(Overlap read) : read(read) {}

    Overlap read;
    std::vector<Vertex*> vertices;
    Vertex* parent = NULL;
};

std::unordered_map<std::string, std::uint32_t> Overlap::name_to_id;

inline bool isSuffix(const std::string& src, const std::string& suffix) {
    return src.size() < suffix.size() ? false :
        src.compare(src.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> annotate(std::vector<Annotation>& dst,
    std::vector<std::unique_ptr<Overlap>>& overlaps);

void reconstruct(std::vector<Annotation>& dst,
    std::vector<std::unique_ptr<Overlap>>& overlaps);

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

    auto ret = annotate(annotations, overlaps);
    reconstruct(annotations, overlaps);

    return 0;
}

struct Mapping {
    std::uint32_t seq_start;
    std::uint32_t seq_end;
    std::uint32_t gen_start;
    std::uint32_t gen_end;
    std::uint32_t gen_length;
    std::uint64_t ref_id;

    Mapping(std::uint32_t seq_s, std::uint32_t seq_e, std::uint32_t gen_s, std::uint32_t gen_e, std::uint32_t gen_l, std::uint64_t ref_n) :
        seq_start(seq_s),
        seq_end(seq_e),
        gen_start(gen_s),
        gen_end(gen_e),
        gen_length(gen_l),
        ref_id(ref_n)
    { }
};

bool sort_paf(const std::unique_ptr<Overlap>& s1, const std::unique_ptr<Overlap>& s2) {
    if (s1->q_id == s2->q_id) {
        if(s1->q_begin == s2->q_begin) {
            return (s1->q_end > s2->q_end);
        }
        return (s1->q_begin < s2->q_begin);
    }
    return (s1->q_id < s2->q_id);
}

bool sort_reference_repeating_sections(std::pair<std::uint64_t, std::uint64_t>& p1, std::pair<std::uint64_t, std::uint64_t>& p2) {
    if (std::get<0>(p1) == std::get<0>(p2)) {
        return std::get<1>(p1) > std::get<1>(p2);
    }
    return std::get<0>(p1) < std::get<0>(p2);
}

bool unique_reference_repeating_sections(std::pair<std::uint64_t, std::uint64_t>& p1, std::pair<std::uint64_t, std::uint64_t>& p2){
    return (std::get<0>(p1) == std::get<0>(p2) && std::get<1>(p1) == std::get<1>(p2));
}

std::vector<std::pair<std::uint64_t, std::uint64_t>> annotate(std::vector<Annotation>& dst,
    std::vector<std::unique_ptr<Overlap>>& overlaps) {

    std::sort(overlaps.begin(), overlaps.end(), sort_paf);

    std::unordered_map<std::uint64_t, std::vector<Mapping>> sequence_mapping_details;

    for(auto &i : overlaps) {
        std::uint64_t key = i->q_id;
        Mapping mapp(i->q_begin, i->q_end, i->t_begin, i->t_end, i->t_length, i->t_id);
        if(sequence_mapping_details.find(key) == sequence_mapping_details.end()) {
            std::vector<Mapping> vec;
            vec.push_back(mapp);
            sequence_mapping_details[key] = vec;
        } else {
            (sequence_mapping_details[key]).push_back(mapp);
        }
    }

    std::unordered_map<std::uint64_t, std::vector<Mapping>> chimeric_reads;
    std::unordered_map<std::uint64_t, std::vector<Mapping>> repeating_reads;
    std::unordered_set<std::uint64_t> chimers;
    std::unordered_set<std::uint64_t> repeatings;
    std::vector<Mapping> all_repeatings;

    for (auto itr = sequence_mapping_details.begin(); itr != sequence_mapping_details.end(); itr++) {
        if ((itr->second).size() > 1) {
            uint32_t seq_start = (itr->second)[0].seq_start;
            uint32_t seq_end = (itr->second)[0].seq_end;
            uint32_t gen_start = (itr->second)[0].gen_start;
            uint32_t gen_end = (itr->second)[0].gen_end;
            for (int i = 1; i < (itr->second).size(); i++) {
                if (!(seq_start <= (itr->second)[i].seq_start && seq_end >= (itr->second)[i].seq_end)) {
                    chimers.emplace(itr->first);
                    int seq_gap = abs((int)(std::max(seq_start, (itr->second)[i].seq_start) - std::min(seq_end, (itr->second)[i].seq_end)));
                    int ref_gap = abs((int)(std::max(gen_start, (itr->second)[i].gen_start) - std::min(seq_end, (itr->second)[i].seq_end)));
                    if (std::min(gen_start, (itr->second)[i].gen_start) > 100 || std::max(gen_end, (itr->second)[i].gen_end) < ((itr->second)[i].gen_length - 100)
                    && (std::min(ref_gap, seq_gap)/std::max(ref_gap, seq_gap)) > 0.12) {
                        chimeric_reads[itr->first] = itr->second;
                        break;
                    }
                }
            } if (chimers.find(itr->first) == chimers.end()) {
                int i = 1;
                if ((itr->second)[0].seq_start == (itr->second)[1].seq_start && (itr->second)[0].seq_end == (itr->second)[1].seq_end) {
                    i = 0;
                }
                all_repeatings.insert(all_repeatings.end(), (itr->second).begin() + i, (itr->second).end());
                repeating_reads[itr->first] = itr->second;
                repeatings.emplace(itr->first);
            }
        }
    }

    std::vector<std::pair<std::uint64_t, std::uint64_t>> reference_repeating_regions;

    for (int i = 0; i < all_repeatings.size(); i++) {
        for(int j = 0; j < all_repeatings.size(); j++) {
            if (all_repeatings[i].ref_id == all_repeatings[j].ref_id &&  i != j) {
                if (abs((int)(all_repeatings[i].gen_start - all_repeatings[j].gen_start)) + abs((int)(all_repeatings[i].gen_end - all_repeatings[j].gen_end)) < 1000) {
                    std::uint64_t position = all_repeatings[i].gen_start;
                    position = position << 32;
                    position += all_repeatings[i].gen_end;
                    reference_repeating_regions.push_back(std::make_pair(all_repeatings[i].ref_id, position));
                }
            }
        }
    }


    std::sort(reference_repeating_regions.begin(), reference_repeating_regions.end(), sort_reference_repeating_sections);
    std::vector<std::pair<std::uint64_t, std::uint64_t>>::iterator vec_it;
    vec_it = std::unique(reference_repeating_regions.begin(), reference_repeating_regions.end(), unique_reference_repeating_sections);
    reference_repeating_regions.resize(std::distance(reference_repeating_regions.begin(), vec_it));

    for(auto itr = chimeric_reads.begin(); itr != chimeric_reads.end(); itr++) {
        Annotation ann = dst[itr->first];
        std::vector<uint64_t> chimeric_annotations;
        for(int i = 0; i < (itr->second).size() - 1; i++) {
            if((itr->second)[i].seq_start != (itr->second)[i+1].seq_start && (itr->second)[i].seq_end != (itr->second)[i+1].seq_end){
                uint64_t positions = (itr->second)[i].seq_end;
                positions = positions << 32;
                positions = positions | (itr->second)[i+1].seq_start;
                chimeric_annotations.push_back(positions);

            }
        }
        ann.chimeric = chimeric_annotations;
        dst[itr->first] = ann;
    }

    for(auto itr = repeating_reads.begin(); itr != repeating_reads.end(); itr++) {
        Annotation ann = dst[itr->first];
        std::vector<uint64_t> repeating_annotations;
        int i = 1;
        if ((itr->second[0]).seq_start == (itr->second[1]).seq_start && (itr->second[0]).seq_end == (itr->second[1]).seq_end) {
            i = 0;
        }
        while(i < (itr->second).size()) {
            uint64_t positions = (itr->second)[i].seq_start;
            positions = positions << 32;
            positions = positions | (itr->second)[i].seq_end;
            if(std::find(repeating_annotations.begin(), repeating_annotations.end(), positions) == repeating_annotations.end()) {
                repeating_annotations.push_back(positions);
            }
            i++;
        }
        ann.repeat = repeating_annotations;
        dst[itr->first] = ann;
    }

    return reference_repeating_regions;
}

std::vector<Vertex*> DepthFirstSearch(std::unordered_map<std::uint64_t, std::vector<Vertex*>> heads);
void clear_contained_reads(std::vector<std::unique_ptr<Overlap>> &overlaps, std::vector<Annotation>& dst);
std::unordered_map<std::uint64_t, std::vector<Vertex*>> create_graph(std::vector<Vertex> &vertices, std::vector<std::unique_ptr<Overlap>> &overlaps);
void statistics(std::vector<Vertex*> ends);

void reconstruct(std::vector<Annotation>& dst,
    std::vector<std::unique_ptr<Overlap>>& overlaps) {
    clear_contained_reads(overlaps, dst);
    //remove_covered_repeats(repeats, paf_objects);
    std::vector<Vertex> vertices;
    std::unordered_map<std::uint64_t, std::vector<Vertex*>> heads = create_graph(vertices, overlaps);
    std::vector<Vertex*> ends = DepthFirstSearch(heads);
    statistics(ends);
}

bool paf_unique(const std::unique_ptr<Overlap>& a, const std::unique_ptr<Overlap>& b) {
    return a->q_id == b->q_id;
}

void clear_contained_reads(std::vector<std::unique_ptr<Overlap>> &overlaps, std::vector<Annotation>& dst) {
    std::vector<std::unique_ptr<Overlap>>::iterator it = overlaps.begin();
    auto long_cmp = [](const std::unique_ptr<Overlap>& a, const std::unique_ptr<Overlap>& b) {
        if (a->t_id == b->t_id) {
            return (a->t_end - a->t_begin > b->t_end - b->t_begin);
        }
        return (a->t_id > b->t_id);
    };
    std::sort(overlaps.begin(), overlaps.end(), long_cmp);
    it = std::unique (overlaps.begin(), overlaps.end(), paf_unique);
    overlaps.resize(std::distance(overlaps.begin(), it));
    auto start_cmp = [](const std::unique_ptr<Overlap>& a, const std::unique_ptr<Overlap>& b) {
        if (a->t_id == b->t_id) {
            if (a->t_begin == b->t_begin) {
                return (a->t_end > b->t_end);
            }
            return (a->t_begin < b->t_begin);
        }
        return (a->t_id > b->t_id);
    };
    std::sort(overlaps.begin(), overlaps.end(), start_cmp);
    it = overlaps.begin();
    std::vector<std::unique_ptr<Overlap>>::iterator next;
    std::vector<std::unique_ptr<Overlap>>::iterator to_remove;
    while (it != --overlaps.end()) {
        do {next = std::next(it);} while ((*next) == NULL);
        while ((*next)->t_end <= (*it)->t_end) {
            to_remove = next;
            do {
                next++;
            } while ((*next) == NULL && next != overlaps.end());
            (*to_remove) = NULL;
            if(next == overlaps.end()) {
                overlaps.erase(std::remove_if(overlaps.begin(), overlaps.end(), [](std::unique_ptr<Overlap> &p) {return p == NULL;}), overlaps.end());
                return;
            }
        }
        it = next;
    }
    overlaps.erase(std::remove_if(overlaps.begin(), overlaps.end(), [](std::unique_ptr<Overlap> &p) {return p == NULL;}), overlaps.end());
}

std::unordered_map<std::uint64_t, std::vector<Vertex*>> create_graph(std::vector<Vertex> &vertices, std::vector<std::unique_ptr<Overlap>> &overlaps) {
    for (auto const& overlap: overlaps) {
        vertices.emplace_back(Vertex(*overlap));
    }
    std::vector<Vertex>::iterator it = vertices.begin();
    std::vector<Vertex>::iterator next;
    std::vector<Vertex*> heads;
    std::unordered_map<std::uint64_t, std::vector<Vertex*>> re;
    heads.emplace_back(&(*it));
    while(it != --vertices.end()) {
        next = std::next(it);
        if ((*next).read.t_id != (*it).read.t_id) {
            re[(*it).read.t_id] = heads;
            heads.clear();
            heads.emplace_back(&(*next));
            it++;
            continue;
        }
        while ((*next).read.t_begin <= (*it).read.t_end && (*next).read.t_id == (*it).read.t_id) {
            it->vertices.emplace_back(&(*next));
            next++;
            if(next == vertices.end()) break;
        }
        if (next == std::next(it) && (*next).read.t_id == (*std::next(it)).read.t_id) heads.emplace_back(&(*next));
        it++;
    }
    re[(*it).read.t_id] = heads;
    heads.clear();
    for (auto const& vertex : vertices) {
        printf("S\t%llu\t%c\tLN:i:%u\n", vertex.read.q_id, '*', vertex.read.q_length);
        for (auto const& next: vertex.vertices) {
            printf("L\t%llu\t%c\t%llu\t%c\t%c\n", vertex.read.q_id, '+', next->read.q_id, '-', '*');
        }
    }
    return re;
}

std::vector<Vertex*> DepthFirstSearch(std::unordered_map<std::uint64_t, std::vector<Vertex*>> heads) {
    Vertex* max;
    Vertex* head;
    Vertex* longest;
    std::vector<Vertex*> ends;
    int begin;
    for (auto const& seq: heads) {
        unsigned int length=0;
        for (unsigned int i = 0; i < seq.second.size(); i++) {
            head = seq.second[i];
            begin = head->read.t_begin;
            while(!head->vertices.empty()) {
                max = head->vertices[0];
                for (auto const& vertex : head->vertices) {
                    if (vertex->read.t_begin > max->read.t_begin) max = vertex;
                }
                max->parent = head;
                head = max;
            }
            if (length < head->read.t_begin - begin) {
                length = head->read.t_begin - begin;
                longest = head;
            }
        }
        ends.emplace_back(longest);
    }
    return ends;
}

void statistics(std::vector<Vertex*> ends) {
    int count, last_index;
    for (auto const& end: ends) {
        count = 1;
        auto curr = end;
        last_index = end->read.t_end;
        while (curr->parent != NULL) {
            count++;
            curr = curr->parent;
        }
        fprintf(stderr, "Genome coverage: %f%%\n", (last_index - curr->read.t_begin) / (float) end->read.t_length * 100);
        fprintf(stderr, "Number of used reads: %d\n", count);
    }
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
