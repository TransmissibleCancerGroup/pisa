//
// Created by Kevin Gori on 11/09/2017.
//

#include <fstream>
#include <sstream>
#include <iostream>
#include "utility.h"

void printContigs(const std::vector<std::string> &contigs, int minSize) {
    unsigned long counter = 0;
    for (auto &ctg : contigs) {
        if (ctg.size() > minSize) {
            std::cout << ">Contig" << counter++ << "_" << ctg.size() << '\n' << ctg << std::endl;
        }
    }
}

void writeContigs(const std::string &filename, const std::vector<std::string> &contigs) {
    std::ofstream f(filename);
    uint counter = 0;
    for (auto &ctg : contigs) {
        f << ">Contig" << counter++ << '\n' << ctg << std::endl;
    }
    f.close();
}

SeqLib::BamRecordVector toBamRecordVector(const std::vector<std::string> &seqs) {
    SeqLib::BamRecordVector brv;
    brv.reserve(seqs.size());
    unsigned long count = 0;
    for (const auto &seq : seqs) {
        std::stringstream ss;
        ss << "Contig_" << count++;
        SeqLib::BamRecord b;
        b.init();
        b.SetSequence(seq);
        b.SetQualities(std::string(seq.size(), 'F'), 30);
        b.SetQname(ss.str());
        brv.push_back(b);
    }
    return brv;
}

int readLength(SeqLib::BamReader &reader) {
    SeqLib::BamRecord br;
    std::vector<int> lengths;
    while (reader.GetNextRecord(br)) {
        lengths.push_back(br.Length());
        if (lengths.size() > 100) break;
    }
    reader.Reset();
    return *std::max_element(lengths.begin(), lengths.end());
}