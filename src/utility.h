//
// Created by Kevin Gori on 11/09/2017.
//

#ifndef PISA_UTILITY_H
#define PISA_UTILITY_H
#include <vector>
#include "SeqLib/BamRecord.h"
#include "SeqLib/BamReader.h"

void printContigs(const std::vector<std::string> &contigs, int minSize = 250);
void writeContigs(const std::string &filename, const std::vector<std::string> &contigs);
SeqLib::BamRecordVector toBamRecordVector(const std::vector<std::string> &seqs);
int readLength(SeqLib::BamReader &reader);

#endif //PISA_UTILITY_H
