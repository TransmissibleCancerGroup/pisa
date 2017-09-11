//
// Created by Kevin Gori on 11/09/2017.
//

#include "pisaTest.h"
#include <iostream>
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"
#include "readFilter.h"
#include "utility.h"
#include "threadedFermiAssembler.h"

using std::cerr;
using std::cout;
using std::endl;

namespace Pisa {
    void test(const pisaOptions options) {
        cout << "Pisa test" << endl;
        cout << "Working on file: " << options.filename << endl;
        cout << "Using " << options.threads << " threads for assembly" << endl;

        SeqLib::BamReader reader;
        SeqLib::BamWriter writer;
        auto pool = SeqLib::ThreadPool(options.threads);

        reader.SetThreadPool(pool);
        writer.SetThreadPool(pool);

        reader.Open(options.filename);
        SeqLib::BamHeader header = reader.Header();
        cout << header.IDtoName(0) << endl;
        cout << reader << endl;
        SeqLib::BamRecord br;
        std::vector<SeqLib::BamRecord> semimapped;
        std::vector<SeqLib::BamRecord> unmapped;
        ReadFilter filter(10, 10);

        writer.SetHeader(header);
        writer.Open("semimapped.bam");
        writer.WriteHeader();

        unsigned long counter = 0;
        while (reader.GetNextRecord(br)) {
            if (++counter % 100000 == 0) cout << "read " << counter << " " << br << endl;
            switch (filter.check(br)) {
                case Pisa::ReadType::HalfMappedPair:
                    semimapped.push_back(br);
                    break;
                case Pisa::ReadType::UnmappedPair:
                    unmapped.push_back(br);
                    break;
                case Pisa::ReadType::Split:
                    semimapped.push_back(br);
                    break;
                case Pisa::ReadType::FilterFail:
                    break;
                case Pisa::ReadType::MappedPair:
                    break;
                default:
                    assert(false);
                    break;
            }
        }

        cout << "Found " << semimapped.size()
             << " semi-mapped reads, and "
             << unmapped.size() << " unmapped reads" << endl;

        // Build contigs from semi-mappable read set
//    cout << "Training read corrector on " << unmapped.size() + semimapped.size() << " reads" << endl;
//    auto bfc = pisa::bfcTrain(semimapped);
//    pisa::bfcUpdate(unmapped, bfc);
//    cout << "Trained." << endl;
//    cout << "Correcting reads" << endl;
//    pisa::bfcCorrect(semimapped, bfc);
//    pisa::bfcCorrect(unmapped, bfc);
//    cout << "Done." << endl;

        for (auto &rec : semimapped) {
            writer.WriteRecord(rec);
        }

        auto semimappedContigs = fermiAssemble(semimapped, true, false, 4);
        cout << "ASSEMBLE SEMIMAPPED" << endl;
        cout << "assembled " << semimappedContigs.size() << " contigs." << endl;
        printContigs(semimappedContigs, 1000);
        writeContigs("round1.fa", semimappedContigs);

        // Build contigs from unmappable read set
        auto unmappedContigs = Pisa::fermiAssemble(unmapped, true, false, 4);
        cout << "ASSEMBLE UNMAPPED" << endl;
        cout << "assembled " << unmappedContigs.size() << " contigs." << endl;
        printContigs(unmappedContigs, 1000);
        writeContigs("round2.fa", unmappedContigs);

        auto contigsAsBam = toBamRecordVector(semimappedContigs);
        for (auto &ctg : toBamRecordVector(unmappedContigs)) {
            contigsAsBam.push_back(ctg);
        }

        auto contigs = Pisa::fermiAssemble(contigsAsBam, true, false, 4);
        cout << "ASSEMBLE SEMIMAPPED AND UNMAPPED ASSEMBLY PRODUCTS" << endl;
        cout << "assembled " << contigs.size() << " contigs." << endl;
        printContigs(contigs, 1000);
        writeContigs("round3.fa", contigs);

        std::vector<SeqLib::BamRecord> all;
        all.reserve(semimapped.size() + unmapped.size());
        for (auto &rec : semimapped) all.push_back(rec);
        semimapped.clear();
        for (auto &rec : unmapped) all.push_back(rec);
        auto allContigs = Pisa::fermiAssemble(all, true, false, 4);
        cout << "ASSEMBLE SEMIMAPPED AND UNMAPPED READS" << endl;
        cout << "assembled " << allContigs.size() << " contigs." << endl;
        printContigs(allContigs, 1000);
        writeContigs("round4.fa", allContigs);

        contigsAsBam = toBamRecordVector(semimappedContigs);
        for (auto &ctg : unmapped) {
            contigsAsBam.push_back(ctg);
        }

        auto mixed = Pisa::fermiAssemble(contigsAsBam, false, false, 4);
        cout << "ASSEMBLE SEMIMAPPED ASSEMBLY PRODUCTS WITH UNMAPPED READS" << endl;
        cout << "assembled " << mixed.size() << " contigs." << endl;
        printContigs(mixed, 1000);
        writeContigs("round5.fa", mixed);

        writer.Close();
        writer.BuildIndex();
    }
}