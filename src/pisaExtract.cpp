//
// Created by Kevin Gori on 11/09/2017.
//

#include <array>
#include "pisaExtract.h"
#include "readFilter.h"
#include "SeqLib/BamWriter.h"
#include "SeqLib/BamReader.h"

using std::cerr;
using std::cout;
using std::endl;

namespace Pisa {
    void extract(const pisaOptions options) {
        const std::string filename = options.filename;
        const int nthreads = options.threads;

        cout << "Pisa extract" << endl;
        cout << "Scanning " << filename << " for insertion-associated reads, using "
             << nthreads << " threads" << endl;

        // Set up a thread pool for faster IO
        auto pool = SeqLib::ThreadPool(nthreads);

        SeqLib::BamReader reader;
        if (nthreads > 1) reader.SetThreadPool(pool);

        // Open the input file
        if (!reader.Open(filename)) {
            cerr << "Failed to open " << filename << endl;
            std::exit(1);
        }

        SeqLib::BamHeader header = reader.Header();
        cout << header.IDtoName(0) << endl;
        cout << reader << endl;
        SeqLib::BamRecord br;
        Pisa::ReadFilter filter(10, 20);

        // Three output files
        std::array<SeqLib::BamWriter, 3> outfiles; // semimapped, unmapped and split
        std::array<std::string, 3> filenames = {"semimapped.bam", "unmapped.bam", "split.bam"};

        for (int i = 0; i < outfiles.size(); i++) {
            if (nthreads > 1) {
                outfiles[i].SetThreadPool(pool);
            }
            outfiles[i].SetHeader(header);
            outfiles[i].Open(filenames[i]);
            outfiles[i].WriteHeader();
        }

        unsigned long counter = 0;
        while (reader.GetNextRecord(br)) {
            if (++counter % 1000000 == 0) cout << "read " << counter << " " << br << endl;
            switch (filter.check(br)) {
                case Pisa::ReadType::HalfMappedPair:
                    outfiles[0].WriteRecord(br);
                    break;
                case Pisa::ReadType::UnmappedPair:
                    outfiles[1].WriteRecord(br);
                    break;
                case Pisa::ReadType::Split:
                    outfiles[2].WriteRecord(br);
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
        reader.Close();
        cout << "File scan ended" << endl;

        // Close files and clean up
        for (int i = 0; i < 3; i++) {
            outfiles[i].Close();
            cout << "Wrote " << filenames[i] << endl;
            outfiles[i].BuildIndex();
        }
    }
}