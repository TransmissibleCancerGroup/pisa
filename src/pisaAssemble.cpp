//
// Created by Kevin Gori on 11/09/2017.
//

#include "pisaAssemble.h"
#include "threadedFermiAssembler.h"
#include "utility.h"
#include "SeqLib/BamWriter.h"
using std::cerr;
using std::cout;
using std::endl;

namespace Pisa {
    void assemble(const pisaOptions options) {

        cout << "Pisa assemble" << endl;
        std::string outfile = options.filename.substr(0, options.filename.find_last_of(".")) + ".contigs.fa";
        cout << "Assembling " << options.filename << endl;
        SeqLib::BamReader reader;
        if (options.threads > 1) reader.SetThreadPool(options.threads);
        reader.Open(options.filename);

        SeqLib::BamRecord br;
        std::vector<SeqLib::BamRecord> brv;
        while (reader.GetNextRecord(br)) {
            brv.push_back(br);
        }

        cout << "Found " << brv.size() << " reads." << endl;
        Pisa::ThreadedFermiAssembler fermi;
        auto contigs = Pisa::fermiAssemble(brv, true, false, options.threads);
        cout << "Writing " << contigs.size() << " contigs to " << outfile << endl;
        printContigs(contigs);
            //writeContigs(outfile, contigs);
    }
}