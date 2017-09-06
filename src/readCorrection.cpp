//
// Created by Kevin Gori on 26/07/2017.
//

#include "SeqLib/BFC.h"  // for more control with read error correction
#include "SeqLib/BWAWrapper.h"  // for realigning contigs
#include "SeqLib/FermiAssembler.h"  // for SeqLib assembly
#include "SeqLib/BamReader.h"
#include "readCorrection.h"

namespace pisa {
// Use a batch of reads to train an error corrector
    SeqLib::BFC bfcTrain(const SeqLib::BamRecordVector &brv) {
        if (brv.empty()) {
            std::cerr << "Can't train with no reads" << std::endl;
            exit(1);
        }
        auto corrector = SeqLib::BFC();
        for (const auto &br: brv) {
            corrector.AddSequence(br.Sequence().c_str(), br.Qualities().c_str(), br.Qname().c_str());
        }
        corrector.Train();
        corrector.clear();
        return corrector;
    }

    void bfcUpdate(const SeqLib::BamRecordVector &brv, SeqLib::BFC &corrector) {
        if (brv.empty()) {
            std::cerr << "Can't train with no reads" << std::endl;
            return;
        }
        for (const auto &br: brv) {
            corrector.AddSequence(br.Sequence().c_str(), br.Qualities().c_str(), br.Qname().c_str());
        }
        corrector.Train();
        corrector.clear();
    }

// Use a trained read corrector to correct a batch of reads
    void bfcCorrect(SeqLib::BamRecordVector &brv, SeqLib::BFC &readCorrector) {
        if (brv.empty()) return;
        unsigned long nCorrected = 0;
        for (auto &read : brv) {
            readCorrector.AddSequence(read.Sequence().c_str(), read.Qualities().c_str(), read.Qname().c_str());
        }

        if (!readCorrector.ErrorCorrect()) {
            std::cerr << "Fatal error correcting sequences" << std::endl;
            readCorrector.clear();
            exit(1);
        }

        for (auto &i : brv) {
            std::string newseq;
            std::string newqual;
            readCorrector.GetSequence(newseq, newqual);
            if (newseq != i.Sequence()) {
                i.SetSequence(newseq);
                nCorrected++;
            }
        }
        std::cout << "BFC Corrected " << nCorrected << " reads" << std::endl;
        readCorrector.clear();
    }

} // namespace pisa