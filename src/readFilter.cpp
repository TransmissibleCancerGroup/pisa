//
// Created by Kevin Gori on 31/07/2017.
//

#include "SeqLib/BFC.h"  // for more control wiht read error correction
#include "SeqLib/BWAWrapper.h"  // for realigning contigs
#include "SeqLib/FermiAssembler.h"  // for pure SeqLib assembly
#include "SeqLib/BamReader.h"
#include <SeqLib/BamWriter.h>
#include "readFilter.h"

namespace Pisa {

// Handles the read-identifying logic
    ReadType ReadFilter::check(const SeqLib::BamRecord &br) {
        if (br.DuplicateFlag() || br.QCFailFlag()) {
            return ReadType::FilterFail;
        }
        // Ignore low quality
        if (br.MeanPhred() < _minMeanPhred) {
            return ReadType::FilterFail;
        }
        // Ignore very short alignments
        if (br.Sequence().size() < 25) {
            return ReadType::FilterFail;
        }
        // Ignore supplementary alignments
        if (br.SecondaryFlag()) {
            return ReadType::FilterFail;
        }

        if (br.CountBWAChimericAlignments() > 1) {
            return ReadType::FilterFail;
        }

        if (br.MappedFlag()) {
            if (br.MapQuality() < _minMapQual) {
//            std::cout << br << std::endl;
                return ReadType::FilterFail;
            }
            if (!br.MateMappedFlag()) {
                return ReadType::HalfMappedPair;
            }

            // collect Split-Mapped / Split-Split / Unmapped-Split
            const auto cig = br.GetCigar();
            const auto frField = cig.front();
            const auto bkField = cig.back();
            if ((frField.Type() == 'S' || frField.Type() == 'H') && frField.Length() > 25) {
                return ReadType::Split;
            }
            if ((bkField.Type() == 'S' || bkField.Type() == 'H') && bkField.Length() > 25) {
                return ReadType::Split;
            }
            return ReadType::MappedPair;
        } else if (br.MateMappedFlag()) {
            return ReadType::HalfMappedPair;
        } else { // collect Unmapped-Unmapped
            return ReadType::UnmappedPair;
        }
    }

} // namespace pisa