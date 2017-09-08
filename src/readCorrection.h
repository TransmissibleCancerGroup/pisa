//
// Created by Kevin Gori on 26/07/2017.
//

#ifndef PISA_READCORRECTION_H
#define PISA_READCORRECTION_H

namespace Pisa {
// Use a batch of reads to train an error corrector
    SeqLib::BFC bfcTrain(const SeqLib::BamRecordVector &brv);

// Use a batch of reads to update a trained BFC error corrector
    void bfcUpdate(const SeqLib::BamRecordVector &brv, SeqLib::BFC &corrector);

// Use a trained read corrector to correct a batch of reads
    void bfcCorrect(SeqLib::BamRecordVector &brv, SeqLib::BFC &readCorrector);

} // namespace pisa
#endif //PISA_READCORRECTION_H
