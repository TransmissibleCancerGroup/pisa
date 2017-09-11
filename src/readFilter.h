//
// Created by Kevin Gori on 31/07/2017.
//

#ifndef PISA_READFILTER_H
#define PISA_READFILTER_H
#include "SeqLib/BamRecord.h"

namespace Pisa {

    enum ReadType {
        HalfMappedPair,
        MappedPair,
        Split,
        UnmappedPair,
        FilterFail
    };

    class ReadFilter {
    public:
        ReadFilter() = default;

        ReadFilter(int minMapQual, int minMeanPhred) : _minMapQual(minMapQual), _minMeanPhred(minMeanPhred) {}

        ReadType check(const SeqLib::BamRecord &br);

    private:
        int _minMapQual{0};
        int _minMeanPhred{0};
    };
} // namespace pisa

#endif //PISA_READFILTER_H
