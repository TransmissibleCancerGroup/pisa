#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <svaba/SeqLib/SeqLib/BamWriter.h>
//#include "svabaAssemblerEngine.h"  // for SvABA assembly
//#include "KmerFilter.h"  // for SvABA assembly
#include "SeqLib/BamReader.h"
#include "SeqLib/FermiAssembler.h"  // for pure SeqLib assembly
#include "SeqLib/BWAWrapper.h"  // for realigning contigs
#include "SeqLib/BFC.h"  // for more control wiht read error correction
#include "ctpl/ctpl_stl.h"  // threadpool
//#include "readCorrection.h"
//#include "readFilter.h"

using namespace SeqLib;
using std::cerr;
using std::cout;
using std::endl;

static const char *TOWER =
        "      PISA \n"
                "    -------\n"
                "   /n n n/ \n"
                "  /n n n/  \n"
                " /n n n/   \n"
                "/ / / /    \n";

static const char *USAGE =
        "Usage: pisa <bamfile> <reference.fa>\n";


//svabaReadVector toSvabaReadVector(const SeqLib::BamRecordVector &brv) {
//    svabaReadVector srv;
//    for (const auto &br : brv) {
//        svabaRead r(br, "t001");
//        r.SetSeq(br.Sequence());
//        srv.push_back(r);
//    }
//    return srv;
//}


// Does read correction on SvabaReadVector and updates reads inplace.
// The corrected sequence is accessed by its .Seq() accessor.
// The uncorrected sequence still exists and is accessed by its inherited
// BamRecord .Sequence() accessor.
//svabaReadVector &kmerCorrect(svabaReadVector &srv) {
//    KmerFilter kmer;
//    int kmer_corrected = 0;
//
//    std::vector<char*> charVec;  // kmer correction needs raw char pointers
//    charVec.reserve(srv.size());
//    for (const auto &r : srv) {
//        charVec.push_back(strdup(r.Seq().c_str()));  // strdup copy needs to be freed
//    }
//    kmer.makeIndex(charVec);
//
//    for (auto &i : charVec) {
//        if (i != nullptr) free(i);  // free memory alloc'ed in strdup
//    }
//
//    kmer_corrected = kmer.correctReads(srv);
//    cout << "Corrected " << kmer_corrected << " reads." << endl;
//
//    return srv;
//}


//std::vector<std::string> svabaAssemble(const SeqLib::BamRecordVector &brv, bool correct = false) {
//    auto srv = toSvabaReadVector(brv);
//    if (correct) kmerCorrect(srv);
//    auto engine = svabaAssemblerEngine("test", 0.025, 40, 125);
//    engine.fillReadTable(srv);
//    engine.performAssembly(3);
//
//    std::vector<std::string> contigs;
//    for (auto &contig : engine.getContigs()) {
//        contigs.push_back(contig.Seq);
//    }
//
//    sort(contigs.begin(), contigs.end(),
//         [](std::string a, std::string b) {
//             return a.size() > b.size();
//         }
//    );
//    return contigs;
//}


// Correct, filter and assemble reads in a BamRecordVector, using Fermi
std::vector<std::string> fermiAssemble(const SeqLib::BamRecordVector &brv,
                                       bool correct = false,
                                       bool aggressiveTrim = false,
                                       uint32_t minOverlap = 33) {
    if (brv.empty()) {
        return std::vector<std::string>();
    }

    FermiAssembler fermi;
    if (aggressiveTrim) fermi.SetAggressiveTrim();
    fermi.SetMinOverlap(minOverlap > 63 ? 63 : minOverlap);
    fermi.AddReads(brv);
    if (correct) fermi.CorrectReads();
    fermi.PerformAssembly();

    std::vector<std::string> contigs = fermi.GetContigs();
    sort(contigs.begin(), contigs.end(),
         [](std::string a, std::string b) {
             return a.size() > b.size();
         }
    );
    return contigs;
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


std::vector<std::string> fermiAssembleStrings(const std::vector<std::string> &seqs,
                                              bool correct = false) {
    SeqLib::BamRecordVector brv = toBamRecordVector(seqs);
    return fermiAssemble(brv, correct);
}


void printContigs(const std::vector<std::string> &contigs, int minSize = 250) {
    unsigned long counter;
    for (auto &ctg : contigs) {
        if (ctg.size() > minSize) {
            cout << ">Contig" << counter++ << "_" << ctg.size() << '\n' << ctg << endl;
        }
    }
}


void writeContigs(const std::string &filename, const std::vector<std::string> &contigs) {
    std::ofstream f(filename);
    uint counter = 0;
    for (auto &ctg : contigs) {
        f << ">Contig" << counter++ << '\n' << ctg << endl;
    }
    f.close();
}


enum Correction {
    NONE,
    BFC,
    FERMI
};


struct LockingBamWriter {
    LockingBamWriter(const std::string &outfile, const SeqLib::BamHeader &header) {
        writer.SetHeader(header);
        writer.Open(outfile);
        writer.WriteHeader();
    }

    LockingBamWriter() = delete;

    LockingBamWriter(LockingBamWriter const &) = delete;

    LockingBamWriter &operator=(const LockingBamWriter &) = delete;

    LockingBamWriter &operator=(LockingBamWriter &&) = delete;

    LockingBamWriter(LockingBamWriter &&) = delete;

    ~LockingBamWriter() {
        Close();
    }

    bool BuildIndex() {
        std::lock_guard<std::mutex> lock(mutex);
        return writer.BuildIndex();
    }

    bool Close() {
        std::lock_guard<std::mutex> lock(mutex);
        return writer.Close();
    }

    bool WriteRecord(const SeqLib::BamRecord &record) {
        std::lock_guard<std::mutex> lock(mutex);
        return writer.WriteRecord(record);
    }

private:
    SeqLib::BamWriter writer;
    std::mutex mutex;
};


void processRegion(int thread_id, std::vector<SeqLib::BamReader> &readers,
                   std::vector<std::mutex> &mutexes,
                   const SeqLib::GenomicRegion &region,
                   const std::string &reference,
                   LockingBamWriter &writer,
                   unsigned int region_id,
                   Correction correction = Correction::BFC) {

    BamRecord br;
    BamRecordVector all_reads;
    BamRecordVector half_mapped;
    BamRecordVector split;

    unsigned long nUnmapped = 0;
    {
        std::lock_guard<std::mutex> lock(mutexes.at(thread_id));
        readers.at(thread_id).SetRegion(region);
        while (readers.at(thread_id).GetNextRecord(br)) {
            // Ignore dupes and failures
            if (br.DuplicateFlag() || br.QCFailFlag()) {
                continue;
            }
            // Ignore low quality
            if (br.MeanPhred() < 0 || br.MapQuality() < 0) {
                continue;
            }
            // Ignore very short alignments
            if (br.Sequence().size() < 25) {
                continue;
            }

            all_reads.push_back(br);

            // Read is half mapped, so keep
            if ((br.MappedFlag() && !br.MateMappedFlag()) ||
                (!br.MappedFlag() && br.MateMappedFlag())) {
                half_mapped.push_back(br);
            }

            // Read is split, so keep
            if (br.MappedFlag()) {
                SeqLib::Cigar cigar = br.GetCigar();
                auto front = cigar.front();
                auto back = cigar.front();
                if ((front.Type() == 'S' || front.Type() == 'H') && front.Length() > 10) split.push_back(br);
                else if ((back.Type() == 'S' || back.Type() == 'H') && back.Length() > 10) split.push_back(br);
            } else {
                nUnmapped++;
            }
        }
    }

    if (nUnmapped == 0) {
        cerr << "No unmapped reads found" << endl;
        return;
    }

    cout << "Processing " << region.ChrName(readers.at(thread_id).Header())
         << ":" << region.pos1 << "-" << region.pos2 << endl;

    // Read correction
    if (correction == Correction::BFC) {
        SeqLib::BFC readCorrector = pisa::bfcTrain(all_reads);
        pisa::bfcCorrect(half_mapped, readCorrector);
        // bfcCorrect(unmapped, readCorrector);
    }

    // Do assembly!
    auto fe_half_mapped_contigs = fermiAssemble(half_mapped, correction == Correction::FERMI);
    // auto fe_unmapped_contigs = svabaAssemble(unmapped); //fermiAssemble(unmapped, correction == Correction::FERMI);

//    printContigs(fe_half_mapped_contigs);
//    printContigs(fe_unmapped_contigs);

    // Align contigs back to reference
    SeqLib::BWAWrapper bwa;
    SeqLib::BamRecordVector aligned;
    bwa.LoadIndex(reference);

    unsigned int contigNumber = 0;
    for (auto &contig : fe_half_mapped_contigs) {
        std::stringstream ss;
        ss << "contig_" << thread_id << "_" << region_id << "_" << contigNumber++;
        bwa.AlignSequence(contig, ss.str(), aligned, true, .75, 1);
    }
    for (auto &alignedBamRec : aligned) {
        writer.WriteRecord(alignedBamRec);
    }
}


struct RegionGenerator {
    RegionGenerator(const SeqLib::BamHeader &header, int blocksize, int pad) : _header(header) {
        for (int i = 0; i < _header.NumSequences(); i++) {
            for (int j = 1, chrlen = _header.GetSequenceLength(i); j < chrlen; j += blocksize) {
                std::stringstream ss;
                ss << _header.IDtoName(i) << ":" << j << "-" << std::min(j + blocksize - 1, chrlen);
                cerr << ss.str() << endl;
                GenomicRegion region{ss.str(), header};
                if (region.pos1 > pad) region.pos1 -= pad;
                if (region.pos2 + pad < chrlen) region.pos2 += pad;
                regions.push_back(region);
            }
        }
    }

    SeqLib::GenomicRegionVector regions;

    std::string regionToString(const SeqLib::GenomicRegion &reg) {
        std::stringstream ss;
        ss << reg.ChrName(_header) << ":" << reg.pos1 << "-" << reg.pos2;
        return ss.str();
    }

private:
    SeqLib::BamHeader _header;
};


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


void run(unsigned int nthread, const std::string &filename, const std::string &reference,
         Correction correctionMethod = Correction::BFC) {

    nthread = nthread > std::thread::hardware_concurrency() ? std::thread::hardware_concurrency() : nthread;
    nthread = nthread < 1 ? 1 : nthread;
    cout << "USING " << nthread << " THREAD" << (nthread > 1 ? "S" : "") << endl;

    // Maintain pool of open BamReaders here (better if thread-local, but
    // not supported in CTPL threadpool and I don't want to write my own)

    SeqLib::BamReader reader;
    reader.Open(filename);
    int rl = readLength(reader);
    SeqLib::BamHeader header = reader.Header();
    LockingBamWriter writer("/Users/kg8/scratch/output.bam", header);

    ctpl::thread_pool pool;

    if (nthread > 1) {
        pool.resize(nthread - 1);
    }

    RegionGenerator gen(header, 25000, rl);
/*
    // Dispatch regions to be processed by threads
    unsigned int regionNumber = 0;
    for (int i = 0; i < header.NumSequences(); i++) {
        for (int j = 1; j < header.GetSequenceLength(i); j += 100000) {
            std::stringstream ss;
            ss << header.IDtoName(i) << ":" << j << "-" << j + 99999;
            cerr << ss.str() << endl;
            GenomicRegion region{ss.str(), header};
            if (nthread > 1) {
                pool.push(
                        processRegion, std::ref(readers), std::ref(mutexes),
                        region, reference, std::ref(writer), regionNumber++,
                        correctionMethod
                );
            }
            else {
                processRegion(0, readers, mutexes, region, reference, writer, regionNumber++, correctionMethod);
            }
        }
    }

    SeqLib::BamReader reader;
    reader.Open(filename);
    SeqLib::BamRecord br;
    SeqLib::BamRecordVector double_unmapped;
    while (reader.GetNextRecord(br)) {
        if (!br.MappedFlag() && !br.MateMappedFlag() && br.MeanPhred() > 30) {
            double_unmapped.push_back(br);
        }
    }

    if (!double_unmapped.empty()) {
        if (correctionMethod == Correction::BFC) {
            SeqLib::BFC readCorrector = bfcTrain(double_unmapped);
            bfcCorrect(double_unmapped, readCorrector);
        }
        auto contigs = fermiAssemble(double_unmapped, correctionMethod == Correction::FERMI);

        // Align contigs back to reference
        SeqLib::BWAWrapper bwa;
        SeqLib::BamRecordVector aligned;
        bwa.LoadIndex(reference);

        unsigned int contigNumber = 0;
        for (auto &contig : contigs) {
            std::stringstream ss;
            ss << "contig_" << 0 << "_" << -1 << "_" << contigNumber++;
            bwa.AlignSequence(contig, ss.str(), aligned, true, .75, 1);
        }
        for (auto &br : aligned) {
            writer.WriteRecord(br);
        }
    }
    if (nthread > 1) pool.stop(true);

    for (auto &reader : readers) reader.Close();
    reader.Close();
    writer.Close();
    // writer.BuildIndex(); // needs to be sorted
    */
}


void extract(const std::string &filename) {
    cout << "Scanning " << filename << " for insertion-associated reads" << endl;

    // Open the input file
    SeqLib::BamReader reader;
    reader.Open(filename);
    SeqLib::BamHeader header = reader.Header();
    cout << header.IDtoName(0) << endl;
    cout << reader << endl;
    SeqLib::BamRecord br;
    pisa::ReadFilter filter(10, 20);

    // Three output files
    SeqLib::BamWriter w_semimapped;
    SeqLib::BamWriter w_unmapped;
    SeqLib::BamWriter w_split;
    w_semimapped.SetHeader(header);
    w_unmapped.SetHeader(header);
    w_split.SetHeader(header);
    w_semimapped.Open("semimapped.bam");
    w_semimapped.WriteHeader();
    w_unmapped.Open("unmapped.bam");
    w_unmapped.WriteHeader();
    w_split.Open("split.bam");
    w_split.WriteHeader();

    unsigned long counter = 0;
    while (reader.GetNextRecord(br)) {
        if (++counter % 1000000 == 0) cout << "read " << counter << " " << br << endl;
        switch (filter.check(br)) {
            case pisa::ReadType::HalfMappedPair:
                w_semimapped.WriteRecord(br);
                break;
            case pisa::ReadType::UnmappedPair:
                w_unmapped.WriteRecord(br);
                break;
            case pisa::ReadType::Split:
                w_split.WriteRecord(br);
                break;
            case pisa::ReadType::FilterFail:
                break;
            case pisa::ReadType::MappedPair:
                break;
            default:
                assert(false);
                break;
        }
    }
    reader.Close();
    cout << "File scan ended" << endl;

    // Close files and clean up
    w_semimapped.Close();
    w_unmapped.Close();
    w_split.Close();

    w_semimapped.BuildIndex();
    w_unmapped.BuildIndex();
    w_split.BuildIndex();
}


void assemble(const std::string &filename) {
    std::string outfile = filename.substr(0, filename.find_last_of(".")) + ".contigs.fa";
    cout << "Assembling " << filename << endl;
    SeqLib::BamReader reader;
    reader.Open(filename);

    SeqLib::BamRecord br;
    SeqLib::BamRecordVector brv;
    while (reader.GetNextRecord(br)) {
        brv.push_back(br);
    }

    cout << "Found " << brv.size() << " reads." << endl;
    SeqLib::FermiAssembler fermi;
    auto contigs = fermiAssemble(brv, true);
    cout << "Writing " << contigs.size() << " contigs to " << outfile << endl;
    writeContigs(outfile, contigs);
}


void test(const std::string &filename) {
    cout << filename << endl;
    SeqLib::BamReader reader;
    reader.Open(filename);
    SeqLib::BamHeader header = reader.Header();
    cout << header.IDtoName(0) << endl;
    cout << reader << endl;
    SeqLib::BamRecord br;
    SeqLib::BamRecordVector semimapped;
    SeqLib::BamRecordVector unmapped;
    pisa::ReadFilter filter(10, 10);
    SeqLib::BamWriter writer;
    writer.SetHeader(header);
    writer.Open("semimapped.bam");
    writer.WriteHeader();

    unsigned long counter = 0;
    while (reader.GetNextRecord(br)) {
        if (++counter % 100000 == 0) cout << "read " << counter << " " << br << endl;
        switch (filter.check(br)) {
            case pisa::ReadType::HalfMappedPair:
                semimapped.push_back(br);
                break;
            case pisa::ReadType::UnmappedPair:
                unmapped.push_back(br);
                break;
            case pisa::ReadType::Split:
                semimapped.push_back(br);
                break;
            case pisa::ReadType::FilterFail:
                break;
            case pisa::ReadType::MappedPair:
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

    auto semimappedContigs = svabaAssemble(semimapped, true);
    cout << "ASSEMBLE SEMIMAPPED" << endl;
    cout << "assembled " << semimappedContigs.size() << " contigs." << endl;
    printContigs(semimappedContigs, 1000);
    writeContigs("round1.fa", semimappedContigs);

    // Build contigs from unmappable read set
    auto unmappedContigs = svabaAssemble(unmapped, true);
    cout << "ASSEMBLE UNMAPPED" << endl;
    cout << "assembled " << unmappedContigs.size() << " contigs." << endl;
    printContigs(unmappedContigs, 1000);
    writeContigs("round2.fa", unmappedContigs);

    auto contigsAsBam = toBamRecordVector(semimappedContigs);
    for (auto &ctg : toBamRecordVector(unmappedContigs)) {
        contigsAsBam.push_back(ctg);
    }

    auto contigs = fermiAssemble(contigsAsBam, true);
    cout << "ASSEMBLE SEMIMAPPED AND UNMAPPED ASSEMBLY PRODUCTS" << endl;
    cout << "assembled " << contigs.size() << " contigs." << endl;
    printContigs(contigs, 1000);
    writeContigs("round3.fa", contigs);

    SeqLib::BamRecordVector all;
    all.reserve(semimapped.size() + unmapped.size());
    for (auto &rec : semimapped) all.push_back(rec);
    semimapped.clear();
    for (auto &rec : unmapped) all.push_back(rec);
    auto allContigs = fermiAssemble(all, true);
    cout << "ASSEMBLE SEMIMAPPED AND UNMAPPED READS" << endl;
    cout << "assembled " << allContigs.size() << " contigs." << endl;
    printContigs(allContigs, 1000);
    writeContigs("round4.fa", allContigs);

    contigsAsBam = toBamRecordVector(semimappedContigs);
    for (auto &ctg : unmapped) {
        contigsAsBam.push_back(ctg);
    }

    auto mixed = fermiAssemble(contigsAsBam);
    cout << "ASSEMBLE SEMIMAPPED ASSEMBLY PRODUCTS WITH UNMAPPED READS" << endl;
    cout << "assembled " << mixed.size() << " contigs." << endl;
    printContigs(mixed, 1000);
    writeContigs("round5.fa", mixed);

    writer.Close();
    writer.BuildIndex();
}


// Handle command line inputs and launch `run`
int main(int argc, char *args[]) {
    if (argc < 3) {
        cerr << TOWER << USAGE << endl;
        exit(1);
    }
    cout << TOWER << endl;
    const std::string program = args[1];
    const std::string filename = args[2];
    std::string reference;
    if (argc > 3) reference = args[3];

    //run(4, filename, reference);
    if (program == "extract") extract(filename);
    if (program == "assemble") assemble(filename);
    if (program == "test") test(filename);
    return 0;
}
