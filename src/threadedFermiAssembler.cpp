//
// Created by Kevin Gori on 07/09/2017.
//

#include "threadedFermiAssembler.h"
#define MAG_MIN_NSR_COEF .1

namespace Pisa {

    ThreadedFermiAssembler::ThreadedFermiAssembler()  : m_seqs(0), m(0), size(0), n_seqs(0), n_utg(0), m_utgs(0)  {
        fml_opt_init(&opt);
    }

    ThreadedFermiAssembler::~ThreadedFermiAssembler() {
        ClearReads();
        ClearContigs();
    }

    // code copied and slightly modified from
    // fermi-lite/misc.c by Heng Li
    void ThreadedFermiAssembler::DirectAssemble(float kcov) {

        rld_t *e = fml_seq2fmi(&opt, n_seqs, m_seqs);
        mag_t *g = fml_fmi2mag(&opt, e);

        opt.mag_opt.min_ensr = opt.mag_opt.min_ensr > kcov * MAG_MIN_NSR_COEF? opt.mag_opt.min_ensr : (int)(kcov * MAG_MIN_NSR_COEF + .499);
        //opt.mag_opt.min_ensr = opt.mag_opt.min_ensr < opt0->max_cnt? opt.mag_opt.min_ensr : opt0->max_cnt;
        //opt.mag_opt.min_ensr = opt.mag_opt.min_ensr > opt0->min_cnt? opt.mag_opt.min_ensr : opt0->min_cnt;
        opt.mag_opt.min_insr = opt.mag_opt.min_ensr - 1;
        fml_mag_clean(&opt, g);
        m_utgs = fml_mag2utg(g, &n_utg);
    }

    void ThreadedFermiAssembler::AddRead(const SeqLib::BamRecord& r) {
        AddRead(SeqLib::UnalignedSequence(r.Qname(), r.Sequence(), r.Qualities())); // probably faster way
    }

    void ThreadedFermiAssembler::AddRead(const SeqLib::UnalignedSequence& r) {

        if (r.Seq.empty())
            return;
        if (r.Name.empty())
            return;

        // dynamically alloc the memory
        if (m <= n_seqs)
            m = m <= 0 ? 32 : (m*2); // if out of mem, double it
        m_seqs = (fseq1_t*)realloc(m_seqs, m * sizeof(fseq1_t));

        // add the name
        m_names.push_back(r.Name);

        // construct the seq
        fseq1_t *s;
        s = &m_seqs[n_seqs];
        s->seq   = strdup(r.Seq.c_str());
        s->qual = r.Qual.empty() ? NULL : strdup(r.Qual.c_str());

        s->l_seq = r.Seq.length();
        size += m_seqs[n_seqs++].l_seq;

    }

    void ThreadedFermiAssembler::AddReads(const SeqLib::UnalignedSequenceVector& v) {

        // alloc the memory
        m = n_seqs + v.size();
        m_seqs = (fseq1_t*)realloc(m_seqs, m * sizeof(fseq1_t));

        for (SeqLib::UnalignedSequenceVector::const_iterator r = v.begin(); r != v.end(); ++r) {
            m_names.push_back(r->Name);
            fseq1_t *s;

            s = &m_seqs[n_seqs];

            s->seq   = strdup(r->Seq.c_str());
            s->qual  = strdup(r->Qual.c_str());

            s->l_seq = r->Seq.length();
            size += m_seqs[n_seqs++].l_seq;
        }


    }
    void ThreadedFermiAssembler::AddReads(const SeqLib::BamRecordVector& brv) {

        // alloc the memory
        m_seqs = (fseq1_t*)realloc(m_seqs, (n_seqs + brv.size()) * sizeof(fseq1_t));

        uint64_t size = 0;
        for (SeqLib::BamRecordVector::const_iterator r = brv.begin(); r != brv.end(); ++r) {
            m_names.push_back(r->Qname());
            fseq1_t *s;

            s = &m_seqs[n_seqs];

            s->seq   = strdup(r->Sequence().c_str());
            s->qual  = strdup(r->Qualities().c_str());

            s->l_seq = r->Sequence().length();
            size += m_seqs[n_seqs++].l_seq;
        }

    }

    void ThreadedFermiAssembler::ClearContigs() {
        fml_utg_destroy(n_utg, m_utgs);
        m_utgs = 0;
        n_utg = 0;
    }

    void ThreadedFermiAssembler::ClearReads() {
        if (!m_seqs)
            return; //already cleared

        for (size_t i = 0; i < n_seqs; ++i) {
            fseq1_t * s = &m_seqs[i];
            if (s->qual)
                free(s->qual);
            s->qual = NULL;
            if (s->seq)
                free(s->seq);
            s->seq = NULL;
        }
        free(m_seqs);
        m_seqs = NULL;

    }

    void ThreadedFermiAssembler::CorrectReads() {
        fml_correct(&opt, n_seqs, m_seqs);
    }

    void ThreadedFermiAssembler::CorrectAndFilterReads() {
        fml_fltuniq(&opt, n_seqs, m_seqs);
    }

    void ThreadedFermiAssembler::PerformAssembly() {
        m_utgs = fml_assemble(&opt, n_seqs, m_seqs, &n_utg); // assemble!
    }

    std::vector<std::string> ThreadedFermiAssembler::GetContigs() const {
        std::vector<std::string> c;
        for (size_t i = 0; i < n_utg; ++i)
            c.push_back(std::string(m_utgs[i].seq));
        return c;
    }

    /*void ThreadedFermiAssembler::count() {

      // initialize BFC options
      uint64_t tot_len = 0;
      for (int i = 0; i < n_seqs; ++i)
        tot_len += m_seqs[i].l_seq; // compute total length
      int l_pre = tot_len - 8 < 20? tot_len - 8 : 20;

      //bfc_ch_t *ch = fml_count(n_seqs, m_seqs, opt.ec_k, 20, l_pre, opt.n_threads);
      //std::cerr << " ch->k " << ch->k << " ch->l_pre " << ch->l_pre << std::endl;

      // directly from fml count
      cnt_step_t cs;
      cs.n_seqs = n_seqs, cs.seqs = m_seqs, cs.k = opt.ec_k, cs.q = 20;
      cs.ch = bfc_ch_init(cs.k, l_pre);
      }*/

    SeqLib::UnalignedSequenceVector ThreadedFermiAssembler::GetSequences() const {

        SeqLib::UnalignedSequenceVector r;
        for (size_t i = 0; i < n_seqs; ++i) {
            fseq1_t * s = &m_seqs[i];
            SeqLib::UnalignedSequence read;
            if (s->seq)
                read.Seq = (std::string(s->seq));
            read.Name = m_names[i];
            r.push_back(read);
        }
        return r;
    }

}
