#include <algorithm>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include "SeqLib/BamWriter.h"
#include "SeqLib/BamReader.h"
#include "SeqLib/BWAWrapper.h"  // for realigning contigs
#include "ctpl_stl.h"  // threadpool
#include "pisaAssemble.h"
#include "pisaExtract.h"
#include "pisaTest.h"
#include "args.hxx"

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

enum Correction {
    NONE,
    BFC,
    FERMI
};


// Handle command line inputs and launch `run`
int main(int argc, char *argv[]) {
    auto start = std::chrono::steady_clock::now();

    args::ArgumentParser parser("This is a test program.", "This goes after the options.");
    args::HelpFlag help(parser, "help", USAGE, {'h', "help"});
    args::Positional<std::string> program_arg(parser, "program", "The subprogram to use");
    args::ValueFlag<int> threads_arg(parser, "threads", "The number of threads to use", {'t'});
    args::ValueFlag<std::string> filename_arg(parser, "filename", "The filename to operate on", {'f'});

    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help&)
    {
        std::cout << parser;
        return 0;
    }
    catch (args::ParseError &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    catch (args::ValidationError &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    cout << TOWER << endl;

    const std::string program = args::get(program_arg); //argv[1];
    const std::string filename = args::get(filename_arg); // argv[2];
    const int threads = args::get(threads_arg);

    pisaOptions opts;
    opts.filename = filename;
    opts.threads =  threads;
    if (program == "extract") Pisa::extract(opts);
    if (program == "assemble") Pisa::assemble(opts);
    if (program == "test") Pisa::test(opts);

    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    cout << std::chrono::duration <double, std::milli> (diff).count() << " ms" << endl;
    return 0;
}
