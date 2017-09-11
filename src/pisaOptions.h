//
// Created by Kevin Gori on 11/09/2017.
//

/*
 * A struct to hold various commandline options to be passed to pisa subprograms
 */

#ifndef PISA_PISAOPTIONS_H
#define PISA_PISAOPTIONS_H
#include <string>

struct pisaOptions {
    std::string filename;
    int threads = 1;
};

#endif //PISA_PISAOPTIONS_H
