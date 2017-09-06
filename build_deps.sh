#!/usr/bin/env bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo build_deps.sh: $DIR
echo build_deps.sh: CPPFLAGS=${CPPFLAGS}
cd $DIR/deps/SeqLib
./autogen.sh
./configure CPPFLAGS=${CPPFLAGS}
make CPPFLAGS=${CPPFLAGS}