# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/home/vunguyen/tin_sinh/simseq_sc/build/_deps/seqan3_fetch_content-src"
  "/home/vunguyen/tin_sinh/simseq_sc/build/_deps/seqan3_fetch_content-build"
  "/home/vunguyen/tin_sinh/simseq_sc/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix"
  "/home/vunguyen/tin_sinh/simseq_sc/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/tmp"
  "/home/vunguyen/tin_sinh/simseq_sc/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src/seqan3_fetch_content-populate-stamp"
  "/home/vunguyen/tin_sinh/simseq_sc/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src"
  "/home/vunguyen/tin_sinh/simseq_sc/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src/seqan3_fetch_content-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/home/vunguyen/tin_sinh/simseq_sc/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src/seqan3_fetch_content-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/home/vunguyen/tin_sinh/simseq_sc/build/_deps/seqan3_fetch_content-subbuild/seqan3_fetch_content-populate-prefix/src/seqan3_fetch_content-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
