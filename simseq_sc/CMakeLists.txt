cmake_minimum_required (VERSION 3.14)
project (bam2eqv LANGUAGES CXX VERSION 1.0.0)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSEQAN_HAS_ZLIB=1 -std=c++11 -O3 -DNDEBUG")

set (seqan3_git_tag "3.3.0") # adapt as needed, e.g. "3.2.0" or "master"

message (STATUS "Fetching SeqAn3 ${seqan3_git_tag}:")
include (FetchContent)
FetchContent_Declare (
      seqan3_fetch_content
      GIT_REPOSITORY "https://github.com/seqan/seqan3.git"
      GIT_TAG "${seqan3_git_tag}"
)
# Download and make SeqAn3 available.
FetchContent_MakeAvailable (seqan3_fetch_content)

add_executable (simseq_sc src/simseq_sc.cpp)

# Set up everything needed to use SeqAn3 with my_app:
target_link_libraries (simseq_sc PUBLIC seqan3::seqan3)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(simseq_sc PUBLIC OpenMP::OpenMP_CXX)
endif()

target_include_directories(simseq_sc
  PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)
