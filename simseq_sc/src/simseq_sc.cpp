#include <fstream>
#include <iostream>
#include <random>
#include <thread>
#include <chrono>
#include <algorithm>
#include <map>
#include <sys/stat.h>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/slice.hpp>
#include <seqan3/alphabet/views/complement.hpp>

#include "clipp.h"

using namespace clipp;

using types = seqan3::type_list<std::vector<seqan3::dna5>, std::string>;
using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
using sequence_record_type = seqan3::sequence_record<types, fields>;

// shorten ref name
// samtools only works when read name is shorter than 256 characters
// get only name of species
std::vector<std::string> strsplit(std::string s, std::string delimiter){
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;

  while ((pos_end = s.find(delimiter, pos_start)) != std::string::npos) {
    token = s.substr(pos_start, pos_end - pos_start);
      pos_start = pos_end + delim_len;
      res.push_back (token);
    }

  res.push_back (s.substr (pos_start));
  return res;
}

/*
std::string shorten_name(const std::string ref_name){
  std::vector<std::string> objs = strsplit(ref_name, "|");
  std::vector<std::string> species_paths = strsplit(objs[1], ".");
  std::string res = objs[0] + "|" + species_paths.back() + "|" +
                    objs[2] + "|" + objs[3] + "|" + objs[4];
  return res;
}
*/

void simulate_read_3prime(int& count_read, const auto& ref_seq, const std::string& ref_name,
                          std::mt19937& norm_gen, std::default_random_engine& uniform_gen,
                          auto& out1,
                          auto& out2,
                          size_t read_len = 100, double depth = 1, int min_n_read = 200, int read_count = 0,
                          double frag_len_mean = 250, double frag_len_std = 25){

    std::normal_distribution<double> d{frag_len_mean, frag_len_std};
    std::uniform_int_distribution<int> true_false(0, 1);

    int seq_len = ref_seq.size();
    int n_seq = 0;
    if (read_count > 0){
      n_seq = read_count;
    }else{
      if(n_seq == 0) n_seq = depth * seq_len;
      if(n_seq < min_n_read){
        n_seq = min_n_read;
      }
    }

    int j = 0;

    while(j < n_seq){

      int frag_len = d(norm_gen); // get random fragment length based on the fragment length distribution
      while(frag_len < 1) frag_len = d(norm_gen);
      int start_pos = 0;
      int end_pos = seq_len;

      // randomly extract the fragment if sequence is longer than fragment length
      if(seq_len > frag_len){
        //std::uniform_int_distribution<int> start(0, seq_len - frag_len);
        //std::uniform_int_distribution<int> end(seq_len - frag_len,seq_len);
        //end_pos = end(uniform_gen);
        end_pos = seq_len; // We fix the second read at the 3 primes, but the first read varies depending on fragment length
        start_pos = end_pos - frag_len;
      }
      auto fragment = ref_seq | seqan3::views::slice(start_pos, end_pos);

      auto reverse_fragment = fragment | std::views::reverse;

      bool longer_than_read = fragment.size() > read_len;

      // view for both reverse and non-reverse
      auto read1_for_fragment = longer_than_read ? (fragment | seqan3::views::slice(0, read_len)) : fragment;

      auto read2_for_fragment_temp = longer_than_read ? (reverse_fragment | seqan3::views::slice(0, read_len)) : reverse_fragment;
      auto read2_for_fragment = read2_for_fragment_temp | std::views::reverse;

      // final read
      count_read++;
      std::string id1 = "read" + std::to_string(count_read) + "/" + ref_name + ";mate1:" + std::to_string(start_pos) + ":" + std::to_string(start_pos + read_len) + ";mate2:" + std::to_string(end_pos - read_len) + ":" + std::to_string(end_pos);
      std::string id2 = "read" + std::to_string(count_read) + "/" + ref_name + ";mate1:" + std::to_string(start_pos) + ":" + std::to_string(start_pos + read_len) + ";mate2:" + std::to_string(end_pos - read_len) + ":" + std::to_string(end_pos);

      std::vector<seqan3::dna5> read1;
      std::vector<seqan3::dna5> read2;


      read1.assign(read1_for_fragment.begin(), read1_for_fragment.end());
      read2.assign(read2_for_fragment.begin(), read2_for_fragment.end());

      // write to file
      sequence_record_type record1{read1, id1};
      sequence_record_type record2{read2, id2};
      out1.push_back(record1);
      out2.push_back(record2);
      j++;

    }
}

void simulate_read_5prime(int& count_read, const auto& ref_seq, const std::string& ref_name,
                          std::mt19937& norm_gen, std::default_random_engine& uniform_gen,
                          auto& out1,
                          auto& out2,
                          size_t read_len = 100, double depth = 1, int min_n_read = 200, int read_count = 0,
                          double frag_len_mean = 250, double frag_len_std = 25){

    std::normal_distribution<double> d{frag_len_mean, frag_len_std};
    std::uniform_int_distribution<int> true_false(0, 1);

    int seq_len = ref_seq.size();
    int n_seq = 0;
    if (read_count > 0){
      n_seq = read_count;
    }else{
      if(n_seq == 0) n_seq = depth * seq_len;
      if(n_seq < min_n_read){
        n_seq = min_n_read;
      }
    }

    int j = 0;

    while(j < n_seq){

      int frag_len = d(norm_gen); // get random fragment length based on the fragment length distribution
      while(frag_len < 1) frag_len = d(norm_gen);
      int start_pos = 0;
      int end_pos = seq_len;

      // randomly extract the fragment if sequence is longer than fragment length
      if(seq_len > frag_len){
        start_pos = 0;// We fix the first read at the 5 primes, but the second read varies depending on fragment length
        end_pos = start_pos + frag_len;
      }
      auto fragment = ref_seq | seqan3::views::slice(start_pos, end_pos);

      auto reverse_fragment = fragment | std::views::reverse;

      bool longer_than_read = fragment.size() > read_len;

      // view for both reverse and non-reverse
      auto read1_for_fragment = longer_than_read ? (fragment | seqan3::views::slice(0, read_len)) : fragment;

      auto read2_for_fragment_temp = longer_than_read ? (reverse_fragment | seqan3::views::slice(0, read_len)) : reverse_fragment;
      auto read2_for_fragment = read2_for_fragment_temp | std::views::reverse;

      // final read
      count_read++;
      std::string id1 = "read" + std::to_string(count_read) + "/" + ref_name + ";mate1:" + std::to_string(start_pos) + ":" + std::to_string(start_pos + read_len) + ";mate2:" + std::to_string(end_pos - read_len) + ":" + std::to_string(end_pos);
      std::string id2 = "read" + std::to_string(count_read) + "/" + ref_name + ";mate1:" + std::to_string(start_pos) + ":" + std::to_string(start_pos + read_len) + ";mate2:" + std::to_string(end_pos - read_len) + ":" + std::to_string(end_pos);

      std::vector<seqan3::dna5> read1;
      std::vector<seqan3::dna5> read2;


      read1.assign(read1_for_fragment.begin(), read1_for_fragment.end());
      read2.assign(read2_for_fragment.begin(), read2_for_fragment.end());

      // write to file
      sequence_record_type record1{read1, id1};
      sequence_record_type record2{read2, id2};
      out1.push_back(record1);
      out2.push_back(record2);
      j++;

    }
}

void simulate_read(int& count_read, const auto& ref_seq, const std::string& ref_name,
                          std::mt19937& norm_gen, std::default_random_engine& uniform_gen,
                          auto& out1,
                          auto& out2,
                          size_t read_len = 100, double depth = 1, int min_n_read = 200, int read_count = 0,
                          double frag_len_mean = 250, double frag_len_std = 25){

    std::normal_distribution<double> d{frag_len_mean, frag_len_std};
    std::uniform_int_distribution<int> true_false(0, 1);

    int seq_len = ref_seq.size();
    int n_seq = 0;
    if (read_count > 0){
      n_seq = read_count;
    }else{
      if(n_seq == 0) n_seq = depth * seq_len;
      if(n_seq < min_n_read){
        n_seq = min_n_read;
      }
    }

    int j = 0;

    while(j < n_seq){

      int frag_len = d(norm_gen);
      while(frag_len < 1) frag_len = d(norm_gen);
      int start_pos = 0;
      int end_pos = seq_len;

      // randomly extract the fragment if sequence is longer than fragment length
      if(seq_len > frag_len){
        std::uniform_int_distribution<int> start(0, seq_len - frag_len);
        start_pos = start(uniform_gen);
        end_pos = start_pos + frag_len;
      }
      auto fragment = ref_seq | seqan3::views::slice(start_pos, end_pos);

      auto reverse_fragment = fragment | std::views::reverse;

      bool longer_than_read = fragment.size() > read_len;

      // view for both reverse and non-reverse
      auto read1_for_fragment = longer_than_read ? (fragment | seqan3::views::slice(0, read_len)) : fragment;

      auto read2_for_fragment_temp = longer_than_read ? (reverse_fragment | seqan3::views::slice(0, read_len)) : reverse_fragment;
      auto read2_for_fragment = read2_for_fragment_temp | std::views::reverse;

      // final read
      count_read++;
      std::string id1 = "read" + std::to_string(count_read) + "/" + ref_name + ";mate1:" + std::to_string(start_pos) + ":" + std::to_string(start_pos + read_len) + ";mate2:" + std::to_string(end_pos - read_len) + ":" + std::to_string(end_pos);
      std::string id2 = "read" + std::to_string(count_read) + "/" + ref_name + ";mate1:" + std::to_string(start_pos) + ":" + std::to_string(start_pos + read_len) + ";mate2:" + std::to_string(end_pos - read_len) + ":" + std::to_string(end_pos);

      std::vector<seqan3::dna5> read1;
      std::vector<seqan3::dna5> read2;


      read1.assign(read1_for_fragment.begin(), read1_for_fragment.end());
      read2.assign(read2_for_fragment.begin(), read2_for_fragment.end());

      // write to file
      sequence_record_type record1{read1, id1};
      sequence_record_type record2{read2, id2};
      out1.push_back(record1);
      out2.push_back(record2);
      j++;

    }
}

int main(int argc, char* argv[]){

  std::string input = "";
  std::string out_path = "out";
  size_t n_threads = 1;
  double depth = 1;
  int min_n_read = 200;
  double frag_len_mean = 250;
  double frag_len_std = 25;
  int read_len = 100;
  int bias_type = 0; // bias type: 0 (uniform), 3 (3 prime) and 5 (5 prime)
  std::string read_count = "";

  auto cli = (
      value("input folder, containing genome files", input),
      option("-t", "--nthreads").doc("Number of threads [default: 1]") & value("times", n_threads),
      required("-o", "--output").doc("Path to output folder") & value("outfile", out_path),
      option("-v", "--version").call([]{std::cout << "version 0.1\n\n";}).doc("show version"),
      option("-d", "--depth").doc("Generate (depth) * (length) for each reference sequence [default: 1]") & value("times", depth),
      option("-b", "--biastype").doc("Bias type of read distribution: 0 (uniform), 3 (3 prime) and 5 (5 prime) [default: 0]") & integer("times", bias_type),
      option("-m", "--min").doc("Generate at least [m] for each reference seq [default: 200]") & value("times", min_n_read),
      option("-fm", "--fragmean").doc("Mean of frangment length [default: 250]") & number("ratio", frag_len_mean),
      option("-fs", "--fragstd").doc("Std of frangment length [default: 25]") & number("ratio", frag_len_std),
      option("-r", "--readlen").doc("Read length [default: 100]") & integer("times", read_len),
      option("-rc", "--readcount").doc("Read count file:") & value("readfile",read_count)
  );

  if(!parse(argc, argv, cli)){
    std::cout << "Program: simseq_sc (simulator for sequencing reads of single-cell RNA-sequencing data)\n";
    std::cout << "Version: 0.1\n";
    std::cout << "Authors: Thinh Trac, Trung Nghia Vu, et al. \n\n";
    std::cout << make_man_page(cli, argv[0]);

    return 0;
  }

  if(std::filesystem::is_directory(input) == false){
    std::cout << "[ERROR] Input " << input << " is not a directory" << std::endl;
    return 0;
  }

  if(std::filesystem::is_directory(out_path) == false){
    if(std::filesystem::create_directory(out_path) == false){
      std::cout << "[ERROR] Can not creaete the output directory " << out_path  << std::endl;
      return 0;
    };
  } else {
    std::cout << "[WARNING] Output directory " << out_path << " already exists" << std::endl;
  }

  std::vector<std::string> input_files;

  for (const auto & entry : std::filesystem::directory_iterator(input))
    input_files.push_back(entry.path());

  // if number of files is smaller than number of thread
  if(input_files.size() < n_threads){
    n_threads = input_files.size();
  }

  // get input per thread
  std::vector<std::string> input_files_per_thread[n_threads];
  int k = 0;
  for(auto& input_file: input_files){
    int id = k % n_threads;
    input_files_per_thread[id].push_back(input_file);
    k++;
  }

  for(int i = 0 ; i < n_threads ; i++){
    std::cout << "Thread " << i << ":" << std::endl;
    for(auto& input_file: input_files_per_thread[i]){
      std::cout << input_file << std::endl;
    }
  }

  std::map<std::string, int> id_map;

  std::ifstream file(read_count);
  if (!file) {
    std::cout << "[ERROR] Read count file " << read_count << " does not exist." << std::endl;
    return 0;
  } else {
    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string id;
        int number;
        
        if (ss >> id >> number) {  // Read ID and Number
            id_map[id] = number;   // Store in the map
        } else {
            std::cout << "Invalid read count file." << std::endl;
            return 0;
        }
    }
  }

  std::vector<std::thread> readers;
  for (size_t i = 0; i < n_threads; ++i) {
    readers.emplace_back([&, i]() {

      // random generator for fragment length and read starting position
      std::random_device rd{};
      std::mt19937 norm_gen{rd()};
      std::default_random_engine uniform_gen;

      for(auto& input_file: input_files_per_thread[i]){
        seqan3::sequence_file_input fin{input_file};
        std::string output_file_1 = out_path + "/" + std::filesystem::path(input_file).stem().string() + "_sim_1.fasta.gz";
        std::string output_file_2 = out_path + "/" + std::filesystem::path(input_file).stem().string() + "_sim_2.fasta.gz";
        seqan3::sequence_file_output fout1{output_file_1};
        seqan3::sequence_file_output fout2{output_file_2};

        int count = 0;
        for(auto& rec: fin){

          auto seq = rec.sequence();
          std::string id = rec.id();

          // shorten ref name
          // std::string short_id = shorten_name(id);
          std::vector<std::string> objs = strsplit(id, " ");
          int myreadcount = (int) id_map[objs[0]];

          // nghiavtr/19Jan2025: add functions for 3 prime and 5 prime bias
          if(bias_type == 3){
            simulate_read_3prime(count, seq, id, norm_gen, uniform_gen, fout1, fout2,
                          read_len = read_len,
                          depth = depth, min_n_read = min_n_read, myreadcount,
                          frag_len_mean = frag_len_mean, frag_len_std = frag_len_std);
          }else{
            if(bias_type == 5){
              simulate_read_5prime(count, seq, id, norm_gen, uniform_gen, fout1, fout2,
                            read_len = read_len,
                            depth = depth, min_n_read = min_n_read, myreadcount,
                            frag_len_mean = frag_len_mean, frag_len_std = frag_len_std);
            }else{
              simulate_read(count, seq, id, norm_gen, uniform_gen, fout1, fout2,
                            read_len = read_len,
                            depth = depth, min_n_read = min_n_read, myreadcount,
                            frag_len_mean = frag_len_mean, frag_len_std = frag_len_std);
                }
          }
        }

        std::cout << "Finish " << std::filesystem::path(input_file).filename() << std::endl;

      }
    });
  }
  for (auto& t : readers) {
      t.join();
  }

  return 0;
}
