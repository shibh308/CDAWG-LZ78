#ifndef CDAWG_LZ78_UTILS_HPP
#define CDAWG_LZ78_UTILS_HPP

#include <string>
#include <bitset>
#include "suffix_array.hpp"
#include "dawg.hpp"
#include "compression_measure.hpp"

std::string load_file(const std::string& data_path, std::size_t length_limit){
  std::clog << "loading: " << data_path << std::endl;
  std::ifstream file(data_path);
  assert(file.is_open());
  std::string text;
  text.resize(length_limit);
  file.read(&text[0], length_limit);
  std::size_t read_bytes = file.gcount();
  if(read_bytes != length_limit){
    text.resize(read_bytes);
  }
  std::clog << "Loading file \"" << data_path << "\" is finished." << std::endl;
  return text;
}

void compute_compression_measures(std::string filename, int length){
  auto text = load_file("./data/" + filename, length);
  std::bitset<256> bs;
  int n = text.size();
  std::clog << "computing sigma..." << std::endl;
  for(int i = 0; i < n; ++i){
    bs.set(static_cast<unsigned char>(text[i]));
  }
  int sigma = bs.count();
  std::clog << "computing e..." << std::endl;
  int e = [&](){
    int e = 0;
    CDAWGBase cdawg(text);
    for(auto& node : cdawg.nodes){
      e += node.ch.size();
    }
    return e;
  }();
  std::clog << "computing r..." << std::endl;
  int r = compute_RLBWT_length(text);
  std::clog << "constructing SA..." << std::endl;
  text += '\1';
  SALCP sa(text);
  std::clog << "constructing ST..." << std::endl;
  NormalSuffixTree st(text);
  std::clog << "computing lz77..." << std::endl;
  int z77 = compute_lz77_length(st, sa.sa);
  std::clog << "computing lz78..." << std::endl;
  int z78 = compute_lz78_length(st);
  std::cout << "filename,n,sigma,e,r,lz77,lz78\n";
  std::cout << filename << "," << n << "," << sigma << "," << e << "," << r << "," << z77 << "," << z78 << std::endl;
}


#endif //CDAWG_LZ78_UTILS_HPP
