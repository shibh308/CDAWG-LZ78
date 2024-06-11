#ifndef CDAWG_LZ78_COMPRESSION_MEASURE_HPP
#define CDAWG_LZ78_COMPRESSION_MEASURE_HPP


#include "sdsl_suffixtree.hpp"
#include "waveletmatrix.hpp"

int compute_lz77_length(NormalSuffixTree& st, std::vector<int>& sa){
  // SA range min
  int n = st.text.size();
  assert(n < 1 << 28);
  WaveletMatrix<int, 28> wm(sa);

  int i = 0;
  int k = 0;
  while(++k, i < n){
    int node = st.leaves[i];
    while(true){
      if(node == 0){
        ++i;
        break;
      }
      auto [l, r] = st.sa_ranges[node];
      int cnt = wm.count_lower(l, r, i);
      if(cnt > 0){
        i += st.depths[node];
        break;
      }
      else{
        node = st.parents[node];
      }
    }
  }
  return k;
}

int compute_lz78_length(NormalSuffixTree& st){
  // SA range min
  int n = st.text.size();
  std::vector<int> lz78_depth(st.parents.size(), 0);
  int i = 0;
  int k = 0;
  while(++k, i < n){
    int node = st.leaves[i];
    while(true){
      if(lz78_depth[node] != 0){
        ++lz78_depth[node];
        i += lz78_depth[node];
        break;
      }
      else{
        if(st.depths[st.parents[node]] == lz78_depth[st.parents[node]]){
          lz78_depth[node] = st.depths[st.parents[node]] + 1;
          i += lz78_depth[node];
          break;
        }
        node = st.parents[node];
      }
    }
  }
  return k;
}

int compute_RLBWT_length(std::string& text){
  sdsl::csa_wt csa;
  construct_im(csa, text, 1);
  int r = 0;
  int p = -1;
  int n = csa.bwt.size();
  for(int i = 0; i < n; ++i){
    int b = csa.bwt[i];
    if(p != b){
      ++r;
      p = b;
    }
  }
  return r;
}

#endif //CDAWG_LZ78_COMPRESSION_MEASURE_HPP
