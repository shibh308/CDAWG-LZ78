#ifndef CDAWG_LZ78_COMPRESSION_MEASURE_HPP
#define CDAWG_LZ78_COMPRESSION_MEASURE_HPP


#include "sdsl_suffixtree.hpp"
#include "waveletmatrix.hpp"
#include "lz78.hpp"

int compute_lz77_length(NormalSuffixTree& st, std::vector<int>& sa){
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

struct LZ78SlinkTrie{
  std::vector<std::unordered_map<unsigned char, int>> lz_trie;
  std::vector<int> slinks;
  LZ78SlinkTrie(std::vector<std::pair<int, unsigned char>>& lz78_phrases) : lz_trie(lz78_phrases.size() + 1), slinks(lz78_phrases.size() + 1){
    for(int i = 0; i < lz78_phrases.size(); ++i){
      lz_trie[lz78_phrases[i].first + 1][lz78_phrases[i].second] = i + 1;
    }
    std::queue<int> que;
    que.emplace(0);
    while(!que.empty()){
      auto x = que.front();
      que.pop();
      for(auto [c, y] : lz_trie[x]){
        if(x == 0){
          slinks[y] = 0;
        }
        else{
          int sl = slinks[x];
          std::stack<int> stack;
          while(sl != 0 && !lz_trie[sl].contains(c)){
            stack.emplace(sl);
            sl = slinks[sl];
          }
          if(sl == 0 && !lz_trie[sl].contains(c)){
            lz_trie[sl][c] = lz_trie.size();
            lz_trie.emplace_back();
            slinks.emplace_back(0);
          }
          while(!stack.empty()){
            int sl_next = stack.top();
            stack.pop();
            lz_trie[sl_next][c] = lz_trie.size();
            lz_trie.emplace_back();
            slinks.emplace_back(lz_trie[sl].at(c));
            sl = sl_next;
          }
          slinks[y] = lz_trie[slinks[x]].at(c);
        }
        que.emplace(y);
      }
    }
  }
  std::vector<int> reinsert_text(std::string_view text){
    std::vector<int> phrases;
    int i = 0;
    while(i != text.size()){
      int node = 0;
      for(;i != text.size() && lz_trie[node].contains(text[i]); ++i){
        node = lz_trie[node][text[i]];
      }
      phrases.emplace_back(node);
    }
    return phrases;
  }
};
std::pair<int,int> compute_lz78_slink_trie_size(NormalSuffixTree& st){
  auto random_access_func = [&](int i){ return st.text[i]; };
  auto sa_range_func = [&](int start_pos, int target_depth){
    return st.sa_range(start_pos, target_depth);
  };
  auto lz78_phrases = compute_lz78(st.text.size(), 0, st.text.size(), random_access_func, sa_range_func).first;
  auto slink_trie = LZ78SlinkTrie(lz78_phrases);
  return std::make_pair(slink_trie.lz_trie.size(), slink_trie.reinsert_text(st.text).size());
}

int compute_lz78_length(NormalSuffixTree& st){
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

std::vector<int> compute_lz78_length_array(NormalSuffixTree& st, int start_pos){
  int n = st.text.size();
  std::vector<int> length;
  length.reserve(n - start_pos);
  std::vector<int> lz78_depth(st.parents.size(), 0);
  int i = start_pos;
  while(i < n){
    int node = st.leaves[i];
    while(true){
      if(lz78_depth[node] != 0){
        ++lz78_depth[node];
        i += lz78_depth[node];
        length.emplace_back(lz78_depth[node]);
        break;
      }else{
        if(st.depths[st.parents[node]] == lz78_depth[st.parents[node]]){
          lz78_depth[node] = st.depths[st.parents[node]] + 1;
          i += lz78_depth[node];
          length.emplace_back(lz78_depth[node]);
          break;
        }
        node = st.parents[node];
      }
    }
  }
  return length;
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
