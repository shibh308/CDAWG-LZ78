#ifndef CDAWG_LZ78_SDSL_SUFFIXTREE_HPP
#define CDAWG_LZ78_SDSL_SUFFIXTREE_HPP

#include "sdsl/suffix_trees.hpp"
#include "sdsl/construct.hpp"

struct SuccinctSuffixTree{
  int n;
  sdsl::cst_sada<> st;
  SuccinctSuffixTree(std::string& s) : n(s.size()){
    construct_im(st, s, 1);
  };
  std::size_t memory_usage(){
    return size_in_bytes(st);
  }
  unsigned char random_access(int i){
    return st.csa.text[i];
  }
  std::pair<int,int> sa_range(int start_pos, int target_depth){
    auto isa_idx = st.csa.isa[start_pos];
    auto node = st.select_leaf(isa_idx + 1);
    while(node != st.root()){
      auto parent = st.parent(node);
      if(st.depth(parent) < target_depth){
        break;
      }
      node = parent;
    }
    if(node == st.root()){
      return {0, n};
    }
    int left = st.lb(node);
    int right = st.rb(node);
    return {left - 1, right};
  }
};

struct NormalSuffixTree{
  int n;
  std::vector<int> depths;
  std::vector<int> parents;
  std::vector<std::pair<int,int>> sa_ranges;
  std::vector<int> leaves;
  std::string text;
  NormalSuffixTree(std::string& s) : n(s.size()), depths(2 * n), parents(2 * n, -1), sa_ranges(2 * n), text(s), leaves(n + 1){
    sdsl::cst_sct3<> st;
    construct_im(st, s, 1);
    std::vector<int> stack;
    std::vector<std::tuple<int, unsigned char, int>> children;
    int i = 0;
    for(auto it = st.begin(); it != st.end(); ++it){
      if(it.visit() == 1){
        auto v = *it;
        int depth = st.depth(v);
        sa_ranges[i] = { st.lb(v) - 1, st.rb(v) };
        depths[i] = depth;
        while(!stack.empty() && sa_ranges[stack.back()].second < sa_ranges[i].second)stack.pop_back();
        if(!stack.empty()){
          parents[i] = stack.back();
        }
        stack.emplace_back(i);
        if(st.is_leaf(v)){
          if(depth != 1){
            leaves[n - (depth - 1)] = i;
          }
        }
        ++i;
      }
    }
    depths.resize(i);
    parents.resize(i);
    sa_ranges.resize(i);
    depths.shrink_to_fit();
    parents.shrink_to_fit();
    sa_ranges.shrink_to_fit();
    leaves.shrink_to_fit();
  };
  std::size_t memory_usage(){
    return sizeof(unsigned char) * text.size()
         + sizeof(int) * depths.size()
         + sizeof(int) * parents.size()
         + sizeof(int) * 2 * sa_ranges.size()
         + sizeof(int) * leaves.size();
  }
  unsigned char random_access(int i){
    return text[i];
  }
  std::pair<int,int> sa_range(int start_pos, int target_depth){
    int node = leaves[start_pos];
    while(node){
      int parent = parents[node];
      if(depths[parent] < target_depth){
        break;
      }
      node = parent;
    }
    if(node == 0){
      return {0, n};
    }
    return sa_ranges[node];
  }
};

#endif //CDAWG_LZ78_SDSL_SUFFIXTREE_HPP
