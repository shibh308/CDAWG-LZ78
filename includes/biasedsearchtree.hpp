
#ifndef CDAWG_LZ78_BIASEDSEARCHTREE_HPP
#define CDAWG_LZ78_BIASEDSEARCHTREE_HPP

#endif //CDAWG_LZ78_BIASEDSEARCHTREE_HPP

#include <cassert>
#include <functional>
#include <vector>
#include <utility>

#define BST_COMPARE_EQUAL 0
#define BST_COMPARE_GREATER 1
#define BST_COMPARE_LESS -1

template<typename weight_type, typename value_type>
struct BiasedSearchTree{
  struct Node{
    Node(value_type value) : value(value), lch(-1), rch(-1){}
    value_type value;
    int lch = -1, rch = -1;
  };
  std::vector<Node> nodes;
  // TODO: linear time construction
  explicit BiasedSearchTree(const std::vector<std::pair<weight_type, value_type>>& vec) : nodes(vec.size()){
    recursive_build(vec, 0, vec.size());
  }
  int recursive_build(const std::vector<std::pair<weight_type, value_type>>& vec, int l, int r){
    if(l == r)return -1;
    weight_type s = 0;
    for(int i = l; i < r; ++i){
      s += vec[i].first;
    }
    int m = l;
    weight_type t = vec[l].first;
    while(t * 2 < s){
      ++m;
      t += vec[m].first;
    }
    int node_id = nodes.size();
    nodes.emplace_back(vec[m].second);
    nodes[node_id].lch = recursive_build(vec, l, m);
    nodes[node_id].rch = recursive_build(vec, m + 1, r);
    return node_id;
  }
  int search(std::function<int(value_type&)> compare_func, int node_id = 0){
    assert(node_id != -1);
    auto& node = nodes[node_id];
    int res = compare_func(node.value);
    if(res == BST_COMPARE_GREATER){
      return search(node.rch);
    }
    else if(res == BST_COMPARE_LESS){
      return search(node.lch);
    }
    return node_id;
  }
};