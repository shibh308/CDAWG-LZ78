#ifndef CDAWG_LZ78_DAWG_HPP
#define CDAWG_LZ78_DAWG_HPP

#include <iostream>
#include <queue>
#include "map.hpp"
#include "vector.hpp"
#include "biasedsearchtree.hpp"


std::uint32_t msb32(std::uint32_t v){
  v |= (v >> 1);
  v |= (v >> 2);
  v |= (v >> 4);
  v |= (v >> 8);
  v |= (v >> 16);
  return v ^ (v >> 1);
}

struct DAWGBase{
  struct Node{
    DynamicHashMap<unsigned char, int> ch;
    int slink, len;
    explicit Node(int len) : slink(-1), len(len){}
  };

  std::vector<Node> nodes;
  int final_node = 0;

  explicit DAWGBase(std::string_view text){
    nodes.emplace_back(0);
    for(int i = 0; i < text.size(); ++i){
      add_node(i, text[i]);
    }
    nodes.shrink_to_fit();
    std::clog << "Construction of DAWGBase is finished." << std::endl;
  }

  void add_node(int i, unsigned char c){
    int new_node = nodes.size();
    int target_node = (nodes.size() == 1 ? 0 : final_node);
    final_node = new_node;
    nodes.emplace_back(i + 1);

    for(; target_node != -1 &&
          !nodes[target_node].ch.find(c).has_value(); target_node = nodes[target_node].slink){
      nodes[target_node].ch.add(c, new_node);
    }
    if(target_node == -1){
      nodes[new_node].slink = 0;
    }else{
      int sp_node = nodes[target_node].ch.find(c).value();
      if(nodes[target_node].len + 1 == nodes[sp_node].len){
        nodes[new_node].slink = sp_node;
      }else{
        int clone_node = nodes.size();
        nodes.emplace_back(nodes[target_node].len + 1);
        nodes[clone_node].ch = nodes[sp_node].ch;
        nodes[clone_node].slink = nodes[sp_node].slink;
        for(; target_node != -1 && nodes[target_node].ch.find(c).has_value() &&
              nodes[target_node].ch.find(c).value() == sp_node; target_node = nodes[target_node].slink){
          nodes[target_node].ch.add(c, clone_node);
        }
        nodes[sp_node].slink = nodes[new_node].slink = clone_node;
      }
    }
  }
};

// TODO: construct CDAWG from SA/LCP/BWT
struct CDAWGBase{

  struct Node{
    Vector<std::tuple<unsigned char, int, int>, std::uint16_t> ch;
    int len;
    explicit Node(int len) : len(len){}
  };

  std::vector<Node> nodes;

  explicit CDAWGBase(std::string_view text){
    DAWGBase dawg = DAWGBase(text);
    assert(dawg.nodes[0].ch.n != 1);
    std::vector<int> skip(dawg.nodes.size(), -1);
    std::vector<int> skip_length(dawg.nodes.size(), 0);
    std::iota(skip.begin(), skip.end(), 0);

    // topological sort
    std::vector<int> in_deg(dawg.nodes.size());
    for(auto& node : dawg.nodes){
      for(auto [label, to] : node.ch.items()){
        ++in_deg[to];
      }
    }
    std::vector<int> tps_order;
    tps_order.reserve(dawg.nodes.size());
    std::queue<int> que;
    que.emplace(0);
    while(!que.empty()){
      int x = que.front();
      tps_order.emplace_back(x);
      que.pop();
      for(auto [k, y] : dawg.nodes[x].ch.items()){
        if(--in_deg[y] == 0){
          que.emplace(y);
        }
      }
    }

    for(int i = tps_order.size() - 1; i >= 0; --i){
      int x = tps_order[i];
      auto& node = dawg.nodes[x];
      if(node.ch.size() == 1){
        int to = node.ch.items().front().second;
        skip[x] = skip[to];
        skip_length[x] = skip_length[to] + 1;
      }
    }
    std::vector<int> node_id_mapping(tps_order.size());
    for(auto x : tps_order){
      auto &node = dawg.nodes[x];
      if(node.ch.size() != 1){
        Node cdawg_node(node.len);
        std::vector<std::tuple<unsigned char, int, int>> children;
        for(auto [label, to]: node.ch.items()){
          children.emplace_back(label, skip[to], skip_length[to] + 1);
        }
        cdawg_node.ch = Vector<std::tuple<unsigned char, int, int>, std::uint16_t>(children);
        std::sort(cdawg_node.ch.begin(), cdawg_node.ch.end());
        node_id_mapping[x] = nodes.size();
        nodes.emplace_back(cdawg_node);
      }
    }
    for(auto& node : nodes){
      for(auto& [label, to, length] : node.ch){
        to = node_id_mapping[to];
      }
    }
    nodes.shrink_to_fit();
    std::clog << "Construction of CDAWGBase is finished." << std::endl;
  }
};

struct SimpleCDAWG{
  int n;
  std::vector<std::tuple<int,int,int>> parents; // parent, length, num_pred_edges
  std::vector<int> edge_offset;
  std::vector<int> node_length;
  std::vector<int> num_path_to_sink;
  std::vector<std::pair<int,int>> root_coming_edges; // from, length
  std::vector<unsigned char> root_coming_edge_labels;
  int num_iter;

  explicit SimpleCDAWG(const CDAWGBase& cdawg) :
    n(cdawg.nodes.size()),
    parents(n, std::make_tuple(-1, -1, -1)),
    edge_offset(n + 1),
    node_length(n),
    num_path_to_sink(n)
  {
    num_path_to_sink[n - 1] = 1;
    for(int x = n - 1; x >= 0; --x){
      for(auto [label, to, length] : cdawg.nodes[x].ch){
        num_path_to_sink[x] += num_path_to_sink[to];
      }
    }
    std::vector<int> in_degs(n + 1);
    for(int x = 0; x < n; ++x){
      for(auto [label, to, length] : cdawg.nodes[x].ch){
        ++in_degs[to + 1];
      }
    }
    for(int i = 0; i < n; ++i){
      in_degs[i + 1] += in_degs[i];
    }
    std::copy(in_degs.begin(), in_degs.end(), edge_offset.begin());

    parents.resize(in_degs[n]);
    for(int x = 0; x < n; ++x){
      node_length[x] = cdawg.nodes[x].len;
      int num_pred_edges = 0;
      for(auto [label, to, length]: cdawg.nodes[x].ch){
        parents[in_degs[to]++] = {x, length, num_pred_edges};
        num_pred_edges += num_path_to_sink[to];
      }
    }
    for(int x = 0; x < n; ++x){
      std::sort(std::next(parents.begin(), edge_offset[x]), std::next(parents.begin(), edge_offset[x + 1]),
                [&](auto x, auto y){
                  auto [x_node, x_len, x_pred] = x;
                  auto [y_node, y_len, y_pred] = y;
                  return node_length[x_node] + x_len < node_length[y_node] + y_len;
                });
    }
    std::vector<std::tuple<int,int,int>> to_and_labels;
    for(auto [label, to, length]: cdawg.nodes[0].ch){
      to_and_labels.emplace_back(to, length, label);
    }
    std::sort(to_and_labels.begin(), to_and_labels.end());
    // to de sort suru no ga mazui
    for(auto [to, length, label] : to_and_labels){
      root_coming_edges.emplace_back(to, length);
      root_coming_edge_labels.emplace_back(label);
    }
  }
  /// find the range in the SA corresponding to T[start_pos, target_depth]
  std::pair<std::pair<int,int>, unsigned char> climb_linear(int start_pos, int target_depth){
    num_iter = 0;
    int x = n - 1;
    int remain_length = node_length[n - 1] - start_pos;
    int left = 0;
    int stop_x = 0;
    unsigned char first_label = 0;
    while(x){
      ++num_iter;
      auto it = std::lower_bound(std::next(parents.begin(), edge_offset[x]), std::next(parents.begin(), edge_offset[x + 1]), remain_length, [&](auto x, int remain_length){
        auto [x_node, x_len, x_pred] = x;
        return remain_length - x_len > node_length[x_node];
      });
      std::pair<int,int> next = {std::get<0>(*it), remain_length - std::get<1>(*it)};
      if(next.second < target_depth){
        if(stop_x == 0){
          stop_x = x;
        }
        left += std::get<2>(*it);
      }
      if(next.first == 0){
        int root_edge_index = std::distance(root_coming_edges.begin(), std::lower_bound(root_coming_edges.begin(), root_coming_edges.end(), std::make_pair(x, std::get<1>(*it))));
        first_label = root_coming_edge_labels[root_edge_index];
      }
      std::tie(x, remain_length) = next;
    }
    return {std::make_pair(left, left + num_path_to_sink[stop_x]), first_label};
  }
  size_t get_memory_usage() const{
    size_t sum = 0;
    sum += parents.size() * sizeof(parents[0]);
    sum += edge_offset.size() * sizeof(edge_offset[0]);
    sum += node_length.size() * sizeof(node_length[0]);
    sum += num_path_to_sink.size() * sizeof(num_path_to_sink[0]);
    return sum;
  }
};

struct HeavyPathDecomposedCDAWG{
  int n;
  std::vector<std::tuple<int,int,int>> heavy_parent, light_parent; // parent, length, num_pred_edges
  std::vector<std::vector<std::tuple<int,int,int>>> heavy_paths; // (node, length to bottom, num_pred_edges_sum)
  std::vector<std::pair<int,int>> heavy_paths_inv; // (heavy path id, node idx in heavy path)
  std::vector<int> light_parent_offset;
  std::vector<int> node_length;
  std::vector<int> num_path_to_sink, num_path_to_source;
  int num_iter;

  explicit HeavyPathDecomposedCDAWG(const CDAWGBase& cdawg) :
    n(cdawg.nodes.size()),
    heavy_parent(n, std::make_tuple(-1, -1, -1)),
    heavy_paths_inv(n, std::make_pair(-1, -1)),
    light_parent_offset(n + 1),
    node_length(n),
    num_path_to_source(n),
    num_path_to_sink(n)
  {
    num_path_to_source[0] = 1;
    for(int x = 0; x < n; ++x){
      for(auto [label, to, length] : cdawg.nodes[x].ch){
        num_path_to_source[to] += num_path_to_source[x];
      }
    }
    num_path_to_sink[n - 1] = 1;
    for(int x = n - 1; x >= 0; --x){
      for(auto [label, to, length] : cdawg.nodes[x].ch){
        num_path_to_sink[x] += num_path_to_sink[to];
      }
    }
    assert(num_path_to_sink[0] == num_path_to_source[n - 1]);

    std::vector<int> light_in_deg(n + 1);
    // Symmetric Centroid Decomposition
    for(int x = 0; x < n; ++x){
      node_length[x] = cdawg.nodes[x].len;
      std::pair<int,int> flag(msb32(num_path_to_source[x]), msb32(num_path_to_sink[x]));
      for(auto [label, to, length] : cdawg.nodes[x].ch){
        std::pair<int,int> identifier(msb32(num_path_to_source[to]), msb32(num_path_to_sink[to]));
        if(flag == identifier){
          heavy_parent[to] = {x, length, -1};
        }
        else{
          ++light_in_deg[to + 1];
        }
      }
    }
    for(int i = 0; i < n; ++i){
      light_in_deg[i + 1] += light_in_deg[i];
    }
    std::copy(light_in_deg.begin(), light_in_deg.end(), light_parent_offset.begin());
    light_parent.resize(light_in_deg[n]);
    for(int x = 0; x < n; ++x){
      node_length[x] = cdawg.nodes[x].len;
      int num_pred_edges = 0;
      for(auto [label, to, length]: cdawg.nodes[x].ch){
        if(std::get<0>(heavy_parent[to]) != x){
          light_parent[light_in_deg[to]] = {x, length, num_pred_edges};
          ++light_in_deg[to];
        }
        else{
          std::get<2>(heavy_parent[to]) = num_pred_edges;
        }
        num_pred_edges += num_path_to_sink[to];
      }
    }

    for(int x = n - 1; x >= 0; --x){
      if(heavy_paths_inv[x].first == -1){
        std::vector<std::tuple<int,int,int>> path;
        heavy_paths_inv[x] = {heavy_paths.size(), path.size()};
        path.emplace_back(x, 0, 0);
        int length = 0;
        int pred_sum = 0;
        int y = x;
        while(true){
          auto [nex, edge_len, num_pred_edges] = heavy_parent[y];
          if(nex == -1)break;
          y = nex;
          length += edge_len;
          pred_sum += num_pred_edges;
          heavy_paths_inv[y] = {heavy_paths.size(), path.size()};
          path.emplace_back(y, length, pred_sum);
        }
        heavy_paths.emplace_back(path);
      }
    }

    for(int x = 0; x < n; ++x){
      std::sort(std::next(light_parent.begin(), light_parent_offset[x]), std::next(light_parent.begin(), light_parent_offset[x + 1]),
                [&](auto x, auto y){
                  auto [x_node, x_len, x_pred] = x;
                  auto [y_node, y_len, y_pred] = y;
                  return node_length[x_node] + x_len < node_length[y_node] + y_len;
                });
    }
  }
  std::vector<std::tuple<int,int,int>> get_parents(int x) const{
    std::vector<std::tuple<int,int,int>> parents;
    if(std::get<0>(heavy_parent[x]) != -1){
      parents.emplace_back(heavy_parent[x]);
    }
    for(int i = light_parent_offset[x]; i < light_parent_offset[x + 1]; ++i){
      parents.emplace_back(light_parent[i]);
    }
    return parents;
  }
  /// find the range in the SA corresponding to T[start_pos, target_depth]
  std::pair<int,int> climb_linear(int start_pos, int target_depth){
    num_iter = 0;
    int x = n - 1;
    int remain_length = node_length[n - 1] - start_pos;
    int left = 0;
    int stop_x = 0;
    while(x){
      ++num_iter;
      auto parents = get_parents(x);
      std::sort(parents.begin(), parents.end(), [&](auto x, auto y){
        auto [x_node, x_len, x_pred] = x;
        auto [y_node, y_len, y_pred] = y;
        return node_length[x_node] + x_len < node_length[y_node] + y_len;
      });
      auto it = std::lower_bound(parents.begin(), parents.end(), remain_length, [&](auto x, int remain_length){
        auto [x_node, x_len, x_pred] = x;
        return remain_length - x_len > node_length[x_node];
      });
      std::pair<int,int> next = {std::get<0>(*it), remain_length - std::get<1>(*it)};
      assert(it != parents.end());
      if(next.second < target_depth){
        if(stop_x == 0){
          stop_x = x;
        }
        left += std::get<2>(*it);
      }
      std::tie(x, remain_length) = next;
    }
    return {left, left + num_path_to_sink[stop_x]};
  }
  /// find the range in the SA corresponding to T[start_pos, target_depth]
  std::pair<int,int> climb_log2(int start_pos, int target_depth){
    num_iter = 0;
    int x = n - 1;
    int remain_length = node_length[n - 1] - start_pos;
    int left = 0;
    int stop_x = 0;
    while(x){
      ++num_iter;
      auto [hp_idx, idx_in_hp] = heavy_paths_inv[x];
      int orig_length = std::get<1>(heavy_paths[hp_idx][idx_in_hp]);
      int orig_pred = std::get<2>(heavy_paths[hp_idx][idx_in_hp]);
      int next_hp_idx = [&](){
        int top = heavy_paths[hp_idx].size();
        int bottom = idx_in_hp;
        while(top - bottom > 1){
          int check_idx = (top + bottom) / 2;
          auto [check_node, length, _] = heavy_paths[hp_idx][check_idx];
          int skip_length = length - orig_length;
          int check_remain_length = remain_length - skip_length;
          int node_length_r = node_length[check_node];
          int node_length_l = node_length[check_node] - num_path_to_source[check_idx];
          bool contains = (target_depth <= check_remain_length) && node_length_l < check_remain_length &&
                          check_remain_length <= node_length_r;
          (contains ? bottom : top) = check_idx;
        }
        return bottom;
      }();
//      std::clog << heavy_paths[hp_idx].size() << " " << idx_in_hp << " _> " << next_hp_idx << std::endl;
      auto [skip_node, skip_length, num_pred_edges] = heavy_paths[hp_idx][next_hp_idx];
      remain_length -= (skip_length - orig_length);
      if(remain_length < 0){
        left += (num_pred_edges - orig_pred);
      }
      x = skip_node;
      if(!x){
        break;
      }
      auto parents = get_parents(x);
      std::sort(parents.begin(), parents.end(), [&](auto x, auto y){
        auto [x_node, x_len, x_pred] = x;
        auto [y_node, y_len, y_pred] = y;
        return node_length[x_node] + x_len < node_length[y_node] + y_len;
      });
      auto it = std::lower_bound(parents.begin(), parents.end(), remain_length, [&](auto x, int remain_length){
        auto [x_node, x_len, x_pred] = x;
        return remain_length - x_len > node_length[x_node];
      });
      std::pair<int, int> next = {std::get<0>(*it), remain_length - std::get<1>(*it)};
      assert(it != parents.end());
      if(next.second < target_depth){
        if(stop_x == 0){
          stop_x = x;
        }
        left += std::get<2>(*it);
      }
      std::tie(x, remain_length) = next;
    }
    return {left, left + num_path_to_sink[stop_x]};
  }
};

#endif //CDAWG_LZ78_DAWG_HPP
