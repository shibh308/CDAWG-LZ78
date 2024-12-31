#ifndef CDAWG_LZ78_BWT_INDEX_HPP
#define CDAWG_LZ78_BWT_INDEX_HPP

#include "sdsl_suffixtree.hpp"


int floor(int numerator, int denominator) {
  assert(denominator != 0);
  int result = numerator / denominator;
  if ((numerator % denominator != 0) && ((numerator < 0) != (denominator < 0))) {
    result -= 1;
  }
  return result;
}

int ceil(int numerator, int denominator) {
  assert(denominator != 0);
  int result = numerator / denominator;
  if ((numerator % denominator != 0) && ((numerator > 0) == (denominator > 0))) {
    result += 1;
  }
  return result;
}

class TextIndex{
  int root_block_size;
  std::vector<sdsl::bit_vector> block_types;
  std::vector<sdsl::rank_support_v<>> block_ranks;
  std::vector<std::vector<std::pair<int, int>>> leaf_links;
  std::vector<uint8_t> bottom_nodes;
public:
  TextIndex(int root_block_size, std::vector<sdsl::bit_vector> block_types_, std::vector<std::vector<std::pair<int, int>>>&& leaf_links, std::vector<uint8_t>&& bottom_nodes) : root_block_size(root_block_size), block_types(block_types_), leaf_links(std::move(leaf_links)), bottom_nodes(std::move(bottom_nodes)){
    for(int i = 0; i < block_types.size(); ++i){
      block_ranks.emplace_back(&block_types[i]);
    }
  }
  std::size_t memory_usage(){
    std::size_t sum = 0;
    for(auto& v : block_types){
      sum += (v.size() + 7) / 8;
    }
    sum += bottom_nodes.size() * sizeof(uint8_t);
    for(auto& v : leaf_links){
      sum += v.size() * 8 + 8;
    }
    return sum;
  }
  uint8_t operator[](int idx) const{
    // T[i] へのrandom accessのテスト
    // 必要なのは block_types, block_ranks, leaf_links, bottom_nodes
    int block_idx = idx / root_block_size;
    int block_ofs = idx % root_block_size;
    int block_size = root_block_size;
    int d = 0;
    while(block_size > 1){
      bool internal = block_types[d][block_idx];
      if(internal){
        int internal_rank = block_ranks[d](block_idx);
        int next_idx = internal_rank * 2 + (block_ofs >= block_size / 2);
        block_idx = next_idx;
        if(block_ofs >= block_size / 2){
          block_ofs -= block_size / 2;
        }
        ++d;
        block_size >>= 1;
      }else{
        int leaf_rank = block_idx - block_ranks[d](block_idx);
        block_idx = leaf_links[d][leaf_rank].first;
        block_ofs += leaf_links[d][leaf_rank].second;
        if(block_ofs >= block_size){
          ++block_idx;
          block_ofs -= block_size;
        }
      }
    }
    return bottom_nodes[block_idx];
  }
};

class ISAIndex{
  int root_block_size;
  std::vector<sdsl::bit_vector> block_types;
  std::vector<sdsl::rank_support_v<>> block_ranks;
  std::vector<std::vector<int>> internal_ofs;
  std::vector<std::vector<std::pair<int, int>>> leaf_links;
  std::vector<std::vector<std::pair<int, int>>> leaf_link_ofs;
  std::vector<int> top_block_isas;
public:
  ISAIndex(int root_block_size, std::vector<sdsl::bit_vector>&& block_types_, std::vector<std::vector<int>>&& internal_ofs, std::vector<std::vector<std::pair<int, int>>>&& leaf_links, std::vector<std::vector<std::pair<int, int>>>&& leaf_link_ofs, std::vector<int>&& top_block_isas) : root_block_size(root_block_size), block_types(block_types_), internal_ofs(std::move(internal_ofs)), leaf_links(std::move(leaf_links)), leaf_link_ofs(std::move(leaf_link_ofs)), top_block_isas(std::move(top_block_isas)){
    for(int i = 0; i < block_types.size(); ++i){
      block_ranks.emplace_back(&block_types[i]);
    }
  }
  std::size_t memory_usage(){
    std::size_t sum = 0;
    for(auto& v : internal_ofs){
      sum += v.size() * sizeof(int) + 4;
    }
    sum += top_block_isas.size() * sizeof(int);
    for(auto& v : block_types){
      sum += (v.size() + 7) / 8;
    }
    for(auto& v : leaf_links){
      sum += v.size() * 16 + 16;
    }
    return sum;
  }
  int operator[](int idx) const{
    int block_idx = idx / root_block_size;
    int block_ofs = idx % root_block_size;
    int block_size = root_block_size;
    int isa_sum = top_block_isas[block_idx];
    int d = 0;
    while(block_size > 1){
      bool internal = block_types[d][block_idx];
      if(internal){
        int internal_rank = block_ranks[d](block_idx);
        int next_idx = internal_rank * 2 + (block_ofs >= block_size / 2);
        block_idx = next_idx;
        if(block_ofs >= block_size / 2){
          isa_sum += internal_ofs[d][internal_rank];
          block_ofs -= block_size / 2;
        }
        ++d;
        block_size >>= 1;
      }else{
        int leaf_rank = block_idx - block_ranks[d](block_idx);
        block_idx = leaf_links[d][leaf_rank].first;
        block_ofs += leaf_links[d][leaf_rank].second;
        if(block_ofs >= block_size){
          isa_sum += leaf_link_ofs[d][leaf_rank].second;
          ++block_idx;
          block_ofs -= block_size;
        }
        else{
          isa_sum -= leaf_link_ofs[d][leaf_rank].first;
        }
      }
    }
////      std::cout << isa[idx] << " = " << isa_sum << " " << bottom_nodes[block_idx] << std::endl;
//    assert(isa[idx] == isa_sum);
//    assert(block_ofs == 0);
////      std::cout << disa[idx] << " " << bottom_nodes[block_idx] << std::endl;
//    assert(disa[idx] == bottom_nodes[block_idx]);
    return isa_sum;
  }
};

class LCPIndex{
  std::vector<std::tuple<uint8_t,int,int>> symbols;
  std::vector<int> node_str_lengths;
  std::vector<int> node_dlcp_sums;
  std::vector<int> node_lcp_mins_from_start;
  int root;
  int lcp_0;
  int f(int node, int idx_in_node, int relative_d){
    auto [type, i, j] = symbols[node];
    if(type == 0){
      assert(idx_in_node == 0 || idx_in_node == -1);
      return relative_d >= 0 ? 0 : -1;
    }else if(type == 1){
      int target_idx = -1;
      if(idx_in_node != -1){
        target_idx = idx_in_node / node_str_lengths[i];
        int target_ofs = idx_in_node % node_str_lengths[i];
        int res = f(i, target_ofs, relative_d - target_idx * node_dlcp_sums[i]);
        if(res != -1){
          return res + node_str_lengths[i] * target_idx;
        }
      }
      int idx, res;
      // target_idx番目以降で探す
      if(node_dlcp_sums[i] > 0){
        idx = floor(relative_d - node_lcp_mins_from_start[i], node_dlcp_sums[i]);
      }
      else if(node_dlcp_sums[i] < 0){
        idx = ceil(relative_d - node_lcp_mins_from_start[i], node_dlcp_sums[i]);
      }
      else{
        if(node_lcp_mins_from_start[i] <= relative_d){
          idx = 0;
        }
        else{
          return -1;
        }
      }
      idx = std::max(idx, target_idx + 1);
      if(idx < j){
        res = f(i, -1, relative_d - node_dlcp_sums[i] * idx);
        if(res == -1){
          return -1;
        }
        return res + node_str_lengths[i] * idx;
      }
      return -1;
    }else{
      if(idx_in_node == -1){
        if(node_lcp_mins_from_start[i] <= relative_d){
          return f(i, -1, relative_d);
        }
        if(node_dlcp_sums[i] + node_lcp_mins_from_start[j] <= relative_d){
          return f(j, -1, relative_d - node_dlcp_sums[i]) + node_str_lengths[i];
        }
        return -1;
      }
      if(idx_in_node < node_str_lengths[i]){
        int res = f(i, idx_in_node, relative_d);
        if(res != -1){
          return res;
        }
        if(node_dlcp_sums[i] + node_lcp_mins_from_start[j] <= relative_d){
          res = f(j, -1, relative_d - node_dlcp_sums[i]);
          if(res == -1){
            return -1;
          }
          return res + node_str_lengths[i];
        }
        return -1;
      }
      else{
        int res = f(j, idx_in_node - node_str_lengths[i], relative_d - node_dlcp_sums[i]);
        if(res != -1){
          return res + node_str_lengths[i];
        }
        return -1;
      }
    }
  }
  int g(int node, int idx_in_node, int relative_d){
    auto [type, i, j] = symbols[node];
    assert(idx_in_node <= node_str_lengths[node]);
    if(type == 0){
      assert(idx_in_node == 0 || idx_in_node == -1);
      return relative_d >= 0 ? 0 : -1;
    }else if(type == 1){
      int target_idx = j;
      if(idx_in_node != -1){
        target_idx = idx_in_node / node_str_lengths[i];
        int target_ofs = idx_in_node % node_str_lengths[i];
        int res = g(i, target_ofs, relative_d - target_idx * node_dlcp_sums[i]);
        if(res != -1){
          return res + node_str_lengths[i] * target_idx;
        }
      }
      int idx, res;
      int d_diff = relative_d - node_lcp_mins_from_start[i];
      if(node_dlcp_sums[i] > 0){
        if(d_diff < 0){
          idx = -1;
        }
        else{
          idx = std::min(floor(d_diff, node_dlcp_sums[i]), target_idx - 1);
        }
      }
      else if(node_dlcp_sums[i] < 0){
        if(d_diff < 0){
          idx = ceil(d_diff, node_dlcp_sums[i]);
          if(idx >= target_idx){
            idx = -1;
          }
        }
        else{
          idx = target_idx - 1;
        }
      }
      else{
        if(node_lcp_mins_from_start[i] <= relative_d){
          idx = target_idx - 1;
        }
        else{
          return -1;
        }
      }
      if(0 <= idx){
        res = g(i, -1, relative_d - node_dlcp_sums[i] * idx);
        if(res == -1){
          return -1;
        }
        return res + node_str_lengths[i] * idx;
      }
      return -1;
    }else{
      if(idx_in_node == -1){
        if(node_dlcp_sums[i] + node_lcp_mins_from_start[j] <= relative_d){
          int res = g(j, -1, relative_d - node_dlcp_sums[i]);
          if(res != -1){
            return res + node_str_lengths[i];
          }
        }
        if(node_lcp_mins_from_start[i] <= relative_d){
//          std::cout << "Jump to lchild" << std::endl;
          return g(i, -1, relative_d);
        }
        return -1;
      }
      else if(idx_in_node >= node_str_lengths[i]){
//        std::cout << "Check: rchild" << std::endl;
        int res = g(j, idx_in_node - node_str_lengths[i], relative_d - node_dlcp_sums[i]);
        if(res != -1){
          return res + node_str_lengths[i];
        }
        if(node_lcp_mins_from_start[i] <= relative_d){
//          std::cout << "Jump to lchild" << std::endl;
          return g(i, -1, relative_d);
        }
        return -1;
      }
      else{
        int res = g(i, idx_in_node, relative_d);
        if(res != -1){
          return res;
        }
        return -1;
      }
    }
  }
public:
  LCPIndex(int root, int lcp_0, std::vector<std::tuple<uint8_t,int,int>>&& symbols, std::vector<int>&& node_str_lengths, std::vector<int>&& node_dlcp_sums, std::vector<int>&& node_lcp_mins_from_start) : root(root), lcp_0(lcp_0), symbols(std::move(symbols)), node_str_lengths(std::move(node_str_lengths)), node_dlcp_sums(std::move(node_dlcp_sums)), node_lcp_mins_from_start(std::move(node_lcp_mins_from_start)){}
  std::size_t memory_usage(){
    std::size_t sum = 0;
    sum += symbols.size() * sizeof(std::tuple<uint8_t,int,int>);
    sum += node_str_lengths.size() * sizeof(int);
    sum += node_dlcp_sums.size() * sizeof(int);
    sum += node_lcp_mins_from_start.size() * sizeof(int);
    return sum;
  }
  int find_next_idx_lcp_less_than_or_equal_to_d(int idx, int d){
    int res = f(root, idx, d - lcp_0);
    if(res != -1){
      return res;
    }
    return -1;
  };
  int find_prev_idx_lcp_less_than_or_equal_to_d(int idx, int d){
    int res = g(root, idx, d - lcp_0);
    if(res != -1){
      return res;
    }
    return -1;
  };
};

struct RLBWTBasedIndex{
  TextIndex* text_index;
  ISAIndex* isa_index;
  LCPIndex* lcp_index;
  RLBWTBasedIndex(TextIndex* text_index, ISAIndex* isa_index, LCPIndex* lcp_index) : text_index(std::move(text_index)), isa_index(std::move(isa_index)), lcp_index(std::move(lcp_index)){}
  RLBWTBasedIndex(std::string& _text);
  std::size_t memory_usage(){
    return text_index->memory_usage() + isa_index->memory_usage() + lcp_index->memory_usage();
  }
};

RLBWTBasedIndex build_rlbwt_index(std::string _text){
  _text += '\x01';
  int n = _text.size();
  sdsl::csa_wt csa;
  construct_im(csa, _text, 1);

//  {
//    sdsl::lcp_wt lcp;
//    // 0, a, aa, aaa, aaaa
//    // 0, 0, 1, 2, 3
//    std::string _t = "aaaa";
//    construct_im(lcp, _t);
//    sdsl::csa_wt csa;
//    construct_im(csa, _t);
//    for(int i = 0; i < csa.size(); ++i){
//      std::cout << csa[i] << " ";
//    }
//    std::cout << std::endl;
//    for(int i = 0; i < lcp.size(); ++i){
//      std::cout << lcp[i] << " ";
//    }
//    std::cout << std::endl;
//  }

  auto& text = csa.text;
  auto& sa = csa;
  auto& bwt = csa.bwt;
  auto& isa = csa.isa;
  auto& psi = csa.psi;
  auto lcp = sdsl::lcp_wt();
  construct_im(lcp, _text, 1);
  std::vector<int> bwt_cnt(256, 0);
  for(int i = 0; i < bwt.size(); ++i){
    if(bwt[i] != 255){
      ++bwt_cnt[bwt[i] + 1];
    }
  }
  for(int i = 1; i < 256; ++i){
    bwt_cnt[i] += bwt_cnt[i - 1];
  }

  std::vector<int> phi(n + 1);
  std::vector<int> phi_inv(n + 1);
  for(int i = 0; i <= n; ++i){
    phi[i] = isa[i] == 0 ? sa[n] : sa[isa[i] - 1];
    phi_inv[i] = sa[i] == n ? isa[0] : isa[sa[i] + 1];
  }

  sdsl::bit_vector run_border(n);
  for(int i = 0; i < n; ++i){
    if(i == 0 || i == n - 1 || bwt[i] != bwt[i - 1] || bwt[i] != bwt[i + 1]){
      run_border[i] = true;
    }
  }sdsl::rank_support_v<> r1(&run_border);
  sdsl::select_support_mcl<> s1(&run_border);
  int num_sampled_poses = std::accumulate(run_border.begin(), run_border.end(), 0);
  sdsl::bit_vector phrase_start(n);
  phrase_start[0] = true;
  for(int p = 0; p < n; ++p){
    if(run_border[p]){
      if(sa[p] > 0){
        int i = sa[p] - 1;
        phrase_start[i] = true;
      }
    }
  }
  sdsl::rank_support_v<> r2(&phrase_start);
  sdsl::select_support_mcl<> s2(&phrase_start);
  sdsl::bit_vector phrase_start_and_second(n + 1);
  phrase_start_and_second[0] = true;
  phrase_start_and_second[1] = true;
  for(int p = 0; p < n; ++p){
    if(p == 0 || bwt[p] != bwt[p - 1]){
      int i = sa[p];
      phrase_start_and_second[i] = true;
      if(i + 1 <= n){
        phrase_start_and_second[i + 1] = true;
      }
    }
  }
  sdsl::rank_support_v<> r3(&phrase_start_and_second);
  sdsl::select_support_mcl<> s3(&phrase_start_and_second);
  int r3_max = r3(n + 1);
  std::vector<int> disa(n + 1);
  disa[0] = isa[0];
  for(int i = 1; i <= n; ++i){
    disa[i] = isa[i] - isa[i - 1];
  }

//  std::vector<int> s;
//  for(int p = 0; p < n; ++p){
//    if(p == 0 || bwt[p] != bwt[p - 1]){
//      s.emplace_back(sa[p]);
//    }

//  }
//  s.emplace_back(n);
//  std::sort(s.begin(), s.end());
//  for(int iter = 0; iter + 1 < s.size(); ++iter){
//    int l = s[iter];
//    int r = s[iter + 1];
//    std::cout << "phrase: [" << l << " " << r << ")" << std::endl;
//    for(int i = l + 1; i < r; ++i){
//      if(disa[i] != disa[phi[i]]){
//        std::cout << i << ": " << disa[i] << " " << disa[phi[i]] << std::endl;
//      }
////      assert(disa[i] == disa[phi[i]]);
//    }
//  }
//  exit(0);

  std::vector<int> dlcp(n + 1);
  dlcp[0] = lcp[0];
  for(int i = 1; i <= n; ++i){
    dlcp[i] = lcp[i] - lcp[i - 1];
  }

  TextIndex* text_index = [&](){
    int root_block_size = 1;
    while(root_block_size < n / num_sampled_poses){
      root_block_size <<= 1;
    }
    // blocks: 層毎のblock一覧
    // 何番目のinternal blockであるかが分かれば次の層で何番目のblockになるかが分かる
    // leafの場合, 該当のinternalへのリンクを用意する
    std::vector<sdsl::bit_vector> block_types;
    std::vector<std::vector<std::pair<int, int>>> leaf_links;
    std::vector<uint8_t> bottom_nodes;
    std::vector<int> block_start_poses((n + root_block_size - 1) / root_block_size);
    for(int i = 0; i < block_start_poses.size(); ++i){
      block_start_poses[i] = i * root_block_size;
    }
    for(int block_size = root_block_size; block_size > 1; block_size >>= 1){
      sdsl::bit_vector block_types_(block_start_poses.size());
      std::vector<std::pair<int, int>> leaf_links_;
      std::vector<int> next_block_start_poses;
      std::vector<int> leaf_source_poses;

      for(int i = 0; i < block_start_poses.size(); ++i){
        int lb = block_start_poses[i];
        int rb = std::min(lb + block_size, n);
        int pred_sampled = r2(rb) == 0 ? -1 : s2(r2(rb)); // prev sampled position
        int succ_sampled = r2(lb) == num_sampled_poses ? -1 : s2(r2(lb) + 1); // next sampled position
        bool internal = false;
        internal |= (succ_sampled != -1 && (succ_sampled - rb + 1) <= block_size);
        internal |= (pred_sampled != -1 && (lb - pred_sampled) <= block_size);
        block_types_[i] = internal;
        if(internal){
          // internal node
          next_block_start_poses.emplace_back(lb);
          if(lb + block_size / 2 < n){
            next_block_start_poses.emplace_back(lb + block_size / 2);
          }
        }else{
          // leaf node
          int src_pos = -1;
          int l_sa = 0;
          int r_sa = n;
          for(int idx = block_size - 1; idx >= 0; --idx){
            int c = text[lb + idx];
            l_sa = csa.bwt.rank(l_sa, c) + bwt_cnt[c];
            r_sa = csa.bwt.rank(r_sa, c) + bwt_cnt[c];
          }
          int p = isa[lb];
          r_sa -= 1;

          assert(l_sa <= p && p <= r_sa);
          for(int k = 0; k < block_size; ++k){
            l_sa = phi_inv[l_sa];
            r_sa = phi_inv[r_sa];
            int next_run_border = s1(r1(l_sa) + 1);
            if(next_run_border <= r_sa){
              src_pos = sa[next_run_border] - k - 1;
              break;
            }
          }
          assert(src_pos != -1);
          leaf_source_poses.emplace_back(src_pos);
          auto it = std::prev(std::upper_bound(block_start_poses.begin(), block_start_poses.end(), src_pos));
          int block_idx = std::distance(block_start_poses.begin(), it);
          int block_ofs = src_pos - block_start_poses[block_idx];
          assert(block_ofs < block_size);
          leaf_links_.emplace_back(block_idx, block_ofs);
        }
      }
      block_start_poses = std::move(next_block_start_poses);
      block_types.emplace_back(std::move(block_types_));
      leaf_links.emplace_back(std::move(leaf_links_));
    }
    bottom_nodes.resize(block_start_poses.size());
    for(int i = 0; i < block_start_poses.size(); ++i){
      bottom_nodes[i] = text[block_start_poses[i]];
    }
//    std::vector<sdsl::rank_support_v<>> block_ranks(block_types.size());
//    for(int i = 0; i < block_types.size(); ++i){
//      block_ranks[i] = sdsl::rank_support_v<>(&block_types[i]);
//    }
    return new TextIndex(root_block_size, std::move(block_types), std::move(leaf_links), std::move(bottom_nodes));
  }();

  ISAIndex* isa_index = [&](){
    int root_block_size = 1;
    while(root_block_size < n / r3_max){
      root_block_size <<= 1;
    }
    std::vector<sdsl::bit_vector> block_types;
    std::vector<std::vector<int>> internal_ofs;
    std::vector<std::vector<std::pair<int, int>>> leaf_links;
    std::vector<std::vector<std::pair<int, int>>> leaf_link_ofs;
    std::vector<int> bottom_nodes;
    std::vector<int> block_start_poses((n + root_block_size - 1) / root_block_size);
    std::vector<int> top_block_isas((n + root_block_size - 1) / root_block_size);
    for(int i = 0; i < block_start_poses.size(); ++i){
      block_start_poses[i] = i * root_block_size;
      top_block_isas[i] = isa[block_start_poses[i]];
    }
    for(int block_size = root_block_size; block_size > 1; block_size >>= 1){
      sdsl::bit_vector block_types_(block_start_poses.size());
      std::vector<std::pair<int, int>> leaf_links_;
      std::vector<std::pair<int, int>> leaf_link_ofs_;
      std::vector<int> internal_ofs_;
      std::vector<int> next_block_start_poses;
      std::vector<int> leaf_source_poses;
      for(int i = 0; i < block_start_poses.size(); ++i){
        int lb = block_start_poses[i];
        int rb = std::min(lb + block_size, n + 1);
        int pred_sampled = r3(rb) == 0 ? -1 : s3(r3(rb)); // prev sampled position
        int succ_sampled = r3(lb) == r3_max ? -1 : s3(r3(lb) + 1); // next sampled position
        bool internal = false;
        internal |= (succ_sampled != -1 && (succ_sampled - rb + 1) <= block_size);
        internal |= (pred_sampled != -1 && (lb - pred_sampled) <= block_size);
        block_types_[i] = internal;
        if(internal){
          // internal node
          int left = lb;
          int mid = lb + block_size / 2;
          next_block_start_poses.emplace_back(left);
          if(mid <= n){
            internal_ofs_.emplace_back(isa[mid] - isa[left]);
            next_block_start_poses.emplace_back(mid);
          }
        }else{
          // leaf node
          int src_pos = -1;
          while(true){
            assert(disa[lb] == disa[phi[lb]]);
            lb = phi[lb];
            rb = std::min(lb + block_size - 1, n);
            int next_phrase_start = r3(lb) == r3_max ? -1 : s3(r3(lb) + 1);
            if(next_phrase_start != -1 && next_phrase_start <= rb){
              src_pos = lb;
              break;
            }
          }
          assert(src_pos != -1);
          leaf_source_poses.emplace_back(src_pos);
          auto it = std::upper_bound(block_start_poses.begin(), block_start_poses.end(), src_pos);
          assert(it != block_start_poses.begin());
          --it;
          int block_idx = std::distance(block_start_poses.begin(), it);
          int block_ofs = src_pos - block_start_poses[block_idx];
          int ofs_l = isa[src_pos] - isa[block_start_poses[block_idx]];
          int ofs_r = block_ofs == 0 ? 0 : isa[block_start_poses[block_idx + 1]] - isa[src_pos];
          assert(block_ofs < block_size);
          leaf_links_.emplace_back(block_idx, block_ofs);
          leaf_link_ofs_.emplace_back(ofs_l, ofs_r);
        }
      }
      block_start_poses = std::move(next_block_start_poses);
      block_types.emplace_back(std::move(block_types_));
      leaf_links.emplace_back(std::move(leaf_links_));
      leaf_link_ofs.emplace_back(std::move(leaf_link_ofs_));
      internal_ofs.emplace_back(std::move(internal_ofs_));
    }
    bottom_nodes.resize(block_start_poses.size());
    for(int i = 0; i < block_start_poses.size(); ++i){
      bottom_nodes[i] = disa[block_start_poses[i]];
    }
//    std::vector<sdsl::rank_support_v<>> block_ranks(block_types.size());
//    for(int i = 0; i < block_types.size(); ++i){
//      block_ranks[i] = sdsl::rank_support_v<>(&block_types[i]);
//    }
    return new ISAIndex(root_block_size, std::move(block_types), std::move(internal_ofs), std::move(leaf_links), std::move(leaf_link_ofs), std::move(top_block_isas));
  }();

//  // generate bidirectional scheme (unused)
//  std::vector<std::tuple<int,int,int>> factors;
//  int last_incorrect = -2;
//  for(int p = 0; p <= n; ++p){
//    if(p == 0 || bwt[p] != bwt[p - 1]){
//      int len = p - last_incorrect - 3;
//      if(len > 0){
//        if(len == 1){
//          factors.emplace_back(p, 1, dlcp[p]);
//        }
//        else{
//          factors.emplace_back(last_incorrect + 3, len, csa.lf[last_incorrect + 3]);
//        }
//        for(int l = 0; l < len; ++l){
//          assert(dlcp[last_incorrect + 3 + l] == dlcp[csa.lf[last_incorrect + 3] + l]);
//        }
//      }
//      last_incorrect = p;
//      factors.emplace_back(p, 1, dlcp[p]);
//    }
//    else if(p - last_incorrect <= 2){
//      factors.emplace_back(p, 1, dlcp[p]);
//    }
//  }
//  int last_len = n + 1 - last_incorrect - 3;
//  if(last_len > 0){
//    if(last_len == 1){
//      factors.emplace_back(n, 1, dlcp[n]);
//    }
//    else{
//      factors.emplace_back(last_incorrect + 3, last_len, csa.lf[last_incorrect + 3]);
//      for(int l = 0; l < last_len; ++l){
//        assert(dlcp[last_incorrect + 3 + l] == dlcp[csa.lf[last_incorrect + 3] + l]);
//      }
//    }
//  }
//
//  std::vector<int> v;
//  int len = 0;
//  for(auto [len_, l, src] : factors){
//    std::cout << len << " " << len_ << std::endl;
//    assert(len == len_);
//    if(l == 1){
//      std::cout << "CHAR: " << src << std::endl;
//      v.emplace_back(src);
//      ++len;
//    }
//    else{
//      std::cout << "COPY: " << l << " " << src << std::endl;
//      for(int j = 0; j < l; ++j){
//        v.emplace_back(dlcp[src + j]);
//      }
//      len += l;
//    }
//  }
//  for(int i = 0; i < dlcp.size(); ++i){
//    std::cout << v[i] << " " << dlcp[i] << std::endl;
//    assert(v[i] == dlcp[i]);
//  }

  LCPIndex* lcp_index = [&](){
    std::vector<int> terminals;
    for(int i = 1; i <= n; ++i){
      terminals.emplace_back(dlcp[i]); // we omit dlcp[0] = lcp[0]
    }
    std::sort(terminals.begin(), terminals.end());
    terminals.erase(std::unique(terminals.begin(), terminals.end()), terminals.end());

    std::vector<std::tuple<uint8_t,int,int>> symbols;
    for(auto x : terminals){
      symbols.emplace_back(0, x, 0);
    }

    std::map<std::pair<int,int>, int> map_rles;
    std::map<std::pair<int,int>, int> map_pairs;

    std::vector<int> now_seq(n);
    for(int i = 0; i < n; ++i){
      int terminal_idx = std::distance(terminals.begin(), std::lower_bound(terminals.begin(), terminals.end(), dlcp[i + 1]));
      now_seq[i] = terminal_idx;
    }

    auto get_rle_id = [&](int i, int len){
      if(map_rles.find({i, len}) == map_rles.end()){
        map_rles[{i, len}] = symbols.size();
        symbols.emplace_back(1, i, len);
      }
  //    std::cout << "RLE  " << map_rles[{i, len}] << ": " << i << " " << len << std::endl;
      return map_rles[{i, len}];
    };
    auto get_pair_id = [&](int i, int j){
      if(map_pairs.find({i, j}) == map_pairs.end()){
        map_pairs[{i, j}] = symbols.size();
        symbols.emplace_back(2, i, j);
      }
  //    std::cout << "PAIR " << map_pairs[{i, j}] << ": " << i << " " << j << std::endl;
      return map_pairs[{i, j}];
    };

    // recompression
    while(num_sampled_poses < now_seq.size()){
      std::vector<std::tuple<int,int,int>> next_factors;
      std::vector<int> next_seq;
      // rle
      for(int i = 0; i < now_seq.size(); ++i){
        if(i + 1 < now_seq.size() && now_seq[i] == now_seq[i + 1]){
          int j = i + 1;
          while(j + 1 < now_seq.size() && now_seq[j] == now_seq[j + 1]){
            ++j;
          }
          int len = j - i + 1;
          next_seq.emplace_back(get_rle_id(now_seq[i], len));
          i = j;
        }
        else{
          next_seq.emplace_back(now_seq[i]);
        }
      }
      std::cout << "rle: " << now_seq.size() << " / " << next_seq.size() << std::endl;
      now_seq = next_seq;
      next_seq.clear();
      std::vector<int> appear_symbols(now_seq.begin(), now_seq.end());
      std::sort(appear_symbols.begin(), appear_symbols.end());
      appear_symbols.erase(std::unique(appear_symbols.begin(), appear_symbols.end()), appear_symbols.end());
      std::vector<std::vector<int>> neighbors(appear_symbols.size());
      std::vector<int> mapped_now_seq(now_seq.size());
      for(int i = 0; i < now_seq.size(); ++i){
        mapped_now_seq[i] = std::distance(appear_symbols.begin(), std::lower_bound(appear_symbols.begin(), appear_symbols.end(), now_seq[i]));
      }
      for(int i = 0; i + 1 < now_seq.size(); ++i){
        int idx1 = mapped_now_seq[i];
        int idx2 = mapped_now_seq[i + 1];
        assert(idx1 != idx2);
        neighbors[idx1].emplace_back(idx2);
        neighbors[idx2].emplace_back(idx1);
      }
      std::vector<int> appear_symbols_lr(appear_symbols.size(), 0);
      for(int c = 0; c < appear_symbols.size(); ++c){
        int sum = 0;
        for(auto d : neighbors[c]){
          sum += appear_symbols_lr[d];
        }
        if(sum <= 0){
          appear_symbols_lr[c] = 1;
        }
        else{
          appear_symbols_lr[c] = -1;
        }
      }

      int sum = 0;
      for(int i = 0; i + 1 < mapped_now_seq.size(); ++i){
        if(appear_symbols_lr[mapped_now_seq[i]] != appear_symbols_lr[mapped_now_seq[i + 1]]){
          sum += appear_symbols_lr[mapped_now_seq[i]];
        }
      }

      if(sum < 0){
        for(auto& x : appear_symbols_lr){
          x *= -1;
        }
      }
      sum = 0;
      for(int i = 0; i + 1 < mapped_now_seq.size(); ++i){
        if(appear_symbols_lr[mapped_now_seq[i]] == 1 && appear_symbols_lr[mapped_now_seq[i + 1]] == -1){
          sum += 1;
        }
      }
  //    std::cout << sum << " / " << now_seq.size() << std::endl;
      assert(sum * 4 >= now_seq.size());
      // pair
      for(int i = 0; i < now_seq.size(); ++i){
        if(i + 1 < now_seq.size() && appear_symbols_lr[mapped_now_seq[i]] == 1 && appear_symbols_lr[mapped_now_seq[i + 1]] == -1){
          next_seq.emplace_back(get_pair_id(now_seq[i], now_seq[i + 1]));
          ++i;
        }
        else{
          next_seq.emplace_back(now_seq[i]);
        }
      }
      std::cout << "pair: " << now_seq.size() << " / " << next_seq.size() << std::endl;
      now_seq = next_seq;
    }

    /*
    std::vector<int> v;
    std::function<void(int)> dfs = [&](int node){
      auto [type, i, j] = symbols[node];
      if(type == 0){
        v.emplace_back(i);
      }
      else if(type == 1){
        for(int k = 0; k < j; ++k){
          dfs(i);
        }
      }
      else{
        dfs(i);
        dfs(j);
      }
    };
    for(auto i : now_seq){
      dfs(i);
    }
    for(int i = 0; i < n; ++i){
      assert(v[i] == dlcp[i + 1]);
    }*/

    while(now_seq.size() > 1){
      std::vector<int> next_seq;
      for(int i = 0; i < now_seq.size(); ++i){
        if(i + 1 < now_seq.size()){
          next_seq.emplace_back(get_pair_id(now_seq[i], now_seq[i + 1]));
          ++i;
        }
        else{
          next_seq.emplace_back(now_seq[i]);
        }
      }
      now_seq = next_seq;
    }
    int root = now_seq[0];

    std::vector<int> node_str_lengths(symbols.size(), 0);
    std::vector<int> node_dlcp_sums(symbols.size(), 0);
    std::vector<int> node_lcp_mins_from_start(symbols.size(), 0);

    std::function<void(int)> nonterminal_calc = [&](int node){
      if(node_str_lengths[node] != 0){
        return;
      }
      auto [type, i, j] = symbols[node];
      if(type == 0){
        node_str_lengths[node] = 1;
        node_dlcp_sums[node] = i;
        node_lcp_mins_from_start[node] = 0;
      }
      else if(type == 1){
        nonterminal_calc(i);
        node_str_lengths[node] = node_str_lengths[i] * j;
        node_dlcp_sums[node] = node_dlcp_sums[i] * j;
        node_lcp_mins_from_start[node] = node_dlcp_sums[node] > 0 ? node_lcp_mins_from_start[i] : node_dlcp_sums[i] * (j - 1) + node_lcp_mins_from_start[i];
      }
      else{
        nonterminal_calc(i);
        nonterminal_calc(j);
        node_str_lengths[node] = node_str_lengths[i] + node_str_lengths[j];
        node_dlcp_sums[node] = node_dlcp_sums[i] + node_dlcp_sums[j];
        node_lcp_mins_from_start[node] = std::min(node_lcp_mins_from_start[i], node_dlcp_sums[i] + node_lcp_mins_from_start[j]);
      }
    };
    nonterminal_calc(root);
    return new LCPIndex(root, dlcp[0], std::move(symbols), std::move(node_str_lengths), std::move(node_dlcp_sums), std::move(node_lcp_mins_from_start));
  }();


//// test code
//  for(int i = 0; i < n; ++i){
//    assert(text_index[i] == _text[i]);
//  }
//  for(int i = 0; i < n; ++i){
//    std::cout << isa[i] << " " << (int)isa_index[i] << std::endl;
//    assert(isa_index[i] == isa[i]);
//  }
//  for(int i = 0; i < 50000; ++i){
//    int j = (rand() + 20) % n;
//    int d = (rand() + 10) % 3000;
//    std::cout << "QUERY: " << j << " " << d << std::endl;
//    int r = lcp_index.find_succ_idx_lcp_less_than_or_equal_to_d(j, d);
//    if(r != -1){
//      std::cout << "query number: " << i << ", (j, d, result) = " << j << " " << d << " " << r << std::endl;
//      for(int k = j; k > r; --k){
//        if(lcp[k] <= d){
//          std::cout << "ERROR: " << k << " " << lcp[k] << std::endl;
//        }
//        assert(lcp[k] > d);
//      }
//      assert(lcp[r] <= d);
//    }
//  }
//
//  for(int i = 0; i < 50000; ++i){
//    int j = (rand() + 20) % n;
//    int d = (rand() + 10) % 3000;
//    std::cout << "QUERY: " << j << " " << d << std::endl;
//    int r = lcp_index.find_prev_idx_lcp_less_than_or_equal_to_d(j, d);
//    if(r != -1){
//      std::cout << "query number: " << i << ", (j, d, result) = " << j << " " << d << " " << r << std::endl;
//      for(int k = j; k < r; ++k){
//        std::cout << lcp[k] << " " << std::flush;
//        assert(lcp[k] > d);
//      }
//      std::cout << lcp[r] << std::endl;
//      assert(lcp[r] <= d);
//    }
//  }
  return RLBWTBasedIndex(text_index, isa_index, lcp_index);
}

RLBWTBasedIndex::RLBWTBasedIndex(std::string& _text) : RLBWTBasedIndex(build_rlbwt_index(_text)){}

#endif //CDAWG_LZ78_BWT_INDEX_HPP

