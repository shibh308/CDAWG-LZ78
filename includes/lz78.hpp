#ifndef CDAWG_LZ78_LZ78_HPP
#define CDAWG_LZ78_LZ78_HPP

#include <functional>
#include <string_view>
#include "markedancestor.hpp"


using random_access_func_type = const std::function<int(int)>;
using sa_range_func_type = const std::function<std::pair<int,int>(int,int)>;
//template<typename random_access_func_type, typename sa_range_func_type>
std::pair<std::vector<std::pair<int,unsigned char>>, size_t> compute_lz78(int n, random_access_func_type random_access_func, sa_range_func_type sa_range_func){

  std::vector<std::pair<int,unsigned char>> phrases;
  MarkedAncestor ma(n);

  for(int i = 0, k = 0; i < n; ++k){
    auto [ls, rs] = sa_range_func(i, n - i);

    assert(rs - ls == 1);
    auto [mark_depth, phrase_id] = ma.get_mark(ls);
    if(i + mark_depth == n){
      phrases.emplace_back(phrase_id, '\0');
      break;
    }
    else{
      phrases.emplace_back(phrase_id, random_access_func(i + mark_depth));
    }
    int phrase_length = mark_depth + 1;
    auto [l, r] = sa_range_func(i, phrase_length);
    ma.mark(l, r, phrase_length, k);
    i += phrase_length;
  }
  return {phrases, ma.get_memory_usage()};
}

#endif //CDAWG_LZ78_LZ78_HPP
