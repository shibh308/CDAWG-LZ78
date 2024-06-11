#ifndef CDAWG_LZ78_WAVELETMATRIX_HPP
#define CDAWG_LZ78_WAVELETMATRIX_HPP

#include <vector>
#include <array>
#include <cstdint>

struct BitVector{
  std::vector<std::uint64_t> v;
  std::vector<int> r;
  BitVector(){}
  void build(){
    r.assign(v.size() + 1, 0);
    for(int i = 0; i < v.size(); ++i)
      r[i + 1] = r[i] + __builtin_popcountll(v[i]);
  }
  bool access(int x){
    return (v[x >> 6] >> (x & 63)) & 1;
  }
  // [0, x)の1の出現回数
  int rank(int x){
    return r[x >> 6] + __builtin_popcountll(v[x >> 6] & ((1uLL << (x & 63)) - 1));
  }
  int rank(int x, bool fl){
    return fl ? rank(x) : x - rank(x);
  }
};

template <typename T, int W>
struct WaveletMatrix{

  std::array<BitVector, W> bv;
  std::array<int, W> zero_cnt;

  WaveletMatrix(std::vector<T>& a){
    int n = a.size();
    std::vector<T> v(a);
    for(int i = W - 1; i >= 0; --i){
      std::vector<uint64_t> b((n >> 6) + 1, 0);
      std::vector<T> v1, v2;
      for(int j = 0; j < n; ++j){
        ((v[j] >> i) & 1 ? v2 : v1).push_back(v[j]);
        b[j >> 6] |= uint64_t((v[j] >> i) & 1) << (j & 63);
      }
      for(int j = 0; j < v.size(); ++j)
        v[j] = (j < v1.size() ? v1[j] : v2[j - v1.size()]);
      bv[i].v = move(b);
      bv[i].build();
      zero_cnt[i] = bv[i].rank(n, 0);
    }
  }

  // [l, r)内で[0, x)を満たす値の数
  int count_lower(int l, int r, T x){
    int cnt = 0;
    for(int i = W - 1; i >= 0; --i){
      bool fl = (x >> i) & 1;
      int st = bv[i].rank(l, fl);
      int en = bv[i].rank(r, fl);
      if(fl){
        st += zero_cnt[i];
        en += zero_cnt[i];
        cnt += (bv[i].rank(r, 0) - bv[i].rank(l, 0));
      }
      l = st, r = en;
    }
    return cnt;
  }
};

#endif //CDAWG_LZ78_WAVELETMATRIX_HPP
