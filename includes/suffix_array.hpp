#ifndef CDAWG_LZ78_SUFFIX_ARRAY_HPP
#define CDAWG_LZ78_SUFFIX_ARRAY_HPP


#include <string>
#include <vector>
#include <numeric>

struct SALCP{
  std::vector<unsigned char> str;
  std::vector<int> sa,rev;
  std::vector<int> lcp;

  template<typename Seq>
  std::vector<int> sa_is(const Seq& s, int upper){
    int n = s.size();
    std::vector<int> ret(n);
    std::vector<int> is_s(n), is_lms(n);
    int m = 0;
    for(int i = n - 2; i >= 0; i--) {
      is_s[i] = (s[i] > s[i + 1]) || (s[i] == s[i + 1] && is_s[i + 1]);
      m += (is_lms[i + 1] = is_s[i] && !is_s[i + 1]);
    }
    auto induced_sort = [&](std::vector<int> &lms) {
      std::vector<int> sum(upper + 1);
      for(auto v: s) ++sum[v + 1];
      std::partial_sum(sum.begin(), sum.end(), sum.begin());
      fill(ret.begin(), ret.end(), -1);
      auto buf = sum;
      for(int i = m - 1; i >= 0; i--) {
        ret[--buf[s[lms[i]] + 1]] = lms[i];
      }
      buf = sum;
      for(auto v: ret) {
        if(v >= 1 && is_s[v - 1]) ret[buf[s[v - 1]]++] = v - 1;
      }
      buf = sum;
      for(int k = n - 1, i = ret[k]; k >= 1; i = ret[--k]) {
        if(i >= 1 && !is_s[i - 1]) {
          ret[--buf[s[i - 1] + 1]] = i - 1;
        }
      }
    };
    std::vector< int > lms;
    for(int i = 1; i < n; i++) {
      if(is_lms[i]) lms.push_back(i);
    }
    induced_sort(lms);
    std::vector<int> new_lms;
    for(int i = 0; i < n; i++) {
      if(!is_s[ret[i]] && ret[i] > 0 && is_s[ret[i] - 1]) {
        new_lms.push_back(ret[i]);
      }
    }
    auto is_same = [&](int a, int b) {
      if(s[a++] != s[b++]) return false;
      for(;; ++a, ++b) {
        if(s[a] != s[b]) return false;
        if(is_lms[a] || is_lms[b]) return is_lms[a] && is_lms[b];
      }
    };
    int rank = 0;
    ret[n - 1] = 0;
    for(int i = 1; i < m; i++) {
      if(not is_same(new_lms[i - 1], new_lms[i])) ++rank;
      ret[new_lms[i]] = rank;
    }
    if(rank + 1 < m) {
      std::vector<int> new_s(m);
      for(int i = 0; i < m; i++) {
        new_s[i] = ret[lms[i]];
      }
      auto lms_sa = sa_is(new_s, rank + 1);
      for(int i = 0; i < m; i++) {
        new_lms[i] = lms[lms_sa[i]];
      }
    }
    induced_sort(new_lms);
    return ret;
  }
  SALCP(const std::string& s, bool construct_isa = true, bool construct_lcp = true) : str(s.begin(), s.end()){
    assert(s.back() == '\1');

    int m = *std::max_element(str.begin(), str.end());
    sa = sa_is(str, m + 1); // O(n + m)
//    std::vector<int> s2(s.size());
//    std::iota(s2.begin(), s2.end(), 0);
//    std::sort(s2.begin(), s2.end(), [&](auto i, auto j){ return s.substr(i) < s.substr(j); });
//    assert(sa == s2);
    if(construct_isa){
      rev.resize(str.size());
      for(int i = 0; i < str.size(); ++i) rev[sa[i]] = i;
    }
    if(construct_lcp){
      build_lcp();
    }
  }
  void build_lcp(){
    int n = str.size();
    lcp.resize(n);
    int t=0;
    for(int i = 0; i < n; ++i) {
      int j = sa[rev[i] - 1];
      if (t > 0) t--;
      for (; j + t < n && i + t < n; ++t)
        if(str[j + t] != str[i + t]) break;
      lcp[rev[i] - 1] = t;
    }
  }
  bool lt_substr(std::string& t, int si, int ti, int tn){
    int sn = str.size();
    while(si < sn and ti < tn){
      if(str[si] < (unsigned char)t[ti]) return 1;
      if(str[si] > (unsigned char)t[ti]) return 0;
      si++, ti++;
    }
    return si == sn && ti < tn;
  }
  template<bool USE_LCP>
  int lower_bound(std::string& t, int st, int en){
    if(USE_LCP){
      // todo: O(en-st + log n) time solution with <O(n),O(1)> RmQ and LCP table
    }
    else{
      int l = -1, r = str.size();
      while(r - l > 1){
        int m = (l + r) >> 1;
        bool flag = lt_substr(t, sa[m], st, en);
        if(flag)l = m;
        else r = m;
      }
      return r;
    }
  }
  template<bool USE_LCP>
  int upper_bound(std::string& t, int st, int en){
    ++t[en - 1];
    int res = lower_bound<USE_LCP>(t, st, en);
    --t[en - 1];
    return res;
  }
};

#endif //CDAWG_LZ78_SUFFIX_ARRAY_HPP
