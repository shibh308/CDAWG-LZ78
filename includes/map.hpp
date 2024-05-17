#ifndef CDAWG_LZ78_MAP_HPP
#define CDAWG_LZ78_MAP_HPP


#include <bit>
#include <cassert>
#include <cstdint>
#include <numeric>
#include <optional>
#include <map>
#include <vector>
#include "vector.hpp"


template<typename K, typename V>
struct Map {
  virtual std::optional<V> find(K key) const = 0;
  virtual int size() const = 0;
  virtual std::uint64_t num_bytes() const = 0;
};

template <typename T, typename U>
struct DynamicHashMap : Map<T, U> {
  static constexpr std::uint64_t z = 65521;
  static constexpr T null = 0;
  // static constexpr std::uint64_t z = 60xf332ac987401cba5;
  std::uint16_t n, d;

  Vector<std::pair<T, U>, std::uint16_t> v;

  DynamicHashMap() : n(0), d(1){
    std::vector<std::pair<T, U>> v_(2, std::make_pair(null, U()));
    v = decltype(v)(v_);
  }
  explicit DynamicHashMap(const std::vector<T>& keys, const std::vector<U>& values) : n(0), d(1), v(2, std::make_pair(null, U())){
    for(int i = 0; i < keys.size(); ++i){
      add(keys[i], values[i]);
    }
  }
  inline std::uint64_t hash(T key) const{return (z * key) & ((1u << d) - 1); }

  int size() const override { return int(n); }

  std::optional<U> find(T x) const override {
    for(std::uint64_t i = hash(x); v[i].first != null; i = (i + 1) & ((1u << d) - 1)){
      if(v[i].first == x){
        return v[i].second;
      }
    }
    return std::nullopt;
  }

  void add(T x, U val){
    assert(x != null);
    if((1u << d) < ((n + 1) * 2u)){
      resize();
    }
    std::uint64_t i = hash(x);
    for(; v[i].first != null && v[i].first != x; i = (i + 1) & ((1u << d) - 1));
    n += v[i].first == null;
    v[i] = {x, val};
  }

  std::vector<std::pair<T, U>> items() const{
    std::vector<std::pair<T, U>> items;
    for(int i = 0; i < v.size(); ++i){
      auto& item = v[i];
      if(item.first != null){
        items.emplace_back(item);
      }
    }
    sort(items.begin(), items.end());
    return items;
  }

  void resize(){
    ++d;
    decltype(v) old_table;
    swap(old_table, v);
    std::vector<std::pair<T, U>> v_(1u << d, std::make_pair(null, U()));
    v = decltype(v)(v_);
    assert(v.size() <= 512);
    n = 0;
    for(int i = 0; i < old_table.size(); ++i){
      auto& item = old_table[i];
      if(item.first != null){
        add(item.first, item.second);
      }
    }
  }
  std::uint64_t num_bytes() const{
    return sizeof(n) + sizeof(d) + v.num_bytes();
  }
};




template <typename T, typename U>
struct HashMap : Map<T, U> {
  static constexpr std::uint64_t z = 65521;
  static constexpr T null = 0;
  // static constexpr std::uint64_t z = 60xf332ac987401cba5;
  std::uint64_t n, d;
  // TOOD: type check

  // std::vector<std::pair<T, U>> v;
  Vector<std::pair<T, U>, std::uint16_t> v;
  HashMap() : n(0), d(1), v(){}

  explicit HashMap(const std::vector<T>& keys, const std::vector<U>& values) : n(0), d(1){
    std::vector<std::pair<T, U>> vec(2, std::make_pair(null, U()));
    for(int i = 0; i < keys.size(); ++i){
      add(vec, keys[i], values[i]);
    }
    v = vec;
  }

  explicit HashMap(const DynamicHashMap<T, U>& hashmap) : n(hashmap.n), d(hashmap.d), v(hashmap.v){
  }
  explicit HashMap(const HashMap<T, U>& other) : n(other.n), d(other.d), v(other.v){
  }

  inline std::uint64_t hash(T key) const{return (z * key) & ((1u << d) - 1); }

  int size() const override { return int(n); }

  std::optional<U> find(T x) const override {
    for(std::uint64_t i = hash(x); v[i].first != null; i = (i + 1) & ((1u << d) - 1)){
      if(v[i].first == x){
        return v[i].second;
      }
    }
    return std::nullopt;
  }

  std::vector<std::pair<T, U>> items(){
    std::vector<std::pair<T, U>> items;
    for(int i = 0; i < v.size(); ++i){
      auto& item = v[i];
      if(item.first != null){
        items.emplace_back(item);
      }
    }
    sort(items.begin(), items.end());
    return items;
  }
  std::uint64_t num_bytes() const{
    return sizeof(n) + sizeof(d) + v.num_bytes();
  }

private:
  void add(std::vector<std::pair<T, U>>& vec, T x, U val){
    assert(x != null);
    if((1u << d) < ((n + 1) * 2u)){
      resize(vec);
    }
    std::uint64_t i = hash(x);
    for(; vec[i].first != null && vec[i].first != x; i = (i + 1) & ((1u << d) - 1));
    n += vec[i].first == null;
    vec[i] = {x, val};
  }

  void resize(std::vector<std::pair<T, U>>& vec){
    ++d;
    std::vector<std::pair<T, U>> old_table;
    swap(old_table, vec);
    vec.assign(1u << d, {null, U()});
    assert(vec.size() <= 512);
    n = 0;
    for(auto item : old_table){
      if(item.first != null){
        add(vec, item.first, item.second);
      }
    }
  }
};


template<typename K, typename V>
struct BinarySearchMap : Map<K, V> {
  Vector<std::pair<K, V>, std::uint8_t> items_;
  explicit BinarySearchMap(){}
  explicit BinarySearchMap(const std::map<K, V>& map){
    std::vector<std::pair<K, V>> items__;
    for(auto [k, v] : map){
      items_.emplace_back(k, v);
    }
    items_ = items__;
  }
  explicit BinarySearchMap(const HashMap<K, V>& map){
    items_ = map.items();
  }
  BinarySearchMap(const std::vector<K>& keys, const std::vector<V>& values){
    assert(keys.size() == values.size());
    for(int i = 0; i < static_cast<int>(keys.size()) - 1; ++i){
      assert(keys[i] < keys[i + 1]);
    }
    std::vector<std::pair<K, V>> items__;
    for(int i = 0; i < keys.size(); ++i){
      items__.emplace_back(keys[i], values[i]);
    }
    items_ = items__;
  }
  explicit BinarySearchMap(const DynamicHashMap<K, V>& hashmap) : items_(hashmap.items()){
  }
  std::optional<V> find(const K key) const override{
    constexpr int linear_search_border = 3;
    unsigned int l = 0, r = items_.size();
    while(r - l > linear_search_border){
      unsigned int mid = (r + l) >> 1u;
      if(items_[mid].first == key){
        return items_[mid].second;
      }else if(items_[mid].first < key){
        l = mid;
      }else{
        r = mid;
      }
    }
    for(unsigned int i = l; i < r; ++i){
      if(items_[i].first == key){
        return items_[i].second;
      }
      else if(key < items_[i].first){
        return std::nullopt;
      }
    }
    return std::nullopt;
  }
  std::vector<std::pair<K, V>> items() const{
    std::vector<std::pair<K, V>> items__(items_.size());
    for(int i = 0; i < items_.size(); ++i){
      items__[i] = items_[i];
    }
    return items__;
  }
  int size() const override{
    return items_.size();
  }
  virtual std::uint64_t num_bytes() const override{
    return items_.num_bytes();
  }
};

/*
template <typename T, typename U>
struct StdMapWrapper : Map<T, U>{
    std::map<T, U> map;
    explicit StdMapWrapper(){}
    explicit StdMapWrapper(const HashMap<T, U>& map_){
        for(auto item : map_.items()){
            map.insert(item);
        }
    }
    StdMapWrapper(const std::vector<T>& keys, const std::vector<U>& values){
        for(int i = 0; i < keys.size(); ++i){
            map.insert({keys[i], values[i]});
        }
    }
    int size() const override{ return map.size(); }
    std::optional<U> find(T key) const override{
        auto iter = map.find(key);
        if(iter == map.end()){
            return std::nullopt;
        }
        else{
            return iter->second;
        }
    }
    bool add(T x, U key){
        map.erase(x);
        return map.insert({x, key}).second;
    }
    std::vector<std::pair<T, U>> items() const{
        return std::vector<std::pair<T, U>>(map.begin(), map.end());
    }
};
*/


#endif //CDAWG_LZ78_MAP_HPP
