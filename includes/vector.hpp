#ifndef CDAWG_LZ78_VECTOR_HPP
#define CDAWG_LZ78_VECTOR_HPP

#ifndef PACKED_DAWG_VECTOR_HPP
#define PACKED_DAWG_VECTOR_HPP

#include <memory>
#include <vector>
#include <cstdint>

template<typename T, typename size_type>
class Vector{
  std::unique_ptr<T[]> pointer;
  size_type _size;
public:
  explicit Vector(){
    _size = 0;
    pointer = std::make_unique<T[]>(0);
  }
  Vector(const std::vector<T>& vector){
    _size = vector.size();
    pointer = std::make_unique<T[]>(_size);
    for(int i = 0; i < _size; ++i){
      pointer[i] = vector[i];
    }
  }
  Vector(const Vector& other){
    _size = other._size;
    pointer = std::make_unique<T[]>(_size);
    for(int i = 0; i < _size; ++i){
      pointer[i] = other.pointer[i];
    }
  }
  ~Vector() = default;
  Vector& operator=(const Vector& other) {
    if (this != &other) {
      Vector tmp(other);
      std::swap(tmp._size, _size);
      std::swap(tmp.pointer, pointer);
    }
    return *this;
  }
  T& operator[](std::size_t index){
    return pointer[index];
  }
  const T& operator[](std::size_t index) const{
    return pointer[index];
  }
  static constexpr std::uint64_t offset_bytes = sizeof(T*) + sizeof(size_type);
  size_type size() const{
    return _size;
  }
  std::uint64_t num_bytes() const{
    return offset_bytes + sizeof(T) * _size;
  }

  // イテレーターの定義
  class iterator {
    T* ptr;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = T*;
    using reference = T&;

    iterator(pointer ptr) : ptr(ptr) {}

    iterator() = default;
    ~iterator() = default;

    iterator(const iterator&) = default;
    iterator& operator=(const iterator&) = default;
    iterator(iterator&&) noexcept = default;
    iterator& operator=(iterator&&) noexcept = default;

    reference operator*() const { return *ptr; }
    pointer operator->() const { return ptr; }
    iterator& operator++() { ++ptr; return *this; }
    iterator operator++(int) { iterator tmp = *this; ++ptr; return tmp; }
    iterator& operator--() { --ptr; return *this; }
    iterator operator--(int) { iterator tmp = *this; --ptr; return tmp; }
    iterator& operator+=(difference_type n) { ptr += n; return *this; }
    iterator operator+(difference_type n) const { return iterator(ptr + n); }
    iterator& operator-=(difference_type n) { ptr -= n; return *this; }
    iterator operator-(difference_type n) const { return iterator(ptr - n); }
    difference_type operator-(const iterator& rhs) const { return ptr - rhs.ptr; }
    reference operator[](difference_type n) const { return *(ptr + n); }

    bool operator==(const iterator& rhs) const { return ptr == rhs.ptr; }
    bool operator!=(const iterator& rhs) const { return ptr != rhs.ptr; }
    bool operator<(const iterator& rhs) const { return ptr < rhs.ptr; }
    bool operator>(const iterator& rhs) const { return ptr > rhs.ptr; }
    bool operator<=(const iterator& rhs) const { return ptr <= rhs.ptr; }
    bool operator>=(const iterator& rhs) const { return ptr >= rhs.ptr; }
  };

  class const_iterator {
    const T* ptr;

  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type = T;
    using difference_type = std::ptrdiff_t;
    using pointer = const T*;
    using reference = const T&;

    // コンストラクター
    const_iterator(pointer ptr) : ptr(ptr) {}

    const_iterator() = default;
    ~const_iterator() = default;

    const_iterator(const const_iterator&) = default;
    const_iterator& operator=(const const_iterator&) = default;
    const_iterator(const_iterator&&) noexcept = default;
    const_iterator& operator=(const_iterator&&) noexcept = default;

    reference operator*() const { return *ptr; }
    pointer operator->() const { return ptr; }
    const_iterator& operator++() { ++ptr; return *this; }
    const_iterator operator++(int) { const_iterator tmp = *this; ++ptr; return tmp; }
    const_iterator& operator--() { --ptr; return *this; }
    const_iterator operator--(int) { const_iterator tmp = *this; --ptr; return tmp; }
    const_iterator& operator+=(difference_type n) { ptr += n; return *this; }
    const_iterator operator+(difference_type n) const { return const_iterator(ptr + n); }
    const_iterator& operator-=(difference_type n) { ptr -= n; return *this; }
    const_iterator operator-(difference_type n) const { return const_iterator(ptr - n); }
    difference_type operator-(const const_iterator& rhs) const { return ptr - rhs.ptr; }
    reference operator[](difference_type n) const { return *(ptr + n); }

    bool operator==(const const_iterator& rhs) const { return ptr == rhs.ptr; }
    bool operator!=(const const_iterator& rhs) const { return ptr != rhs.ptr; }
    bool operator<(const const_iterator& rhs) const { return ptr < rhs.ptr; }
    bool operator>(const const_iterator& rhs) const { return ptr > rhs.ptr; }
    bool operator<=(const const_iterator& rhs) const { return ptr <= rhs.ptr; }
    bool operator>=(const const_iterator& rhs) const { return ptr >= rhs.ptr; }
  };

  iterator begin() { return iterator(pointer.get()); }
  iterator end() { return iterator(pointer.get() + _size); }
  const_iterator begin() const { return const_iterator(pointer.get()); }
  const_iterator end() const { return const_iterator(pointer.get() + _size); }
};


#endif //PACKED_DAWG_VECTOR_HPP

#endif //CDAWG_LZ78_VECTOR_HPP