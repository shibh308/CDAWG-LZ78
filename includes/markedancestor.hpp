#ifndef CDAWG_LZ78_MARKEDANCESTOR_HPP
#define CDAWG_LZ78_MARKEDANCESTOR_HPP

#endif //CDAWG_LZ78_MARKEDANCESTOR_HPP

template <typename T, typename U>
struct SplayTree{
  struct Node{
    int size;
    T key;
    U value, lazy;
    Node* par;
    Node* c[2];
    Node(){}
    Node(T key, U val, Node* nil) : key(key), value(val), lazy(std::make_pair(0, -1)), size(1), par(nil){ c[0] = nil; c[1] = nil;}
  };
  using NodePtr = Node*;
  NodePtr nil;

  SplayTree(){
    nil = new Node(T(), std::make_pair(0, -1), nullptr);
    nil->size = 0;
    nil->par = nil->c[0] = nil->c[1] = nil;
  }

  void delete_tree(NodePtr root){
    delete_node(root);
  }

  void delete_node(NodePtr node){
    if(node->c[0] != nil)delete_node(node->c[0]);
    if(node->c[1] != nil)delete_node(node->c[1]);
    delete node;
  }

  NodePtr make(T key, U val){
    return new Node(key, val, nil);
  }

  void _update(NodePtr x){
    if(x == nil)
      return;
    assert(x != x->c[0]);
    assert(x != x->c[1]);
    assert(x->c[0] == nil || x->c[0] != x->c[1]);
    if(x->c[0] != nil)x->c[0]->lazy = std::max(x->c[0]->lazy, x->lazy);
    if(x->c[1] != nil)x->c[1]->lazy = std::max(x->c[1]->lazy, x->lazy);
    x->value = std::max(x->value, x->lazy);
    x->lazy = std::make_pair(0, -1);
    x->size = x->c[0]->size + x->c[1]->size + 1;
    assert(x->size > 0);
  }

  void rotate(NodePtr p, bool p_right){
    NodePtr x = p->c[p_right];
    NodePtr q = p->par;
    assert(p->c[p_right] == x);
    p->c[p_right] = x->c[p_right ^ 1];
    if(x->c[p_right ^ 1] != nil){
      p->c[p_right]->par = p;
      assert(p != p->c[p_right]);
    }
    // xp間の辺の張り替え
    p->par = x;
    x->c[p_right ^ 1] = p;
    // pq間の辺の張り替え
    x->par = q;
    if(q != nil){
      bool q_right = (q->c[1] == p);
      assert(q->c[q_right] == p);
      q->c[q_right] = x;
    }
    _update(p), _update(x), _update(q);
  }

  void reroot(NodePtr x){
    while(x->par != nil){
      NodePtr p = x->par;
      NodePtr q = p->par;
      bool q_right = (q->c[1] == p);
      bool p_right = (p->c[1] == x);
      if(p->par == nil){
        rotate(p, p_right);
        break;
      }
      // 同じ向きの二回回転
      if(q_right == p_right){
        rotate(q, q_right), rotate(p, p_right);
      }
      else{
        rotate(p, p_right), rotate(q, q_right);
      }
    }
  }

  std::pair<NodePtr, bool> _lower_bound(NodePtr p, T key){
    if(p == nil)
      return std::make_pair(nil, false);
    if(p->key < key){
      auto res = _lower_bound(p->c[1], key);
      return res.second ? res : std::make_pair(p, false);
    }
    else{
      auto res = _lower_bound(p->c[0], key);
      return res.second ? res : std::make_pair(p, true);
    }
  }

  std::pair<NodePtr, bool> lower_bound(NodePtr p, T key){
    if(p == nil)
      return std::make_pair(p, false);
    auto res = _lower_bound(p, key);
    reroot(res.first);
    assert(res.first != nil);
    return res;
  }

  NodePtr access(NodePtr p, int idx){
    if(p == nil)
      return nil;
    while(p->c[0]->size != idx){
      if(p->c[0]->size < idx)
        idx -= p->c[0]->size + 1, p = p->c[1];
      else
        p = p->c[0];
      if(p == nil)
        return nil;
    }
    reroot(p);
    return p;
  }

  std::pair<NodePtr, bool> insert(NodePtr root, T key, U val){
    if(root == nil)
      return std::make_pair(make(key, val), true);
    NodePtr l, r, np;
    bool exist;
    // lower_boundの結果からsplitする時、lower_boundの結果がnilだとバグるので注意
    std::tie(np, exist) = lower_bound(root, key);
    if(exist){
      if(np->key == key)
        return std::make_pair(np, false);
      std::tie(l, r) = split(np);
      return std::make_pair(merge(merge(l, make(key, val)), r), true);
    }
    else{
      return std::make_pair(merge(np, make(key, val)), true);
    }
  }
  // [0, p), [p, n)でsplit
  std::pair<NodePtr, NodePtr> split(NodePtr p){
    if(p == nil)
      return std::make_pair(nil, nil);
    reroot(p);
    NodePtr l = p->c[0];
    l->par = nil;
    p->c[0] = nil;
    _update(p);
    return std::make_pair(l, p);
  }

  NodePtr merge(NodePtr p, NodePtr q){
    reroot(p);
    reroot(q);
    if(q == nil)
      return p;
    if(p == nil)
      return q;
    while(p->c[1] != nil)
      p = p->c[1];
    reroot(p);
    assert(p->c[1] == nil);
    p->c[1] = q;
    q->par = p;
    _update(p);
    assert(p != nil);
    return p;
  }
};

template<typename T>
struct Segtree{
  int n;
  T op;
  std::vector<T> elm;

  Segtree(int n, T init, T op) :
    n(n),
    op(op),
    elm(2 * n, op){
    for(int i = 0; i < n; ++i)
      elm[i + n] = init;
  }
  void update(int x, int y, T val){
    for(x += n, y += n - 1; x <= y; x >>= 1, y >>= 1){
      if(x & 1){
        elm[x] = std::max(elm[x], val);
        x++;
      }
      if(!(y & 1)){
        elm[y] = std::max(elm[y], val);
        y--;
      }
    }
  }
  T get(int x) const{
    x += n;
    T val = elm[x];
    while(x >>= 1)
      val = std::max(val, elm[x]);
    return val;
  }
};

struct MarkedAncestorBySegmentTree{
  Segtree<std::pair<int,int>> seg;
  explicit MarkedAncestorBySegmentTree(int n) : seg(n, std::make_pair(0, -1), std::make_pair(0, -1)){}
  void mark(int l, int r, int mark_depth, int phrase_id){
    seg.update(l, r, {mark_depth, phrase_id});
  }
  /// return (mark_depth, phrase_id)
  /// if there is no mark on path, return (0, -1)
  std::pair<int,int> get_mark(int x){
    return seg.get(x);
  }
};

struct MarkedAncestor{
  int n;
  SplayTree<int, std::pair<int,int>> splay;
  using SplayNode = SplayTree<int, std::pair<int,int>>::NodePtr;
  SplayNode root;
  explicit MarkedAncestor(int n) : n(n){
    root = splay.nil;
    root = splay.insert(root, n, std::make_pair(0, -1)).first;
  }
  ~MarkedAncestor(){
    splay.delete_tree(root);
    delete splay.nil;
  }
  void mark(int l, int r, int mark_depth, int phrase_id){
    root = splay.lower_bound(root, l - 1).first;

    auto left_value = root->value;
    root = splay.insert(root, l - 1, left_value).first;

    SplayNode left, right;
    root = splay.lower_bound(root, l).first;
    std::tie(left, root) = splay.split(root);

    root = splay.lower_bound(root, r - 1).first;
    auto rleft_value = std::max(root->value, std::make_pair(mark_depth, phrase_id));
    root = splay.lower_bound(root, r).first;
    std::tie(root, right) = splay.split(root);
    if(root != splay.nil){
      root->lazy = std::max(root->lazy, std::make_pair(mark_depth, phrase_id));
      splay._update(root);
    }
    root = splay.insert(root, r - 1, rleft_value).first;
    root = splay.merge(left, root);
    root = splay.merge(root, right);
  }
  /// return (mark_depth, phrase_id)
  /// if there is no mark on path, return (0, -1)
  std::pair<int,int> get_mark(int x){
    root = splay.lower_bound(root, x).first;
    auto result = root->value;
    return result;
  }
  std::size_t get_memory_usage(){
    std::size_t sum = 0;
    sum += root->size * sizeof(*root);
    return sum;
  }
};


//struct MarkedAncestor{
//
//#define NO_ELEMENT -1
//#define MULTIPLE_ELEMENT -2
//  struct Node{
//    int value, phrase_id;
//    int l_range;
//    int r_range;
//    int lc, rc;
//    Node(int value, int phrase_id, int l_range, int r_range) : value(value), phrase_id(phrase_id), l_range(l_range), r_range(r_range), lc(-1), rc(-1){}
//  };
//
//  int height;
//  std::vector<Node> nodes;
//
//  explicit MarkedAncestor(int n){
//    for(height = 0; (1 << height) < n; ++height);
//    add_node(0, -1, 0, 1 << height);
//  }
//  int add_node(int value, int phrase_id, int l_range, int r_range){
//    nodes.emplace_back(value, phrase_id, l_range, r_range);
//    return nodes.size() - 1;
//  }
//  void mark(int l, int r, int mark_depth, int phrase_id, int node_id = 0, int lb = 0, int rb = 0){
//    if(node_id == 0){
//      rb = 1 << height;
//    }
//    auto& node = nodes[node_id];
//    if(node.leaf_idx == NO_ELEMENT){
//
//    }
//    else if(node.leaf_idx == MULTIPLE_ELEMENT){
//
//    }
//    if(l == lb && r - l == length){
//
//    }
//  }
//  /// return (mark_depth, phrase_id)
//  /// if there is no mark on path, return (0, -1)
//  std::pair<int,int> get_mark(int x) const{
//
//  }
//};
