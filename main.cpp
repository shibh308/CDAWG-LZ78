
#include <iostream>
#include <fstream>
#include <filesystem>
#include <chrono>
#include "includes/dawg.hpp"
#include "includes/suffix_array.hpp"
#include "includes/lz78.hpp"

std::string load_file(const std::string& data_path, std::size_t length_limit){
  std::clog << "loading: " << data_path << std::endl;
  std::ifstream file(data_path);
  assert(file.is_open());
  std::string text;
  text.resize(length_limit);
  file.read(&text[0], length_limit);
  std::size_t read_bytes = file.gcount();
  if(read_bytes != length_limit){
    text.resize(read_bytes);
  }
  std::clog << "Loading file \"" << data_path << "\" is finished." << std::endl;
  return text;
}

bool is_valid_text(std::string_view text){
  for(auto c : text){
    if(c == '\0' || c == '\1'){
      return false;
    }
  }
  return true;
}

void correctness_check(HeavyPathDecomposedCDAWG& hp_cdawg, const std::string& text){
  SALCP sa(text);
  bool fl = true;
  for(int i = 0; i < text.size(); ++i){
    for(int j = 0; j <= text.size() - i; ++j){
      auto [l, r] = hp_cdawg.climb_linear(i, j);
      std::string substr = text.substr(i, j);
      if(!substr.empty()){
        int lb = sa.lower_bound<false>(substr, 0, substr.size());
        int rb = sa.upper_bound<false>(substr, 0, substr.size());
        if(lb != l || rb != r){
          std::cerr << "Error: The output is incorrect." << std::endl;
          std::cerr << "range : [" << i << ", " << i + j << ")" << std::endl;
          std::cerr << "string: " << text << std::endl;
          std::cerr << "substr: " << text.substr(i, j) << std::endl;
          std::cerr << "SA range (CDAWG)   : [" << l << ", " << r << ")" << std::endl;
          std::cerr << "SA range (correct) : [" << lb << ", " << rb << ")" << std::endl;
          fl = false;
        }
      }
    }
  }
  if(fl){
    std::clog << "Result: correct" << std::endl;
  }
  else{
    std::clog << "Result: incorrect" << std::endl;
  }
}

std::vector<std::pair<int,unsigned char>> compute_lz78_by_suffix_array(std::string& text){
  std::clog << "Computing LZ78 by suffix array..." << std::endl;
  auto random_access_func = [&](int i){ return text[i]; };
  SALCP sa(text);
  auto sa_range_func = [&](int start_pos, int target_depth){

    // i=58, k =36, 3275, 3279

    int lb = sa.lower_bound<false>(text, start_pos, start_pos + target_depth);
    int rb = sa.upper_bound<false>(text, start_pos, start_pos + target_depth);
    return std::make_pair(lb, rb);
  };
  return compute_lz78(text.length(), random_access_func, sa_range_func).first;
}

std::vector<std::pair<int,unsigned char>> compute_lz78_by_CDAWG(std::string& text){
  std::clog << "Computing LZ78 by CDAWG..." << std::endl;
  auto scdawg = [&](){
    CDAWGBase cdawg(text);
    return SimpleCDAWG(cdawg);
  }();
  auto random_access_func = [&](int i){
    auto c1 = text[i];
    auto c2 = scdawg.climb_linear(i, 1).second;
    assert((unsigned char)c1 == c2);
    return c2;
  };
  auto sa_range_func = [&](int start_pos, int target_depth){
    return scdawg.climb_linear(start_pos, target_depth).first;
  };
  return compute_lz78(text.length(), random_access_func, sa_range_func).first;
}

void check_correctness(std::string& text){
  CDAWGBase cdawg(text);
  HeavyPathDecomposedCDAWG hp_cdawg(cdawg);
  SimpleCDAWG scdawg(cdawg);
  for(int i = 0; i < text.length(); ++i){
    std::clog << "start: " << i << std::endl;
    for(int j = 0; i + j <= text.length(); ++j){
      auto a = hp_cdawg.climb_linear(i, j);
      auto b = scdawg.climb_linear(i, j).first;
      auto c = hp_cdawg.climb_log2(i, j);
      assert(a == b);
      assert(a == c);
    }
  }
}


struct BenchMarkResult{
  std::string filename;
  int text_length, num_vertices, num_edges;
  size_t memory_usage_raw_text, memory_usage_cdawg, memory_usage_ma;
  double elapsed_time;
  std::vector<int> num_iter_bins;
  void output_clog(){
    std::clog << std::endl;
    std::clog << "filename                : " << filename << std::endl;
    std::clog << "text length             : " << text_length << std::endl;
    std::clog << "number of vertices      : " << num_vertices << std::endl;
    std::clog << "number of edges         : " << num_edges << std::endl;
    std::clog << "memory usage (raw text) : " << memory_usage_raw_text / 1024.0 << " [kB]" << std::endl;
    std::clog << "memory usage (CDAWG)    : " << memory_usage_cdawg / 1024.0 << " [kB]" << std::endl;
    std::clog << "                          " << memory_usage_cdawg / text_length << " [bytes/character]" << std::endl;
    std::clog << "memory usage (LZ78)     : " << memory_usage_ma / 1024.0 << " [kB]" << std::endl;
    std::clog << "elapsed time (LZ78)     : " << elapsed_time << " [ms]" << std::endl;
  }
  static void output_csv_header(std::ofstream& of){
    of << "filename" << ",";
    of << "text_length" << ",";
    of << "num_vertices" << ",";
    of << "num_edges" << ",";
    of << "memory_usage_text" << ",";
    of << "memory_usage_cdawg" << ",";
    of << "memory_usage_ma" << ",";
    of << "elapsed_time_lz78" << ",";
    of << "num_iter_bins" << std::endl;
  }
  void output_csv(std::ofstream& of){
    of << std::fixed;
    of << filename << ",";
    of << text_length << ",";
    of << num_vertices << ",";
    of << num_edges << ",";
    of << memory_usage_raw_text << ",";
    of << memory_usage_cdawg << ",";
    of << memory_usage_ma << ",";
    of << elapsed_time << ",";
    of << "\"[" << num_iter_bins[0];
    for(int i = 1; i < num_iter_bins.size(); ++i){
      of << "," << num_iter_bins[i];
    }
    of << "]\"" << std::endl;
  }
};



void benchmark(std::string filename, size_t length, std::ofstream& of){

  auto text = load_file("../data/" + filename, length);

  int n = text.length();

  auto scdawg = [&](){
    CDAWGBase cdawg(text);
    text.clear();
    return SimpleCDAWG(cdawg);
  }();
  auto random_access_func = [&](int i){ return scdawg.climb_linear(i, 1).second; };
  auto sa_range_func = [&](int start_pos, int target_depth){
    return scdawg.climb_linear(start_pos, target_depth).first;
  };

  auto time_st = std::chrono::high_resolution_clock::now();
  size_t memory_usage = compute_lz78(n, random_access_func, sa_range_func).second;
  auto time_en = std::chrono::high_resolution_clock::now();
  auto computation_time = std::chrono::duration_cast<std::chrono::milliseconds>(time_en - time_st).count();
  std::vector<int> num_iter_vec(n);
  for(int i = 0; i < n; ++i){
    scdawg.climb_linear(i, n - i);
    num_iter_vec[i] = scdawg.num_iter;
  }

  BenchMarkResult res;
  res.filename = filename;
  res.text_length = n;
  res.num_vertices = scdawg.node_length.size();
  res.num_edges = scdawg.parents.size();
  res.memory_usage_cdawg = scdawg.get_memory_usage();
  res.memory_usage_raw_text = n;
  res.memory_usage_ma = memory_usage;
  res.elapsed_time = computation_time;
  std::vector<int> num_iter_bins(*std::max_element(num_iter_vec.begin(), num_iter_vec.end()) + 1);
  for(int i = 0; i < n; ++i){
    ++num_iter_bins[num_iter_vec[i]];
  }
  res.num_iter_bins = std::move(num_iter_bins);

  res.output_clog();
  res.output_csv(of);
}

int main(int argc, char** argv) {

  if(argc == 3){
    std::string output_file = "../results/output.csv";
    bool exists = std::filesystem::exists(output_file);
    std::string filename = argv[2];
    int n = atoi(argv[1]);
    std::ofstream of(output_file, std::ios_base::out | std::ios_base::app);
    if(!exists){
      BenchMarkResult::output_csv_header(of);
    }
    benchmark(filename, n, of);
//    for(int n = (1 << 3); n <= (1 << 25); n <<= 1){
//      benchmark("sources", n, of);
//      benchmark("dna", n, of);
//      benchmark("english", n, of);
//    }
    std::exit(0);
  }
  else{
    std::clog << "expected args: {filename} {n}" << std::endl;
    std::exit(0);
  }



  auto text = load_file("../data/sources", 10000000); // 1024 * 1024 * 10);

//  text = "baabaab" + std::string(1000, 'a') + "b";
  std::clog << "text length: " << text.length() << std::endl;
//  text = "aabaab";
  assert(is_valid_text(text));
  text += '\1';

  auto lz78_2 = compute_lz78_by_CDAWG(text);
  auto lz78 = compute_lz78_by_suffix_array(text);

  assert(lz78 == lz78_2);

  std::vector<std::string> v;
  for(auto [k, pre] : lz78){
    std::string s;
    if(k != -1){
      s += v[k];
    }
    if(pre){
      s += pre;
    }
    v.emplace_back(s);
  }
//  for(auto x : v)std::cout << x << " ";
//  std::cout << std::endl;

  return 0;
}
