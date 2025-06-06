/*
PRINTING, DEBUGGING AND TESTING
*/
#include "rsqf.hpp" 
#include "additional_methods.hpp" 
#include "bqf_ec.hpp" 
#include "bqf_oom.hpp" 

#include <random>
#include <ctime>
#include <getopt.h>
#include <chrono>



#define MEM_UNIT 64ULL
#define MET_UNIT 3ULL
#define OFF_POS 0ULL
#define OCC_POS 1ULL
#define RUN_POS 2ULL

using namespace std;

void print_pair(std::pair<uint64_t,uint64_t> bound){
  std::cout << bound.first << "-" << bound.second << " ";
}

void show(uint64_t value, std::string name){
  std::cout << name << std::endl;
  print_bits(value);
}


std::string generateRandomKMer(int k) {
    static const char alphabet[] = "ACGT";
    static const int alphabetSize = sizeof(alphabet) - 1;
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<> dis(0, alphabetSize - 1);

    std::string randomKMer;
    for (int i = 0; i < k; ++i) {
        randomKMer += alphabet[dis(gen)];
    }

    return randomKMer;
}



void test_lots_of_full_cqf(){
  uint64_t seed = time(NULL);
  default_random_engine generator;
  generator.seed(seed);
  uniform_int_distribution<uint64_t> distribution;


  int qsize = 19; 
  Rsqf small_qf(qsize, 64-qsize, false);
  Rsqf usual_qf(4); //524288, qsize=19

  for (size_t j=0 ; j<1000 ; j++) {
    uint64_t seedTMP = distribution(generator);
    std::cout << "\nseed " << seedTMP << endl;

    
    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    std::cout << "j " << j << endl;
    Rsqf small_qf(qsize, 64-qsize, false);

    for (size_t i=0 ; i<(1ULL<<qsize) ; i++) { //fill to 2^qsize elements (100%-1)
      //cout << "llooooopppp" << endl;
      //cout << i << " ";
      uint64_t val = distributionTMP(generatorTMP);      
      small_qf.insert(val);
    }     
  }
}


void test_lots_of_full_cqf_enumerate() {
  uint64_t seed = time(NULL);
  default_random_engine generator;
  generator.seed(seed);
  uniform_int_distribution<uint64_t> distribution;

  int qsize = 24; 

  for (size_t j=0 ; j<1000 ; j++) {
    
    uint64_t seedTMP = distribution(generator);
    std::cout << "\nseed " << seedTMP << endl;

    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    std::cout << "j " << j << endl;
    Rsqf small_qf(qsize, 64-qsize, false);
    std::unordered_set<uint64_t> verif;

    for (size_t i=0 ; i<(1ULL<<qsize)-1 ; i++) { //fill to 2^qsize elements (100%-1)
      uint64_t val = distributionTMP(generatorTMP);      
      small_qf.insert(val);
      verif.insert(val);
    }   

    if (verif != small_qf.enumerate()) {
      std::cout << "error verif != enum" << endl;
      exit(0);
    }
  }  
}



void test_lots_of_full_cqf_remove() {
  uint64_t seed = time(NULL);
  default_random_engine generator;
  generator.seed(seed);
  uniform_int_distribution<uint64_t> distribution;

  int qsize = 19; 
  uint64_t val;
  std::unordered_set<uint64_t> enu;

  for (size_t j=0 ; j<10000 ; j++) {
    
    uint64_t seedTMP = distribution(generator);
    std::cout << "\nseed " << seedTMP << endl;

    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    std::cout << "j " << j << endl;
    Rsqf small_qf(qsize, 64-qsize, false);
    std::vector<uint64_t> verif;

    //INSERT
    for (size_t i=0 ; i<(1ULL<<qsize)-1 ; i++) { //fill to 2^qsize elements (100%-1) (1ULL<<qsize)-1
      val = distributionTMP(generatorTMP);      
      small_qf.insert(val);
      verif.push_back(val);
    }   

    //REMOVE ELEMS
    for (size_t i=0 ; i<(1ULL<<qsize)-1 ; i++) { //(1ULL<<qsize)/2
      val = verif.back();
      verif.pop_back();
      small_qf.remove(val);
    } 

    //CHECK ENUMERATE
    enu = small_qf.enumerate();
    if (verif.size() != enu.size()) {
      std::cout << "error verif != enum" << endl;
      exit(0);
    }
    for ( auto it = verif.begin(); it != verif.end(); ++it ){
      if (enu.find(*it) == enu.end()){
        std::cout << "error verif != enum" << endl;
        exit(0);
      }
    }
  }  
}




void test_one_cqf(){
  uint64_t seed = 9915022754138594861ULL; 
  default_random_engine generator;
  generator.seed(seed);
  uniform_int_distribution<uint64_t> distribution;

  std::map<uint64_t, uint64_t> enu;
  std::map<uint64_t, uint64_t> verif;

  int qsize = 8;
  int hashsize = 56;
  Bqf_ec cqf(qsize, 5, 32, 32-(hashsize/2), true);

  

  //INSERT ELEMS
  for (size_t i=0 ; i<3 ; i++) {
    //std::cout << "\ni " << i << endl;
    uint64_t val = distribution(generator);
    //val &= mask_right(qsize);
    //if (val == 0) { val += (1ULL << 45); }
    //else { val += (val<<qsize); }

    val &= mask_right(hashsize); 
    val &= mask_left(64-qsize);
    val |= 64; 

    print_bits(val);
    std::cout << "inserting " << val << " => " << (val%63) << endl; 
    cqf.insert(val, val%63);
    verif.insert({ val, val%63 });
  }

  

  //REMOVE ELEMS

  
  //CHECK ENUMERATE
  enu = cqf.enumerate();
  std::cout << "done inserting, verif size " << verif.size() << " enum size " << enu.size() << endl;

  std::cout << cqf.block2string(0) << "\n" << cqf.block2string(1) << "\n" << cqf.block2string(2) << "\n" << cqf.block2string(3);

  if (verif.size() != enu.size()) {
    std::cout << "error verif != enum" << endl;
    exit(0);
  }
  for ( auto it = verif.begin(); it != verif.end(); ++it ){
    std::cout << "checking if " << (*it).first << " => " << verif[(*it).first] << " is in enu " << endl;
    if (enu.find((*it).first) == enu.end()){
      std::cout << "error verif != enum (diff)" << endl;
      /* std::cout << "enu atm:" << endl;
      for ( auto ite = enu.begin(); ite != enu.end(); ++ite ){
        std::cout << (*ite).first << " => " << (*ite).second << endl;
      } */
      exit(0);
    }
  }


   //REMOVE ELEMS
  for (std::map<uint64_t,uint64_t>::iterator it = verif.begin(); it != verif.end(); it++){
    std::cout << "removing " << (*it).first << " => " << (*it).second << endl; 
    cqf.remove((*it).first, (*it).second);
  }
  verif.clear();


  //std::cout << cqf.block2string(0) << "\n" << cqf.block2string(1);

  //CHECK ENUMERATE
  enu = cqf.enumerate();
  if (verif.size() != enu.size()) {
    std::cout << "error verif != enum (post remove, size) verif:" << verif.size() << "  " << enu.size() << endl;
    exit(0);
  }
  for ( auto it = verif.begin(); it != verif.end(); ++it ){
    if (enu.find((*it).first) == enu.end()){
      std::cout << "error verif != enum (post remove, diff)" << endl;
      exit(0);
    }
  }

  
}


void test_time_fill_cqf(int q, int n){
  uint64_t seed = time(NULL);
  default_random_engine generator;
  generator.seed(seed);
  uniform_int_distribution<uint64_t> distribution;

  ofstream myfile;
  myfile.open ("/tmp/tmp");
  
  auto ttot = std::chrono::high_resolution_clock::now();

  

  for (int j=0 ; j<n ; j++) {

    uint64_t seedTMP = distribution(generator);
    default_random_engine generatorTMP;
    generatorTMP.seed(seedTMP);
    uniform_int_distribution<uint64_t> distributionTMP;

    Rsqf small_qf(q, 64-q, false);

    for (size_t i=0 ; i<(1ULL<<q)-1 ; i++) { //fill to 2^qsize elements (100%-1)

      auto t1 = std::chrono::high_resolution_clock::now();

      uint64_t val = distributionTMP(generatorTMP);      
      small_qf.insert(val);

      myfile << to_string( std::chrono::duration<double, std::milli>(
        std::chrono::high_resolution_clock::now() - t1
      ).count()) << "\n";
    }     
  }

  myfile.close();

  std::cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms\n";
}


void nb_false_positive(){
  std::cout << "START EXPERIMENT\n";

  std::cout << "gut_32 ========================================================================\n";
  Bqf_ec bqf = Bqf_ec::load_from_disk("/scratch/vlevallois/bqf_gut_32");
  auto ttot = std::chrono::high_resolution_clock::now(); //timer build structure + inserts 
  //verif taux de fp
  std::ifstream infile("/scratch/vlevallois/data/neg_queries.fasta");
  std::string kmer;
  uint64_t nb_fp = 0;
  while (infile >> kmer) {
      if (bqf.query(kmer).minimum > 0){
        nb_fp ++;
      }
  }
  std::cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms to check all 1 billion random kmers\n";
  std::cout << "nb de FP sur 1billion neg (supposément) queries : " << nb_fp << "\n";

  std::cout << "gut_19 ========================================================================\n";
  bqf = Bqf_ec::load_from_disk("/scratch/vlevallois/bqf_gut_19");
  ttot = std::chrono::high_resolution_clock::now(); //timer build structure + inserts 
  //verif taux de fp
  infile.clear();
  infile.seekg(0, infile.beg);
  nb_fp = 0;
  while (infile >> kmer) {
      if (bqf.query(kmer).minimum > 0){
        nb_fp ++;
      }
  }
  std::cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms to check all 1 billion random kmers\n";
  std::cout << "nb de FP sur 1billion neg (supposément) queries : " << nb_fp << "\n";

  
}

void experiments(){
  std::cout << "START EXPERIMENTs\n";


  Bqf_ec bqf = Bqf_ec::load_from_disk("/scratch/vlevallois/bqf_gut_19");
  std::cout << "/scratch/vlevallois/bqf_gut\n";

  auto ttot = std::chrono::high_resolution_clock::now(); //pos queries 
    
  std::ifstream infile("/scratch/vlevallois/data/gut_reads.fasta");
  std::string read;

  while (infile >> read) {
      bqf.query(read);
  }

  std::cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms for pos queries\n";

  ttot = std::chrono::high_resolution_clock::now(); //pos queries 
    
  std::ifstream infile2("/scratch/vlevallois/data/neg_reads.fasta");

  while (infile2 >> read) {
      bqf.query(read);
  }

  std::cout << to_string( std::chrono::duration<double, std::milli>( std::chrono::high_resolution_clock::now() - ttot ).count()) << " ms for neg queries\n";
}

int main(int argc, char** argv) {
    //test_one_cqf();

    //test_lots_of_full_cqf();

    //test_lots_of_full_cqf_enumerate();

    //test_lots_of_full_cqf_remove();

    //test_time_fill_cqf(22, 1);

    //nb_false_positive();

    experiments();
}
