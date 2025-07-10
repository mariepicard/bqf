#include "choix_signature.hpp"

/*********************************
 * TOOLS
 ********************************/

uint32_t compute_revcomp(uint32_t mmer)
{
    uint32_t revcomp = 0;

    mmer = ~mmer;
    for (int i = 0; i < M; i++)
    {
        revcomp <<= 2;
        revcomp |= mmer & 0b11;
        mmer >>= 2;
    }

    return revcomp;
}

uint32_t get_new_nucl(uint64_t encoded)
{
    return (static_cast<uint32_t>(encoded) >> (2 * M - 2)) & 0b11;
}

uint32_t min(uint32_t a, uint32_t b)
{
    if (a < b)
    {
        return a;
    }
    return b;
}

/*********************************
 * DISPLAY
 ********************************/

void display_MMer(uint64_t encoded, MMer minimizer, uint pos) {
    uint start = 2*minimizer.get_position() + 2;
    uint end = start - 2*M;
    for (uint i = 0; i < 2*K; i++) {
        uint j = 2*K - 1 - i;
        if ((j >= end && j < start) || j == pos) {
            std::cout << ((encoded >> j) & 0b1);
        } else {
            std::cout << ".";
        }
    }
    std::cout << std::endl;
}

void display_mmer(uint32_t mmer) {
    std::string mmer_display = "";
    for (int i = 0; i < M; i++) {
        char conversion[4] = {'A', 'C', 'G', 'T'};
        char c = conversion[(mmer >> (2*(M - 1))) & 0b11];
        mmer_display.push_back(c);
        mmer <<= 2;
    }
    std::cout << mmer_display << std::endl;
}

void display_kmer(uint64_t kmer) {
    std::string kmer_display = "";
    for (int i = 0; i < K; i++) {
        char conversion[4] = {'A', 'C', 'G', 'T'};
        char c = conversion[(kmer >> (2*(K - 1))) & 0b11];
        kmer_display.push_back(c);
        kmer <<= 2;
    }
    std::cout << kmer_display << std::endl;
}

/*********************************
 * MINIMAL ENCODING OF MMERS
 ********************************/

uint32_t minimal_encoding(uint32_t canon) {
    uint32_t encoding = 0;
    uint i = 0;
    uint palindromic = 0; //0 if palindromic, 1 after a specifying pair is found
    //while no specifying pair is found
    for (; 2*i+1 < M; i++){
        uint32_t r = (~canon >> (2*i)) & 0b11;
        uint32_t l = (canon >> (2*(M - i - 1))) & 0b11;

        if (r == l) {
            encoding |= r << (2*i);
        }
        else {
            uint32_t value = l*4 - l*(l+1)/2 + r - l - 1;
            value += 0b0100; //making sure the encoding of the specifying pair does not start with 0
            encoding |= (value << (2*(M - i - 2)));
            palindromic = 1;
            i++;
            break;
        }
    }
    encoding |= (mask(2*i, 2*M - 4*i) & canon)>>(2*palindromic);
    uint32_t missing = palindromic*(1<<(2*(M - i) - 1));
    return encoding - missing;
}

/*********************************
 * SLIDING WINDOW ALGORITHM
 ********************************/

std::vector<uint32_t> Signature::sliding_window_minimum(std::vector<uint32_t>& mmers, uint start) {
/*  basically, the minimum in the sliding window is the minimum of two "offseted" vectors (min_left and min_right)
    but:
        min_left is computed on the fly (we only store its required element)
        min_right is stored in the vector passed as parameter
    the response is computed in the mmers vector directly to avoid memory allocation.
 */
    uint w = K - M + 1;
    uint size_of_array = mmers.size() - start;
    if (size_of_array < w) {
        for (; mmers.size() > start; mmers.pop_back()) {}
        return mmers;
    }

    uint nbWin = size_of_array / w;
    int nb_elem_last_window = size_of_array % w;

    int min_left = mmers[start];
    for (uint i = 1; i < w; i++) {
        min_left = min_mmer(min_left, mmers[i + start]); 
    }

    for (uint i = 0; i < nbWin - 1; i++) {
        int start_window = i * w + start;
        for (int indice = start_window + w - 2; indice >= start_window; indice--) {
            // we compute "min_right" here, directly in mmers vector
            mmers[indice] = min_mmer(mmers[indice + 1], mmers[indice]); 
        }

        for (uint j = 0; j < w; j++) {
            mmers[start_window + j] = min_mmer(mmers[start_window + j], min_left); 
            min_left = (j == 0) ? mmers[start_window + w + j] : min_mmer(min_left, mmers[start_window + w + j]);
        }
    }

    // last window
    // compute min_right for last window
    int start_window = (nbWin - 1) * w + start;
    for (int indice = start_window + w - 2; indice >= start_window; indice--) {
        mmers[indice] = min_mmer(mmers[indice + 1], mmers[indice]); 
    }

    // compute the min for the last window
    for (int j = 0; j < nb_elem_last_window; j++) {
        
        mmers[start_window + j] = min_mmer(mmers[start_window + j], min_left);
        min_left = (j == 0) ? mmers[start_window + w + j] : min_mmer(min_left, mmers[start_window + w + j]);
    }
    
    mmers[start_window + nb_elem_last_window] = min_mmer(mmers[start_window + nb_elem_last_window], min_left);
    for (uint i = 0; i < w - 1; i++) {
        mmers.pop_back();
    }
    return mmers;
}

/*********************************
 * SIGN
 ********************************/

Signature::Signature(std::string filename) : filename(filename), kmers_per_min(1<<(2*M), 0)
{
    begin_build = std::chrono::high_resolution_clock::now();
    end_build = std::chrono::high_resolution_clock::now();
}

uint32_t Signature::last_canonical_mmer_encoding(uint64_t new_nucl, uint32_t& last_mmer, uint32_t& last_revcomp){
    last_revcomp = (last_revcomp >> 2) | ((~new_nucl & 0b11) << (2*(M - 1)));
    last_mmer = ((last_mmer << 2) | new_nucl) & max_value;
    return min(last_mmer, last_revcomp);
}

std::vector<uint32_t> Signature::all_signatures(std::string seq) {
    uint32_t current_mmer = 0;
    uint32_t current_revcomp = 0;
    uint lgth = 0; //lgth of the sequence already computed, starting from the last N character or the beginning
    /*position starting from which all_mmers is a sliding window of all canonical mmers read
    before start position, all_mmers contains all the signatures previously computed*/
    uint start = 0; 
    //up to start : signatures - after start : canonical mmers read
    std::vector<uint32_t> all_mmers = std::vector<uint32_t>(); 
    for (size_t i = 0; i < seq.length(); i++)
    {
        uint32_t new_nucl = 0;
        try {
            new_nucl = nucl_encode(seq[i]);
            lgth++;
        }
        catch (const std::exception &e) {
            lgth = 0;
            current_mmer = 0;
            current_revcomp = 0;
            all_mmers = sliding_window_minimum(all_mmers, start);
            start = all_mmers.size();
        }
        uint32_t canon = last_canonical_mmer_encoding(new_nucl, current_mmer, current_revcomp);//updates mmer & revcomp
        if (lgth >= M) {
            all_mmers.push_back(canon);
        }
    }
    all_mmers = sliding_window_minimum(all_mmers, start);
    return all_mmers;
}

std::vector<uint64_t> Signature::all_canonical_kmers(std::string seq) {
    uint64_t current_kmer = 0;
    uint64_t current_revcomp = 0;
    uint64_t mask_K = mask(0, 2*K);
    uint lgth = 0;
    std::vector<uint64_t> all_kmers = std::vector<uint64_t>();
    for (size_t i = 0; i < seq.length(); i++)
    {
        current_kmer <<= 2;
        current_kmer &= mask_K;
        current_revcomp >>= 2;
        uint64_t new_nucl = 0;
        try {
            new_nucl = nucl_encode(seq[i]);
            lgth++;
            current_kmer |= new_nucl;
            current_revcomp |= (~new_nucl & 0b11) << (2*(K - 1));
        }
        catch (const std::exception &e) {
            lgth = 0;
            current_kmer = 0;
            current_revcomp = 0;
        }
        
        if (lgth >= K) {
            uint64_t canon = min(current_kmer, current_revcomp);
            all_kmers.push_back(canon);
        }
    }
    return all_kmers;
}

std::pair<uint,uint> Signature::add_minimizers(std::string seq, std::vector<uint32_t>& occurences, Bqf_ec& kmers, uint64_t& superkmerlgth)
{
    //std::cout << "Size of mmer array : " << (seq.size() - M) << "\n Size of kmer array : " << (seq.size() - K) << std::endl;
    std::vector<uint32_t> signs = all_signatures(seq);
    std::vector<uint64_t> all_kmers = all_canonical_kmers(seq);
    uint arr_lgth = signs.size();
    if (arr_lgth != all_kmers.size()) {
        std::cerr << "Signature array of size " << signs.size() << std::endl;
        std::cerr << "K-mers array of size " << all_kmers.size() << std::endl;
        throw std::runtime_error("different array lgth");
    }

    uint nb_distinct_min = 0;
    uint nb_skipped = 0;
    for (uint i = 0; i < arr_lgth; i++) {
        if (signs[i] == max_value) {
            nb_skipped ++;
        } else {
            uint32_t position_in_bqf = minimal_encoding(signs[i]);//sign_with_extra_bit(all_kmers[i], minimizer);
            kmers_per_min[position_in_bqf]++;
            if (!kmers.insert(kmer_to_hash(all_kmers[i], K))) {
                occurences[position_in_bqf]++;
                nb_distinct_min++;
            }
            
            if (i > 0 && signs[i - 1] == signs[i])
                superkmerlgth++;
        }
    }
    return std::make_pair(nb_distinct_min, nb_skipped);
}

void Signature::output()
{
    std::string path = "../../";
    std::string distname = path + name() + "_dist.csv";
    std::string timename = path + name() + "_time.txt";
    std::ofstream timefile(timename);
    if (!timefile.is_open())
    {
        throw std::runtime_error("Could not open file " + timename);
        exit(EXIT_FAILURE);
    }
    timefile << "build time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end_build - begin_build).count() << "ms" << std::endl;
    std::vector<std::string> files;
    files.push_back(filename);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(files, 1, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    std::vector<uint32_t> occurences(1 << (2 * M), 0);
    Bqf_ec kmers = Bqf_ec(18,1,K,0,false);
    auto begin = std::chrono::high_resolution_clock::now();
    uint nb_distinct_min = 0;
    uint nb_skipped = 0;
    uint64_t superkmerlgth = 0;
    while (parser.refill(rg))
    {
        for (auto &rp : rg)
        {
            std::pair<uint,uint> min_skipped = add_minimizers(rp.seq, occurences, kmers, superkmerlgth);
            nb_distinct_min += min_skipped.first;
            nb_skipped += min_skipped.second;
        }
    }
    parser.stop();

    auto stop = std::chrono::high_resolution_clock::now();
    auto time_per_minimizer = std::chrono::duration_cast<std::chrono::nanoseconds>((stop - begin) / (nb_distinct_min + nb_skipped));
    timefile << "computing minimizers time: " << time_per_minimizer.count() << "ns" << std::endl;
    timefile << "Found minimizers: " << nb_distinct_min << std::endl;
    timefile << "Skipped minimizers: " << nb_skipped << std::endl;
    timefile << "Number of follow-ups: " << superkmerlgth << std::endl;
    timefile.close();

    std::ofstream repartitionfile(distname);
    if (!repartitionfile.is_open())
    {
        throw std::runtime_error("Could not open file " + distname);
        exit(EXIT_FAILURE);
    }
    std::string output_in_file = "mmer,nb_distinct_occs,nb_total_occs\n";
    for (uint32_t mmer = 0; mmer < max_value/2 + 1; mmer++)
    {
        if (occurences[mmer] > 0) {
            output_in_file += std::to_string(mmer) + "," + std::to_string(occurences[mmer]) + "," + std::to_string(kmers_per_min[mmer]) + "\n";
        }
    }
    repartitionfile << output_in_file;
    repartitionfile << max_value << "," << nb_skipped << std::endl;
    repartitionfile.close();
}

uint32_t Signature::sign_with_extra_bit(uint64_t encoded, MMer minimizer) {
    uint idx = (2*minimizer.get_position() + (K - M + 1))%(2*K);
    return (minimizer.get_min_encoding() << 1) | ((encoded >> idx) & 0b1);
}

/********************************
 * KMC signature
 *******************************/

bool KMC_sign::is_allowed(uint32_t mmer)
{
    if ((mmer & 0x3f) == 0x3f) // TTT suffix
        return false;
    if ((mmer & 0x3f) == 0x3b) // TGT suffix
        return false;
    if ((mmer & 0x3c) == 0x3c) // TG* suffix !!!! consider issue #152
        return false;

    for (uint32_t j = 0; j < M - 3; ++j)
        if ((mmer & 0xf) == 0) // AA inside
            return false;
        else
            mmer >>= 2;

    if (mmer == 0) // AAA prefix
        return false;
    if (mmer == 0x04) // ACA prefix
        return false;
    if ((mmer & 0xf) == 0) // *AA prefix
        return false;
    return true;
}

KMC_sign::KMC_sign(std::string filename) : Signature(filename)
{
    end_build = std::chrono::high_resolution_clock::now();
}

std::string KMC_sign::name()
{
    return "kmc_signature";
}

uint32_t KMC_sign::last_canonical_mmer_encoding(uint64_t new_nucl, uint32_t& last_mmer, uint32_t& last_revcomp){
    last_revcomp = (last_revcomp >> 2) | ((~new_nucl & 0b11) << (2*(M - 1)));
    last_mmer = ((last_mmer << 2) | new_nucl) & max_value;
    uint32_t canon = min(last_revcomp, last_mmer);
    if (is_allowed(canon)) {
        return canon;
    }
    else {
        return max_value;
    }
}


/********************************
 * Roberts signature : mapping
 *******************************/

Roberts_sign::Roberts_sign(std::string filename) : Signature(filename)
{
    end_build = std::chrono::high_resolution_clock::now();
}

std::string Roberts_sign::name()
{
    return "alternate_mappings_signature";
}

uint32_t Roberts_sign::last_canonical_mmer_encoding(uint64_t new_nucl, uint32_t& last_mmer, uint32_t& last_revcomp){
    last_revcomp = (~last_revcomp) & max_value;
    new_nucl = (new_nucl^0b01) & 0b11;
    last_revcomp = (last_revcomp >> 2) | ((~new_nucl & 0b11) << (2*(M - 1)));
    last_mmer = (((~last_mmer) << 2) | new_nucl) & max_value;
    return min(last_revcomp, last_mmer);
}


/********************************
 * Wood signature : xor
 *******************************/

Wood_sign::Wood_sign(std::string filename) : Signature(filename)
{
    end_build = std::chrono::high_resolution_clock::now();
}

std::string Wood_sign::name()
{
    return "xor_signature";
}


/********************************
 * Frequency-based signature
 *******************************/

void add_all_mmers(std::string sequence, std::vector<uint32_t> &occ)
{
    if (sequence.length() < M)
        return;
    uint32_t mask_2M = mask(0,2*M);
    uint32_t mmer = 0;
    uint32_t revcomp = 0;
    uint32_t canon;
    uint lgth = 0;
    uint32_t new_nucl = 0;
    for (size_t i = 0; i < sequence.length(); i++)
    {
        try {
            new_nucl = nucl_encode(sequence[i]);
            lgth++;
        }
        catch (const std::exception &e) {
            lgth = 0;
            mmer = 0;
        }
        mmer <<= 2;
        mmer &= mask_2M;
        mmer |= new_nucl;
        revcomp >>= 2;
        revcomp |= (~new_nucl & 0b11) << (2 * M - 2);
        if (lgth >= M) {
            canon = min(mmer, revcomp);
            occ[canon]++;
        }
    }
}

std::vector<uint32_t> argsort(std::vector<uint32_t> &v)
{
    std::vector<uint32_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(), [&v](uint32_t i, uint32_t j)
                     { return v[i] < v[j]; });
    return idx;
}

Frequency_sign::Frequency_sign(std::string filename) : Signature(filename), order(1 << (2 * M))
{
    std::vector<std::string> files;
    files.push_back(filename);
    fastx_parser::FastxParser<fastx_parser::ReadSeq> parser(files, 1, 1);
    parser.start();
    auto rg = parser.getReadGroup();
    std::vector<uint32_t> occurences(1 << (2 * M), 0);
    while (parser.refill(rg))
    {
        for (auto &rp : rg)
        {
            add_all_mmers(rp.seq, occurences);
        }
    }
    parser.stop();

    std::vector<uint32_t> idx = argsort(occurences);
    most_occ = idx.back();
    
    for (int i = 0; i < (1 << (2 * M)); i++)
    {
        order[idx[i]] = i;
    }
    end_build = std::chrono::high_resolution_clock::now();
};

std::string Frequency_sign::name()
{
    return "frequency_based_signature";
}


/********************************
 * main
 *******************************/

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::cerr << "Please provide filename and signature mode" << std::endl;
        return EXIT_FAILURE;
    }
    std::string filename = std::string(argv[1]);
    if (argc == 2) {
        Signature signature = Signature(filename);
        signature.output();
    }
    else {
        std::string mode = std::string(argv[2]);
        if (mode == "KMC") {
            KMC_sign signature = KMC_sign(filename);
            signature.output();
        } else if (mode == "XOR") {
            Wood_sign signature = Wood_sign(filename);
            signature.output();
        } else if (mode == "MAP") {
            Roberts_sign signature = Roberts_sign(filename);
            signature.output();
        } else if (mode == "FRE") {
            Frequency_sign signature = Frequency_sign(filename);
            signature.output();
        } else {
            std::cerr << "Invalid mode" << std::endl;
            return EXIT_FAILURE;
        }
    }
    return 0;
}
