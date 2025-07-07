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

uint32_t minimal_encoding(uint32_t canon, uint m) {
    uint32_t encoding = 0;
    uint i = 0;
    uint palindromic = 0; //0 if palindromic, 1 after a specifying pair is found
    //while no specifying pair is found
    for (; 2*i+1 < m; i++){
        uint32_t r = (~canon >> (2*i)) & 0b11;
        uint32_t l = (canon >> (2*(m - i - 1))) & 0b11;

        if (r == l) {
            encoding |= r << (2*i);
        }
        else {
            uint32_t value = l*4 - l*(l+1)/2 + r - l - 1;
            value += 0b0100; //making sure the encoding of the specifying pair does not start with 0
            encoding |= (value << (2*(m - i - 2)));
            palindromic = 1;
            i++;
            break;
        }
    }
    encoding |= (mask(2*i, 2*m - 4*i) & canon)>>(2*palindromic);
    uint32_t missing = palindromic*(1<<(2*(m - i) - 1));
    return encoding - missing;
}

/*********************************
 * SLIDING WINDOW ALGORITHM
 ********************************/

 SlidingWindowMinimum::addElement(uint32_t element, std::function<bool(uint32_t, uint32_t)> compare) {
    while (!minima.empty() && compare(minima.back(), element)) {
        minima.pop_back();
    }
    minima.push_back(element);
 }

/*********************************
 * SIGN
 ********************************/

Signature::Signature(std::string filename) : filename(filename), kmers_per_min(1<<(2*M), 0)
{
    begin_build = std::chrono::high_resolution_clock::now();
    end_build = std::chrono::high_resolution_clock::now();
}

std::pair<uint,uint> Signature::add_minimizers(std::string seq, std::vector<uint32_t>& occurences, Bqf_ec& kmers, uint64_t& superkmerlgth)
{
    uint nb_distinct_min = 0;
    uint nb_skipped = 0;
    MMer prev = MMer(max_value, K);
    uint64_t encoded = 0;
    uint lgth = 0;
    for (size_t i = 0; i < seq.length(); i++)
    {
        encoded <<= 2;
        try {
            encoded |= nucl_encode(seq[i]);
            lgth++;
        }
        catch (const std::exception &e) {
            lgth = 0;
            encoded = 0;
        }
        if (lgth >= K) {
            MMer minimizer = sign(encoded, prev);
            uint32_t position_in_bqf = sign_with_extra_bit(encoded, minimizer);
            kmers_per_min[position_in_bqf]++;
            if (kmers.insert(kmer_to_hash(encoded, K))) {
                if (minimizer.get_position() < K) {
                    //std::cout << minimizer.get_mmer() << std::endl;
                    occurences[position_in_bqf]++;
                    nb_distinct_min++;
                }
                else {
                    nb_skipped ++;
                }
            }
            
            if (prev.get_mmer() == minimizer.get_mmer() && lgth > K)
                superkmerlgth++;
            prev = minimizer;
        }
        
    }
    return std::make_pair(nb_distinct_min, nb_skipped);
}

void Signature::output()
{
    std::string timename = "../../" + name() + "_time.txt";
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
    Bqf_ec kmers = Bqf_ec(20,1,K,0,false);
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
    auto time_per_minimizer = std::chrono::duration_cast<std::chrono::nanoseconds>((stop - begin) / nb_distinct_min);
    timefile << "computing minimizers time: " << time_per_minimizer.count() << "ns" << std::endl;
    timefile << "Found minimizers: " << nb_distinct_min << std::endl;
    timefile << "Skipped minimizers: " << nb_skipped << std::endl;
    timefile << "Number of follow-ups: " << superkmerlgth << std::endl;
    timefile.close();

    std::string distname = "../../" + name() + "_dist.csv";
    std::ofstream repartitionfile(distname);
    if (!repartitionfile.is_open())
    {
        throw std::runtime_error("Could not open file " + distname);
        exit(EXIT_FAILURE);
    }
    repartitionfile << "mmer,nb_distinct_occs,nb_total_occs" << std::endl;
    for (uint32_t mmer = 0; mmer < max_value; mmer++)
    {
        uint32_t canon = min(mmer, compute_revcomp(mmer));
        repartitionfile << mmer << "," << occurences[mmer] << "," << kmers_per_min[mmer] << std::endl;

    }
    repartitionfile << max_value << "," << nb_skipped << std::endl;
    repartitionfile.close();
}

MMer Signature::sign(uint64_t encoded){
    uint32_t mask = max_value;
    uint32_t min_mmer = max_value;
    uint pos = K;
    uint32_t current_mmer = encoded & mask;
    uint32_t current_revcomp = compute_revcomp(current_mmer);

    for (uint64_t i = M - 1; i < K; i++)
    {
        uint32_t canon = min(current_mmer, current_revcomp);
        if (canon < min_mmer)
        {
            min_mmer = canon;
            pos = i;
        }

        encoded >>= 2;
        current_mmer = encoded & mask;
        current_revcomp <<= 2;
        current_revcomp |= get_new_nucl(~encoded);
        current_revcomp &= mask;
    }
    return MMer(min_mmer, pos);
}

MMer Signature::sign(uint64_t encoded, MMer prev)
{
    prev.update_pos();
    if (prev.get_position() >= K)
    {
        return sign(encoded);
    }
    uint32_t last_mmer = static_cast<uint32_t>(encoded) & max_value;
    uint32_t canon = min(last_mmer, compute_revcomp(last_mmer));
    if (canon <= prev.get_mmer())
    {
        return MMer(canon, M - 1);
    }
    return prev;
}

uint32_t Signature::sign_with_extra_bit(uint64_t encoded, MMer minimizer) {
    uint idx = (2*minimizer.get_position() + (K - M + 1))%(2*K);
    //display_MMer(encoded, minimizer, idx);
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

MMer KMC_sign::sign(uint64_t encoded)
{
    uint32_t min_mmer = max_value;
    uint pos = K;
    uint32_t mask = max_value;
    uint32_t current_mmer = static_cast<uint32_t>(encoded) & max_value;
    uint32_t current_revcomp = compute_revcomp(current_mmer);
    for (int i = M - 1; i < K; i++)
    {
        uint32_t canon = min(current_mmer, current_revcomp);
        bool allowed = is_allowed(canon);
        if (canon < min_mmer && allowed)
        {
            min_mmer = canon;
            pos = i;
        }
        encoded >>= 2;
        current_mmer = encoded & mask;
        current_revcomp <<= 2;
        current_revcomp |= get_new_nucl(~encoded);
        current_revcomp &= mask;
    }
    return MMer(min_mmer, pos);
};

MMer KMC_sign::sign(uint64_t encoded, MMer prev)
{
    prev.update_pos();
    if (prev.get_position() >= K)
    {
        return sign(encoded);
    }
    uint32_t last_mmer = static_cast<uint32_t>(encoded) & max_value;
    uint32_t revcomp = compute_revcomp(last_mmer);
    uint32_t canon = min(last_mmer, revcomp);
    uint32_t min_mmer = prev.get_mmer();

    if (canon <= min_mmer && is_allowed(last_mmer))
    {
        return MMer(canon, M - 1);
    }
    return prev;
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

MMer Roberts_sign::sign(uint64_t encoded)
{
    uint32_t min_mmer = max_value;
    uint pos = K;
    uint32_t current_mmer = 0;
    uint32_t current_revcomp = 0;
    // init revcomp / mmer
    for (int i = 0; i < M - 1; i++)
    {
        uint32_t nucl = (encoded & 0b11) ^ 0b01;
        current_mmer |= nucl << (2 * i);
        current_revcomp <<= 2;
        current_revcomp |= ~nucl & 0b11;
        encoded >>= 2;
        encoded = ~encoded;
    }
    // loop over all existing mmers
    for (int i = M - 1; i < K; i++)
    {
        uint32_t nucl = (encoded & 0b11) ^ 0b01;
        current_mmer |= nucl << (2 * M - 2);
        current_revcomp <<= 2;
        current_revcomp |= ~nucl & 0b11;
        current_revcomp &= max_value;
        uint32_t canon = min(current_mmer, current_revcomp);
        if (canon < min_mmer)
        {
            min_mmer = canon;
            pos = i;
        }
        encoded >>= 2;
        encoded = ~encoded;
        current_mmer >>= 2;
    }
    return MMer(min_mmer, pos);
}

MMer Roberts_sign::sign(uint64_t encoded, MMer prev)
{
    prev.update_pos();
    if (prev.get_position() >= K)
    {
        return sign(encoded);
    }
    uint32_t last_mmer = static_cast<uint32_t>(encoded) & max_value;
    uint32_t revcomp = compute_revcomp(last_mmer);
    uint32_t canon = min(last_mmer, revcomp);
    uint32_t min_mmer = prev.get_mmer();
    if (canon <= min_mmer)
    {
        return MMer(canon, M - 1);
    }
    return prev;
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

MMer Wood_sign::sign(uint64_t encoded)
{
    uint32_t mask = max_value;

    uint32_t min_bin_key = max_value;
    uint32_t min_mmer = max_value;
    uint pos = K;
    uint32_t current_mmer = encoded & mask;
    uint32_t current_revcomp = compute_revcomp(current_mmer);

    for (uint64_t i = M - 1; i < K; i++)
    {
        uint32_t canon = min(current_mmer, current_revcomp);
        uint32_t temp_bin_key = xor_mask ^ canon;
        if (temp_bin_key < min_bin_key)
        {
            min_bin_key = temp_bin_key;
            min_mmer = canon;
            pos = i;
        }

        encoded >>= 2;
        current_mmer = encoded & mask;
        current_revcomp <<= 2;
        current_revcomp |= get_new_nucl(~encoded);
        current_revcomp &= mask;
    }
    return MMer(min_mmer, pos);
}

MMer Wood_sign::sign(uint64_t encoded, MMer prev)
{
    prev.update_pos();
    if (prev.get_position() >= K)
    {
        return sign(encoded);
    }
    uint32_t last_mmer = static_cast<uint32_t>(encoded) & max_value;
    uint32_t revcomp = compute_revcomp(last_mmer);
    uint32_t canon = min(last_mmer, revcomp);
    uint32_t min_mmer = prev.get_mmer();
    if ((canon ^ xor_mask) <= (min_mmer ^ xor_mask))
    {
        return MMer(canon, M - 1);
    }
    return prev;
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

MMer Frequency_sign::sign(uint64_t encoded)
{
    uint32_t mmer = static_cast<uint32_t>(encoded) & max_value;
    uint32_t revcomp = compute_revcomp(mmer);
    uint32_t mask = max_value;

    uint pos = K;
    uint32_t min_value = most_occ;

    for (int i = M - 1; i < K; i++)
    {
        uint32_t canon = min(mmer, revcomp);
        if (order[canon] < order[min_value])
        {
            min_value = canon;
            pos = i;
        }
        encoded >>= 2;
        mmer = encoded & mask;
        revcomp <<= 2;
        revcomp |= get_new_nucl(~encoded);
        revcomp &= mask;
    }
    return MMer(min_value, pos);
};

MMer Frequency_sign::sign(uint64_t encoded, MMer prev)
{
    prev.update_pos();
    if (prev.get_position() >= K)
    {
        return sign(encoded);
    }
    uint32_t last_mmer = static_cast<uint32_t>(encoded) & max_value;
    uint32_t canon = min(last_mmer, compute_revcomp(last_mmer));
    uint32_t min_mmer = prev.get_mmer();
    if (order[canon] <= order[min_mmer])
    {
        return MMer(canon, M - 1);
    }
    return prev;
}

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
