#ifndef SIGNATURE_HPP
#define SIGNATURE_HPP

#include <vector>
#include <numeric>
#include <stdint.h>
#include <string>
#include <iostream>
#include <fstream>
#include "FastxParser.hpp"
#include "bqf_ec.hpp"
#include "additional_methods.hpp"
#include <algorithm>
#include <deque>
#include <chrono>
#include <functional>
static const uint M = 10;
static const uint K = 32;
static const uint64_t INDEX2_XOR_MASK = 0xe37e28c4271b5a2dULL;

uint32_t compute_revcomp(uint32_t mmer);

void display_bits(uint32_t mmer, std::ostream& str) {
    for (int i = 0; i < 32; i++){
        str << ((mmer>>(32 - 1))&1);
        if (((32 - 2 - i)&0b111) == 0b111) {
            std::cout << " ";
        }
        mmer <<= 1;
    }
}

uint32_t mask(int number_bits_right, int length) {
    uint32_t left = ((~0) << number_bits_right);
    uint32_t right = ~((~0) << (number_bits_right + length));
    return left & right;
}

/**
 * @brief encoding of canonical mmers using 2*m - 1 bits (odd m)
 * it relies on the algorithm proposed by Roland Wittler in "General  encoding  of  canonical k-mers",
 * DOI = https://doi.org/10.1101/2023.03.09.531845
 * Note that since it is only used on signatures, it is only the initial version implemented 
 * and not the sliding window one.
 * \param canon : canonical mmer to encode, encoded with a standard encoding
 * \return an encoding using one less bit
 */

uint32_t minimal_encoding(uint32_t canon);

class MMer
{
    uint32_t mmer;
    uint32_t min_encoding;
    uint pos;

public:
    MMer(uint32_t mmer, uint pos) : mmer(mmer), min_encoding(minimal_encoding(mmer)), pos(pos) {}
    uint32_t get_mmer() { return mmer; }
    uint32_t get_min_encoding() { return min_encoding; }
    uint get_position() { return pos; }
};



/*********************************
 * SIGN
 ********************************/
class Signature
{
    std::string filename;

public:
    //maximum value of a kmer
    const uint32_t max_value = (1 << (2 * M)) - 1;
    //start and end of building the signature (useful if pre-computing is done)
    std::chrono::high_resolution_clock::time_point begin_build;
    std::chrono::high_resolution_clock::time_point end_build;
    //number of kmers for each quotient value
    std::vector<uint32_t> kmers_per_min;

    /********** CONSTRUCTORS **********/

    Signature() {};
    Signature(std::string filename);
    /**
     * @brief name of the signature, used to distinguish signatures when building files
     * */
    virtual std::string name() {return "default_signature";}

    /****** COMPARISON FUNCTIONS ******/
    /** 
     * @brief order of mmers defined by the signature. Default is lexicographic order
     * \returns u < v
     */
    virtual bool compare(uint32_t u, uint32_t v) { return u < v;}
    /** 
     * @brief minimum of two kmers using signature comparison function
     */
    uint32_t min_mmer(uint32_t u, uint32_t v) { 
        if (compare(u,v)) {
            return u;
        } return v;
    }

    /****** COMPUTING SIGNATURES ******/

    /** @brief In place sliding window minimum algorithm
     * using https://github.com/lrobidou/sliding-minimum-windows/blob/main/min.cpp implementation
     * under license GNU Affero General Public License v3.0
     * the algorithm can be found in https://doi.org/10.1093/bioinformatics/btad305 
     * the code was adapted to 
     *      - use uint32_t integers
     *      - use different comparison functions
     *      - use a starting position in the array
     * \param mmers : array of all canonical mmers to use. Up to the position start, 
     * it contains the signatures of each kmer in a sequence. Starting from position start, it contains
     * all the m-factors of the following kmers. 
     * \param start : the position from which to compute the minimizers
     * The position start is necessary because of N characters that breaks kmers consecutivity
     * \returns an array of all minimizers
     * */
    std::vector<uint32_t> sliding_window_minimum(std::vector<uint32_t>& mmers, uint start);
    /**
     * @brief lists all signatures from a sequence, using a sliding window minimum algorithm
     * \param seq is the sequence from which to compute all signatures
     * \returns a vector v containing all signatures, such that for each i, v[i] is the signature of the ith kmer
     */
    std::vector<uint32_t> all_signatures(std::string seq);
    /**
     * @brief lists all kmers from a sequence
     * \param seq is the sequence from which to compute all kmers
     * \returns a vector containing all kmers
     */
    std::vector<uint64_t> all_canonical_kmers(std::string seq);

    /** 
     * @brief knowing the previous mmer and the previous revcomp, and the new nucleotide read, 
     * it updates the mmer and the revcomp so that they are the new mmer and revcomp being read
     * and returns the canonical representation of said mmer
     * \param new_nucl the last nucleotide being read, encoded with traditional encoding 
     * (could be worth changing it to quick encoding)
     * \param last_mmer the previous mmer from the sequence (updated)
     * \param last_revcomp the previous revcomp from the sequence (updated)
     * \returns the canonical representation of the new mmer being read
     */
    virtual uint32_t last_canonical_mmer_encoding(uint64_t new_nucl, uint32_t& last_mmer, uint32_t& last_revcomp);
    /**
     * @brief add an extra bit of padding to the signature, chosen as far away from the signature as possible
     * \param encoded : the encoded kmer
     * \param minimizer : the signature of the kmer
     * \returns a quotient value with a deterministic extra bit
     */
    uint32_t sign_with_extra_bit(uint64_t encoded, MMer minimizer);

    std::pair<uint,uint> add_minimizers(std::string seq, std::vector<uint32_t>& occurences, Bqf_ec& kmers, uint64_t& superkmerlgth);

    void output();
};

class KMC_sign : public Signature
{
    static bool is_allowed(uint32_t mmer);

public:
    
    KMC_sign(std::string filename);
    std::string name() override;
    uint32_t last_canonical_mmer_encoding(uint64_t new_nucl, uint32_t& last_mmer, uint32_t& last_revcomp) override;
};

class Roberts_sign : public Signature
{
public:

    Roberts_sign(std::string filename);

    std::string name() override;
    uint32_t last_canonical_mmer_encoding(uint64_t new_nucl, uint32_t& last_mmer, uint32_t& last_revcomp) override;
};

class Wood_sign : public Signature
{
    uint32_t xor_mask = INDEX2_XOR_MASK & max_value;

public:
    Wood_sign(std::string filename);

    std::string name() override;
    bool compare(uint32_t u, uint32_t v) override {
        return (u^xor_mask) < (v^xor_mask);
    }
};

class Frequency_sign : public Signature
{
    std::vector<uint32_t> order;
    uint32_t most_occ;

public:
    Frequency_sign(std::string filename);

    std::string name() override;
    bool compare(uint32_t u, uint32_t v) override {
        return order[u] < order[v];
    }
};

#endif
