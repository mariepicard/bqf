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
#define M 15
#define K 32
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

uint32_t minimal_encoding(uint32_t canon, uint m);

class MMer
{
    uint32_t mmer;
    uint32_t min_encoding;
    uint pos;

public:
    MMer(uint32_t mmer, uint pos) : mmer(mmer), min_encoding(minimal_encoding(mmer, M)), pos(pos) {}
    uint32_t get_mmer() { return mmer; }
    uint32_t get_min_encoding() { return min_encoding; }
    uint get_position() { return pos; }
    void update_pos() { pos += 1; }
};


class SlidingWindowMinimum {
    std::deque<uint32_t> minima = std::deque<uint32_t>();

public :
    uint32_t getMinimum() {return minima.front();}
    void clear () {minima.clear();}
    void addElement(uint32_t element, std::function<bool(uint32_t, uint32_t)> compare);
    void deleteElement(uint32_t element) {
        if (minima.front() == element) {
            minima.pop_front();
        }
    };
}

/*********************************
 * SIGN
 ********************************/
class Signature
{
    std::string filename;
    

public:
    const uint32_t max_value = (1 << (2 * M)) - 1;
    std::chrono::high_resolution_clock::time_point begin_build;
    std::chrono::high_resolution_clock::time_point end_build;
    std::vector<uint32_t> kmers_per_min;

    Signature() {};
    Signature(std::string filename);

    virtual std::function<bool(uint32_t, uint32_t)> compare(){return [](uint32_t a, uint_32_t b){return a < b;}};
    virtual uint32_t last_canonical_mmer_encoding(uint64_t encoded);
    virtual MMer sign(uint64_t encoded);
    virtual MMer sign(uint64_t encoded, MMer prev);
    uint32_t sign_with_extra_bit(uint64_t encoded, MMer minimizer);
    virtual std::string name() {return "default_signature";}

    std::pair<uint,uint> add_minimizers(std::string seq, std::vector<uint32_t>& occurences, Bqf_ec& kmers, uint64_t& superkmerlgth);

    void output();
};

class KMC_sign : public Signature
{
    static bool is_allowed(uint32_t mmer);

public:

    KMC_sign(std::string filename);
    std::string name() override;
    MMer sign(uint64_t encoded) override;
    MMer sign(uint64_t encoded, MMer prev) override;
};

class Roberts_sign : public Signature
{
public:

    Roberts_sign(std::string filename);

    std::string name() override;
    MMer sign(uint64_t encoded) override;
    MMer sign(uint64_t encoded, MMer prev) override;
};

class Wood_sign : public Signature
{
    uint32_t xor_mask = INDEX2_XOR_MASK & max_value;

public:

    Wood_sign(std::string filename);

    std::string name() override;
    MMer sign(uint64_t encoded) override;
    MMer sign(uint64_t encoded, MMer prev) override;
};

class Frequency_sign : public Signature
{
    std::vector<uint32_t> order;
    uint32_t most_occ;

public:
    Frequency_sign(std::string filename);

    std::string name() override;
    MMer sign(uint64_t encoded) override;
    MMer sign(uint64_t encoded, MMer prev) override;
};

#endif
