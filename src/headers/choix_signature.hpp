#ifndef SIGNATURE_HPP
#define SIGNATURE_HPP

#include <vector>
#include <numeric>
#include <stdint.h>
#include <string>
#include <iostream>
#include <fstream>
#include "FastxParser.hpp"
#include <algorithm>
#include <chrono>
#define M 10
#define K 32
static const uint64_t INDEX2_XOR_MASK = 0xe37e28c4271b5a2dULL;



class MMer
{
    uint32_t mmer;
    uint pos;

public:
    MMer(uint32_t mmer, uint pos) : mmer(mmer), pos(pos) {}
    uint32_t get_mmer() { return mmer; }
    uint get_position() { return pos; }
    void update_pos() { pos += 1; }
};


/*********************************
 * SIGN
 ********************************/
class Signature
{
    std::string filename;

public:
    uint32_t max_value = (1 << (2 * M)) - 1;
    std::chrono::high_resolution_clock::time_point begin_build;
    std::chrono::high_resolution_clock::time_point end_build;

    Signature() {};
    Signature(std::string filename);

    virtual MMer sign(uint64_t encoded);
    virtual MMer sign(uint64_t encoded, MMer prev);
    virtual std::string name() {return "default_signature";}

    std::pair<uint,uint> add_minimizers(std::string seq, std::vector<uint32_t>& occurences);

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