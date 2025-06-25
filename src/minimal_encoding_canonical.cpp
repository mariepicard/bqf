#include <iostream>
#include <fstream>
#include <numeric>
#include <cstdint>
#include <string>
#include <bit>

static const int memory_size_kmer = 16;

/* WARNING : MINIMIZER SIZE NEEDS TO BE < 16 
    '-> IF MEMORY INTEGERS ARE 32 BITS LONG, SHIFT DOES NOT WORK  WITH VALUES OVER 32 

    TO HAVE ONE LESS BIT OF ENCODING : USE ODD VALUES OF M
 */

uint32_t encode(std::string kmer){
    uint32_t encoded = 0;
    for(char& c : kmer) {
        encoded <<= 2;
        //encoded |= ((c >> 1) & 0b11);
        switch (c) {
            case 'G':
                encoded |= 2;
                break;
            case 'T' :
                encoded |= 3;
                break;
            case 'C' :
                encoded |= 1;
                break;
            case 'A' :
                break;
            default :
                throw std::invalid_argument( "received non nucleotidic value");
                break;
        }
    }
    return encoded;
}

void display_bits(uint32_t kmer, uint k) {
    for (uint i = 0; i < 2*k; i++){
        std::cout << ((kmer>>(2*k - 1))&1);
        kmer <<= 1;
        if (((2*k - 2 - i)&0b111) == 0b111) {
            std::cout << " ";
        }
    }
    std::cout << std::endl;
}

uint32_t mask(int number_bits_right, int length) {
    std::cout << "Nb bits right : " << number_bits_right << " - Length : " << length << std::endl;
    uint32_t left = ((~0) << number_bits_right);
    uint32_t right = ~((~0) << (number_bits_right + length));
    return left & right;
}

uint32_t minimal_encoding(uint32_t canon, uint k) {
    uint32_t encoding = 0;
    uint i = 0;
    uint palindromic = 0; //0 if palindromic, 1 after a specifying pair is found
    //while no specifying pair is found
    for (; 2*i+1 < k; i++){
        uint32_t r = (~canon >> (2*i)) & 0b11;
        uint32_t l = (canon >> (2*(k - i - 1))) & 0b11;

        if (r == l) {
            encoding |= r << (2*i);
        }
        else {
            uint32_t value = (l<<2) - (l*(l+1)>>1) + r - l - 1;
            value += 0b0100;
            encoding |= (value << (2*(k - i - 2)));
            palindromic = 1;
            i++;
            break;
        }
    }
    encoding |= (mask(2*i, 2*k - 4*i) & canon)>>(2*palindromic);
    uint32_t missing = (palindromic*((1<<(2*(k - i))) - (1 << (k + (k&0b1)))))>>1;
    return encoding - missing;
}

int main(int argc, char* argv[]){
    
    if (argc < 2) {
        std::cout << "Format : ./" << argv[0] << " k [to_encode1] [to_encode2] ..." << std::endl;
        return 1;
    }
    uint k = std::stoi(argv[1]);
    for (int i = 2; i < argc; i++) {
        uint32_t encoding = minimal_encoding(encode(std::string(argv[i])), k);
        std::cout << "Final encoding : \n\t";
        display_bits(encoding, k);
        std::cout << std::endl;
    }
    return 0;
}