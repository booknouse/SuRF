#ifndef RANK_H_
#define RANK_H_

#include "bitvector.hpp"

#include <assert.h>

#include <vector>

#include "popcount.h"

namespace surf {

class BitvectorRank : public Bitvector {
public:
    BitvectorRank() : basic_block_size_(0), rank_lut_(nullptr) {};
    BitvectorRank(const BitvectorRank& other): Bitvector(other), basic_block_size_(other.basic_block_size_){
        position_t num_blocks = num_bits_ / basic_block_size_ + 1;
        rank_lut_ = new position_t[num_blocks];
        memmove(rank_lut_, other.rank_lut_, num_blocks*sizeof(position_t));
    }
    BitvectorRank(const position_t basic_block_size, 
		  const std::vector<std::vector<word_t> >& bitvector_per_level, 
		  const std::vector<position_t>& num_bits_per_level,
		  const level_t start_level = 0,
		  const level_t end_level = 0/* non-inclusive */) 
	: Bitvector(bitvector_per_level, num_bits_per_level, start_level, end_level) {
	basic_block_size_ = basic_block_size;
	initRankLut();
    }

    ~BitvectorRank() {}

    // Counts the number of 1's in the bitvector up to position pos.
    // pos is zero-based; count is one-based.
    // E.g., for bitvector: 100101000, rank(3) = 2
    position_t rank(position_t pos) const {
        assert(pos <= num_bits_);
        position_t word_per_basic_block = basic_block_size_ / kWordSize;
        position_t block_id = pos / basic_block_size_;
        position_t offset = pos & (basic_block_size_ - 1);
        return (rank_lut_[block_id] 
		+ popcountLinear(bits_, block_id * word_per_basic_block, offset + 1));
    }

    position_t rankLutSize() const {
	return ((num_bits_ / basic_block_size_ + 1) * sizeof(position_t));
    }

    position_t serializedSize() const {
	position_t size = sizeof(num_bits_) + sizeof(basic_block_size_) 
	    + bitsSize() + rankLutSize();
	//sizeAlign(size);
	return size;
    }

    position_t size() const {
	return (sizeof(BitvectorRank) + bitsSize() + rankLutSize());
    }

    void prefetch(position_t pos) const {
	__builtin_prefetch(bits_ + (pos / kWordSize));
	__builtin_prefetch(rank_lut_ + (pos / basic_block_size_));
    }

    void serialize(char*& dst) const {
        *reinterpret_cast<uint32_t*>(dst) = htobe32(num_bits_);
	dst += sizeof(num_bits_);
    *reinterpret_cast<uint32_t*>(dst) = htobe32(basic_block_size_);
	dst += sizeof(basic_block_size_);
    auto num_words = numWords();
    for(int i=0;i<num_words;i++) {
        *reinterpret_cast<uint64_t*>(dst) = htobe64(bits_[i]);
		dst += sizeof(uint64_t);
    }
    position_t num_blocks = num_bits_ / basic_block_size_ + 1;
    for(int i=0;i<num_blocks;i++) {
        *reinterpret_cast<uint32_t*>(dst) = htobe32(rank_lut_[i]);
        dst += sizeof(uint32_t);
    }
	//align(dst);
    }

    int deSerialize(const char*& src) {
        num_bits_ = be32toh(*reinterpret_cast<const uint32_t*>(src));
	src += sizeof(num_bits_);
    basic_block_size_=be32toh(*reinterpret_cast<const uint32_t*>(src));
	src += sizeof(basic_block_size_);
    auto num_words = numWords();
	bits_ = new word_t[num_words];
    for(int i=0;i<num_words;i++) {
        bits_[i] = be64toh(*reinterpret_cast<const uint64_t*>(src));
		src += sizeof(uint64_t);
    }
	position_t num_blocks = num_bits_ / basic_block_size_ + 1;
	rank_lut_ = new position_t[num_blocks];
    for(int i=0;i<num_blocks;i++) {
        rank_lut_[i] = be32toh(*reinterpret_cast<const uint32_t*>(src));
        src += sizeof(uint32_t);
    }
	//align(src);
	return 0;
    }

    void destroy() {
	delete[] bits_;
	delete[] rank_lut_;
    }

private:
    void initRankLut() {
        position_t word_per_basic_block = basic_block_size_ / kWordSize;
        position_t num_blocks = num_bits_ / basic_block_size_ + 1;
	rank_lut_ = new position_t[num_blocks];

        position_t cumu_rank = 0;
        for (position_t i = 0; i < num_blocks - 1; i++) {
            rank_lut_[i] = cumu_rank;
            cumu_rank += popcountLinear(bits_, i * word_per_basic_block, basic_block_size_);
        }
	rank_lut_[num_blocks - 1] = cumu_rank;
    }

    position_t basic_block_size_;
    position_t* rank_lut_; //rank look-up table
};

} // namespace surf

#endif // RANK_H_
