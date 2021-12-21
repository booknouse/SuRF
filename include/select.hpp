#ifndef SELECT_H_
#define SELECT_H_

#include "bitvector.hpp"

#include <assert.h>

#include <vector>

#include "config.hpp"
#include "popcount.h"

namespace surf {

class BitvectorSelect : public Bitvector {
public:
    BitvectorSelect() : sample_interval_(0), num_ones_(0), select_lut_(nullptr) {};
    BitvectorSelect(const BitvectorSelect& other):Bitvector(other), sample_interval_(other.sample_interval_), num_ones_(other.num_ones_) {
        auto select_lut_sz = selectLutSize();
        auto num_samples = select_lut_sz/sizeof(position_t);
        select_lut_ = new position_t[num_samples];
        memmove(select_lut_, other.select_lut_, num_samples*sizeof(position_t));
    }
    BitvectorSelect(const position_t sample_interval, 
		    const std::vector<std::vector<word_t> >& bitvector_per_level, 
		    const std::vector<position_t>& num_bits_per_level,
		    const level_t start_level = 0,
		    const level_t end_level = 0/* non-inclusive */) 
	: Bitvector(bitvector_per_level, num_bits_per_level, start_level, end_level) {
	sample_interval_ = sample_interval;
	initSelectLut();
    }

    ~BitvectorSelect() {}

    // Returns the postion of the rank-th 1 bit.
    // posistion is zero-based; rank is one-based.
    // E.g., for bitvector: 100101000, select(3) = 5
    position_t select(position_t rank) const {
	assert(rank > 0);
	assert(rank <= num_ones_ + 1);
	position_t lut_idx = rank / sample_interval_;
	position_t rank_left = rank % sample_interval_;
	// The first slot in select_lut_ stores the position of the first 1 bit.
	// Slot i > 0 stores the position of (i * sample_interval_)-th 1 bit
	if (lut_idx == 0)
	    rank_left--;

	position_t pos = select_lut_[lut_idx];

	if (rank_left == 0)
	    return pos;

	position_t word_id = pos / kWordSize;
	position_t offset = pos % kWordSize;
	if (offset == kWordSize - 1) {
	    word_id++;
	    offset = 0;
	} else {
	    offset++;
	}
	word_t word = bits_[word_id] << offset >> offset; //zero-out most significant bits
	position_t ones_count_in_word = popcount(word);
	while (ones_count_in_word < rank_left) {
	    word_id++;
	    word = bits_[word_id];
	    rank_left -= ones_count_in_word;
	    ones_count_in_word = popcount(word);
	}
	return (word_id * kWordSize + select64_popcount_search(word, rank_left));
    }

    position_t selectLutSize() const {
	return ((num_ones_ / sample_interval_ + 1) * sizeof(position_t));
    }

    position_t serializedSize() const {
	position_t size = sizeof(num_bits_) + sizeof(sample_interval_) + sizeof(num_ones_)
	    + bitsSize() + selectLutSize();
	//sizeAlign(size);
	return size;
    }

    position_t size() const {
	return (sizeof(BitvectorSelect) + bitsSize() + selectLutSize());
    }

    position_t numOnes() const {
	return num_ones_;
    }

    void serialize(char*& dst) const {
		*reinterpret_cast<uint32_t*>(dst) = htobe32(num_bits_);
	dst += sizeof(num_bits_);
	*reinterpret_cast<uint32_t*>(dst) = htobe32(sample_interval_);
	dst += sizeof(sample_interval_);
	*reinterpret_cast<uint32_t*>(dst) = htobe32(num_ones_);
	dst += sizeof(num_ones_);
	for(int i=0;i<numWords();i++) {
		*reinterpret_cast<uint64_t*>(dst) = htobe64(bits_[i]);
		dst += sizeof(uint64_t);
	}
	auto num_samples = num_ones_ / sample_interval_ + 1;
	for(int i=0;i<num_samples;i++){
		*reinterpret_cast<uint32_t*>(dst) = htobe32(select_lut_[i]);
		dst += sizeof(uint32_t);
	}
	//align(dst);
    }

    int deSerialize(const char*& src) {
		num_bits_ = be32toh(*reinterpret_cast<const uint32_t*>(src));
	src += sizeof(num_bits_);
	sample_interval_ = be32toh(*reinterpret_cast<const uint32_t*>(src));
	src += sizeof(sample_interval_);
	num_ones_ = be32toh(*reinterpret_cast<const uint32_t*>(src));
	src += sizeof(num_ones_);
	auto num_words = numWords();
	bits_ = new word_t[num_words];
	for(int i=0;i<num_words;i++){
		bits_[i] = be64toh(*reinterpret_cast<const uint64_t*>(src));
		src += sizeof(uint64_t);
	}
	auto num_samples = num_ones_ / sample_interval_ + 1;
	select_lut_ = new position_t[num_samples];
	for(int i=0;i<num_samples;i++){
		select_lut_[i] = be32toh(*reinterpret_cast<const uint32_t*>(src));
		src += sizeof(uint32_t);
	}
	//align(src);
	return 0;
    }

    void destroy() {
	delete[] bits_;
	delete[] select_lut_;
    }

private:
    // This function currently assumes that the first bit in the
    // bitvector is one.
    void initSelectLut() {
	position_t num_words = num_bits_ / kWordSize;
	if (num_bits_ % kWordSize != 0)
	    num_words++;

	std::vector<position_t> select_lut_vector;
	select_lut_vector.push_back(0); //ASSERT: first bit is 1
	position_t sampling_ones = sample_interval_;
	position_t cumu_ones_upto_word = 0;
	for (position_t i = 0; i < num_words; i++) {
	    position_t num_ones_in_word = popcount(bits_[i]);
	    while (sampling_ones <= (cumu_ones_upto_word + num_ones_in_word)) {
		int diff = sampling_ones - cumu_ones_upto_word;
		position_t result_pos = i * kWordSize + select64_popcount_search(bits_[i], diff);
		select_lut_vector.push_back(result_pos);
		sampling_ones += sample_interval_;
	    }
	    cumu_ones_upto_word += popcount(bits_[i]);
	}

	num_ones_ = cumu_ones_upto_word;
	position_t num_samples = select_lut_vector.size();
	select_lut_ = new position_t[num_samples];
	for (position_t i = 0; i < num_samples; i++)
	    select_lut_[i] = select_lut_vector[i];
    }


private:
    position_t sample_interval_;
    position_t num_ones_;
    position_t* select_lut_; //select look-up table
};

} // namespace surf

#endif // SELECT_H_
