/**
 * @file bgen.cpp
 * @brief Expected allele counts from BGEN layout-2 genotype probability data.
 *
 * The R side (R/genotypeIo.R) reads one variant's record from the .bgen file,
 * decompresses its probability data (zlib or zstd) and hands the
 * uncompressed bytes here. The layout is BGEN 1.2 / 1.3's
 * (https://www.well.ox.ac.uk/~gav/bgen_format/spec/latest.html):
 *
 *   uint32 N, uint16 K, uint8 Pmin, uint8 Pmax,
 *   N bytes of ploidy (low 6 bits) and missingness (top bit),
 *   uint8 phased, uint8 B,
 *   the probabilities, B bits each, least significant bit first.
 *
 * For each sample, a phased record stores K - 1 probabilities per haplotype
 * (the K-th is one minus their sum); an unphased one stores one per possible
 * genotype but the last, the genotypes enumerated as allele-count vectors in
 * colex order. A stored value v is the probability v / (2^B - 1). Samples
 * flagged missing store zeros, which are skipped. A sample's ploidy is any
 * the spec allows (6 bits, so up to 63), and the expected count of an allele
 * runs from 0 to the ploidy.
 */

#include <cpp11.hpp>
#include <cstdint>
#include <vector>

using namespace cpp11;

namespace {

// Reads B-bit unsigned values from a byte buffer, least significant bit
// first, as the BGEN spec stores them.
class BitReader {
public:
	BitReader(const unsigned char* data, size_t size, int bits)
	    : data_(data), size_(size), bits_(bits), pos_(0), buffer_(0), held_(0) {}

	double next() {
		while (held_ < bits_) {
			if (pos_ >= size_) cpp11::stop("BGEN probability data ends early.");
			buffer_ |= static_cast<uint64_t>(data_[pos_++]) << held_;
			held_ += 8;
		}
		uint64_t mask = (bits_ == 64) ? ~uint64_t(0) : ((uint64_t(1) << bits_) - 1);
		uint64_t value = buffer_ & mask;
		buffer_ >>= bits_;
		held_ -= bits_;
		return static_cast<double>(value) / static_cast<double>(mask);
	}

	// Moves past `count` values without decoding them: a sample nobody asked
	// for costs some arithmetic rather than a division per probability.
	void skip(uint64_t count) {
		uint64_t skipped = count * static_cast<uint64_t>(bits_);
		if (skipped <= static_cast<uint64_t>(held_)) {
			buffer_ >>= skipped;
			held_ -= static_cast<int>(skipped);
			return;
		}
		skipped -= static_cast<uint64_t>(held_);
		buffer_ = 0;
		held_ = 0;
		pos_ += static_cast<size_t>(skipped / 8);
		const int partial = static_cast<int>(skipped % 8);
		if (pos_ > size_ || (partial > 0 && pos_ == size_)) {
			cpp11::stop("BGEN probability data ends early.");
		}
		if (partial > 0) {
			buffer_ = static_cast<uint64_t>(data_[pos_++]) >> partial;
			held_ = 8 - partial;
		}
	}

private:
	const unsigned char* data_;
	size_t size_;
	int bits_;
	size_t pos_;
	uint64_t buffer_;
	int held_;
};

uint32_t read_u32(const unsigned char* p) {
	return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
}

uint16_t read_u16(const unsigned char* p) {
	return uint16_t(p[0]) | (uint16_t(uint16_t(p[1]) << 8));
}

// The allele-count vectors of the genotypes of a sample of ploidy `ploidy`
// with `K` alleles, in the colex order the spec stores their probabilities
// in: the last allele's count varies slowest. counts[g][k] is allele k's
// count in genotype g.
void colex_genotypes(int ploidy, int K, std::vector<int>& current, int k,
                     int remaining, std::vector<std::vector<int>>& out) {
	if (k == 0) {
		current[0] = remaining;
		out.push_back(current);
		return;
	}
	for (int c = 0; c <= remaining; ++c) {
		current[k] = c;
		colex_genotypes(ploidy, K, current, k - 1, remaining - c, out);
	}
}

std::vector<std::vector<int>> genotype_counts(int ploidy, int K) {
	std::vector<std::vector<int>> out;
	std::vector<int> current(K, 0);
	colex_genotypes(ploidy, K, current, K - 1, ploidy, out);
	return out;
}

// How many genotypes a sample of ploidy `ploidy` with `K` alleles has -- the
// multisets of `ploidy` alleles, C(ploidy + K - 1, ploidy) -- counted only up
// to just past `cap`, so a corrupt record claiming an enormous ploidy or
// allele count is caught before its genotypes are enumerated.
double genotype_number(int ploidy, int K, double cap) {
	double n = 1.0;
	for (int j = 1; j <= ploidy && n <= cap; ++j) {
		n = n * (K - 1 + j) / j;
	}
	return n;
}

}  // namespace

/**
 * @brief The expected counts of the requested alleles in the requested samples.
 *
 * A sample's probabilities sit at a position only the samples before it can
 * tell, so those are walked; but only the samples asked for are decoded, and
 * the walk stops at the last of them.
 *
 * @param data One variant's uncompressed layout-2 probability data.
 * @param nSamples The file's sample count, which the data must repeat.
 * @param alleles 1-based indices of the alleles whose expected counts to
 *   return, one output column each.
 * @param samples 1-based, strictly ascending positions of the samples to
 *   return, one output row each; empty for every sample.
 * @return A samples x length(alleles) matrix; NA for a sample flagged missing
 *   or of ploidy 0.
 */
[[cpp11::register]]
cpp11::writable::doubles_matrix<> bgenDosageCpp(cpp11::raws data, int nSamples,
                                                 cpp11::integers alleles,
                                                 cpp11::integers samples) {
	const size_t size = data.size();
	const unsigned char* p = RAW(static_cast<SEXP>(data));
	if (size < 10) cpp11::stop("BGEN probability data is truncated.");
	uint32_t N = read_u32(p);
	if (static_cast<int>(N) != nSamples)
		cpp11::stop("BGEN variant stores %d samples but the file declares %d.",
		            static_cast<int>(N), nSamples);
	int K = read_u16(p + 4);
	if (size < 10 + N) cpp11::stop("BGEN probability data is truncated.");
	const unsigned char* ploidyBytes = p + 8;
	int phased = p[8 + N];
	int bits = p[9 + N];
	if (phased != 0 && phased != 1) cpp11::stop("BGEN phased flag must be 0 or 1.");
	if (bits < 1 || bits > 32) cpp11::stop("BGEN bit depth must be between 1 and 32.");

	const int nOut = alleles.size();
	for (int a = 0; a < nOut; ++a) {
		if (alleles[a] < 1 || alleles[a] > K)
			cpp11::stop("Allele %d is not one of the variant's %d alleles.", alleles[a], K);
	}
	const int nAsked = samples.size();
	const int nRows = (nAsked > 0) ? nAsked : nSamples;
	for (int s = 0; s < nAsked; ++s) {
		if (samples[s] < 1 || samples[s] > nSamples || (s > 0 && samples[s] <= samples[s - 1]))
			cpp11::stop("Samples must be ascending positions among the file's %d samples.", nSamples);
	}
	// Each ploidy's genotypes, enumerated when a sample of it is first met.
	// The data cannot hold more probabilities than it has bits, which bounds
	// the enumeration of a record that is not what it claims to be.
	std::vector<std::vector<std::vector<int>>> countsOf(64);
	const double cap = 8.0 * static_cast<double>(size) / bits + 1.0;

	cpp11::writable::doubles_matrix<> out(nRows, nOut);
	BitReader reader(p + 10 + N, size - 10 - N, bits);
	std::vector<double> stored;
	int row = 0;
	for (int i = 0; i < nSamples && row < nRows; ++i) {
		int ploidy = ploidyBytes[i] & 63;
		bool missing = (ploidyBytes[i] & 128) != 0;
		if (!phased && ploidy > 0 && countsOf[ploidy].empty()) {
			if (genotype_number(ploidy, K, cap) > cap) cpp11::stop("BGEN probability data is truncated.");
			countsOf[ploidy] = genotype_counts(ploidy, K);
		}
		// How many probabilities this sample stores, and what they are.
		int nStored = 0;
		if (ploidy > 0) {
			nStored = phased ? ploidy * (K - 1) : static_cast<int>(countsOf[ploidy].size()) - 1;
		}
		if (nAsked > 0 && samples[row] != i + 1) {
			reader.skip(static_cast<uint64_t>(nStored));
			continue;
		}
		stored.assign(nStored, 0.0);
		for (int s = 0; s < nStored; ++s) stored[s] = reader.next();
		if (missing || ploidy == 0) {
			for (int a = 0; a < nOut; ++a) out(row, a) = NA_REAL;
			++row;
			continue;
		}
		for (int a = 0; a < nOut; ++a) {
			int allele = alleles[a] - 1;
			double expected = 0.0;
			if (phased) {
				// Each haplotype carries the allele with its probability; the
				// last allele's is one minus the stored ones'.
				for (int h = 0; h < ploidy; ++h) {
					const double* hap = stored.data() + h * (K - 1);
					if (allele < K - 1) {
						expected += hap[allele];
					} else {
						double rest = 0.0;
						for (int k = 0; k < K - 1; ++k) rest += hap[k];
						expected += 1.0 - rest;
					}
				}
			} else {
				const std::vector<std::vector<int>>& g = countsOf[ploidy];
				double last = 1.0;
				for (int s = 0; s < nStored; ++s) {
					expected += g[s][allele] * stored[s];
					last -= stored[s];
				}
				expected += g[nStored][allele] * last;
			}
			out(row, a) = expected;
		}
		++row;
	}
	return out;
}

/**
 * @brief The sample ids of a BGEN sample identifier block.
 *
 * @param ids The block after its two leading counts: each id as a 2-byte
 *   length and that many bytes.
 * @param nSamples How many ids the block holds.
 */
[[cpp11::register]]
cpp11::writable::strings bgenSampleIdsCpp(cpp11::raws ids, int nSamples) {
	const unsigned char* p = RAW(static_cast<SEXP>(ids));
	const size_t size = ids.size();
	size_t at = 0;
	cpp11::writable::strings out(nSamples);
	for (int i = 0; i < nSamples; ++i) {
		if (at + 2 > size) cpp11::stop("The BGEN sample block is truncated.");
		size_t len = read_u16(p + at);
		at += 2;
		if (at + len > size) cpp11::stop("The BGEN sample block is truncated.");
		out[i] = std::string(reinterpret_cast<const char*>(p + at), len);
		at += len;
	}
	return out;
}
