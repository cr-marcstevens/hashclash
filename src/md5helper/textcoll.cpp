/**************************************************************************\
|
|    Copyright (C) 2009 Marc Stevens
|
|    This program is free software: you can redistribute it and/or modify
|    it under the terms of the GNU General Public License as published by
|    the Free Software Foundation, either version 3 of the License, or
|    (at your option) any later version.
|
|    This program is distributed in the hope that it will be useful,
|    but WITHOUT ANY WARRANTY; without even the implied warranty of
|    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
|    GNU General Public License for more details.
|
|    You should have received a copy of the GNU General Public License
|    along with this program.  If not, see <http://www.gnu.org/licenses/>.
|
\**************************************************************************/

#include <cmath>
#include <algorithm>
#include <map>
#include <stdexcept>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <random>
#include <atomic>

#define MD5DETAIL_INLINE_IMPL
#include <hashclash/saveload_gz.hpp>
#include <hashclash/md5detail.hpp>
#include <hashclash/differentialpath.hpp>
#include <hashclash/booleanfunction.hpp>
#include <hashclash/rng.hpp>
#include <hashclash/timer.hpp>
#include <hashclash/progress_display.hpp>

#include <boost/lexical_cast.hpp>

#include "main.hpp"

//#define CPUPERFORMANCE
#ifdef CPUPERFORMANCE
#include <hashclash/cpuperformance.hpp>
//uint64 cpu_step_t[64];
#define UPDATE(s) update_performance_counter __update(cpu_step_t[s]);
#ifdef __GNUC__
#include <sched.h>
#endif
#else
#define UPDATE(s)
#endif

#ifdef DOTESTCOUNTS
#define TESTCOUNT(s) ++testcounts[s];
#else
#define TESTCOUNT(s)
#endif

static inline size_t hammingweight(uint32_t x) { return __builtin_popcountl(x); }
static inline size_t hammingweight(uint64_t x) { return __builtin_popcountll(x); }

struct counter_exponential_print {
	counter_exponential_print(const std::string& _name = "counter", uint64_t _c = 0): name(_name), c(_c) {}
	const uint64_t operator()() const { return c; }
	uint64_t operator++() {
		if (hammingweight(++c)==1)
			std::cout << " " << name << "=" << c << std::flush;
		return c;
	}
	std::atomic_uint64_t c;
	std::string name;
};

std::vector< std::pair<size_t, size_t> > split_workload(size_t N, size_t jobs)
{
	std::vector< std::pair<size_t, size_t> > ret(jobs);
	if (jobs == 0)
		return ret;
	size_t frac = N/jobs, rem = N%jobs;
	for (size_t i = 0; i < jobs; ++i)
	{
		if (i == 0)
			ret[i].first = 0;
		else
			ret[i].first = ret[i-1].second;
		ret[i].second = ret[i].first + frac + (i<rem?1:0);
	}
	if (ret.back().second != N)
		throw std::runtime_error("split_wordload-inconsistency");
	return ret;
}

template<typename F>
void run_workload(size_t jobs, F&& f)
{
	std::vector<std::thread> threadgroup;
	for (size_t i = 0; i < jobs; ++i)
		threadgroup.emplace_back(std::forward<F>(f), i, jobs);
	for (size_t i = 0; i < jobs; ++i)
		threadgroup[i].join();
}

template<typename T>
void randomize_vector(std::vector<T>& in)
{
	for (size_t i = 0; i < in.size(); ++i)
	{
		size_t j = xrng64() % in.size();
		if (i == j)
			continue;
		std::swap(in[i], in[j]);
	}
}

struct localxrng {
	uint64_t x,y;
	
	localxrng()
	{
		std::random_device rd;
		x = uint64_t(rd()) + uint64_t(rd())<<32;
		y = uint64_t(rd()) + uint64_t(rd())<<32;
		x |= 1;
	}
	
	uint32_t operator()()
	{
		uint64_t t = x, s = y;
		x = y;
		t ^= t << 23;
		t ^= t >> 18;
		t ^= s ^ (s >> 5);
		y = t;
		return uint32_t((t+s)>>32);
	}
};


template<typename UInt = size_t>
UInt Binomial(UInt N, UInt K)
{
	if (K == 0 || K >= N)
		return 1;
	if (K > N/2)
		return Binomial<UInt>(N, N-K);
	return Binomial<UInt>(N-1, K-1) * N / K;
}

// range-expression for bit mask:
//    input: UInt mask
//    let S = { UInt x | (x & ~mask) == 0 }
//    then |S| = 2^hammingweight(mask)
//    range-expression: enumerate x in S in lexicographical order
// usage:
//    for (auto&& x : bit_mask_range(mask)) { loop-statement; }
template<typename UInt>
class bit_mask_range_t
{
	UInt _mask;
public:
	explicit bit_mask_range_t(UInt mask): _mask(mask) {}
	size_t count(UInt mask) const { return size_t(1)<<hammingweight(mask); }

	struct iterator {
		UInt _cur, _mask;
		iterator(UInt mask = 0) : _cur(0), _mask(mask) {}
		UInt operator*() const { return _cur; }
		// only implement correct comparison between non-end (_mask != 0) and end (_mask = 0)
		bool operator!=(const iterator& other) const { return _mask != other._mask; }
		iterator& operator++()
		{
			// compute next value in lexicographical order
			_cur = ((_cur | ~_mask) + 1) & _mask;
			// on wrap-around to 0 set _mask to 0 to indicate end
			if (_cur == 0)
				_mask = 0;
			return *this;
		}
	};	
	iterator begin() const { return iterator(_mask); }
	iterator end() const { return iterator(); }
};
template<typename UInt = uint32_t>
bit_mask_range_t<UInt> bit_mask_range(UInt mask) { return bit_mask_range_t<UInt>(mask); }

// range-expression for bit mask with minimum weight
//    input: UInt mask, unsigned minweight
//    let S = { UInt x | (x & ~mask) == 0 && hammingweight(x) >= minweight }
//    then |S| = sum_{k=minweight}^{hammingweight(mask)} Binomial(hammingweight(mask), k)
//    range-expression: enumerate x in S in lexicographical order
// usage:
//    for (auto&& x : bit_mask_minweight_range(mask, minweight)) { loop-statement; }

template<typename UInt>
class bit_mask_minweight_range_t
{
	UInt _mask;
	unsigned _minweight;
	UInt _LSBs[sizeof(UInt)*8];
public:
	explicit bit_mask_minweight_range_t(UInt mask, unsigned minweight)
		: _mask(mask), _minweight(minweight)
	{
		if (hammingweight(mask) < minweight)
		{
			_mask = _minweight = 0;
		}
		UInt LSB = 0;
		for (unsigned k = 0; ; ++k)
		{
			LSB |= UInt(1) << k;
			LSB &= _mask;
			_LSBs[hammingweight(LSB)] = LSB;
			if (LSB == _mask)
				break;
		}
	}
	size_t count() const {
		size_t c = 0, n = hammingweight(_mask);
		for (unsigned k = _minweight; k <= n; ++k)
			c += Binomial<size_t>(n, k);
		return c;
	}
	
	struct iterator {
		const bit_mask_minweight_range_t* _parent;
		UInt _cur, _mask;
		
		iterator(const bit_mask_minweight_range_t<UInt>* parent = nullptr, UInt cur = 0, UInt mask = 0) : _parent(parent), _cur(cur), _mask(mask) {}
		
		UInt operator*() const { return _cur; }
		// only implement correct comparison between non-end (_mask != 0) and end (_mask = 0)
		bool operator != (const iterator& other) const { return _mask != other._mask; }
		iterator& operator++()
		{
			// compute next value in lexicographical order
			_cur = ((_cur | ~_mask) + 1) & _mask;
			// on wrap-around to 0 set _mask to 0 to indicate end
			if (_cur == 0)
			{
				_mask = 0;
			}
			else if (hammingweight(_cur) < _parent->_minweight)
			{
				// set the least amount of LSBs to 1 to reach minweight
				for (unsigned k = _parent->_minweight - hammingweight(_cur); hammingweight(_cur) < _parent->_minweight; ++k)
					_cur |= _parent->_LSBs[k];
			}
			return *this;
		}
	};
	iterator begin() const { return iterator(this, _LSBs[_minweight], _mask); }
	iterator end() const { return iterator(); }
};
template<typename UInt = uint32_t>
bit_mask_minweight_range_t<UInt> bit_mask_minweight_range(UInt mask, unsigned minweight) { return bit_mask_minweight_range_t<UInt>(mask, minweight); }


void printword(uint32_t mt)
{
	for (unsigned b = 0; b < 4; ++b)
		std::cout << char( (mt>>(8*b))&0xFF );
}
std::string word2str(uint32_t mt)
{
	std::string ret("0123");
	for (unsigned b = 0; b < 4; ++b)
		ret[b] = char( (mt>>(8*b))&0xFF );
	return ret;
}



struct byte_alphabet {
	uint64_t byte_ok[4];
	uint64_t byte_size;
	uint8_t byte_val[256];

	byte_alphabet()
		: byte_size(0), byte_ok{0,0,0,0}
	{
	}
		
	byte_alphabet(const std::string& alphabet)
		: byte_size(0), byte_ok{0,0,0,0}
	{
		for (size_t byte : alphabet)
			byte_ok[byte/64] |= uint64_t(1)<<(byte%64);
		for (unsigned i = 0; i < 256; ++i)
		{
			if (!check(uint8_t(i)))
				continue;
			byte_val[byte_size] = uint8_t(i);
			++byte_size;
		}
	}

	uint64_t size() const { return byte_size; }

	bool check(uint8_t byte) const
	{
		return 1ULL == (1ULL & ( byte_ok[byte/64] >> (byte%64) ));
	}

	const uint8_t* begin() const { return byte_val+0; }
	const uint8_t* end() const { return byte_val+byte_size; }
};

struct message_alphabet {
	byte_alphabet bytes[64];
	uint64_t word_size[16];

	message_alphabet(const std::string& alphabet = "")
	{
		for (unsigned b = 0; b < 64; ++b)
			bytes[b] = byte_alphabet(alphabet);
		for (unsigned w = 0; w < 16; ++w)
		{
			word_size[w] = 1;
			for (unsigned b = 0; b < 4; ++b)
				word_size[w] *= bytes[w*4+b].size();
		}
	}
	
	void set_byte_alphabet(int b, const std::string& alphabet)
	{
		bytes[b] = byte_alphabet(alphabet);
		// update corresponding word_size
		unsigned w = b/4;
		word_size[w] = 1;
		for (unsigned i = 0; i < 4; ++i)
			word_size[w] *= bytes[w*4+i].size();
	}

	
	bool checkbyte(int b, uint8_t v) const
	{
		return bytes[b].check(v);
	}
	
	bool checkword(int t, uint32_t mt) const
	{
		for (unsigned b = 0; b < 4; ++b)
		{
			uint8_t v = uint8_t((mt>>(b*8)) & 0xFF);
			if (!bytes[4*t+b].check(v))
				return false;
		}
		return true;
	}
	
	uint8_t samplebyte(int b) const
	{
		return bytes[b].byte_val[ xrng64() % bytes[b].size() ];
	}

	uint8_t samplebyte(int b, localxrng& rnd) const
	{
		return bytes[b].byte_val[ rnd() % bytes[b].size() ];
	}
	
	uint32_t sampleword(int t) const
	{
		uint32 w = 0, i = 0, rnd = xrng64();

		i = rnd % bytes[4*t+0].size(); 
		w += uint32_t(bytes[4*t+0].byte_val[i]) << 0;
		rnd /= bytes[4*t+0].size();

		i = rnd % bytes[4*t+1].size();
		w += uint32_t(bytes[4*t+1].byte_val[i]) << 8;
		rnd = xrng64();

		i = rnd % bytes[4*t+2].size();
		w += uint32_t(bytes[4*t+2].byte_val[i]) << 16;
		rnd /= bytes[4*t+2].size();

		i = rnd % bytes[4*t+3].size();
		w += uint32_t(bytes[4*t+3].byte_val[i]) << 24;
		
		return w;
	}
	
	uint32_t sampleword(int t, localxrng& rng) const
	{
		uint32 w = 0, i = 0, rnd = rng();

		i = rnd % bytes[4*t+0].size(); 
		w += uint32_t(bytes[4*t+0].byte_val[i]) << 0;
		rnd /= bytes[4*t+0].size();

		i = rnd % bytes[4*t+1].size();
		w += uint32_t(bytes[4*t+1].byte_val[i]) << 8;
		rnd = rng();

		i = rnd % bytes[4*t+2].size();
		w += uint32_t(bytes[4*t+2].byte_val[i]) << 16;
		rnd /= bytes[4*t+2].size();

		i = rnd % bytes[4*t+3].size();
		w += uint32_t(bytes[4*t+3].byte_val[i]) << 24;
		
		return w;
	}
	
	struct word_range_t {
		const byte_alphabet* _bytes;
		
		explicit word_range_t(const byte_alphabet* bytes)
			: _bytes(bytes)
		{}
		
		uint64_t count() const {
			uint64_t ret = 1;
			for (unsigned b = 0; b < 4; ++b)
				ret *= _bytes[b].size();
			return ret;
		}
		
		struct iterator {
			const byte_alphabet* _bytes;
			size_t _index[4];
			
			iterator(const byte_alphabet* bytes = nullptr)
				: _bytes(bytes), _index{0,0,0,0}
			{
			}
			bool operator != (const iterator& other) const
			{
				return _bytes != other._bytes;
			}
			iterator& operator++()
			{
				int b = 0;
				while (b < 4 && ++_index[b] >= _bytes[b].size())
				{
					_index[b] = 0;
					++b;
				}
				if (b == 4)
					_bytes = nullptr;
				return *this;
			}
			uint32 operator*() const
			{
				uint32 w = 0;
				w += uint32(_bytes[0].byte_val[_index[0]]) << 0;
				w += uint32(_bytes[1].byte_val[_index[1]]) << 8;
				w += uint32(_bytes[2].byte_val[_index[2]]) << 16;
				w += uint32(_bytes[3].byte_val[_index[3]]) << 24;
				return w;
			}
		};
		
		iterator begin() const { return iterator(_bytes); }
		iterator end() const { return iterator(); }
	};
	word_range_t word_range(int t) const { return word_range_t(bytes+(4*t)); }
};

struct masked_value {
	uint32_t mask;
	uint32_t value;

	masked_value(uint32_t m = 0, uint32_t v = 0) : mask(m), value(v&m) {}
	
	bool check(uint32_t c) const { return 0 == ((c^value)&mask); }
	uint32_t sample() const { return (uint32_t(xrng64()) & ~mask) ^ value; }
	uint32_t sample(localxrng& rng) const { return (uint32_t(rng()) & ~mask) ^ value; }

	uint64_t count() const { return uint64_t(1) << hammingweight(~mask); }
	
	struct range_t {
		uint32_t _loopvalue;
		uint32_t _loopmask;
		range_t(uint32_t loopvalue = 0, uint32_t loopmask = 0)
			: _loopvalue(loopvalue), _loopmask(loopmask)
		{}
		struct iterator {
			bit_mask_range_t<uint32_t>::iterator _iterator;
			uint32_t _loopvalue;
			iterator(uint32_t loopvalue = 0, uint32_t loopmask = 0) : _iterator(loopmask), _loopvalue(loopvalue) {}
			bool operator != (const iterator& other) const { return _iterator != other._iterator; }
			iterator& operator++() { ++_iterator; return *this; }
			uint32_t operator*() const { return *_iterator ^ _loopvalue; }
		};
		iterator begin() const { return iterator(_loopvalue, _loopmask); }
		iterator end() const { return iterator(); }
	};
	range_t range()     const { return range_t(value, ~mask); }
	range_t range_rnd() const { return range_t(value ^ (xrng64() & ~mask), ~mask); }
	range_t range_rnd(localxrng& rng) const { return range_t(value ^ (rng() & ~mask), ~mask); }
};


template<size_t N>
struct md5state_t {
	static const int offset = 3;
	uint32 Q[N];
	uint32 m[16];
	
	md5state_t() = default;
	md5state_t(const md5state_t<N>&) = default;
	
	template<size_t N2>
	md5state_t(const md5state_t<N2>& s)
	{
		for (unsigned i = 0; i < 16; ++i)
			m[i] = s.m[i];
		for (unsigned i = 0; i < N && i < N2; ++i)
			Q[i] = s.Q[i];
	}

	uint32& Qt(int t) { return Q[offset+t]; }
	const uint32& Qt(int t) const { return Q[offset+t]; }
	uint32& Wt(int t) { return m[md5_wt[t]]; }
	const uint32& Wt(int t) const { return m[md5_wt[t]]; }

	void computeQtp1(int t)
	{
		Qt(t+1) = md5_step(t, Qt(t), Qt(t-1), Qt(t-2), Qt(t-3), Wt(t));
	}
	void computeQtm3(int t)
	{
		Qt(t-3) = md5_step_bw(t, Qt(t+1), Qt(t), Qt(t-1), Qt(t-2), Wt(t));
	}
	void computeWt(int t)
	{
		Wt(t) = md5_step_bw(t, Qt(t+1), Qt(t), Qt(t-1), Qt(t-2), Qt(t-3));
	}

};

#ifndef NOSERIALIZATION
namespace boost {
        namespace serialization {

                template<class Archive, size_t N>
                void serialize(Archive& ar, ::md5state_t<N>& p, const unsigned int file_version)
                {
                	size_t n = N;
                        ar & make_nvp("N", n);
                        if (n != N) throw std::runtime_error("md5state_t: N mismatch in loading");
                        for (size_t i = 0; i < 16; ++i)
	                        ar & make_nvp("m", p.m[i]);
                        for (size_t i = 0; i < N; ++i)
	                        ar & make_nvp("Q", p.Q[i]);
                }

        }
}
#endif // NOSERIALIZATION


struct textcoll_solver_t {
	textcoll_solver_t(): testcounts(1<<20,0), threads(40) {}
	
	uint64 cpu_step_t[64];
	vector<uint64> testcounts;
	differentialpath diffpath;
	uint32 m_diff[16];
	unsigned threads;

	static const int offset = 3;

	uint32 dQ[68];
	uint32 dT[68];
	uint32 dR[68];

	uint32 dQt(int t) const { return dQ[offset+t]; }
	uint32 dTt(int t) const { return dT[offset+t]; }
	uint32 dRt(int t) const { return dR[offset+t]; }

	uint32 Qvaluemask[68];
	uint32 Qvalue[68];
	uint32 Qprev[68];
	uint32 Qprev2[68];
	
	size_t prefixblocks;
	uint32 ihv1[4];
	uint32 ihv2[4];
	
	uint32 Qtvaluemask(int t) const { return Qvaluemask[offset+t]; }
	uint32 Qtvalue(int t) const { return Qvalue[offset+t]; }
	uint32 Qtprev(int t) const { return Qprev[offset+t]; }
	uint32 Qtprev2(int t) const { return Qprev2[offset+t]; }

	message_alphabet MA;
	
	void filltables();
	void fillalphabet(const std::string& alphabet)
	{
		MA = message_alphabet(alphabet);
	}
	void findcollision(const vector<differentialpath>& diffpaths, const std::string& alphabet)
	{
		fillalphabet(alphabet);

		diffpath = diffpaths.front();
		filltables();

		if (ihv2[0] == ihv1[0])
		{
			std::cout << "Start block 1" << std::endl;
			start();
		} else
		{
			std::cout << "Start block 2" << std::endl;
			start_block2();
		}
	}

	typedef md5state_t<32> halfstate_t;
	typedef md5state_t<68> fullstate_t;
	typedef vector< halfstate_t > vec_halfstate_t;
	
	void start();
	void start_block2();

	void generateQ13Q17(vec_halfstate_t& out, uint64_t N);
	void extend_step_fw(int t, vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N);
	void extend_step_bw(int t, vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N);
	void extend_step_m11(vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N);
	void extend_step_m10(vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N);


	template<size_t N>
	unsigned Q9m9tunnel(const md5state_t<N>& S) const
	{
		return hammingweight(~Qvaluemask[offset+9] & ~Qprev[offset+10] & S.Qt(11) & ~S.Qt(10));
	}
	
	template<size_t N>
	bool checkrotationQtQtp1(int t, const md5state_t<N>& s) const
	{
		uint32_t R1 = s.Qt(t+1) - s.Qt(t);
		uint32_t R2 = R1 + dRt(t);
		uint32_t T1 = rotate_right(R1, md5_rc[t]);
		uint32_t T2 = rotate_right(R2, md5_rc[t]);
		return T2-T1 == dTt(t);
	}

	template<size_t N>
	masked_value masked_value_Qt(int t, const md5state_t<N>& s) const
	{
		// ignore Qprev(t) (included by default) and Qprev(t+1) (not included by default)
		uint32_t mask = Qtvaluemask(t) & ~Qtprev(t);
		return masked_value(mask, Qtvalue(t) & mask);
	}
	template<size_t N>
	masked_value masked_value_QtQtm1(int t, const md5state_t<N>& s) const
	{
		// use Qprev(t) included by default
		return masked_value(Qtvaluemask(t), Qtvalue(t) ^ (s.Qt(t-1) & Qtprev(t)));
	}
	template<size_t N>
	masked_value masked_value_QtQtp1(int t, const md5state_t<N>& s) const
	{
		// ignore Qprev(t) included by default
		uint32_t mask = Qtvaluemask(t) & ~Qtprev(t);
		uint32_t value = Qtvalue(t) & mask;
		// use Qprev(t+1) not included by default
		mask |= Qtprev(t+1);
		value ^= (s.Qt(t+1) ^ Qtvalue(t+1)) & Qtprev(t+1);
		return masked_value(mask, value);
	}

	template<size_t N>
	bool checkQt(int t, const md5state_t<N>& s) const
	{
		return masked_value_Qt(t, s).check(s.Qt(t));
	}
	template<size_t N>
	bool checkQtQtm1(int t, const md5state_t<N>& s) const
	{
		return masked_value_QtQtm1(t, s).check(s.Qt(t));
	}
	template<size_t N>
	bool checkQtQtp1(int t, const md5state_t<N>& s) const
	{
		return masked_value_QtQtp1(t, s).check(s.Qt(t));
	}

	void completeQ7Q24(const halfstate_t& Q7Q24state);
	void compute_good_m10(const fullstate_t& Sin);

	void check_solution(const fullstate_t& sol);

//	std::vector< std::pair<uint32,uint32> > m14Q11;
//	std::vector< std::pair<uint32,uint32> > m14m13;
//	std::vector< std::tuple<uint32,uint32,uint32,uint32> > Q7m12m13m14;

	std::vector<uint32_t> m4rndrange, good_m10;
	
	std::vector< std::array<uint32_t,4> > vecQ7m10m12m13;
	std::unordered_map<uint32_t,size_t> Q7ptr;

	std::mutex mut;
	typedef std::lock_guard<std::mutex> lock_t;
	
};


void textcoll_solver_t::start()
{
	static const uint64_t N = 1ULL<<22;
	vector< halfstate_t > Q7Q24states;
	
	std::cout << "Trying to load 'Q7Q24.bin.gz'..." << std::flush;
	try {
		load_gz(Q7Q24states, "Q7Q24", binary_archive);
		std::cout << "done: " << Q7Q24states.size() << "." << std::endl;
	} catch (std::exception& e)
	{
		Q7Q24states.clear();
		std::cout << "failed." << std::endl;
	}

	// generate partial solutions over Q7 ... Q24
	if (Q7Q24states.empty())
	{
		vector< halfstate_t > in;
		vector< halfstate_t >& out = Q7Q24states;

/*
		MA.set_byte_alphabet(0, "B");
		MA.set_byte_alphabet(1, "A");
		MA.set_byte_alphabet(2, "S");
		MA.set_byte_alphabet(3, "E");
		MA.set_byte_alphabet(4, "X");
		MA.set_byte_alphabet(5, "6");
		MA.set_byte_alphabet(6, "4");
		MA.set_byte_alphabet(7, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ+/");
/*
//		MA.set_byte_alphabet(16, "0123456789");
//		MA.set_byte_alphabet(17, "0123456789");
//		MA.set_byte_alphabet(18, "0123456789");
//		MA.set_byte_alphabet(19, "0123456789");
*/
/*		MA.set_byte_alphabet(20, "bB");
		MA.set_byte_alphabet(21, "eEiI"); // diff: +4
		MA.set_byte_alphabet(22, "eE");
		MA.set_byte_alphabet(23, "rR");

		MA.set_byte_alphabet(60, "wW");
		MA.set_byte_alphabet(61, "hH");
		MA.set_byte_alphabet(62, "aA");
		MA.set_byte_alphabet(63, "tT");
*/

		// Q17-Q13,m1
		generateQ13Q17(out, N);

		// m15 Q12 ... m12 Q9
		for (int t = 15; t >= 12; --t)
		{
			std::swap(in, out);
			randomize_vector(in);
			extend_step_bw(t, out, in, N);
		}

		std::swap(in, out);
		randomize_vector(in);
		extend_step_fw(17, out, in, N); // m6

		std::swap(in, out);
		randomize_vector(in);
		extend_step_m11(out, in, N);    // m11

		std::swap(in, out);
		randomize_vector(in);
		extend_step_fw(19, out, in, N); // m0

		std::swap(in, out);
		randomize_vector(in);
		extend_step_fw(20, out, in, N); // m5
		for (size_t i = 0; i < out.size(); )
		{
			if (MA.checkword(5, out[i].m[5]) && MA.checkword(5, out[i].m[5] + m_diff[5]))
			{
				++i;
				continue;
			}
			std::swap( out[i], out.back() );
			out.pop_back();
		}
		std::cout << "out.size() after m5+m_diff filtering: " << out.size() << std::endl;

		std::swap(in, out);
		randomize_vector(in);
		extend_step_m10(out, in, N); // m10, m15

		std::swap(in, out);
		randomize_vector(in);
		extend_step_fw(23, out, in, N); // m4
		
		std::cout << "Saving 'Q7Q24.bin.gz'..." << std::flush;
		save_gz(out, "Q7Q24", binary_archive);
		
		std::map<char, size_t> charcount;
		for (int wt : { 0, 1, 4, 5, 6, 10, 11, 12, 13, 14, 15 })
		{
			for (unsigned b = 0; b < 4; ++b)
			{
				charcount.clear();
				for (auto& S: out)
					++charcount[ char((S.m[wt]>>(8*b))&0xFF) ];
				std::cout << "m" << wt << "[" << b << "]=msg[" << (wt*4+b) << "]: ";
				std::vector< std::pair<char,size_t> > charcountsorted(charcount.begin(), charcount.end());
				std::sort(charcountsorted.begin(), charcountsorted.end(), [](auto& l, auto& r){ return l.second > r.second; });
				for (auto& cc : charcountsorted)
					std::cout << cc.first;
				std::cout << std::endl;
				std::cout << "m" << wt << "[" << b << "]=msg[" << (wt*4+b) << "]: ";
				std::cout.precision(4);
				for (auto& cc : charcountsorted)
					std::cout << " " << cc.first << ":" << double(cc.second)/double(out.size());
				std::cout << std::endl << std::endl;
			}
		}
		exit(0);
	}

	// sort on maximum Q9m9 tunnel	
	std::sort(Q7Q24states.begin(), Q7Q24states.end(), 
		[this](const halfstate_t& l, const halfstate_t& r)
		{
			return Q9m9tunnel(l) > Q9m9tunnel(r);
		});
	std::cout << "Q9m9 tunnel best state: " << Q9m9tunnel(Q7Q24states.front()) << std::endl;

	size_t skipstates = 0;
	for (auto& Q7Q24state : Q7Q24states)
	{
		if (skipstates > 0)
		{
			--skipstates;
			continue;
		}
		completeQ7Q24(Q7Q24state);
	}
}



void textcoll_solver_t::generateQ13Q17(vec_halfstate_t& out, uint64_t N)
{
	out.resize(N);
	auto job_range = split_workload(N, threads);
	
	std::cout << "\n======\nGenerate Q13-Q17,m1:";
	// to avoid sync we just let progress_display be maintained by the first thread
	progress_display pd(job_range.front().second);
	
	uint64_t attempts = 0;
	run_workload(threads, 
		[this,&job_range,&attempts,&out,&pd](size_t ji, size_t jn)
		{
			auto begin = job_range[ji].first, end = job_range[ji].second;
			uint64_t myattempts = 0;
			localxrng rng;

			auto m1range = MA.word_range(1);
			const size_t m1rangecount = m1range.count(), Q17count = masked_value_QtQtm1(17, halfstate_t()).count();
			for (auto k = begin; k < end; )
			{
				++myattempts;
				auto& S = out[k];
				S.Qt(13) = masked_value_Qt(13,S).sample(rng);
				S.Qt(14) = masked_value_QtQtm1(14,S).sample(rng);
				S.Qt(15) = masked_value_QtQtm1(15,S).sample(rng);
				S.Qt(16) = masked_value_QtQtm1(16,S).sample(rng);
				if (!checkrotationQtQtp1(13, S)) continue;
				if (!checkrotationQtQtp1(14, S)) continue;
				if (!checkrotationQtQtp1(15, S)) continue;

				auto Q17 = masked_value_QtQtm1(17,S);
				if (Q17count < m1rangecount)
				{
					S.Qt(17) = Q17.sample(rng);
					S.computeWt(16);
					if (!MA.checkword(1, S.m[1]))
						continue;
				} else {
					S.m[1] = MA.sampleword(1, rng);
					S.computeQtp1(16);
					if (!Q17.check(S.Qt(17)))
					continue;
				}
				if (!checkrotationQtQtp1(16, S)) continue;

				// S is good
				++k;
				if (ji == 0)
					++pd;
			}
			{ lock_t lock(mut); attempts += myattempts; }
		});

	std::cout << "Generate Q13-Q17,m1: attempts: 2^" << log(double(attempts))/log(2.0) << " success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}



void textcoll_solver_t::extend_step_fw(int t, vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N)
{
	out.resize(N);
	auto job_range = split_workload(N, threads);
	
	unsigned wt = md5_wt[t];
	std::cout << "\n======\nStep " << t << ": (Q" << (t+1) << ",m" << wt << "):";
	// to avoid sync we just let progress_display be maintained by the first thread
	progress_display pd(job_range.front().second);
	
	uint64_t attempts = 0;
	run_workload(threads, 
		[this,&job_range,&attempts,&out,&in,&pd,t,wt](size_t ji, size_t jn)
		{
			auto begin = job_range[ji].first, end = job_range[ji].second;
			uint64_t myattempts = 0;
			localxrng rng;

			auto wtrange = MA.word_range(wt);
			const size_t wtrangecount = wtrange.count(), Qtp1count = masked_value_QtQtm1(t+1, halfstate_t()).count();
			auto inptr = begin % in.size(); // process in elements in round-robin, starting at same offset as out
			for (auto k = begin; k < end; )
			{
				++myattempts;
				auto& S = out[k];
				S = in[inptr];
				if (++inptr >= in.size()) inptr -= in.size();

				auto Qtp1mv = masked_value_QtQtm1(t+1,S);
				if (Qtp1count < wtrangecount)
				{
					S.Qt(t+1) = Qtp1mv.sample(rng);
					S.computeWt(t);
					if (!MA.checkword(wt, S.m[wt]))
						continue;
				} else {
					S.m[wt] = MA.sampleword(wt, rng);
					S.computeQtp1(t);
					if (!Qtp1mv.check(S.Qt(t+1)))
						continue;
				}
				if (!checkrotationQtQtp1(t, S))
					continue;

				// S is good
				++k;
				if (ji == 0)
					++pd;
			}
			{ lock_t lock(mut); attempts += myattempts; }
		});
	std::cout << "Step " << t << ": attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}



void textcoll_solver_t::extend_step_bw(int t, vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N)
{
	out.resize(N);
	auto job_range = split_workload(N, threads);
	
	unsigned wt = md5_wt[t];
	std::cout << "\n======\nStep " << t << ": (Q" << (t-3) << ",m" << wt << "):";
	// to avoid sync we just let progress_display be maintained by the first thread
	progress_display pd(job_range.front().second);
	
	uint64_t attempts = 0;
	run_workload(threads, 
		[this,&job_range,&attempts,&out,&in,&pd,t,wt](size_t ji, size_t jn)
		{
			auto begin = job_range[ji].first, end = job_range[ji].second;
			uint64_t myattempts = 0;
			localxrng rng;

			auto wtrange = MA.word_range(wt);
			const size_t wtrangecount = wtrange.count(), Qtm3count = masked_value_QtQtp1(t-3, halfstate_t()).count();
			auto inptr = begin % in.size(); // process in elements in round-robin, starting at same offset as out
			for (auto k = begin; k < end; )
			{
				++myattempts;
				auto& S = out[k];
				S = in[inptr];
				if (++inptr >= in.size()) inptr -= in.size();

				auto Qtm3mv = masked_value_QtQtp1(t-3,S);
				if (Qtm3count < wtrangecount)
				{
					S.Qt(t-3) = Qtm3mv.sample(rng);
					S.computeWt(t);
					if (!MA.checkword(wt, S.m[wt]))
						continue;
				} else {
					S.m[wt] = MA.sampleword(wt, rng);
					S.computeQtm3(t);
					if (!Qtm3mv.check(S.Qt(t-3)))
						continue;
				}
				if (!checkrotationQtQtp1(t-3, S))
					continue;

				// S is good
				++k;
				if (ji == 0)
					++pd;
			}
			{ lock_t lock(mut); attempts += myattempts; }
		});
	std::cout << "Step " << t << ": attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}



void textcoll_solver_t::extend_step_m11(vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N)
{
	out.resize(N);
	auto job_range = split_workload(N, threads);
	
	std::cout << "\n======\nStep: (m11,Q8,Q19):";
	// to avoid sync we just let progress_display be maintained by the first thread
	progress_display pd(job_range.front().second);
	
	uint64_t attempts = 0;
	run_workload(threads, 
		[this,&job_range,&attempts,&out,&in,&pd](size_t ji, size_t jn)
		{
			auto begin = job_range[ji].first, end = job_range[ji].second;
			uint64_t myattempts = 0;
			localxrng rng;

			auto m11range = MA.word_range(11);
			auto inptr = begin % in.size(); // process in elements in round-robin, starting at same offset as out
			for (auto k = begin; k < end; )
			{
				++myattempts;
				auto& S = out[k];
				S = in[inptr];
				if (++inptr >= in.size()) inptr -= in.size();

				S.m[11] = MA.sampleword(11, rng);
				S.computeQtp1(18);
				auto Q19mv = masked_value_QtQtm1(19, S);
				if (!Q19mv.check(S.Qt(19)))
					continue;
				S.computeQtm3(11);
				auto Q8mv = masked_value_QtQtp1(8, S);
				if (!Q8mv.check(S.Qt(8)))
					continue;
				if (!checkrotationQtQtp1(8, S))
					continue;
				if (!checkrotationQtQtp1(18, S))
					continue;

				// S is good
				++k;
				if (ji == 0)
					++pd;
			}
			{ lock_t lock(mut); attempts += myattempts; }
		});
	std::cout << "Step m11: attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}



void textcoll_solver_t::extend_step_m10(vec_halfstate_t& out, const vec_halfstate_t& in, uint64_t N)
{
	out.resize(N);
	auto job_range = split_workload(N, threads);
	
	std::cout << "\n======\nStep: (m10,Q7,Q22,Q23):";
	// to avoid sync we just let progress_display be maintained by the first thread
	progress_display pd(job_range.front().second);
	
	uint64_t attempts = 0;
	run_workload(threads, 
		[this,&job_range,&attempts,&out,&in,&pd](size_t ji, size_t jn)
		{
			auto begin = job_range[ji].first, end = job_range[ji].second;
			uint64_t myattempts = 0;
			localxrng rng;

			auto m11range = MA.word_range(11);
			auto inptr = begin % in.size(); // process in elements in round-robin, starting at same offset as out
			for (auto k = begin; k < end; )
			{
				++myattempts;
				auto& S = out[k];
				S = in[inptr];
				if (++inptr >= in.size()) inptr -= in.size();

				S.m[10] = MA.sampleword(10, rng);
				
				S.computeQtm3(10);
				auto Q7mv = masked_value_QtQtp1(7,S);
				if (!Q7mv.check(S.Qt(7)))
					continue;

				S.computeQtp1(21); // m10
				auto Q22mv = masked_value_QtQtm1(22,S);
				if (!Q22mv.check(S.Qt(22)))
					continue;

				S.computeQtp1(22); // m15
				auto Q23mv = masked_value_QtQtm1(23,S);
				if (!Q23mv.check(S.Qt(23)))
					continue;
				if (!checkrotationQtQtp1(7,S))
					continue;
				if (!checkrotationQtQtp1(21,S))
					continue;
				if (!checkrotationQtQtp1(22,S))
					continue;

				// S is good
				++k;
				if (ji == 0)
					++pd;
			}
			{ lock_t lock(mut); attempts += myattempts; }
		});
	std::cout << "Step m10: attempts: 2^" << log(double(attempts))/log(2.0) << ", success rate: 2^" << log(double(N)/double(attempts))/log(2.0) << std::endl;
}

void textcoll_solver_t::completeQ7Q24(const halfstate_t& Q7Q24state)
{
	fullstate_t S = Q7Q24state;

	// sanity check S
	for (int t = 10; t <= 23; ++t)
	{
		unsigned wt = md5_wt[t];
		uint32 Wtorg = S.m[wt];
		S.computeWt(t);
		if (S.m[wt] != Wtorg)
			throw std::runtime_error("Q7Q24state: incorrect Wt");
		if (!MA.checkword(wt, S.m[wt]))
			throw std::runtime_error("Q7Q24state: Wt alphabet violation");
	}

	// iterate over all m4 values in random order
	if (m4rndrange.empty())
	{
		auto m4range = MA.word_range(4);
		m4rndrange.reserve(m4range.count());
		for (uint32_t m4 : m4range)
			m4rndrange.emplace_back(m4);
	}
	randomize_vector(m4rndrange);

	size_t m4attempts = 0, m4ok = 0;
	for (auto& m4 :	m4rndrange)
	{
		if (++m4attempts > 32 && m4ok == 0)
			return;
		S.m[4] = m4;
		compute_good_m10(S);
		std::cout << "m4=" << word2str(m4) << ": good_m10 size: " << good_m10.size() << std::endl;
		if (good_m10.size() < 6000) //(1ULL<<14))
			continue;
		++m4ok;
		
		std::cout << "MSG: ";
		for (unsigned t = 0; t < 16; ++t)
			for (unsigned b = 0; b < 4; ++b)
				std::cout << char( (S.m[t]>>(b*8))&0xFF);
		std::cout << std::endl;

		uint32_t m13org = S.m[13], Q10org = S.Qt(10);

		std::vector< std::pair<uint32_t,uint32_t> > m13Q10good;
		for (uint32_t m13 : MA.word_range(13))
		{
			S.m[13] = m13;
			S.computeQtm3(13);
			if (!masked_value_QtQtp1(10,S).check(S.Qt(10)))
				continue;
			if (!checkrotationQtQtp1(10,S))
				continue;
			m13Q10good.emplace_back(S.m[13], S.Qt(10));
		}
		std::cout << "m13Q10good size: " << m13Q10good.size() << std::endl;
		uint32_t Q9m9 = ~Qvaluemask[offset+9] & ~Qprev[offset+10] & S.Qt(11);
		std::sort(m13Q10good.begin(), m13Q10good.end(), 
			[Q9m9](auto& l, auto& r)
			{
				return hammingweight(Q9m9 & ~l.second)
					> hammingweight(Q9m9 & ~r.second);
			});

		// find m12 that still satisfies Q9 and Q8 conditions
		// forget m12 for which there is already another m12 with same Q8 and Q7
		// i.e. the difference is part of the Q9m9 tunnel

		std::map< std::array<uint32_t, 3>, std::array< uint32_t, 2 > > Q7810m1213;

		progress_display pdQ7810(m13Q10good.size());
		std::atomic_size_t m13readptr(0), m13writeptr(0);

		run_workload(threads, 
			[this,&m13Q10good,&pdQ7810,&m13readptr,&m13writeptr,&S,&Q7810m1213]
			(size_t ji, size_t jn)
			{
				fullstate_t myS = S;
				auto m12range = MA.word_range(12);
				std::map< std::array<uint32_t, 3>, std::array< uint32_t, 2 > > tmp;
				bool stop = false;
				while (!stop)
				{
					
					size_t m13ptr = m13readptr.fetch_add(1);
					if (m13ptr >= m13Q10good.size())
						break;

					auto m13Q10 = m13Q10good[m13ptr];
					myS.m[13] = m13Q10.first;
					myS.Qt(10) = m13Q10.second;

					// compute all results for m12range and store in tmp					
					for (uint32_t m12 : m12range)
					{
						myS.m[12] = m12;
						myS.m[10] = 0; // fixed dummy value

						myS.computeQtm3(12);
						if (!masked_value_QtQtp1(9,myS).check(myS.Qt(9)))
							continue;

						myS.computeQtm3(11);
						if (!masked_value_QtQtp1(8,myS).check(myS.Qt(8)))
							continue;
					
						myS.computeQtm3(10); // virtual Q7-under-fixed-dummy-m10

						if (!checkrotationQtQtp1(9,myS))
							continue;
						if (!checkrotationQtQtp1(8,myS))
							continue;

						tmp.insert( std::make_pair( std::array<uint32_t,3>({myS.Qt(7),myS.Qt(8),myS.Qt(10)}), std::array<uint32_t,2>({ myS.m[12], myS.m[13] })));
					}

					// wait till ready to write
					while (m13ptr != m13writeptr)
						std::this_thread::yield();
					// write results in tmp to Q7810m1213
					{
						lock_t lock(mut);
						++pdQ7810;
						for (auto& kv : tmp)
						{
							auto itf = Q7810m1213.insert(kv);
							if (itf.second && hammingweight(Q7810m1213.size())==1) std::cout << Q7810m1213.size() << " " << std::flush;

						}
						if (Q7810m1213.size() >= 1ULL<<20)
							stop = true;
					}
					// finished writing, increase write ptr for next thread
					++m13writeptr;
					
					tmp.clear();
				}	
			});
		std::cout << "\nQ7810m1213size: " << Q7810m1213.size() << std::endl;
		

		// iterate over those m12, and find all (m10, m12) satisfying Q7-Q9, Q22-Q23
		static const size_t maxLUT = 1ULL<<30;
		
		progress_display pdm12m13(Q7810m1213.size());


		vecQ7m10m12m13.clear();
		vecQ7m10m12m13.reserve(maxLUT);
		
		auto m1213it = Q7810m1213.begin();
		run_workload(threads, 
			[this,&S,&Q7810m1213,&m1213it,&pdm12m13]
			(size_t ji, size_t jn)
			{
				vector< std::array<uint32_t,4> > tmp;
				tmp.reserve(1ULL<<20);
				fullstate_t myS = S;
				while (true)
				{
					auto it = Q7810m1213.begin();
					// grab next m1213it
					{
						lock_t lock(mut);
						it = m1213it;
						++m1213it;
						++pdm12m13;
						if (vecQ7m10m12m13.size() == vecQ7m10m12m13.capacity())
							break;
					}
					if (it == Q7810m1213.end())
						break;

					myS.m[12] = it->second[0];
					myS.m[13] = it->second[1];
					myS.computeQtm3(13);
					myS.computeQtm3(12);
					myS.computeQtm3(11);
					for (uint32_t m10 : good_m10)
					{
						myS.m[10] = m10;

						myS.computeQtm3(10); // m10
						if (!masked_value_QtQtp1(7,myS).check(myS.Qt(7)))
							continue;

						myS.computeQtp1(21); // m10
						if (!masked_value_QtQtm1(22,myS).check(myS.Qt(22)))
							continue;

						myS.computeQtp1(22); // m15
						if (!masked_value_QtQtm1(23,myS).check(myS.Qt(23)))
							continue;

						myS.computeQtp1(23); // m4
						if (!masked_value_QtQtm1(24,myS).check(myS.Qt(24)))
							continue;

						if (!checkrotationQtQtp1(7,myS))
							continue;
						if (!checkrotationQtQtp1(21,myS))
							continue;
						if (!checkrotationQtQtp1(22,myS))
							continue;
						if (!checkrotationQtQtp1(23,myS))
							continue;
					
						tmp.push_back(std::array<uint32_t,4>({myS.Qt(7), myS.m[10], myS.m[12], myS.m[13]}));
						if (tmp.size() >= (1ULL<<20))
						{
							lock_t lock(mut);
							for (auto& v : tmp)
							{
								if (vecQ7m10m12m13.size() == vecQ7m10m12m13.capacity())
									break;
								vecQ7m10m12m13.emplace_back(v);
							}
							tmp.clear();
						}
					}
				}
				lock_t lock(mut);
				for (auto& v : tmp)
				{
					if (vecQ7m10m12m13.size() == vecQ7m10m12m13.capacity())
						break;
					vecQ7m10m12m13.emplace_back(v);
				}
			});
		std::cout << "Converting into look-up table..." << std::flush;
#if 1
		std::sort(vecQ7m10m12m13.begin(), vecQ7m10m12m13.end(), [](auto& l, auto& r) { return l[0] < r[0]; });
		for (size_t i = 0; i < vecQ7m10m12m13.size(); )
		{
			Q7ptr[ vecQ7m10m12m13[i][0] ] = i;
			size_t j = i+1;
			while (vecQ7m10m12m13[j][0] == vecQ7m10m12m13[i][0])
				++j;
			i = j;
		}
		std::cout << "vecQ7m10m12m13size: " << vecQ7m10m12m13.size() << " Q7ptrsize: " << Q7ptr.size() << std::endl;
#else
		// Q7LUT: Q7 => (beginidx,endidx)
		std::unordered_map<uint32_t, std::pair<size_t,size_t> > Q7LUT;
		// first count each Q7
		std::cout << "counting..." << std::flush;
		progress_display pd(vecQ7m10m12m13.size());
		for (auto& v : vecQ7m10m12m13)
			++Q7LUT[v[0]].second, ++pd;
		std::cout << "indexing..." << std::flush;
		size_t lastidx = 0;
		for (auto it = Q7LUT.begin(), itend = Q7LUT.end(); it != itend; ++it)
		{
			// first set proper indices
			it->second.first = lastidx;
			it->second.second += lastidx;
			lastidx = it->second.second;
			// now set end to begin for reorder process
			it->second.second = it->second.first;
		}
		std::cout << "reordering..." << std::flush;
		progress_display pd2(vecQ7m10m12m13.size());
		for (size_t i = 0; i < vecQ7m10m12m13.size(); )
		{
			auto& Q7entry = Q7LUT[ vecQ7m10m12m13[i][0] ];
			// check if i has already been processed
			if (Q7entry.first <= i && i < Q7entry.second)
			{
				++i;
				continue;
			}
			++pd2;
			size_t j = Q7entry.second;
			if (j != i)
				std::swap( vecQ7m10m12m13[i], vecQ7m10m12m13[j] );
			++Q7entry.second;
		}
		std::cout << "done." << std::endl;
//		if (Q7LUT.end()->second != vecQ7m10m12m13.size())
//			throw std::runtime_error("error");		
		std::cout << "vecQ7m10m12m13size: " << vecQ7m10m12m13.size() << " Q7LUTsize: " << Q7LUT.size() << std::endl;
#endif		
		
//		exit(0);




		/* MAIN LOOP */
		S.Qt(-3) = ihv1[0]; S.Qt(-2) = ihv1[3]; S.Qt(-1) = ihv1[2]; S.Qt(0) = ihv1[1];
		
		S.computeQtp1(0); // m0 is fixed
		if (!masked_value_QtQtm1(1,S).check(S.Qt(1)))
			throw std::runtime_error("Q1-inconsistency");
		if (!checkrotationQtQtp1(0, S))
			throw std::runtime_error("Q1-inconsistency");
		S.computeQtp1(1); // m1 is fixed
		if (!masked_value_QtQtm1(2,S).check(S.Qt(2)))
			throw std::runtime_error("Q2-inconsistency");
		if (!checkrotationQtQtp1(1, S))
			throw std::runtime_error("Q2-inconsistency");
		
		counter_exponential_print m2m3cnt("m2m3cnt"), Q7attempts("Q7attempts"), Q7match("Q7match"), Q7success("Q7success"), 
					  m7ok("m7ok"), m8ok("m8ok"), m9ok("m9ok"), m12ok("m12ok"), Q24ok("Q24ok");

		std::vector<uint32_t> m2rndrange;
		for (uint32_t m2 : MA.word_range(2))
			m2rndrange.emplace_back(m2);
		randomize_vector(m2rndrange);
	
		std::atomic_size_t m2index(0);

		fullstate_t Scopy = S;

run_workload(threads, 
	[this,&Scopy,&m2index,&m2rndrange,&Q24ok]
	(size_t ji, size_t jn)
	{
		auto S = Scopy;
		auto m3range = MA.word_range(3);
		while (true)
		{
			size_t m2i = m2index.fetch_add(1);
			if (m2i >= m2rndrange.size())
				break;

			S.m[2] = m2rndrange[m2i];

			S.computeQtp1(2);
			if (!masked_value_QtQtm1(3,S).check(S.Qt(3)))
				continue;
			if (!checkrotationQtQtp1(2,S))
				continue;

			for (uint32_t m3 : m3range)
			{
//				++m2m3cnt;

				S.m[3] = m3;
				S.computeQtp1(3);
				if (!masked_value_QtQtm1(4,S).check(S.Qt(4)))
					continue;
				if (!checkrotationQtQtp1(3,S))
					continue;

				S.computeQtp1(4); // m4 fixed
				if (!masked_value_QtQtm1(5,S).check(S.Qt(5)))
					continue;

				S.computeQtp1(5); // m5 fixed
				if (!masked_value_QtQtm1(6,S).check(S.Qt(6)))
					continue;

				S.computeQtp1(6); // m6 fixed
				if (!masked_value_QtQtm1(7,S).check(S.Qt(7)))
					continue;

				if (!checkrotationQtQtp1(4,S))
					continue;
				if (!checkrotationQtQtp1(5,S))
					continue;
				if (!checkrotationQtQtp1(6,S))
					continue;

//				++Q7attempts;
				auto it = Q7ptr.find(S.Qt(7));
				if (it == Q7ptr.end())
					continue;

				for (size_t q7i = it->second; q7i < vecQ7m10m12m13.size() && vecQ7m10m12m13[q7i][0] == S.Qt(7); ++q7i)
				{
//					++Q7match;

					S.m[10] = vecQ7m10m12m13[q7i][1];
					S.m[12] = vecQ7m10m12m13[q7i][2];
					S.m[13] = vecQ7m10m12m13[q7i][3];

					S.computeQtm3(13); // Q10
					S.computeQtm3(12); // Q9
					S.computeQtm3(11); // Q8
					S.computeWt(7);
					if (!MA.checkword(7, S.m[7]))
						continue;
					if (!checkrotationQtQtp1(7,S))
						continue;

//					++m7ok;


#if 0
					S.computeWt(8);
					if (!MA.checkword(8, S.m[8]))
						continue;
#endif
#if 0
					S.computeWt(9);
					if (!MA.checkword(9, S.m[9]))
						continue;
#endif

					
					// now iterate over Q9m9 tunnel
					S.computeQtp1(21); // m10 => Q22 already checked
					S.computeQtp1(22); // m15 => Q23 already checked
					S.computeQtp1(23); // m4 => Q24 already checked

					uint32_t Q9m9tunnel = ~Qvaluemask[offset+9] & ~Qprev[offset+10] & S.Qt(11) & ~S.Qt(10);
					uint32_t Q9org = S.Qt(9);
					uint32_t Q9cur = 0;
					do {
						Q9cur -= 1; Q9cur &= Q9m9tunnel;
						S.Qt(9) = Q9org ^ Q9cur;
						
						S.computeWt(9);
						if (!MA.checkword(9, S.m[9]))
							continue;
//						++m9ok;
						S.computeWt(8);
						if (!MA.checkword(8, S.m[8]))
							continue;
//						++m8ok;
						S.computeWt(12);
						if (!MA.checkword(12, S.m[12]))
							continue;
//						++m12ok;

						S.computeWt(10);
						if (!MA.checkword(10, S.m[10]))
							throw std::runtime_error("Q9m9tunnel: m10-invalid");
						S.computeWt(11);
						if (!MA.checkword(11, S.m[11]))
							throw std::runtime_error("Q9m9tunnel: m11-invalid");

						int tQtok = 24;
						for (int t = 24; t < 26; ++t)
						{
							S.computeQtp1(t);
							if (!masked_value_QtQtm1(t+1,S).check(S.Qt(t+1)))
								break;
							if (!checkrotationQtQtp1(t,S))
								break;
							tQtok = t+1;
						}
						if (tQtok < 24)
							continue;
						++Q24ok;
						
						check_solution(S);
						
					} while (Q9cur != 0);
				}
			}
		}
	}); // run_workload
	
	} // m4-loop
}

void textcoll_solver_t::check_solution(const fullstate_t& sol)
{
	fullstate_t S = sol;
#if 0
	// check consistency of solution
	for (int t = 0; t <= 23; ++t)
	{
		int wt = md5_wt[t];
		S.computeWt(t);
		if (S.m[wt] != sol.m[wt])
			std::cerr << "sol:W" << t << "-inconsistency!" << std::flush;
	}
#endif

	uint32 block1[16], block2[16];
	for (int k = 0; k < 16; ++k)
	{
		block1[k] = sol.m[k];
		block2[k] = block1[k] + m_diff[k];
	}
	uint32 solihv1[4], solihv2[4];
	for (unsigned i = 0; i < 4; ++i)
	{
		solihv1[i] = solihv2[i] = ihv1[i];
	}
	md5compress(solihv1, block1);
	md5compress(solihv2, block2);
	
	if (solihv2[0]-solihv1[0] != dQ[0] + dQ[offset+61])
		return;
		
	std::cout << "\n ================================ \n Partial dIHV match !!!" << std::endl;
	for (unsigned t = 0; t < 16; ++t)
	{
		std::cout << "block" << t << "=";
		for (unsigned b = 0; b < 4; ++b)
			std::cout << char((block1[t]>>(8*b))&0xFF);
		std::cout << " ";
	}
	std::cout << std::endl;
	
	if (solihv2[1]-solihv1[1] != dQ[3] + dQ[offset+64])
		return;
	if (solihv2[2]-solihv1[2] != dQ[2] + dQ[offset+63])
		return;
	if (solihv2[3]-solihv1[3] != dQ[1] + dQ[offset+62])
		return;

	// necessary conditions for second block of the attack
	if ( (solihv1[1]>>31) != (solihv1[2]>>31) || (solihv1[2]>>31) != (solihv1[3]>>31) )
		return;
		
	std::cout << "\n ================================ \n FULL SOLUTION FOUND \n ================================" << std::endl;
        uint32 x = xrng128();
        std::string filename1 = "textcoll1_" + boost::lexical_cast<string>(x) + ".txt";
        std::string filename2 = "textcoll2_" + boost::lexical_cast<string>(x) + ".txt";
        ofstream of1(filename1.c_str(), ios::binary);
        ofstream of2(filename2.c_str(), ios::binary);
        if ((!of1) || (!of2)) {
                cerr << "Cannot open output files!" << endl;
        }
        save_block(of1, block1);
        save_block(of2, block2);
        of1.close();
        of2.close();
        exit(0);
}

void textcoll_solver_t::compute_good_m10(const fullstate_t& Sin)
{
	fullstate_t S = Sin;
	good_m10.clear();

	auto m10range = MA.word_range(10);
	for (uint32_t m10 : m10range)
	{
		S.m[10] = m10;

		S.computeQtp1(21); // m10
		if (!masked_value_QtQtm1(22,S).check(S.Qt(22)))
			continue;

		S.computeQtp1(22); // m15
		if (!masked_value_QtQtm1(23,S).check(S.Qt(23)))
			continue;

		S.computeQtp1(23); // m4
		if (!masked_value_QtQtm1(24,S).check(S.Qt(24)))
			continue;

		if (!checkrotationQtQtp1(21,S))
			continue;
		if (!checkrotationQtQtp1(22,S))
			continue;
		if (!checkrotationQtQtp1(23,S))
			continue;

		good_m10.emplace_back(m10);
	}
	
}



void textcoll_solver_t::filltables()
{
	for (int t = -3; t <= 64; ++t)
	{
		dQ[offset+t] = dT[offset+t] = dR[offset+t] = 0;
		Qvaluemask[offset+t] = Qvalue[offset+t] = 0;
		Qprev[offset+t] = Qprev2[offset+t] = 0;
	}

	// build tables
	for (int t = diffpath.tbegin(); t < diffpath.tend(); ++t)
	{
		dQ[offset+t] = diffpath[t].diff();
		Qprev[offset+t] = diffpath[t].prev() | diffpath[t].prevn();
		Qprev2[offset+t] = diffpath[t].prev2() | diffpath[t].prev2n();
		Qvaluemask[offset+t] = (~diffpath[t].set0()) | diffpath[t].set1()
							| diffpath[t].prev() | diffpath[t].prevn()
							| diffpath[t].prev2() | diffpath[t].prev2n();
		Qvalue[offset+t] = diffpath[t].set1() | diffpath[t].prevn() | diffpath[t].prev2n();
	}
	for (int t = 0; t < 64; ++t)
	{
		booleanfunction* F = 0;
		if (t < 16) F = & MD5_F_data;
		else if (t < 32) F = & MD5_G_data;
		else if (t < 48) F = & MD5_H_data;
		else F = & MD5_I_data;
		uint32 dF = 0;
		for (unsigned b = 0; b < 32; ++b)
		{
			bf_outcome outcome = F->outcome(diffpath(t,b), diffpath(t-1,b), diffpath(t-2,b));
			if (outcome.size()) {
				if (outcome[0] == bc_plus) 			dF += 1<<b;
				else if (outcome[0] == bc_minus)	dF -= 1<<b;
			} else throw std::runtime_error("ambiguous path!!");
		}
		dT[offset+t] = dQ[offset+t-3] + dF + m_diff[md5_wt[t]];
		dR[offset+t] = dQ[offset+t+1] - dQ[offset+t];
	}

#if 0
	counter_exponential_print attempt("attempt"), ok("ok"), partial("partial");
	fullstate_t S1, S2;
	while (true)
	{
		++attempt;
		S1.Qt(21) = masked_value_Qt(21,S1).sample();
		S1.Qt(22) = masked_value_QtQtm1(22,S1).sample();
		S1.Qt(23) = masked_value_QtQtm1(23,S1).sample();
		S1.Qt(24) = masked_value_QtQtm1(24,S1).sample();
		for (int k = 0; k < 16; ++k)
		{
			S1.m[k] = xrng64();
			S2.m[k] = S1.m[k] + m_diff[k];
		}
		for (int t = 21; t <= 24; ++t)
			S2.Qt(t) = S1.Qt(t) + dQt(t);
		for (int t = 24; t < 64; ++t)
		{
			S1.computeQtp1(t);
			S2.computeQtp1(t);
		}
		if (S2.Qt(61)-S1.Qt(61) != uint32_t(1)<<31)
			continue;
		++partial;
		if (S2.Qt(62)-S1.Qt(62) != uint32_t(1)<<31)
			continue;
		if (S2.Qt(63)-S1.Qt(63) != uint32_t(1)<<31)
			continue;
		if (S2.Qt(64)-S1.Qt(64) != uint32_t(1)<<31)
			continue;
		if ( (S1.Qt(64)>>31) != (S1.Qt(63)>>31) )
			continue;
		if ( (S1.Qt(64)>>31) != (S1.Qt(62)>>31) )
			continue;
		++ok;
		std::cout << "p=" << log(double(ok())/double(attempt()))/log(2.0) << " partial";
		std::cout << "p=" << log(double(ok())/double(partial()))/log(2.0) << std::endl;
	}
#endif
//	Q[0] = Qvalue[0]; Q2[0] = Q[0] + dQ[0];
//	Q[1] = Qvalue[1]; Q2[1] = Q[1] + dQ[1];
//	Q[2] = Qvalue[2]; Q2[2] = Q[2] + dQ[2];
//	Q[3] = Qvalue[3]; Q2[3] = Q[3] + dQ[3];
}


bool test_collision(const uint32_t ihv1[4], const uint32_t ihv2[4], const uint32_t m[16])
{
	uint32_t tihv1[4];
	uint32_t tihv2[4];
	for (unsigned k = 0; k < 4; ++k)
	{
		tihv1[k] = ihv1[k];
		tihv2[k] = ihv2[k];
	}
	md5compress(tihv1, m);
	md5compress(tihv2, m);
	if (tihv1[0] == tihv2[0])
		std::cout << "partial solution found" << std::endl;
	bool collision = true;
	for (unsigned k = 0; k < 4; ++k)
		if (tihv1[k] != tihv2[k])
			collision = false;
	return collision;
	

}

void textcoll_solver_t::start_block2()
{
	fullstate_t S1, S2;
	S1.Qt(-3) = ihv1[0]; S1.Qt(-2) = ihv1[3]; S1.Qt(-1) = ihv1[2]; S1.Qt(0) = ihv1[1];
	S2.Qt(-3) = ihv2[0]; S2.Qt(-2) = ihv2[3]; S2.Qt(-1) = ihv2[2]; S2.Qt(0) = ihv2[1];

	std::cout << "MSBs: ";
	for (int t = 0; t >= -3; --t)
		std::cout << (S1.Qt(t)>>31);
	std::cout << std::endl;	
	if ( (S1.Qt(0)>>31) != (S1.Qt(-1)>>31) || (S1.Qt(-1)>>31) != (S1.Qt(-2)>>31) )
		throw std::runtime_error("IHV1 and IHV2 do not fulfill the necessary conditions such that dF0,dF1 = 1<<31.");

#if 0
	// do padding and strengthening in second block		
	run_workload(threads,[this](size_t ji, size_t jn)
	{
		fullstate_t S1;
		S1.Qt(-3) = ihv1[0]; S1.Qt(-2) = ihv1[3]; S1.Qt(-1) = ihv1[2]; S1.Qt(0) = ihv1[1];
		
		uint64_t totalbits = (64+8) * 8;
		S1.m[2] = 0x80;
		for (unsigned int t = 3; t < 14; ++t)
			S1.m[t] = 0;
		S1.m[14] = uint32_t(totalbits & 0xFFFFFFFF);
		S1.m[15] = uint32_t(totalbits >> 32);
		
		size_t cntr = 0;
		for (uint32_t m0 : MA.word_range(0))
		{
			++cntr; if (cntr > jn) cntr -= jn;
			if (cntr != ji)
				continue;
			
			S1.m[0] = m0;
			S1.computeQtp1(0);
			if ( (S1.Qt(1)>>31) != (S1.Qt(0)>>31) )
				continue;
			
			for (uint32_t m1 : MA.word_range(1))
			{
				S1.m[1] = m1;
				S1.computeQtp1(1);
				if ( (S1.Qt(2)>>31) != (S1.Qt(0)>>31) )
					continue;

				S1.computeQtp1(2);
				if ( (S1.Qt(3)>>31) != (S1.Qt(0)>>31) )
					continue;

				S1.computeQtp1(3);
				if ( (S1.Qt(4)>>31) != (S1.Qt(0)>>31) )
					continue;

				S1.computeQtp1(4);
				if ( (S1.Qt(5)>>31) != (S1.Qt(0)>>31) )
					continue;

				S1.computeQtp1(5);
				if ( (S1.Qt(6)>>31) != (S1.Qt(0)>>31) )
					continue;

				S1.computeQtp1(6);
				if ( (S1.Qt(7)>>31) != (S1.Qt(0)>>31) )
					continue;
					
				int t = 7;
				for (; t < 16; ++t)
				{
					S1.computeQtp1(t);
					if ( (S1.Qt(t+1)>>31) != (S1.Qt(0)>>31) )
						break;
				}
				if (t < 16)
					continue;

				if (test_collision(ihv1, ihv2, S1.m))
				{
					for (unsigned t = 0; t < 16; ++t)
					{
						std::cout << "block" << t << "=";
						for (unsigned b = 0; b < 4; ++b)
							std::cout << char((S1.m[t]>>(8*b))&0xFF);
						std::cout << " ";
					}
					std::cout << std::endl;
					std::cout << "\n ================================ \n FULL SOLUTION FOUND \n ================================" << std::endl;
				        uint32 x = xrng128();
				        std::string filename1 = "textcoll1_block2_" + boost::lexical_cast<string>(x) + ".txt";
				        std::string filename2 = "textcoll2_block2_" + boost::lexical_cast<string>(x) + ".txt";
				        ofstream of1(filename1.c_str(), ios::binary);
				        ofstream of2(filename2.c_str(), ios::binary);
				        if ((!of1) || (!of2)) {
				                cerr << "Cannot open output files!" << endl;
					}
				        save_block(of1, S1.m);
				        save_block(of2, S1.m);
				        of1.close();
				        of2.close();
				        exit(0);
				}

			}
		}
	});
#else
/*
	MA = message_alphabet("()");
	std::vector<uint32_t> words;
	for (uint32_t m0 : MA.word_range(0))
		words.emplace_back(m0);
		
	MA.set_byte_alphabet(0, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	MA.set_byte_alphabet(1, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	MA.set_byte_alphabet(2, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	MA.set_byte_alphabet(3, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	MA.set_byte_alphabet(4, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	MA.set_byte_alphabet(5, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	MA.set_byte_alphabet(6, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
	//MA.set_byte_alphabet(7, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ");
*/
	run_workload(threads,
	   [this]
	   (size_t ji, size_t jn)
	{
		fullstate_t S1;
		counter_exponential_print t15ok("t15ok");
		counter_exponential_print t0ok("t0ok");
		S1.Qt(-3) = ihv1[0]; S1.Qt(-2) = ihv1[3]; S1.Qt(-1) = ihv1[2]; S1.Qt(0) = ihv1[1];
		while (true)
		{
		/*
			S1.m[0] = MA.sampleword(0);
			
			if (S1.m[0] % 41 != ji) // split the search space
				continue;
				
			S1.computeQtp1(0);
			if ( (S1.Qt(1)>>31) != (S1.Qt(0)>>31) )
				continue;
			if (ji == 0) ++t0ok;

		*/	
			int t = 0;
			for (; t < 15; ++t)
			{
				while (true) //for (uint32_t mt : words)
				{
					S1.m[t] = //mt; //words[xrng64() % words.size()];
					MA.sampleword(t);
					S1.computeQtp1(t);
					if ( (S1.Qt(t+1)>>31) == (S1.Qt(t)>>31) )
						break;
				}
/*				if ( (S1.Qt(t+1)>>31) != (S1.Qt(t)>>31) )
				{
					t -= 2; continue;
					break;
				}*/
			}
//			if (t < 12) continue;

/*
			for (uint32_t m1 : MA.word_range(1))
			{
				S1.m[1] = m1; S1.computeQtp1(1); if ( (S1.Qt(2)>>31) != (S1.Qt(0)>>31) ) continue;
				
			for (uint32_t m2 : words)//MA.word_range(13))
			{
				S1.m[2] = m2; S1.computeQtp1(2); if ( (S1.Qt(3)>>31) != (S1.Qt(0)>>31) ) continue;
			for (uint32_t m3 : words)//MA.word_range(13))
			{
				S1.m[3] = m3; S1.computeQtp1(3); if ( (S1.Qt(4)>>31) != (S1.Qt(0)>>31) ) continue;
			for (uint32_t m4 : words)//MA.word_range(13))
			{
				S1.m[4] = m4; S1.computeQtp1(4); if ( (S1.Qt(5)>>31) != (S1.Qt(0)>>31) ) continue;
			for (uint32_t m5 : words)//MA.word_range(13))
			{
				S1.m[5] = m5; S1.computeQtp1(5); if ( (S1.Qt(6)>>31) != (S1.Qt(0)>>31) ) continue;
			for (uint32_t m6 : words)//MA.word_range(13))
			{
				S1.m[6] = m6; S1.computeQtp1(6); if ( (S1.Qt(7)>>31) != (S1.Qt(0)>>31) ) continue;
			for (uint32_t m7 : words)//MA.word_range(13))
			{
				S1.m[7] = m7; S1.computeQtp1(7); if ( (S1.Qt(8)>>31) != (S1.Qt(0)>>31) ) continue;
			for (uint32_t m8 : words)//MA.word_range(13))
			{
				S1.m[8] = m8; S1.computeQtp1(8); if ( (S1.Qt(9)>>31) != (S1.Qt(0)>>31) ) continue;
			for (uint32_t m9 : words)//MA.word_range(13))
			{
				S1.m[9] = m9; S1.computeQtp1(9); if ( (S1.Qt(10)>>31) != (S1.Qt(0)>>31) ) continue;
			for (uint32_t m10 : words)//MA.word_range(13))
			{
				S1.m[10] = m10; S1.computeQtp1(10); if ( (S1.Qt(11)>>31) != (S1.Qt(0)>>31) ) continue;

			for (uint32_t m11 : words)//MA.word_range(13))
			{
				S1.m[11] = m11; S1.computeQtp1(11); if ( (S1.Qt(12)>>31) != (S1.Qt(0)>>31) ) continue;

			for (uint32_t m12 : words)//MA.word_range(13))
			{
				S1.m[12] = m12; S1.computeQtp1(12); if ( (S1.Qt(13)>>31) != (S1.Qt(0)>>31) ) continue;
			for (uint32_t m13 : words)//MA.word_range(13))
			{
				S1.m[13] = m13; S1.computeQtp1(13); if ( (S1.Qt(14)>>31) != (S1.Qt(0)>>31) ) continue;
			for (uint32_t m14 : words)//MA.word_range(14))
			{
				S1.m[14] = m14; S1.computeQtp1(14); if ( (S1.Qt(15)>>31) != (S1.Qt(0)>>31) ) continue;
*/
			for (uint32_t m15 : MA.word_range(15))
			{
				S1.m[15] = m15; //S1.computeQtp1(15); if ( (S1.Qt(16)>>31) != (S1.Qt(0)>>31) ) continue;
				if (ji == 0) ++t15ok;
				if (test_collision(ihv1, ihv2, S1.m))
				{
					for (unsigned t = 0; t < 16; ++t)
					{
						std::cout << "block" << t << "=";
						for (unsigned b = 0; b < 4; ++b)
							std::cout << char((S1.m[t]>>(8*b))&0xFF);
						std::cout << " ";
					}
					std::cout << std::endl;
					std::cout << "\n ================================ \n FULL SOLUTION FOUND \n ================================" << std::endl;
				        uint32 x = xrng128();
				        std::string filename1 = "textcoll1_block2_" + boost::lexical_cast<string>(x) + ".txt";
				        std::string filename2 = "textcoll2_block2_" + boost::lexical_cast<string>(x) + ".txt";
				        ofstream of1(filename1.c_str(), ios::binary);
				        ofstream of2(filename2.c_str(), ios::binary);
				        if ((!of1) || (!of2)) {
				                cerr << "Cannot open output files!" << endl;
					}
				        save_block(of1, S1.m);
				        save_block(of2, S1.m);
				        of1.close();
				        of2.close();
				        exit(0);
				}
			} // m15
/*			} // m14
			} // m13
			} // m12
			} // m11
			} // m10
			} // m9
			} // m8
			} // m7
			} // m6
			} // m5
			} // m4
			} // m3
			} // m2
			} // m1*/
		}
	});
#endif
}

int textcoll(parameters_type& parameters)
{
	parameters.show_mdiffs();

	differentialpath diffpath;
	vector<differentialpath> vecpath;
	bool failed = true;
	try {
		load_gz(vecpath, binary_archive, parameters.infile1);
		failed = false;
	} catch (...) {}
	if (failed)
	{
		vecpath.clear();
		try {
			load_gz(diffpath, binary_archive, parameters.infile1);
			vecpath.push_back(diffpath);
			failed = false;
		} catch (...) {}
	}
	if (failed || vecpath.size() == 0) {
		cerr << "Error: could not load path(s) in '" << parameters.infile1 << "'!" << endl;
		return 1;
	}
	show_path(vecpath[0], parameters.m_diff);
	cout << "Starting..." << endl;

	textcoll_solver_t solver;
	for (unsigned i = 0; i < 16; ++i)
		solver.m_diff[i] = parameters.m_diff[i];
	solver.threads = parameters.threads;

        uint32 ihv1[4] = { md5_iv[0], md5_iv[1], md5_iv[2], md5_iv[3] };
        uint32 ihv2[4] = { md5_iv[0], md5_iv[1], md5_iv[2], md5_iv[3] };
        size_t prefixblocks = 0;

	if (parameters.infile2.size() > 0)
	{
	        uint32 msg1[16];
	        uint32 msg2[16];

                ifstream if2(parameters.infile2.c_str(), ios::binary);
                if (!if2)
                {
                        cerr << "Error: cannot open inputfile 2 '" << parameters.infile2 << "'!" << endl;
                        return 1;
                }
                while (load_block(if2, msg1) > 0)
                {
                	++prefixblocks;
                	for (unsigned k = 0; k < 4; ++k)
                		ihv2[k] = ihv1[k];
			for (unsigned k = 0; k < 16; ++k)
				msg2[k] = msg1[k] + parameters.m_diff[k];
                        md5compress(ihv1, msg1);
                        md5compress(ihv2, msg2);
                }
                bool do_block2 = true;
                for (unsigned k = 0; k < 4; ++k)
                	if (ihv2[k]-ihv1[k] != (uint32_t(1)<<31))
                		do_block2 = false;

	        if (!do_block2)
	        {
		        for (unsigned k = 0; k < 4; ++k)
		        	ihv2[k] = ihv1[k];
		}
	}

	solver.prefixblocks = prefixblocks;
	for (unsigned k = 0; k < 4; ++k)
	{
		solver.ihv1[k] = ihv1[k];
		solver.ihv2[k] = ihv2[k];
	}

        cout << "IHV1   = " << hex;
        for (unsigned k = 0; k < 4; ++k)
                for (unsigned c = 0; c < 4; ++c)
                {
                        cout.width(2); cout.fill('0');
                        cout << ((ihv1[k]>>(c*8))&0xFF);
                }
        cout << dec << endl;
        
        if (ihv1[0] != ihv2[0])
        {
	        cout << "IHV2   = " << hex;
	        for (unsigned k = 0; k < 4; ++k)
	                for (unsigned c = 0; c < 4; ++c)
	                {
	                        cout.width(2); cout.fill('0');
	                        cout << ((ihv2[k]>>(c*8))&0xFF);
	                }
	        cout << dec << endl;
	        uint32 dihv[4] = { ihv2[0] - ihv1[0], ihv2[1] - ihv1[1], ihv2[2] - ihv1[2], ihv2[3] - ihv1[3] };
	        cout << "dIHV   = {" << naf(dihv[0]) << "," << naf(dihv[1]) << "," << naf(dihv[2]) << "," << naf(dihv[3]) << "}" << endl;
	}
	solver.findcollision(vecpath, parameters.alphabet);
	return 0;
}

