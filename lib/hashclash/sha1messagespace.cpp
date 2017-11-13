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

#include "sha1messagespace.hpp"
#include "progress_display.hpp"
#include <cmath>

using namespace std;

namespace hashclash {

	void sweep_matrix(vector< vector<uint32> >& matrix, int columns)
	{
		int rows = matrix.size();
		if (rows == 0) return;
		if (columns == -1) columns = matrix[0].size()*32;

		int fixed = 0;
		for (int col = columns - 1; col >= 0; --col) {
			int firstrow = rows;
			for (int row = fixed; row < rows; ++row)
				if (matrix[row][col>>5] & (1<<(col&31))) {
					firstrow = row;
					break;
				}
			if (firstrow == rows) continue;
			matrix[fixed].swap(matrix[firstrow]);
			for (int row = 0; row < rows; ++row)
				if (row != fixed && matrix[row][col>>5] & (1<<(col&31)))
					xorvec(matrix[row], matrix[fixed]);
			++fixed;
			if (fixed == matrix.size())
				return;
		}
	}

	void complement_basis(const vector< vector<uint32> >& basisin, vector< vector<uint32> >& basisout)
	{
		if (basisin.size() == 0) throw std::runtime_error("complement_basis: trivial basis given");
		int columns = basisin[0].size()*32;
		vector< vector<uint32> > b = basisin;
		reduce_basis(b);
		vector< unsigned > blastcol(b.size());
		for (unsigned row = 0; row < b.size(); ++row)
			for (int col = columns-1; col >= 0; --col)
				if (b[row][col>>5]&(1<<(col&31))) {
					blastcol[row] = col;
					break;
				}
		basisout.clear();
		vector< uint32 > tmp;
		for (int col = 0; col < columns; ++col) {
			bool incomplement = true;
			for (unsigned row = 0; row < blastcol.size(); ++row)
				if (blastcol[row] == col)
					incomplement = false;
			if (!incomplement) continue;
			tmp.clear();
			tmp.resize(basisin[0].size(),0);
			tmp[col>>5] |= 1<<(col&31);
			for (unsigned row = 0; row < b.size(); ++row)
				if (b[row][col>>5]&(1<<(col&31))) {
					unsigned c = blastcol[row];
					tmp[c>>5] |= 1<<(c&31);
				}
			basisout.push_back(tmp);
		}
	}

	void sha1messagespace::tobitrelations_80(vector< vector<uint32> >& bitrelations) {
		// convert basis to bitrelations
		reduce();
		complement_basis(basis, bitrelations);
		for (unsigned i = 0; i < bitrelations.size(); ++i) {
			unsigned v = 0;
			for (unsigned j = 0; j < 80; ++j)
				v += hw(offset[j] & bitrelations[i][j]);
			if (v&1)
				bitrelations[i].push_back(~uint32(0));
			else
				bitrelations[i].push_back(0);
		}
	}
	void sha1messagespace::tobitrelations_16(vector< vector<uint32> >& bitrelations) {
		// convert basis to bitrelations
		tobitrelations_80(bitrelations);
		// compute msg-mask for each message-expansion-bit
		vector< vector< vector<uint32> > > bittovec16;
		uint32 me[80];
		uint32 msg[16];
		vector<uint32> bitrelation(16);
		if (bittovec16.size() == 0) {
			bittovec16.resize(80);
			for (int t = 0; t < 80; ++t) {
				bittovec16[t].resize(32);
				for (int b = 0; b < 32; ++b) {
					// clear current bitrelation
					bitrelation.assign(16,0);
					for (int tt = 0; tt < 16; ++tt)
						for (int bb = 0; bb < 32; ++bb)
						{
							memset(msg, 0, sizeof(msg));
							msg[tt] ^= 1<<bb;
							sha1_me_simple(me, msg);
							if (me[t]&(1<<b))
								bitrelation[tt] ^= 1<<bb;
						}
					bittovec16[t][b] = bitrelation;
				}
			}
		}
		// shorten each bitrelations to message mask with value
		for (unsigned i = 0; i < bitrelations.size(); ++i) {
			bitrelation.assign(16,0);
			for (unsigned t = 0; t < 80; ++t)
				for (unsigned b = 0; b < 32; ++b)
					if (bitrelations[i][t]&(1<<b))
						xorvec(bitrelation, bittovec16[t][b]);
			bitrelation.push_back(bitrelations[i][80]);
			bitrelations[i] = bitrelation;
		}
		// check for contradiction and remove trivial vectors
		sweep_matrix(bitrelations,16*32);
		vector<uint32> nulvector(17,0);
		nulvector[16] = ~uint32(0);
		for (unsigned i = 0; i < bitrelations.size(); ++i)
			if (bitrelations[i] == nulvector)
				throw std::runtime_error("sha1messagespace::tobitrelations_16(): contradicting bitrelations");
		nulvector[16] = 0;
		while (bitrelations.size() > 0 && bitrelations[bitrelations.size()-1] == nulvector)
			bitrelations.pop_back();
	}
	void sha1messagespace::frombitrelations_80(const vector< vector<uint32> >& bitrelations) {
		vector< vector<uint32> > tmp = bitrelations;
		sweep_matrix(tmp, 80*32);
		// check for contradictions and remove trivial vectors
		vector<uint32> nulvector(81,0);
		nulvector[80] = ~uint32(0);
		for (unsigned i = 0; i < tmp.size(); ++i)
			if (tmp[i] == nulvector)
				throw std::runtime_error("sha1messagespace::frombitrelations_80(): contradicting bitrelations");
		nulvector[80] = 0;
		while (tmp.size() > 0 && tmp[tmp.size()-1] == nulvector)
			tmp.pop_back();
		// convert bitrelations to basis
		vector< vector<uint32> > tmp2 = tmp;
		delete_lastcol(tmp2);
		complement_basis(tmp2, basis);
		offset.clear();
		offset.resize(80,0);
		for (unsigned i = 0; i < tmp.size(); ++i)
			if (tmp[i][80] != 0) {
				int lastcol = -1;
				for (int col = 0; col < 80*32; ++col)
					if (tmp[i][col>>5]&(1<<(col&31)))
						lastcol = col;
				offset[lastcol>>5] ^= 1<<(lastcol&31);
			}
	}
	void sha1messagespace::addbitrelations(const vector< vector<uint32> >& bitrelations) {
		// convert basis to bitrelations
		vector< vector<uint32> > basis_complement;
		tobitrelations_80(basis_complement);
		// add bitrelations
		for (unsigned i = 0; i < bitrelations.size(); ++i) {
			basis_complement.push_back(bitrelations[i]);
			basis_complement[basis_complement.size()-1].resize(81,0);
		}
		frombitrelations_80(basis_complement);
	}
/*	void sha1messagespace::reduce_vector(std::vector<uint32>& vec)
	{
		vec.resize(80,0);
		vector<int> lastcol(basis.size(),-1);
		for (unsigned i = 0; i < basis.size(); ++i) {
			for (int col = 80*32-1; col >= 0; --col)
				if (basis[i][col>>5]&(1<<(col&31))) {
					lastcol[i] = col;
					break;
				}
			if (lastcol[i] == -1) continue;
			for (unsigned j = 0; j < basis.size(); ++j)
				if (j != i && (basis[j][lastcol[i]>>5]&(1<<(lastcol[i]&31))))
				{
					reduce();
					reduce_vector(vec);
					return;
				}
		}
		for (unsigned i = 0; i < basis.size(); ++i)
			if (lastcol[i] >= 0 && (vec[lastcol[i]>>5]&(1<<(lastcol[i]&31))))
				xorvec(vec, basis[i]);
	}
	*/

	template<class iterator>
	class iterator_wrapper {
	public:
		typedef iterator iterator_type;
		typedef typename iterator_type::value_type value_type;

		iterator_wrapper() {}
		iterator_wrapper(const iterator_wrapper& r): it(r.it) {}
		iterator_wrapper(const iterator_type& ptr): it(ptr) {}

		iterator_type it;
		bool operator<(const iterator_wrapper& r) const
		{ return **this < *r; }
		bool operator==(const iterator_wrapper& r) const
		{ return **this < *r; }
		bool operator!=(const iterator_wrapper& r) const
		{ return **this != *r; }
		const value_type& operator*() const { return *it; }
	};
	typedef iterator_wrapper< vector< vector<sdr> >::const_iterator > dme_wrapper;

	void derive_sha1messagespace(sha1messagespace& space, unsigned tbegin, unsigned tend, vector< vector<sdr> >& dmes
		, std::vector<double>& dmes_prob //= vector<double>()
		, double minavgprob // = 0
		)
	{
		if (dmes.size() == 0)
			return;
//		vector< vector<sdr> > dmes(_dmes);
//		vector<double> dmes_prob(_dmes_prob);
		dmes_prob.resize(dmes.size(),0);
		vector<bool> dmes_inuse(dmes.size(), false);
		// clean up dmes
		for (unsigned i = 0; i < dmes.size(); ++i) {
			dmes[i].resize(80);
			for (int t = 0; t < tbegin; ++t)
				dmes[i][t].clear();
			for (int t = tend; t < 80; ++t)
				dmes[i][t].clear();
		}
		// sort dmes
		vector< dme_wrapper > dmes_sorted;
		for (dme_wrapper::iterator_type cit = dmes.begin(); cit != dmes.end(); ++cit)
			dmes_sorted.push_back( dme_wrapper(cit) );
		sort(dmes_sorted.begin(), dmes_sorted.end());
		for (unsigned i = 0; i < dmes_sorted.size(); ++i) {
			if (i > 0 && dmes_sorted[i] == dmes_sorted[i-1]) continue;
			unsigned j = i+1;
			while (j < dmes_sorted.size() && dmes_sorted[j] == dmes_sorted[i])
				++j;
			unsigned k = i;
			for (unsigned m = i; m < j; ++m)
				if (dmes_prob[ dmes_sorted[m].it - dmes.begin() ] > dmes_prob[ dmes_sorted[k].it - dmes.begin() ])
					k = m;
			swap(dmes_sorted[k].it, dmes_sorted[i].it);
		}
/*		dmes_sorted.erase( unique(dmes_sorted.begin(), dmes_sorted.end()), dmes_sorted.end() );
		if (dmes_sorted.size() != dmes.size())
			cerr << "derive_sha1messagespace(): warning duplicate dmes, inconsistent behaviour possible!" << endl;
*/

		// sanity check on dmes
 		vector<uint32> bitmask(80);
		for (unsigned t = 0; t < 80; ++t)
			bitmask[t] = dmes[0][t].mask;
		for (unsigned i = 1; i < dmes.size(); ++i)
			for (unsigned t = 0; t < 80; ++t)
				if (dmes[i][t].mask != bitmask[t]) throw std::runtime_error("derive_sha1messagespace(): dmes not sane!");
		// first pick best dme and build static newspace
		sha1messagespace newspace;
		for (int t = 0; t < 80; ++t)
			for (unsigned b = 0; b < 32; ++b)
				if (t < tbegin || t >= tend || (bitmask[t]&(1<<b))==0 || b == 31)
					newspace.buildbasis_addfreebit(t,b);
		unsigned bestindex = 0;
		for (unsigned i = 0; i < dmes.size(); ++i)
			if (dmes_prob[i] > dmes_prob[bestindex])
				bestindex = i;
		for (int t = 0; t < 80; ++t)
			newspace.offset[t] = dmes[bestindex][t].set1conditions();
		double avgprob = dmes_prob[bestindex];
		dmes_inuse[bestindex] = true;

		bool equiprobs = true;
		for (unsigned kk = 1; kk < dmes_prob.size(); ++kk)
			if (dmes_prob[kk] != dmes_prob[kk-1]) {
				equiprobs = false;
				break;
			}

		// iteratively pick best basisextension vector so that newspace remains completely contained in dmes
		vector< vector<sdr> > tmpdme; // single vector<sdr>(80) encapsulated in vector to facilitate dme_wrapper binary search
		tmpdme.push_back(vector<sdr>(80));
		vector<uint32> tmpxor(80), tmpxor2(80);
		vector< vector<uint32> > tmpxoradded;
		set< vector<uint32> > tmpxortried;

		newspace.tobitrelations_80(tmpxoradded);
		cout << "Current number of bitrestrictions: " << tmpxoradded.size() << endl;
		tmpxoradded.clear();

		unsigned freedomgained = 0;
		while (true) {
			int bestaddindex = -1;
			double newavgprob = -1;
			double newuseprob = -1;
			tmpxortried.clear();
			hashclash::progress_display pd(dmes.size());
			for (unsigned i = 0; i < dmes.size(); ++i,++pd) {
				if (!dmes_inuse[i]) {
					// take difference with dmes[bestindex]
					for (int t = tbegin; t < tend; ++t)
						tmpxor[t] = dmes[i][t].sign ^ dmes[bestindex][t].sign;
					tmpxor2 = tmpxor;
					reduce_vector(tmpxor2, tmpxoradded);
					if (tmpxortried.find(tmpxor2) != tmpxortried.end())
						continue;
					// test if forall inuse dmes: dmes^diff in dmes
					bool i_is_ok = true;
					double tmpavgprob = 0;
					unsigned tmpavgcnt = 0;
					for (unsigned j = 0; j < dmes.size() && i_is_ok; ++j) {
						if (dmes_inuse[j]) {
							for (int t = tbegin; t < tend; ++t) {
								tmpdme[0][t] = dmes[j][t];
								tmpdme[0][t].sign ^= tmpxor[t];
							}
							vector< dme_wrapper >::const_iterator cit = lower_bound( dmes_sorted.begin(), dmes_sorted.end(), dme_wrapper(tmpdme.begin()) );
							if (cit == dmes_sorted.end() || *(cit->it) != tmpdme[0])
								i_is_ok = false;
							else {
								unsigned index = cit->it - dmes.begin();
								tmpavgprob += dmes_prob[index];
								++tmpavgcnt;
							}
						}
					}
					if (!i_is_ok) continue;
					// determine newavgprob
					tmpavgprob += avgprob * double(tmpavgcnt);
					tmpavgprob /= double(2*tmpavgcnt);
					if (tmpavgprob >= newavgprob) {
						if (bestaddindex == -1) cout << "." << flush;
						double useprob = -1;
#if 1
						if (equiprobs) {
							unsigned usecnt = 0;
							for (unsigned j = 0; j < dmes.size(); ++j) {
								for (int t = tbegin; t < tend; ++t) {
									tmpdme[0][t] = dmes[j][t];
									tmpdme[0][t].sign ^= tmpxor[t];
								}
								vector< dme_wrapper >::const_iterator cit = lower_bound( dmes_sorted.begin(), dmes_sorted.end(), dme_wrapper(tmpdme.begin()) );
								if (cit != dmes_sorted.end() && *(cit->it) == tmpdme[0])
									++usecnt;
							}
							useprob = double(usecnt)/double(dmes.size());
						}
#endif
						if (useprob >= newuseprob) {
							if (bestaddindex == -1) cout << "." << flush;
							newavgprob = tmpavgprob;
							bestaddindex = i;
							newuseprob = useprob;
						}
						if (equiprobs && (useprob == -1 || useprob == 1)) {
							unsigned iplus = dmes.size()-i;
							i += iplus;
							pd += iplus;
						}
					}
					reduce_vector(tmpxor, tmpxoradded);
					tmpxortried.insert(tmpxor2);
				}
			}
			// add basis vector if any and break otherwise
			if (newavgprob >= 0 && newavgprob >= minavgprob) {
				++freedomgained;
				cout << "added new basis vector: avgprob(" << log(avgprob)/log(2.0) << "=>" << log(newavgprob)/log(2.0) << ") freedom gained(" << freedomgained << ")" << endl;
				for (int t = tbegin; t < tend; ++t)
					tmpxor[t] = dmes[bestaddindex][t].sign ^ dmes[bestindex][t].sign;
				newspace.add2basis(tmpxor);
				avgprob = newavgprob;
				tmpxoradded.push_back(tmpxor);
				reduce_basis(tmpxoradded);
				// update inuse vector taking into account possible double values
				vector<bool> tmpinuse(dmes_inuse);
				for (unsigned j = 0; j < dmes.size(); ++j)
					if (dmes_inuse[j]) {
						tmpinuse[j] = true;
						for (int t = tbegin; t < tend; ++t) {
							tmpdme[0][t] = dmes[j][t];
							tmpdme[0][t].sign ^= tmpxor[t];
						}
						pair< vector< dme_wrapper >::const_iterator, vector< dme_wrapper >::const_iterator>
							citpair = equal_range( dmes_sorted.begin(), dmes_sorted.end(), dme_wrapper(tmpdme.begin()) );
						vector< dme_wrapper >::const_iterator cit = citpair.first; //lower_bound( dmes_sorted.begin(), dmes_sorted.end(), dme_wrapper(tmpdme.begin()) );
						vector< dme_wrapper >::const_iterator citend = citpair.second;
						for (; cit != citend; ++cit) {
							unsigned index = cit->it - dmes.begin();
							tmpinuse[index] = true;
						}
					}
				dmes_inuse = tmpinuse;
				unsigned totinuse = 0;
				for (unsigned j = 0; j < dmes.size(); ++j)
					if (dmes_inuse[j])
						++totinuse;
				cout << ((double(totinuse)*100.0)/double(dmes.size())) << "% in use: " << totinuse << endl;
			} else {
				break;
				cout << "no new basis vector: avgprob(" << log(avgprob)/log(2.0) << ") bestnewavgprob(" << log(newavgprob)/log(2.0) << ") freedom gained(" << freedomgained << ")" << endl;
			}
		}
		space.basis.swap(newspace.basis);
		space.offset.swap(newspace.offset);
		space.tobitrelations_80(tmpxoradded);
		cout << "Total number of bitrestrictions: " << tmpxoradded.size() << endl;
	}

	void getline(std::istream& ifs)
	{
		cout << "{";
	        char c = 0;
	        while (ifs && c != '\n' && c != '\r') {
	                ifs.get(c);
	                cout << c;
		}
	        while (ifs && (c == '\n' || c == '\r')) {
	                ifs.get(c);
	                cout << c;
		}
		cout << "}" << flush;
	        if (ifs)
	                ifs.putback(c);
	}

	bool read_until_char(std::istream& is, const std::string& chars, bool putback = true) {
		char c;
		if (!(is >> c)) return !!is;
		cout << "<" << c;
		while (chars.find(c) == std::string::npos) {
			if (!(is >> c)) return !!is;
			cout << c;
		}
		cout << ">" << flush;
		if (putback)
			is.putback(c);
		return !!is;
	}
	void read_message_bitconditions(std::istream& is, vector< vector<uint32> >& bitrelations80)
	{
		// format: mt[!sdr!] + mt[!sdr!] = 0/1
		// example: m1[!0,5!] + m20[!31!] = 0
		bitrelations80.clear();
		vector<uint32> bitrel;
		char c;
		cout << "Parsing message bitconditions:" << endl;
		while (is) {
			bitrel.assign(81,0);
			if (!read_until_char(is, "#m", true)) return;
			is >> c;
			if (c =='#') {
				cout << "#" << flush;
				getline(is);
				continue;
			}
			cout << "m" << flush;
			int t = -1;
			if (!(is >> t)) return;
			cout << t << flush;
			if (t < 0 || t >= 80) return;
			if (!read_until_char(is, "[")) return;
			sdr mask;
			if (!(is >> mask)) return;
			cout << mask << " " << flush;
			bitrel[t] = mask.mask;
			// first mt[!sdr!], now keep checking for "+mt[!sdr!]" or "=0/1"
			while (true) {
				if (!read_until_char(is, "+=")) return;
				is >> c;
				if (c == '+') {
					cout << "+ " << flush;
					if (!read_until_char(is, "m", false)) return;
					cout << "m" << flush;
					int t2 = -1;
					if (!(is >> t2)) { cout << t2; return; }
					cout << t2 << flush;
					if (t2 < 0 || t2 >= 80) return;
					if (!read_until_char(is, "[")) return;
					if (!(is >> mask)) return;
					cout << mask << " " << flush;
					bitrel[t2] = mask.mask;
				} else {
					if (!read_until_char(is, "01")) return;
					is >> c;
					if (c == '1')
						bitrel[80] = ~uint32(0);
					cout << "=" << (bitrel[80]&1) << endl;
					bitrelations80.push_back(bitrel);
					getline(is);
					break;
				}
			}
		} 
		return;
	}
	
} // namespace hashclash
