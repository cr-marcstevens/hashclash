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

#ifndef HASHCLASH_SHA1MESSAGESPACE_HPP
#define HASHCLASH_SHA1MESSAGESPACE_HPP

#ifndef NOSERIALIZATION
#include <boost/serialization/serialization.hpp>
#endif // NOSERIALIZATION

#include "types.hpp"
#include "sdr.hpp"
#include "sha1detail.hpp"
#include "saveload_bz2.hpp"

namespace hashclash {

	void sweep_matrix(std::vector< std::vector<uint32> >& matrix, int columns = -1);
	void complement_basis(const std::vector< std::vector<uint32> >& basisin, std::vector< std::vector<uint32> >& basisout);

	inline void append_col(std::vector< std::vector<uint32> >& matrix1, const std::vector< uint32 >& vector) {
		if (matrix1.size() != vector.size()) throw;
		for (unsigned i = 0; i < matrix1.size(); ++i)
			matrix1[i].push_back(vector[i]);
	}
	inline void delete_lastcol(std::vector< std::vector<uint32> >& matrix1) {
		for (unsigned i = 0; i < matrix1.size(); ++i)
			matrix1[i].pop_back();
	}
	inline void xorvec(std::vector<uint32>& srcdest, const std::vector<uint32>& toxorwith) {
		if (srcdest.size() != toxorwith.size()) throw std::runtime_error("xorvec()::different sizes!");
		for (unsigned k = 0; k < srcdest.size(); ++k)
			srcdest[k] ^= toxorwith[k];
	}	
	inline void andvec(std::vector<uint32>& srcdest, const std::vector<uint32>& toandwith) {
		if (srcdest.size() != toandwith.size()) throw std::runtime_error("andvec()::different sizes!");
		for (unsigned k = 0; k < srcdest.size(); ++k)
			srcdest[k] &= toandwith[k];
	}


	inline void reduce_basis(std::vector< std::vector<uint32> >& matrix)
	{
		if (matrix.size() == 0) return;
		sweep_matrix(matrix);
		std::vector<uint32> nulvector(matrix[0].size(),0);
		while (matrix.size() > 0 && matrix[matrix.size()-1] == nulvector)
			matrix.pop_back();
	}
	inline void reduce_vector(std::vector< uint32 >& vec, const std::vector< std::vector<uint32> >& basis)
	{
		if (basis.size() == 0) return;
		if (vec.size() != basis[0].size()) throw std::runtime_error("reduce_vector()::different sizes!");
		int col = int(basis[0].size()*32) - 1;
		for (unsigned k = 0; k < basis.size(); ++k)
			for (; col >= 0; --col)
				if (basis[k][col>>5]&(1<<(col&31))) {
					if (vec[col>>5]&(1<<(col&31)))
						xorvec(vec, basis[k]);
					break;
				}
	}

	class sha1messagespace;
	void derive_sha1messagespace(sha1messagespace& space, unsigned tbegin, unsigned tend, std::vector< std::vector<sdr> >& dmes, std::vector<double>& dmes_prob, double minavgprob = 0);
	inline void derive_sha1messagespace(sha1messagespace& space, unsigned tbegin, unsigned tend, std::vector< std::vector<sdr> >& dmes) {
		std::vector<double> tmp;
		derive_sha1messagespace(space, tbegin, tend, dmes, tmp);
	}

	class sha1messagespace {
	public:
		sha1messagespace(): offset(80) {}

		bool operator==(const sha1messagespace& r) const {
			return r.basis == basis && r.offset == offset;
		}
		bool operator!=(const sha1messagespace& r) const {
			return r.basis != basis || r.offset != offset;
		}
		void clear() {
			basis.clear();
			offset.assign(80,0);
		}

		unsigned columns() const { return 80*32; }
		unsigned rows() const { return unsigned(basis.size()); }

		void buildbasis_setbit(unsigned t, unsigned b, bool isone) {
			if (isone)
				offset[t] |= 1<<b;
			else
				offset[t] &= ~uint32(1<<b);
		}
		void buildbasis_addfreebit(unsigned t, unsigned b) {
			std::vector<uint32> tmp(80,0);
			tmp[t] |= 1<<b;
			basis.push_back(tmp);
		}
		void add2basis(const std::vector<uint32>& vec) {
			basis.push_back(vec);
			basis[basis.size()-1].resize(80,0);
		}
		void add2basis(const std::vector< std::vector<uint32> >& basis2) {
			for (unsigned i = 0; i < basis2.size(); ++i)
				add2basis(basis2[i]);
		}

		void reduce() { reduce_basis(basis); reduce_vector(offset, basis); }
		// reduces sha1messagespace if necessary !!
		//void reduce_vector(std::vector<uint32>& vec);

		bool isinmessagespace(const std::vector<uint32>& expandedmessage) {
			std::vector<uint32> tmp = expandedmessage;
			xorvec(tmp, offset);
			int col = columns()-1;
			for (unsigned k = 0; k < basis.size(); ++k)
				for (; col >= 0; --col)
					if (basis[k][col>>5]&(1<<(col&31))) {
						if (tmp[col>>5]&(1<<(col&31)))
							xorvec(tmp, basis[k]);
						break;
					}
			for (unsigned i = 0; i < tmp.size(); ++i)
				if (tmp[i]) return false;
			return true;
		}
		void tobitrelations_80(std::vector< std::vector<uint32> >& bitrelations);
		void tobitrelations_16(std::vector< std::vector<uint32> >& bitrelations);
		void frombitrelations_80(const std::vector< std::vector<uint32> >& bitrelations);
		void addbitrelations(const std::vector< std::vector<uint32> >& bitrelations);
		void addbitrelation(const std::vector<uint32>& bitrelation) {
			std::vector< std::vector<uint32> > tmp;
			tmp.push_back(bitrelation);
			addbitrelations(tmp);
		}

		void swap(sha1messagespace& r) {
			basis.swap(r.basis);
			offset.swap(r.offset);
		}

		std::vector< std::vector<uint32> > basis;
		std::vector<uint32> offset;
	};

	void read_message_bitconditions(std::istream& is, std::vector< std::vector<uint32> >& bitrelations80);
	
} // namespace hashclash

#ifndef NOSERIALIZATION

namespace boost {
	namespace serialization {

		template<class Archive>
		void serialize(Archive& ar, hashclash::sha1messagespace& s, const unsigned int file_version)
		{			
			ar & make_nvp("offset", s.offset);
			ar & make_nvp("basis", s.basis);
		}

	}
}
#endif // NOSERIALIZATION

#endif //HASHCLASH_SHA1MESSAGESPACE_HPP
