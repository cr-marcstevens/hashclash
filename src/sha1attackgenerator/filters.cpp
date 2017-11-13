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

#include "main.hpp"
#include <hashclash/rng.hpp>

void mespace_to_pathbitrelationsmatrix(const sha1messagespace& mespace, vector< vector< vector<uint32> > >& pathbitrelationsmatrix)
{
	sha1messagespace tmpspace = mespace;
	vector< vector<uint32> > bitrels;
	tmpspace.tobitrelations_16(bitrels);
	pathbitrelationsmatrix.assign(16, vector< vector<uint32> >(32));
	for (unsigned i = 0; i < bitrels.size(); ++i) {
		for (int col = (16*32)-1; true; --col) {
			if (col < 0) throw std::runtime_error("mespace_to_pathbitrelationsmatrix(): col < 0");
			unsigned t = col>>5;
			unsigned b = col&31;
			if (bitrels[i][t] & (1<<b)) {
				pathbitrelationsmatrix[t][b] = bitrels[i];
				break;
			}
		}
	}
}

void random_me(uint32 me[80], const vector< vector< vector<uint32> > >& pathbitrelationsmatrix)
{
	for (int t = 0; t < 16; ++t) {
		uint32 metmask = 0;
		uint32 metset1 = 0;
		for (unsigned b = 0; b < 32; ++b) {
			if (pathbitrelationsmatrix[t][b].size()) {
				metmask |= pathbitrelationsmatrix[t][b][t];
				uint32 v = pathbitrelationsmatrix[t][b][16]&1;
				for (int t1 = 0; t1 < t; ++t1)
					v ^= me[t1]&pathbitrelationsmatrix[t][b][t1];
				if (hw(v)&1)
					metset1 |= 1<<b;
			}
		}
		me[t] = (xrng128() & ~metmask) | metset1;
	}
	for (int i = 16; i < 80; ++i)
		me[i]=rotate_left(me[i-3] ^ me[i-8] ^ me[i-14] ^ me[i-16], 1);
}

bool add_bitcondition(sha1differentialpath& path, int t, unsigned b, bitcondition bc)
{
	bitcondition bctb = path(t,b);
	if (bc == bc_constant) {
			if (bctb == bc_plus || bctb == bc_minus)
				return false;
			return true;
	}
	if (bctb == bc)
		return true;
	if (bctb == bc_constant) {
		path.setbitcondition(t,b,bc);
		return true;
	}
	if (bctb == bc_plus || bctb == bc_minus || bc == bc_plus || bc == bc_minus) return false;
	// since bctb!= bc these ifs test for a contradiction
	if ((bctb == bc_one || bctb == bc_zero) && (bc == bc_one || bc == bc_zero)) return false;
	if ((bctb == bc_prev || bctb == bc_prevn) && (bc == bc_prev || bc == bc_prevn)) return false;
	if ((bctb == bc_prev2 || bctb == bc_prev2n) && (bc == bc_prev2 || bc == bc_prev2n)) return false;
	if ((bctb == bc_next || bctb == bc_nextn) && (bc == bc_next || bc == bc_nextn)) return false;
	if ((bctb == bc_next2 || bctb == bc_next2n) && (bc == bc_prev2 || bc == bc_next2n)) return false;

	if (bctb == bc_one || bctb == bc_zero) {
		if (bc == bc_prev) {
			return add_bitcondition(path, t-1, b, bctb);
		} else if (bc == bc_prevn) {
			if (bctb == bc_one)
				return add_bitcondition(path, t-1, b, bc_zero);
			else
				return add_bitcondition(path, t-1, b, bc_one);
		}
		if (bc == bc_next) {
			return add_bitcondition(path, t+1, b, bctb);
		} else if (bc == bc_nextn) {
			if (bctb == bc_one)
				return add_bitcondition(path, t+1, b, bc_zero);
			else
				return add_bitcondition(path, t+1, b, bc_one);
		}			
		throw std::runtime_error("add_bitcondition(): unsupported bitcondition (as of yet) (code 1)!");
	}
	if (bctb == bc_prev) {
		if (bc == bc_one || bc == bc_zero)
			return add_bitcondition(path, t-1, b, bc);
		throw std::runtime_error("add_bitcondition(): unsupported bitcondition (as of yet) (code 2)!");
	}
	if (bctb == bc_prevn) {
		if (bc == bc_one)
			return add_bitcondition(path, t-1, b, bc_zero);
		if (bc == bc_zero)
			return add_bitcondition(path, t-1, b, bc_one);
		throw std::runtime_error("add_bitcondition(): unsupported bitcondition (as of yet) (code 3)!");
	}
	throw std::runtime_error("add_bitcondition(): unsupported bitcondition (as of yet) (code 4)!");
}

void filter_tunnels_bitconditions(vector< sha1differentialpath >& tunnels, const sha1differentialpath& path)
{
	sha1differentialpath tmp2 = path;
	cleanup_path(tmp2);
	// tunnel is orthogonal on message 1, so + => 0, - => 1
	for (int t = tmp2.tbegin(); t < tmp2.tend(); ++t) {
		for (unsigned b = 0; b < 32; ++b) {
			bitcondition bc = tmp2(t,b);
			if (bc == bc_plus)
				tmp2.setbitcondition(t,b,bc_zero);
			else if (bc == bc_minus)
				tmp2.setbitcondition(t,b,bc_one);
		}
	}
	int pathsize = tmp2.offset+20;
	tmp2.path.resize(pathsize);
	tmp2.me.resize(pathsize);

	unsigned idx = 0;
	while (idx < tunnels.size()) {
		bool ok = true;
		const sha1differentialpath& tunnel = tunnels[idx];
		unsigned nc = tunnel.nrcond();
		sha1differentialpath tmp = tmp2;
		for (int t = tunnel.tbegin(); ok && t < tunnel.tend(); ++t) {
			bool changed = false;
			for (unsigned b = 0; ok && b < 32; ++b) {
				bitcondition bc = tunnel(t,b);
				if (bc == bc_constant)
					continue;
				if (bc == bc_plus || bc == bc_minus) {
					if (tmp(t,b) != bc_constant || tmp(t+1,b) == bc_prev || tmp(t+1,b) == bc_prevn)
						ok = false;
					continue;
				}
				if (!add_bitcondition(tmp, t, b, bc))
					ok = false;
				else
					changed = true;
			}
			if (ok && changed)
				cleanup_path(tmp);
		}
		if (ok) {
			++idx;
		} else {
			swap(tunnels[idx], tunnels[tunnels.size()-1]);
			tunnels.pop_back();
		}
	}
}

bool filter_tunnel_bitrelations(const sha1differentialpath& tunnel, const vector< vector<uint32> >& bitrels)
{
	vector<uint32> memask(16,0);
	bool ok = true;
	for (int t = tunnel.tbegin()+4; t < tunnel.tend()-1; ++t)
		if (t >= 0 && t < 16)
			memask[t] = tunnel.getme(t).mask;
	for (unsigned i = 0; i < bitrels.size(); ++i) {
		uint32 cummask = 0;
		for (unsigned k = 0; k < 16; ++k)
			cummask ^= (memask[k] & bitrels[i][k]);
		if (hw(cummask)&1) {
			ok = false;
			break;
		}
	}
	return ok;
}

void filter_tunnels_bitrelations(vector< sha1differentialpath >& tunnels, const vector< vector<uint32> >& bitrels)
{
	unsigned idx = 0;
	vector<uint32> memask(16);
	while (idx < tunnels.size()) {
		bool ok = true;
		const sha1differentialpath& tunnel = tunnels[idx];
		memask.assign(16, 0);
		for (int t = tunnel.tbegin()+4; t < tunnel.tend()-1; ++t)
			if (t >= 0 && t < 16)
				memask[t] = tunnel.getme(t).mask;
		for (unsigned i = 0; i < bitrels.size(); ++i) {
			uint32 cummask = 0;
			for (unsigned k = 0; k < 16; ++k)
				cummask ^= (memask[k] & bitrels[i][k]);
			if (hw(cummask)&1) {
				ok = false;
				break;
			}
		}
		if (ok) {
			++idx;
		} else {
			swap(tunnels[idx], tunnels[tunnels.size()-1]);
			tunnels.pop_back();
		}
	}
}