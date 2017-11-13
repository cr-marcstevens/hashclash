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
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <locale>

#define MAX_LOOP_SIZE 22

template<class _Ty>    
string array_to_cppstring(const _Ty vals[], unsigned cnt)
{
	if (cnt == 0) return string("{}");
	string o = "{ " + boost::lexical_cast<string>(vals[0]);
	for (unsigned i = 1; i < cnt; ++i)
		o += ", " + boost::lexical_cast<string>(vals[i]);
	o += "}";
	return o;
}
template<class _Ty>
string vector_to_cppstring(const vector<_Ty>& vals)
{	return array_to_cppstring(&vals[0], vals.size()); }

#define LCSTRING(s) boost::lexical_cast<string>(s)
string mask2bit_to_string(const string& variable, const uint32 mask, const unsigned bit) {
	if (hw(mask)==0) return "0";
	if (hw(mask)==1) {
		for (unsigned b = 0; b < 32; ++b)
			if (mask & (1<<b)) {
				if (b == bit) return "(" + variable + "&" + LCSTRING(1<<b) + ")";
				if (b < bit) return "((" + variable + "&" + LCSTRING(1<<b) + ")<<" + LCSTRING(bit-b) + ")";
				if (b > bit) return "((" + variable + "&" + LCSTRING(1<<b) + ")>>" + LCSTRING(b-bit) + ")";
			}
	}
	if (hw(mask)<=5) {
		string res;
		for (unsigned b = 0; b < 32; ++b)
			if (mask & (1<<b)) {
				if (res.size())
					res += "+";
				res += "(" + variable + ">>" + LCSTRING(b) + ")";
			}
		return "(((" + res + ")&1)<<" + LCSTRING(bit) + ")";
	}
	return "((hw("+variable+"&"+LCSTRING(mask)+")&1)<<"+LCSTRING(bit)+")";
}

uint32 rotationmask(const vector<uint32>& rels, unsigned rot)
{
	uint32 mask = 0;
	for (unsigned b = 0; b < 32; ++b)
		if (rotate_left(rels[b], rot)&(1<<b))
			mask |= 1<<b;
	return mask;
}
uint32 fixedmask(const vector<uint32>& rels, unsigned bit)
{
	uint32 mask = 0;
	for (unsigned b = 0; b < 32; ++b)
		if (rels[b] & (1<<bit))
			mask |= 1<<b;
	return mask;
}
// rels contains for each bit 0,...,31 a mask that specifies which bits of 'variable' have to be summed
string optimized_mecond_to_string(const string& variable, const vector<uint32>& _rels)
{
	vector<uint32> rels = _rels;
	rels.resize(32,0); // sanitize
	vector< pair<uint32, unsigned> > fixed_mask_bit;
	vector< pair<uint32, unsigned> > rotation_mask_rot;
	while (rels != vector<uint32>(32,0)) {
		int bestgain = -1;
		int bestoption = -1;
		for (unsigned bit = 0; bit < 32; ++bit) {
			uint32 mask = fixedmask(rels, bit);
			int gain = int(hw(mask))*3 - 5;              //^((0-((var>>bit)&1))&mask)
			if (0 == (mask & ~uint32(0-uint32(1<<bit))))
				gain += 1; // don't need to shift 'var' by 'bit' bits
			if (gain > bestgain) {
				bestgain = gain;
				bestoption = bit;
			}
		}
		for (unsigned rot = 0; rot < 32; ++rot) {
			uint32 mask = rotationmask(rels, rot);
			int gain = int(hw(mask))*3 - 3;
			if (gain > bestgain) {
				bestgain = gain;
				bestoption = rot+32;
			}
		}
		if (bestoption == -1) {
			cerr << "don't know what to do anymore ?!" << endl;
			throw;
		}
		if (bestoption < 32) {
			unsigned bit = bestoption;
			uint32 mask = fixedmask(rels, bit);
			fixed_mask_bit.push_back( pair<uint32,unsigned>(mask, bit) );
			for (unsigned b = 0; b < 32; ++b)
				if (mask & (1<<b))
					rels[b] &= ~uint32(1<<bit);
		} else {
			unsigned rot = bestoption-32;
			uint32 mask = rotationmask(rels, rot);
			rotation_mask_rot.push_back( pair<uint32,unsigned>(mask, rot) );
			for (unsigned b = 0; b < 32; ++b)
				if (mask & (1<<b))
					rels[b] &= rotate_right(~uint32(1<<b), rot);
		}
	}
	string res = "0";
	
	for (unsigned i = 0; i < fixed_mask_bit.size(); ++i) {
		uint32 mask = fixed_mask_bit[i].first;
		unsigned bit = fixed_mask_bit[i].second;		
		if (0 == (mask & ~uint32(0-uint32(1<<bit)))) {
			//^((0-(var&(1<<bit)))&mask)
			res += "\n\t\t\t\t^((0-(" + variable + "&(1<<" + boost::lexical_cast<string>(bit) + ")))&" + boost::lexical_cast<string>(mask) + ")";
			res += " //" + boost::lexical_cast<string>(mask) + "=[!";
			bool firstone = true;
			for (unsigned b = 0; b < 32; ++b)
				if (mask & (1<<b)) {
					if (firstone) firstone = false; else res += ",";
					res += boost::lexical_cast<string>(b);
				}
			res += "!]";
		} else {
			//^((0-((var>>bit)&1))&mask)
			res += "\n\t\t\t\t^((0-((" + variable + ">>" + boost::lexical_cast<string>(bit) + ")&1))&" + boost::lexical_cast<string>(mask) + ")";
			res += " //" + boost::lexical_cast<string>(mask) + "=[!";
			bool firstone = true;
			for (unsigned b = 0; b < 32; ++b)
				if (mask & (1<<b)) {
					if (firstone) firstone = false; else res += ",";
					res += boost::lexical_cast<string>(b);
				}
			res += "!]";
		}
	}
	for (unsigned i = 0; i < rotation_mask_rot.size(); ++i) {
		uint32 mask = rotation_mask_rot[i].first;
		unsigned rot = rotation_mask_rot[i].second;
		// ^(rotate_left(var,rot)&mask)
		res += "\n\t\t\t\t^(rotate_left(" + variable + "," + boost::lexical_cast<string>(rot) + ")&" + boost::lexical_cast<string>(mask) + ")";
		res += " //" + boost::lexical_cast<string>(mask) + "=[!";
		bool firstone = true;
		for (unsigned b = 0; b < 32; ++b)
			if (mask & (1<<b)) {
				if (firstone) firstone = false; else res += ",";
				res += boost::lexical_cast<string>(b);
			}
		res += "!]";
		mask = rotate_right(mask, rot);
		res += ", " + boost::lexical_cast<string>(mask) + "=[!";
		firstone = true;
		for (unsigned b = 0; b < 32; ++b)
			if (mask & (1<<b)) {
				if (firstone) firstone = false; else res += ",";
				res += boost::lexical_cast<string>(b);
			}
		res += "!]";
	}
	res += "\n\t\t\t\t";
	return res;
}

// tunnel_analysis.cpp
void fill_tables(const sha1differentialpath& diffpath);
extern uint32 Q[85];
extern uint32 Q2[85];
extern uint32 dQ[85];
extern uint32 dT[85];
extern uint32 dR[85];
extern uint32 dF[80];

extern uint32 Qvaluemask[85];
extern uint32 Qvalue[85];
extern uint32 Qprev[85];
extern uint32 Qprev2[85];
uint32 Qvaluemask_adj[85];
vector< vector< vector<uint32> > > mainpathbitrelationsmatrix;
vector<uint32> mainmask_nextstepcorrection(16,0), mainmask_delayfreedom(16,0), mainmask_medependencies(16,0);

extern const char program_header[];
extern const char program_statistics[];
string generate_program_step(unsigned t)
{
	uint32 metmask = 0, metset1 = 0;
	vector<uint32> met_dependency_metm1(32,0);
	for (unsigned b = 0; b < 32; ++b)
		if (mainpathbitrelationsmatrix[t][b].size()) {
			/***** WARNING: DON'T FORGET TO ADD EXTRA CHECKS FOR inside m[t] bitrelations *****/
			metmask |= 1<<b; // or use mainpathbitrelationsmatrix[t][b][t] if you are lazy!!: costs freedom
			metset1 |= mainpathbitrelationsmatrix[t][b][16]&(1<<b);
			if (t>0)
				met_dependency_metm1[b] = mainpathbitrelationsmatrix[t][b][t-1];
		}
	
	
	std::ostringstream of;
	of << 
		"\nbool step_" << boost::lexical_cast<string>(t) << "()"
		"\n{"
		"\n	" << (t>=16?"//":"") << "static bool firsttime = true;"
		"\n	" << (t>=16?"//":"") << "if (firsttime) { cout << \"t" << t << " \" << endl; firsttime = false; }"
		"\n	const unsigned t = " << t << ";"
		"\n	" << (t>=16?"//":"") << "UPDATE(t);"
		"\n	" << (t>=16?"//":"") << "show_performance_data(t);"
		"\n	bool step" << boost::lexical_cast<string>(t) << "ok = false;"
		"\n	Qr30[offset+t-2] = rotate_left(Q[offset+t-2],30);"
		"\n	const uint32 mpartial2 = 0 - sha1_f1(Q[offset+t-1], Qr30[offset+t-2], Qr30[offset+t-3]) - sha1_ac[0] - Qr30[offset+t-4];"
		"\n	const uint32 metmask = " << metmask << ";"
		;
	if (t == 0) {
		// alternative beginning for t = 0
		of <<
		"\n	{"
		"\n		const uint32 metset1 = " << metset1 << ";"
		;
	} else {
		of <<
		"\n	uint32* Qtcheckptr = Qtp1valsvec[20];"
		"\n	for (uint32* Qtvalsptr = Qtp1valsvecptr[t-1]-1; Qtvalsptr >= Qtp1valsvec[t-1]; --Qtvalsptr) {"
		"\n		Q[offset+t] = *Qtvalsptr;"
		"\n		m[t-1] = *Qtvalsptr + metpartialvec[t-1];"
		"\n		const uint32 metset1 = metset1precomp[" << t-1 << "][" << t << "] ^ " << optimized_mecond_to_string("m["+boost::lexical_cast<string>(t-1)+"]", met_dependency_metm1) << ";"
		;
	}
	if (hw(Qvaluemask_adj[offset+t+1])>=hw(metmask)) {
		uint32 Qtp1mask = ~Qvaluemask_adj[offset+t+1];
		uint32 metmaskfix = metmask & Qtp1mask;
		uint32 metmasknofix = metmask ^ metmaskfix;
		Qtp1mask ^= metmaskfix;
		if (hw(Qtp1mask & (~metmask)) <= MAX_LOOP_SIZE || (metmask&~Qtp1mask)==0) {
			uint32 metmasknofixrange = 0;
			for (unsigned rr = 0; rr < 6; ++rr)
				metmasknofixrange |= metmasknofix>>rr;
			uint32 curmask = Qtp1mask & metmasknofixrange;
			of << 
				"\n		const uint32 Qtp1val = Qvalue[offset+t+1] ^ (Qprev[offset+t+1]&Q[offset+t]) ^ ((~metset1) & " << metmaskfix << ") ^ (xrng64()&" << Qtp1mask << ");"
				"\n		const uint32 mpartial = mpartial2 - rotate_left(Q[offset+t],5);"
				"\n		metpartialvec[t] = mpartial;"
				"\n		uint32* Qtp1valsptr = Qtp1valsvec[t];"
				"\n		uint32 cur = 0; const uint32 curmask = " << curmask << ";"
				"\n		do { cur -= 1; cur &= curmask;"
				"\n		//uint32* ptrbu = Qtp1valsptr;"
				"\n		uint32 Qtp1cur = 0;"
				"\n		do {"
				"\n			Qtp1cur -= 1; Qtp1cur &= " << Qtp1mask << "^curmask;"
				"\n			Q[offset+t+1] = Qtp1cur ^ cur ^ Qtp1val;"
				"\n			m[t] = Q[offset+t+1] + mpartial;"
				"\n			if (0 != ((m[t]^metset1)&" << metmasknofix << ")) break;"
				"\n			uint32 metxorfix = (m[t]^metset1) & " << metmaskfix << ";"
				"\n			*(Qtp1valsptr++) = Q[offset+t+1]^metxorfix;"
				"\n		} while (Qtp1cur != 0);"
				"\n		// loop_data_vec[40+t].add(Qtp1valsptr-ptrbu, 1<<hw(uint32(" << Qtp1mask << "^curmask)));"
				"\n		} while (cur != 0);"
				"\n		time_avg[t] += Qtp1valsptr-Qtp1valsvec[t];"
				"\n		loop_data_vec[t].add(Qtp1valsptr-Qtp1valsvec[t],1<<" << hw(Qtp1mask) << ");"
				"\n		if (Qtp1valsptr-Qtp1valsvec[t]>0) loop_data_vec[t+20].add(Qtp1valsptr-Qtp1valsvec[t],1<<" << hw(Qtp1mask) << ");"
				"\n		const uint32 checkmask_mask = " << Qtp1mask << ";"
				"\n		const uint32 checkmask_add = 0;"
				;
		} else {
			of <<
				"\n		const uint32 Qtp1val = Qvalue[offset+t+1] ^ (Qprev[offset+t+1]&Q[offset+t]) ^ ((~metset1) & " << metmaskfix << ");"
				"\n		const uint32 mpartial = mpartial2 - rotate_left(Q[offset+t],5);"
				"\n		metpartialvec[t] = mpartial;"
				"\n		uint32* Qtp1valsptr = Qtp1valsvec[t];"
				"\n		for (unsigned j = 0; j < (1<<MAX_LOOP_SIZE); ++j) {"
				"\n			Q[offset+t+1] = (xrng64()&" << Qtp1mask << ") ^ Qtp1val;"
				"\n			m[t] = Q[offset+t+1] + mpartial;"
				"\n			if (0 != ((m[t]^metset1)&" << metmasknofix << ")) " << (t>=13||t==4||t==6||t==7?"break;":"continue;") <<
				"\n			uint32 metxorfix = (m[t]^metset1) & " << metmaskfix << ";"
				"\n			*(Qtp1valsptr++) = Q[offset+t+1]^metxorfix;"
				"\n		}"
				"\n		time_avg[t] += Qtp1valsptr-Qtp1valsvec[t];"
				"\n		loop_data_vec[t].add(Qtp1valsptr-Qtp1valsvec[t],1<<MAX_LOOP_SIZE);"
				"\n		if (Qtp1valsptr-Qtp1valsvec[t]>0) loop_data_vec[20+t].add(Qtp1valsptr-Qtp1valsvec[t],1<<MAX_LOOP_SIZE);"
				"\n		const uint32 checkmask_mask = " << Qtp1mask << ";"
				"\n		const uint32 checkmask_add = 0;"
				;
		}
	} else {
		uint32 metmaskinv = ~metmask;
		uint32 Qtp1maskinv = Qvaluemask_adj[offset+t+1];
		uint32 Qtp1maskfix = Qtp1maskinv & metmaskinv;
		uint32 Qtp1masknofix = Qtp1maskinv ^ Qtp1maskfix;
		metmaskinv ^= Qtp1maskfix;

		if (hw(metmaskinv & (~Qtp1maskinv)) <= MAX_LOOP_SIZE || (Qtp1maskinv&~metmaskinv)==0) {
			uint32 Qtp1masknofixrange = 0;
			for (unsigned rr = 0; rr < 6; ++rr)
				Qtp1masknofixrange |= Qtp1masknofix>>rr;
			uint32 curmask = metmaskinv & Qtp1masknofixrange;
			of << 
				"\n		const uint32 Qtp1val = Qvalue[offset+t+1] ^ (Qprev[offset+t+1]&Q[offset+t]);"
				"\n		const uint32 metset1adj = metset1 ^ ((~Qtp1val) & " << Qtp1maskfix << ") ^ (xrng64()&" << metmaskinv << ");"
				"\n		const uint32 mpartial = mpartial2 - rotate_left(Q[offset+t],5);"
				"\n		metpartialvec[t] = mpartial;"
				"\n		uint32* Qtp1valsptr = Qtp1valsvec[t];"
				"\n		uint32 cur = 0; const uint32 curmask = " << curmask << ";"
				"\n		do { cur -= 1; cur &= curmask;"
				"\n		uint32 metcur = 0;"
				"\n		//uint32* ptrbu = Qtp1valsptr;"
				"\n		do {"
				"\n			metcur -= 1; metcur &= " << metmaskinv << "^curmask;"
				"\n			m[t] = metcur ^ cur ^ metset1adj;"
				"\n			Q[offset+t+1] = m[t] - mpartial;"
				"\n			if (0 != ((Q[offset+t+1]^Qtp1val)&" << Qtp1masknofix << ")) break;"
				"\n			uint32 Qtp1xorfix = (Q[offset+t+1]^Qtp1val) & " << Qtp1maskfix << ";"
				"\n			*(Qtp1valsptr++) = Q[offset+t+1]^Qtp1xorfix;"
				"\n		} while (metcur != 0);"
				"\n		// loop_data_vec[40+t].add(Qtp1valsptr-ptrbu, 1<<hw(uint32(" << metmaskinv << "^curmask)));"
				"\n		} while (cur != 0);"
				"\n		time_avg[t] += Qtp1valsptr-Qtp1valsvec[t];"
				"\n		loop_data_vec[t].add(Qtp1valsptr-Qtp1valsvec[t],1<<" << hw(metmaskinv) << ");"
				"\n		if (Qtp1valsptr-Qtp1valsvec[t]>0) loop_data_vec[20+t].add(Qtp1valsptr-Qtp1valsvec[t],1<<" << hw(metmaskinv) << ");"
				"\n		const uint32 checkmask_mask = " << metmaskinv << ";"
				"\n		const uint32 checkmask_add = mpartial;"
				;
		} else {
			of << 
				"\n		const uint32 Qtp1val = Qvalue[offset+t+1] ^ (Qprev[offset+t+1]&Q[offset+t]);"
				"\n		const uint32 metset1adj = metset1 ^ ((~Qtp1val) & " << Qtp1maskfix << ");"
				"\n		const uint32 mpartial = mpartial2 - rotate_left(Q[offset+t],5);"
				"\n		metpartialvec[t] = mpartial;"
				"\n		uint32* Qtp1valsptr = Qtp1valsvec[t];"
				"\n		for (unsigned j = 0; j < (1<<MAX_LOOP_SIZE); ++j) {"
				"\n			m[t] = (xrng64() & " << metmaskinv << ") ^ metset1adj;"
				"\n			Q[offset+t+1] = m[t] - mpartial;"
				"\n			if (0 != ((Q[offset+t+1]^Qtp1val)&" << Qtp1masknofix << ")) " << (t>=13||t==4||t==6||t==7?"break;":"continue;") <<
				"\n			uint32 Qtp1xorfix = (Q[offset+t+1]^Qtp1val) & " << Qtp1maskfix << ";"
				"\n			*(Qtp1valsptr++) = Q[offset+t+1]^Qtp1xorfix;"
				"\n		}"
				"\n		time_avg[t] += Qtp1valsptr-Qtp1valsvec[t];"
				"\n		loop_data_vec[t].add(Qtp1valsptr-Qtp1valsvec[t],1<<MAX_LOOP_SIZE);"
				"\n		if (Qtp1valsptr-Qtp1valsvec[t]>0) loop_data_vec[20+t].add(Qtp1valsptr-Qtp1valsvec[t],1<<MAX_LOOP_SIZE);"
				"\n		const uint32 checkmask_mask = " << metmaskinv << ";"
				"\n		const uint32 checkmask_add = mpartial;"
				;
		}
	}
	of << 
		"\n		if (Qtp1valsptr == Qtp1valsvec[t]) " << (t==0?"return step0ok;":"continue;") <<
		"\n		step" << boost::lexical_cast<string>(t) << "ok = true;"
		"\n		//process_checkmask(Qtp1valsvec[t], Qtp1valsptr, checkmask_mask, checkmask_add);"
		"\n		//*(Qtcheckptr++) = *Qtvalsptr;"
		"\n		//continue;"
		;
	if (t>0) {
		for (unsigned i = t+1; i < 16; ++i) {
			vector<uint32> dependency_metm1(32,0);
			for (unsigned b = 0; b < 32; ++b)
				if (mainpathbitrelationsmatrix[i][b].size())
					dependency_metm1[b] = mainpathbitrelationsmatrix[i][b][t-1];
			of << "\n		metset1precomp[" << t << "][" << i << "] = metset1precomp[" << t-1 << "][" << i << "] ^ " << optimized_mecond_to_string("m["+boost::lexical_cast<string>(t-1)+"]", dependency_metm1) << ";";
		}
	}
	of <<
		"\n		Qtp1valsvecptr[t] = Qtp1valsptr;"
		"\n		step_" << t+1 << "();"
		"\n	}"
		;
	if (t>0) {
		of << 
			"\n	if (Qtcheckptr-Qtp1valsvec[20] > 0) {"
			"\n		loop_data_vec[16+t].add(Qtcheckptr-Qtp1valsvec[20], Qtp1valsvecptr[t-1]-Qtp1valsvec[t-1]);"
			"\n		process_checkmask(Qtp1valsvec[20], Qtcheckptr);"
			"\n	}"
			;
	}
	of <<
		"\n	return step" << boost::lexical_cast<string>(t) << "ok;"
		"\n}"
		;
	return of.str();
}

void generate_program()
{
	uint32 mxordiff[80];
	for (unsigned t = 0; t < 80; ++t) {
		mxordiff[t] = maindiffpath.getme(t).mask;
		if (t >= 16 && mxordiff[t] != rotate_left(mxordiff[t-3] ^ mxordiff[t-8] ^ mxordiff[t-14] ^ mxordiff[t-16], 1)) {
			cerr << "mxordiff inconsistent @ " << t << endl;
			throw;
		}
	}
	fill_tables(maindiffpath);
	vector< vector< vector<uint32> > > pathbitrelationsmatrix;
	mespace_to_pathbitrelationsmatrix(mainmespace, pathbitrelationsmatrix);
	mainpathbitrelationsmatrix = pathbitrelationsmatrix;
	fill_tables(maindiffpath);
	sha1differentialpath maindiffpath_adj = maindiffpath;
	
	vector<uint32> mask_nextstepcorrection(16,0), mask_delayfreedom(16,0), mask_medependencies(16,0);
	for (unsigned t = 0; t < 16; ++t)
		for (unsigned b = 0; b < 32; ++b)
			if (t+1 < 16 && pathbitrelationsmatrix[t+1][b].size())
				mask_medependencies[t] |= pathbitrelationsmatrix[t+1][b][t];
	mainmask_nextstepcorrection = mask_nextstepcorrection;
	mainmask_delayfreedom = mask_delayfreedom;
	mainmask_medependencies = mask_medependencies;

	{
		// fill tables with adjusted bitconditions
		// however! keep the adjusted Qvaluemask in a seperate variable: Qvaluemask_adj
		uint32 Qtmp[85];
		memcpy(Qtmp, Qvaluemask, sizeof(uint32)*85);
		fill_tables(maindiffpath_adj);
		memcpy(Qvaluemask_adj, Qvaluemask, sizeof(uint32)*85);
		memcpy(Qvaluemask, Qtmp, sizeof(uint32)*85);
	}

	std::string pathbitrelationsmatrixstring;
	{
		std::stringstream tmp;
		{
			boost::archive::text_oarchive oa(tmp);
			oa << boost::serialization::make_nvp("pathbitrelationsmatrix", pathbitrelationsmatrix);
		} // oa destructor
		pathbitrelationsmatrixstring = tmp.str();
		string::size_type pos = 0;
		while (pos < pathbitrelationsmatrixstring.size()) {
			if (pathbitrelationsmatrixstring[pos] == '\n' || pathbitrelationsmatrixstring[pos] == '\r')
				pathbitrelationsmatrixstring.erase(pos,1);
			else
				++pos;
		}
	}

	ofstream of("collfind.cpp");
	of
		<< program_header <<
		"\n#include <hashclash/saveload.hpp>"
		"\n"
		"\n#define MAX_LOOP_SIZE " << MAX_LOOP_SIZE <<
		"\n"
		"\nuint64 tendcount = 0;"
		"\nuint64 testcounts[1<<16];"
		"\nuint32 m_diff[80] = " << array_to_cppstring(mxordiff, 80) << ";"
		"\nconst uint32 dQ[85] = " << array_to_cppstring(dQ, 85) << ";"
		"\nconst uint32 Qvaluemask[85] = " << array_to_cppstring(Qvaluemask, 85) << ";"
		"\nconst uint32 Qvaluemask_adj[85] = " << array_to_cppstring(Qvaluemask_adj, 85) << ";"
		"\nconst uint32 Qvalue[85] = " << array_to_cppstring(Qvalue, 85) << ";"
		"\nconst uint32 Qprev[85] = " << array_to_cppstring(Qprev, 85) << ";"
		"\nconst uint32 Qprev2[85] = " << array_to_cppstring(Qprev2, 85) << ";"
		"\nconst uint32 dF[80] = " << array_to_cppstring(dF, 80) << ";"
		"\nvector< vector< vector<uint32> > > pathbitrelationsmatrix;"
		"\nchar pathbitrelationsmatrixstring[] = \"" << pathbitrelationsmatrixstring << "\";"
		"\nuint32 Qtp1valsvec[30][1<<22];"
		"\nuint32* Qtp1valsvecptr[30];"
		"\nuint32 metpartialvec[30];"
		"\n"
		<< program_statistics <<
		"\nvoid step0();"
		"\nint main(int argc, char** argv) {"
		"\n	{ stringstream tmp; tmp.str(string(pathbitrelationsmatrixstring)); boost::archive::text_iarchive ia(tmp); ia >> boost::serialization::make_nvp(\"pathbitrelationsmatrix\", pathbitrelationsmatrix); }"
		"\n	for (unsigned t = 0; t < 16; ++t)"
		"\n		for (unsigned b = 0; b < 32; ++b)"
		"\n			if (pathbitrelationsmatrix[t][b].size()!=0 && pathbitrelationsmatrix[t][b][t] != 1<<b) {"
		"\n				cout << \"Inside message relation @ t=\" << t << \" b=\" << b << \": \";"
		"\n				for (unsigned b2 = 0; b2 < 32; ++b2)"
		"\n					if (pathbitrelationsmatrix[t][b][t]&(1<<b2))"
		"\n						cout << \"(\" << t << \",\" << b2 << \")\";"
		"\n				cout << endl;"
		"\n			}"
		"\n	cout << \"Starting collision search...\" << endl;"
		"\n	for (unsigned i = 0; i < (1<<16); ++i) testcounts[i] = 0;"
		"\n	while (true) {"
		"\n		try {"
		"\n			step0();"
		"\n		} catch (std::exception&e) { cerr << e.what() << endl; } catch (...) {}"
		"\n	}"
		"\n}"
		"\nbool test_first_block(const uint32 block[16]) {"
		"\n	static uint32 ihv[5];"
		"\n	static uint32 me[80];"
		"\n	memcpy(ihv, sha1_iv, 5*4);"
		"\n	memcpy(me, block, 4*16);"
		"\n	for (unsigned i = 16; i < 80; ++i)"
		"\n		me[i]=rotate_left(me[i-3] ^ me[i-8] ^ me[i-14] ^ me[i-16], 1);"
		"\n	sha1compress_me(ihv, me);"
		"\n	Q2[0] = Q[0] = rotate_right(ihv[4], 30);"
		"\n	Q2[1] = Q[1] = rotate_right(ihv[3], 30);"
		"\n	Q2[2] = Q[2] = rotate_right(ihv[2], 30);"
		"\n	Q2[3] = Q[3] = ihv[1];"
		"\n	Q2[4] = Q[4] = ihv[0];"
		"\n	if ((Q[0]&Qvaluemask_adj[0])!=Qvalue[0]) return false;"
		"\n	if ((Q[1]&Qvaluemask_adj[1])!=(Qvalue[1]^(Q[0]&Qprev[1]))) return false;"
		"\n	if ((Q[2]&Qvaluemask_adj[2])!=(Qvalue[2]^(Q[1]&Qprev[2])^(Q[0]&Qprev2[2]))) return false;"
		"\n	if ((Q[3]&Qvaluemask_adj[3])!=(Qvalue[3]^(Q[2]&Qprev[3])^(Q[1]&Qprev2[3]))) return false;"
		"\n	if ((Q[4]&Qvaluemask_adj[4])!=(Qvalue[4]^(Q[3]&Qprev[4])^(Q[2]&Qprev2[4]))) return false;"
		"\n	return true;"
		"\n}"
		"\n"
		"\nbool step_0();"
		"\nvoid step0() {"
		"\n	{UPDATE(79);"
		"\n	static uint64 fbcnt = 0;"
		"\n	while (true) {"
		"\n		for (unsigned i = 0; i < 16; ++i)"
		"\n			firstmsg[i] = xrng128();"
		"\n		if (hw(++fbcnt)==1) cout << \"(fb:\" << fbcnt << \")\" << flush;"
		"\n		if (test_first_block(firstmsg))"
		"\n			break;"
		"\n	}"
		"\n	for (unsigned i = 0; i < 5; ++i) {"
		"\n		Qr30[i] = rotate_left(Q[i], 30);"
		"\n		Q2r30[i] = rotate_left(Q2[i], 30);"
		"\n	}"
		"\n	}//UPDATE(79);"
		"\n	//cout << \".\" << flush;"
		"\n"
		;
	for (unsigned j = 0; j < 16; ++j) {
		uint32 metset1j = 0;
	        for (unsigned b = 0; b < 32; ++b)
        	        if (pathbitrelationsmatrix[j][b].size()) {
                	        metset1j |= pathbitrelationsmatrix[j][b][16]&(1<<b);
              		}
		of << "\n	metset1precomp[0]["<<j<<"]=" << metset1j << ";";
	}
	of <<
		"\n 	step_0();"
		"\n}" << flush;
	;
	of << 
		"\nvoid step_16()"
		"\n{"
		"\n	UPDATE(16);"
		"\n	for (int t = 16-5; t < 16; ++t) {"
		"\n		Q2[offset+t] = Q[offset+t] + dQ[offset+t];"
		"\n		//m2[t] = m[t] ^ m_diff[t];"
		"\n	}"
		"\n	for (uint32* Qtvalsptr = Qtp1valsvecptr[15]-1; Qtvalsptr >= Qtp1valsvec[15]; --Qtvalsptr) "
		"\n	{"
		"\n		avg_data_vec[40].add_cnt(1);"
		"\n		Q[offset+16] = *Qtvalsptr;"
		"\n		Q2[offset+16] = Q[offset+16] + dQ[offset+16];"
		"\n		m[15] = *Qtvalsptr + metpartialvec[15];"
		"\n		m2[15] = m[15] ^ m_diff[15];"
		"\n		int t = 16;"
		"\n		for (; t < 20; ++t) {"
		"\n			m[t]=rotate_left(m[t-3] ^ m[t-8] ^ m[t-14] ^ m[t-16], 1);"
		"\n			m2[t] = m[t] ^ m_diff[t];"
		"\n			Q[offset+t+1] = m[t] + sha1_f1(Q[offset+t-1], rotate_left(Q[offset+t-2],30), rotate_left(Q[offset+t-3],30))"
		"\n				+ sha1_ac[0] + rotate_left(Q[offset+t],5) + rotate_left(Q[offset+t-4],30);"
		"\n			Q2[offset+t+1] = m2[t] + sha1_f1(Q2[offset+t-1], rotate_left(Q2[offset+t-2],30), rotate_left(Q2[offset+t-3],30))"
		"\n				+ sha1_ac[0] + rotate_left(Q2[offset+t],5) + rotate_left(Q2[offset+t-4],30);"
		"\n			if (naf(Q2[offset+t+1]-Q[offset+t+1]).mask != naf(dQ[offset+t+1]).mask) break;"
		"\n		}"
		"\n		if (t < 20) continue;"
		"\n		for (; t < 33; ++t) {"
		"\n			m[t]=rotate_left(m[t-3] ^ m[t-8] ^ m[t-14] ^ m[t-16], 1);"
		"\n			m2[t] = m[t] ^ m_diff[t];"
		"\n			Q[offset+t+1] = m[t] + sha1_f2(Q[offset+t-1], rotate_left(Q[offset+t-2],30), rotate_left(Q[offset+t-3],30))"
		"\n				+ sha1_ac[1] + rotate_left(Q[offset+t],5) + rotate_left(Q[offset+t-4],30);"
		"\n			Q2[offset+t+1] = m2[t] + sha1_f2(Q2[offset+t-1], rotate_left(Q2[offset+t-2],30), rotate_left(Q2[offset+t-3],30))"
		"\n				+ sha1_ac[1] + rotate_left(Q2[offset+t],5) + rotate_left(Q2[offset+t-4],30);"
		"\n		}"
		"\n		for (int t2 = 1; t2 <= t+1; ++t2) {"
		"\n			avg_data_vec[t2].add_cnt(1);"
		"\n			if (Q2[offset+t2]-Q[offset+t2]==dQ[offset+t2])"
		"\n				avg_data_vec[t2].add_sum(1);"
		"\n		}"
		"\n		if (t==33 && Q2[offset+33]==Q[offset+33] && Q2[offset+32]==Q[offset+32] && Q2[offset+31]==Q[offset+31] && Q2[offset+30]==Q[offset+30] && Q2[offset+29]==Q[offset+29]) {"
		"\n			cout << \"!\" << endl;"
		"\n			avg_data_vec[40].add_sum(1);"
		"\n			++time_avg[40];"
		"\n			for (int j = 15; j <= 33; ++j)"
		"\n				cout << \"dQ\" << j << \"=\" << sdr(Q[offset+j],Q2[offset+j]) << \"\\tdM\" << j << \"=\" << sdr(m[j],m2[j]) << endl;"
		"\n//			throw std::runtime_error(\".\");"
		"\n		}"
		"\n	}"
		"\n	static uint64 throwcnt = 0;"
		"\n//	if ((++throwcnt & 0xFFFFFF)==0) throw std::runtime_error(\".\");"
		
		"\n}" << endl;
	for (int t = 15; t >= 0; --t) 
		of << generate_program_step(t) << flush;
	of.close();
	system("./sha1_collfind.sh collfind.cpp");
	exit(0);
}

const char program_header[] = 
	"\n#include <cmath>"
	"\n#include <iomanip>"
	"\n#include <algorithm>"
	"\n#include <map>"
	"\n#include <fstream>"
	"\n#include <stdexcept>"
	"\n#define SHA1DETAIL_INLINE_IMPL"
	"\n#include <hashclash/sha1detail.hpp>"
	"\n#include <hashclash/sha1differentialpath.hpp>"
	"\n#include <hashclash/booleanfunction.hpp>"
	"\n#include <hashclash/rng.hpp>"
	"\n#include <hashclash/timer.hpp>"
	"\n#include <hashclash/bestof.hpp>"
	"\n#include <hashclash/progress_display.hpp>"
	"\n#include <hashclash/sha1messagespace.hpp>"
	"\n#include <boost/lexical_cast.hpp>"
	"\nusing namespace std;"
	"\nusing namespace hashclash;"
	"\n#define CPUPERFORMANCE"
	"\n#ifdef CPUPERFORMANCE"
	"\n#include <hashclash/cpuperformance.hpp>"
	"\nuint64 cpu_step_t[80];"
	"\n#define UPDATE(s) update_performance_counter __update(cpu_step_t[s]);"
	"\n#ifdef __GNUC__"
	"\n#include <sched.h>"
	"\n#endif"
	"\n#else"
	"\n#define UPDATE(s)"
	"\n#endif"
	"\nconst int offset = 4;"
	"\nuint32 firstmsg[16];"
	"\nuint32 m[80];"
	"\nuint32 m2[80];"
	"\nuint32 Q[85];"
	"\nuint32 Q2[85];"
	"\nuint32 Qr30[85];"
	"\nuint32 Q2r30[85];"
	"\nuint32 metset1precomp[16][16];"
	""
	"\ntemplate<int tstart, int tend>"
	"\ninline void compute_mt() {"
	"\n        for (int t = tstart; t < tend; ++t)"
	"\n                m[t] = rotate_left(m[t-3] ^ m[t-8] ^ m[t-14] ^ m[t-16], 1);"
	"\n}"
	"\ntemplate<int r> inline uint32 sha1_f(uint32 b, uint32 c, uint32 d);"
	"\ntemplate<> inline uint32 sha1_f<0>(uint32 b, uint32 c, uint32 d) { return sha1_f1(b,c,d) + 0x5A827999; }"
	"\ntemplate<> inline uint32 sha1_f<1>(uint32 b, uint32 c, uint32 d) { return sha1_f2(b,c,d) + 0x6ED9EBA1; }"
	"\ntemplate<> inline uint32 sha1_f<2>(uint32 b, uint32 c, uint32 d) { return sha1_f3(b,c,d) + 0x8F1BBCDC; }"
	"\ntemplate<> inline uint32 sha1_f<3>(uint32 b, uint32 c, uint32 d) { return sha1_f4(b,c,d) + 0xCA62C1D6; }"
	"\n"
	"\ntemplate<int tstart, int tend>"
	"\ninline void compute_qt() {"
	"\n        if (tstart >= tend) return;"
	"\n        uint32 a = Q[offset+tstart], b = Q[offset+tstart-1], c = Qr30[offset+tstart-2], d = Qr30[offset+tstart-3], e = Qr30[offset+tstart-4];"
	"\n        e += rotate_left(a, 5) + sha1_f<(tstart+0)/20>(b,c,d) + m[tstart+0]; Q[offset+tstart+1] = e; Qr30[offset+tstart-1] = b = rotate_left(b, 30);"
	"\n        if (tstart+1 == tend) return;"
	"\n        d += rotate_left(e, 5) + sha1_f<(tstart+1)/20>(a,b,c) + m[tstart+1]; Q[offset+tstart+2] = d; Qr30[offset+tstart+0] = a = rotate_left(a, 30);"
	"\n        if (tstart+2 == tend) return;"
	"\n        c += rotate_left(d, 5) + sha1_f<(tstart+2)/20>(e,a,b) + m[tstart+2]; Q[offset+tstart+3] = c; Qr30[offset+tstart+1] = e = rotate_left(e, 30);"
	"\n        if (tstart+3 == tend) return;"
	"\n        b += rotate_left(c, 5) + sha1_f<(tstart+3)/20>(d,e,a) + m[tstart+3]; Q[offset+tstart+4] = b; Qr30[offset+tstart+2] = d = rotate_left(d, 30);"
	"\n        if (tstart+4 == tend) return;"
	"\n        a += rotate_left(b, 5) + sha1_f<(tstart+4)/20>(c,d,e) + m[tstart+4]; Q[offset+tstart+5] = a; Qr30[offset+tstart+3] = c = rotate_left(c, 30);"
	"\n        if (tstart+5 == tend) return;"
	"\n        e += rotate_left(a, 5) + sha1_f<(tstart+5)/20>(b,c,d) + m[tstart+5]; Q[offset+tstart+6] = e; Qr30[offset+tstart+4] = b = rotate_left(b, 30);"
	"\n        if (tstart+6 == tend) return;"
	"\n        d += rotate_left(e, 5) + sha1_f<(tstart+6)/20>(a,b,c) + m[tstart+6]; Q[offset+tstart+7] = d; Qr30[offset+tstart+5] = a = rotate_left(a, 30);"
	"\n        if (tstart+7 == tend) return;"
	"\n        c += rotate_left(d, 5) + sha1_f<(tstart+7)/20>(e,a,b) + m[tstart+7]; Q[offset+tstart+8] = c; Qr30[offset+tstart+6] = e = rotate_left(e, 30);"
	"\n        if (tstart+8 == tend) return;"
	"\n        b += rotate_left(c, 5) + sha1_f<(tstart+8)/20>(d,e,a) + m[tstart+8]; Q[offset+tstart+9] = b; Qr30[offset+tstart+7] = d = rotate_left(d, 30);"
	"\n        if (tstart+9 == tend) return;"
	"\n        a += rotate_left(b, 5) + sha1_f<(tstart+9)/20>(c,d,e) + m[tstart+9]; Q[offset+tstart+10] = a; Qr30[offset+tstart+8] = c = rotate_left(c, 30);"
	"\n        if (tstart+10 == tend) return;"
	"\n        e += rotate_left(a, 5) + sha1_f<(tstart+10)/20>(b,c,d) + m[tstart+10]; Q[offset+tstart+11] = e; Qr30[offset+tstart+9] = b = rotate_left(b, 30);"
	"\n        if (tstart+11 == tend) return;"
	"\n        d += rotate_left(e, 5) + sha1_f<(tstart+11)/20>(a,b,c) + m[tstart+11]; Q[offset+tstart+12] = d; Qr30[offset+tstart+10] = a = rotate_left(a, 30);"
	"\n        if (tstart+12 == tend) return;"
	"\n        c += rotate_left(d, 5) + sha1_f<(tstart+12)/20>(e,a,b) + m[tstart+12]; Q[offset+tstart+13] = c; Qr30[offset+tstart+11] = e = rotate_left(e, 30);"
	"\n        if (tstart+13 == tend) return;"
	"\n        b += rotate_left(c, 5) + sha1_f<(tstart+13)/20>(d,e,a) + m[tstart+13]; Q[offset+tstart+14] = b; Qr30[offset+tstart+12] = d = rotate_left(d, 30);"
	"\n        if (tstart+14 == tend) return;"
	"\n        a += rotate_left(b, 5) + sha1_f<(tstart+14)/20>(c,d,e) + m[tstart+14]; Q[offset+tstart+15] = a; Qr30[offset+tstart+13] = c = rotate_left(c, 30);"
	"\n        if (tstart+15 < tend) throw; // prevent stupid mistake"
	"\n}"

	;
const char program_statistics[] = 
	"\nstruct avg_data {"
	"\n	uint64 avgcnt, avgsum;"
	"\n	avg_data(): avgcnt(0), avgsum(0) {}"
	"\n	void add_ok() {"
	"\n		++avgcnt;"
	"\n		++avgsum;"
	"\n	}"
	"\n	void add_bad() {"
	"\n		++avgcnt;"
	"\n	}"
	"\n	void add_cnt(uint64 toadd = 1) {"
	"\n		avgcnt += toadd;"
	"\n	}"
	"\n	void add_sum(uint64 toadd = 1) {"
	"\n		avgsum += toadd;"
	"\n	}"
	"\n	void show(const std::string& name) {"
	"\n		double avg = double(avgsum)/double(avgcnt);"
	"\n		double logavg = log(avg)/log(2.0);"
	"\n		if (logavg >= 20 || logavg < -2.0) "
	"\n	                cout << name << \":\t2^(\" << logavg << \")\";"
	"\n		else"
	"\n			cout << name << \":\t\" << avg;"
	"\n		cout << \"\t(#=\" << avgcnt << \")\" << endl;"
	"\n	}"
	"\n};"
	"\nvector<avg_data> avg_data_vec(80);"
	"\n"
	"\nstruct loop_data {"
	"\n        uint64 proc_cnt[101];"
	"\n        uint64 totcnt;"
	"\n        bool fulldetail;"
	"\n        loop_data() { "
	"\n        	fulldetail = false;"
	"\n                for (unsigned i = 0; i < 101; ++i)"
	"\n                        proc_cnt[i] = 0;"
	"\n                totcnt = 0;"
	"\n        }"
	"\n        void add(uint32 ok, uint32 cnt) {"
	"\n                unsigned proc = round(float(ok*100)/float(cnt));"
	"\n                if (proc > 100) throw;"
	"\n                ++proc_cnt[proc];"
	"\n                ++totcnt;"
	"\n        }"
	"\n        void add(unsigned proc) {"
	"\n                if (proc > 100) throw;"
	"\n                ++proc_cnt[proc];"
	"\n                ++totcnt;"
	"\n        }"
	"\n        void show(const std::string& name) {"
	"\n                cout << name << \": \";"
	"\n                uint64 mincnt = totcnt >> 3;"
	"\n                uint64 remcnt = totcnt;"
	"\n                if (mincnt == 0) mincnt = 1;"
	"\n                for (unsigned i = 0; i < 101; ++i)"
	"\n                        if (proc_cnt[i] >= mincnt || i == 0 || i == 100 || (proc_cnt[i]>0 && fulldetail)) {"
	"\n                                cout << \"(\" << i << \"% occurs with p=\" << double(proc_cnt[i])/double(totcnt) << \") \";"
	"\n                                remcnt -= proc_cnt[i];"
	"\n                        }"
	"\n                cout << \"(other p=\" << double(remcnt)/double(totcnt) << \")\" << endl;"
	"\n        }"
	"\n};"
	"\nvector<loop_data> loop_data_vec(80);"
	"\n"
	"\nstruct checkmaskcnt_type {"
	"\n        void process(uint32* Qtbegin, uint32* Qtend, uint32 mask = 0xFFFFFFFF) {"
	"\n                if (Qtend-Qtbegin < 2) return;"
	"\n                bits.resize(32);"
	"\n                for (unsigned b = 0; b < 32; ++b) {"
	"\n			   if ((mask & (1<<b))==0) continue;"
	"\n                        uint32 cnt = 0;"
	"\n                        for (uint32* p = Qtbegin; p != Qtend; ++p)"
	"\n                                if (*p & (1<<b))"
	"\n                                        ++cnt;"
	"\n//                      if (cnt > ((Qtend-Qtbegin+1)>>1))"
	"\n//                              cnt = (Qtend-Qtbegin)-cnt;"
	"\n                        bits[b].add(cnt, Qtend-Qtbegin);"
	"\n                }"
	"\n        }"
	"\n        void show(const std::string& name) {"
	"\n                if (bits.empty()) return;"
	"\n                cout << name << \":\" << endl;"
	"\n                for (unsigned b = 0; b < 32; ++b) {"
	"\n			   if (bits[b].totcnt == 0) continue;"
	"\n                        bits[b].show(\"checkmaskbit \" + boost::lexical_cast<string>(b));"
	"\n                }"
	"\n        }"
	"\n        vector<loop_data> bits;"
	"\n};"
	"\ncheckmaskcnt_type checkmaskcnt;"
	"\nvoid process_checkmask(uint32* Qtbegin, uint32* Qtend, uint32 mask = 0xFFFFFFFF, uint32 add = 0)"
	"\n{"
	"\n	   if (add) { for (uint32* ptr = Qtbegin; ptr != Qtend; ++ptr) *ptr += add; }"
	"\n        checkmaskcnt.process(Qtbegin, Qtend,mask);"
	"\n	   if (add) { for (uint32* ptr = Qtbegin; ptr != Qtend; ++ptr) *ptr -= add; }"
	"\n}"
	"\n"
	"\nvector<uint64> time_avg(1024);"
	"\nvector<uint64> step_cnt(20,0);"
	"\ntimer spd_sw(true), runtime_sw(true), restart_sw(true);"
	"\n#ifdef CPUPERFORMANCE"
	"\nuint64 cpu_timestamp_begin = cpu_timestamp();"
	"\n#endif"
	"\nvoid show_performance_data(const int correct_cpu_step_t = -1)"
	"\n{"
	"\n	if (spd_sw.time() < 15) return;"
	"\n	if (time_avg[40] == 0 && restart_sw.time() > 60) {"
	"\n		restart_sw.start();"
	"\n		throw std::runtime_error(\"restart\");"
	"\n	}"
	"\n	spd_sw.start();"
	"\n	cout << endl;"
	"\n#ifdef CPUPERFORMANCE"
	"\n	uint64 cputime = cpu_timestamp();"
	"\n	uint64 cpuspend = cputime-cpu_timestamp_begin;"
	"\n	if (correct_cpu_step_t != -1) {"
	"\n		for (int t = 0; t <= correct_cpu_step_t; ++t)"
	"\n			cpu_step_t[t] += cputime;"
	"\n	}"
	"\n	for (unsigned i = 0; i < 16; ++i)"
	"\n		if (cpu_step_t[i]-cpu_step_t[i+1])"
	"\n			cout << \"t=\" << i << \"\t\" << round(double(cpu_step_t[i]-cpu_step_t[i+1])*1000.0/double(cpuspend)) << \"\t\" << cpu_step_t[i]-cpu_step_t[i+1] << endl;"
	"\n     for (unsigned i = 16; i < 80; ++i) {"
	"\n             if (cpu_step_t[i] > (cpuspend<<1))"
	"\n                     cout << \"i=\" << i << \"\t\" << round(double(cpu_step_t[i]+cputime)*1000.0/double(cpuspend)) << \"\t\" << (cpu_step_t[i]+cputime) << endl;"
	"\n             else if (cpu_step_t[i])"
	"\n                     cout << \"i=\" << i << \"\t\" << round(double(cpu_step_t[i])*1000.0/double(cpuspend)) << \"\t\" << (cpu_step_t[i]) << endl;"
	"\n     }"
	"\n	if (correct_cpu_step_t != -1) {"
	"\n		for (int t = 0; t <= correct_cpu_step_t; ++t)"
	"\n			cpu_step_t[t] -= cputime;"
	"\n	}"
	"\n#endif"
	"\n	for (unsigned i = 0; i < 80; ++i)"
	"\n		if (loop_data_vec[i].totcnt)"
	"\n			loop_data_vec[i].show(\"loop \" + boost::lexical_cast<string>(i) + \" stats\");"
	"\n	for (unsigned i = 0; i < 80; ++i)"
	"\n		if (avg_data_vec[i].avgcnt)"
	"\n			avg_data_vec[i].show(\"avg \" + boost::lexical_cast<string>(i) + \" stats\");"
	"\n	for (unsigned i = 0; i < time_avg.size(); ++i)"
	"\n		if (time_avg[i])"
	"\n			cout << \"timeavg \" << i << \": \t2^\" << log(double(time_avg[i])/double(runtime_sw.time()))/log(2.0) << \"#/s\" << endl;"
	"\n	for (unsigned i = 0; i < 20; ++i)"
	"\n		if (step_cnt[i])"
	"\n			cout << \"(\" << i << \":\" << step_cnt[i] << \")\t\";"
	"\n	cout << endl;"
	"\n	checkmaskcnt.show(\"checkmask\");"
	"\n}"
	"\n"
	"\nvoid verify_block(unsigned tend = 16)"
	"\n{"
	"\n	static uint32 cnt = 0;"
	"\n	if (hw(++cnt)!=1) return;"
	"\n	UPDATE(76);"
	"\n	for (unsigned t = 16; t < tend; ++t)"
	"\n		m[t]=rotate_left(m[t-3] ^ m[t-8] ^ m[t-14] ^ m[t-16], 1);"
	"\n	for (unsigned t = 0; t < tend && t < 20; ++t)"
	"\n		if (rotate_left(Q[offset+t-2],30)!=Qr30[offset+t-2]) {"
	"\n			cout << \"Precomp Qr30 error @ t=\" << t << endl;"
	"\n			exit(0);"
	"\n		}"
	"\n	for (unsigned t = 0; t < tend && t < 20; ++t) {"
	"\n		bool terror = false;"
	"\n		uint32 Qtp1 = m[t] + sha1_f1(Q[offset+t-1], rotate_left(Q[offset+t-2],30), rotate_left(Q[offset+t-3],30)) + sha1_ac[0] + rotate_left(Q[offset+t],5) + rotate_left(Q[offset+t-4],30);"
	"\n		if (Qtp1 != Q[offset+t+1]) {"
	"\n			cout << \"Hash error @ t=\" << t << \": \" << hex << (Qtp1^Q[offset+t+1]) << dec << endl;"
	"\n			terror = true;"
	"\n		}   "
	"\n		if ((Q[offset+t+1]^Qvalue[offset+t+1]^(Qprev[offset+t+1]&Q[offset+t])) & Qvaluemask[offset+t+1]) {"
	"\n			cout << \"Verify error @ t=\" << t << \": \";"
	"\n			cout << hex << ((Q[offset+t+1]^Qvalue[offset+t+1]^(Qprev[offset+t+1]&Q[offset+t])) & Qvaluemask[offset+t+1]) << dec << endl;"
	"\n			terror = true;"
	"\n		}"
	"\n		if (t < 16) {"
	"\n			for (unsigned b = 0; b < 32; ++b)"
	"\n				if (pathbitrelationsmatrix[t][b].size()) {"
	"\n					uint32 r = pathbitrelationsmatrix[t][b][16]&1;"
	"\n					for (unsigned j = 0; j <= t; ++j)"
	"\n						r ^= m[j]&pathbitrelationsmatrix[t][b][j];"
	"\n					if (hw(r)&1) {"
	"\n 						cout << \"Message relation error @ t=\" << t << \" b=\" << b << \": \";"
	"\n						for (unsigned j = 0; j <= t; ++j)"
	"\n							for (unsigned b2 = 0; b2 < 32; ++b2)"
	"\n								if (pathbitrelationsmatrix[t][b][j]&(1<<b2))"
	"\n									cout << \"(\" << j << \",\" << b2 << \")\";"
	"\n						cout << endl;"
	"\n						terror = true;"
	"\n					}"
	"\n				}"
	"\n		}"
	"\n		if (terror) exit(1);"
	"\n	}"
	"\n}"
	;
	