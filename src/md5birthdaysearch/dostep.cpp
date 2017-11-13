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

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <boost/program_options.hpp>

#include "main.hpp"

#include <hashclash/md5detail.hpp>
#include <hashclash/timer.hpp>
#include <hashclash/rng.hpp>

using namespace hashclash;
using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

unsigned load_block_fillrandom(istream& i, uint32 block[])
{
	unsigned bytes = load_block(i, block);
	unsigned len = 0;
	for (unsigned k = 0; k < 16; ++k)
		for (unsigned c = 0; c < 4; ++c)
		{
			if (len >= bytes)
				block[k] ^= (xrng128()&0xFF) << (c*8);
			++len;
		}
	return bytes;
}

unsigned load_block(istream& i, uint32 block[])
{	
	for (unsigned k = 0; k < 16; ++k)
		block[k] = 0;

	unsigned len = 0;
	char uc;
	for (unsigned k = 0; k < 16; ++k)
		for (unsigned c = 0; c < 4; ++c)
		{
			i.get(uc);
			if (i) {
				++len;
				block[k] += uint32((unsigned char)(uc)) << (c*8);
			} else {
				i.putback(uc);
				i.setstate(ios::failbit);
				return len;
			}
		}	
	return len;
}

void save_block(ostream& o, uint32 block[])
{
	for (unsigned k = 0; k < 16; ++k)
		for (unsigned c = 0; c < 4; ++c)
			o << (unsigned char)((block[k]>>(c*8))&0xFF);
}

int dostep(birthday_parameters& parameters)
{
	uint32 ihv1[4];
	uint32 ihv2[4];
	uint32 msg1[16];
	uint32 msg2[16];
	{
		ifstream if1(parameters.inputfile1.c_str(), ios::binary);
		if (!if1) {
			cerr << "Error: cannot open inputfile 1 '" << parameters.inputfile1 << "'!" << endl;
			return 1;
		}
		ifstream if2(parameters.inputfile2.c_str(), ios::binary);
		if (!if2) {
			cerr << "Error: cannot open inputfile 2 '" << parameters.inputfile2 << "'!" << endl;
			return 1;
		}
		ofstream* of1 = 0;
		ofstream* of2 = 0;
		if (parameters.modi == 0) {
			of1 = new ofstream(parameters.outputfile1.c_str(), ios::binary);
			if (!(*of1)) {
				cerr << "Error: cannot open outputfile 1 '" << parameters.outputfile1 << "'!" << endl;
				return 1;
			}
			of2 = new ofstream(parameters.outputfile2.c_str(), ios::binary);
			if (!(*of2)) {
				delete of1;
				cerr << "Error: cannot open outputfile 2 '" << parameters.outputfile1 << "'!" << endl;
				return 1;
			}
		}
		
		for (unsigned k = 0; k < 4; ++k)
			ihv1[k] = ihv2[k] = md5_iv[k];

		seed(0x12345678);
		// load, md5 and save inputfile1
		unsigned file1blocks = 0;
		while (load_block_fillrandom(if1, msg1) > 64-12) { // stop when at least 12 bytes are not used
			md5compress(ihv1, msg1);
			if (of1) save_block(*of1, msg1);
			++file1blocks;
			addseed(ihv1[0]);
		}
		// load, md5 and save inputfile2
		unsigned file2blocks = 0;
		while (load_block_fillrandom(if2, msg2) > 64-12) { // stop when at least 12 bytes are not used
			md5compress(ihv2, msg2);
			if (of2) save_block(*of2, msg2);
			++file2blocks;
			addseed(ihv2[0]);
		}
		// flush last partial message block when padding blocks are required
		if (file1blocks < file2blocks) {
			md5compress(ihv1, msg1);
			if (of1) save_block(*of1, msg1);
			++file1blocks;
		}
		if (file2blocks < file1blocks) {
			md5compress(ihv2, msg2);
			if (of2) save_block(*of2, msg2);
			++file2blocks;
		}
		// expand length outputfile1 to length outputfile2, by saving and md5'ing
		while (file1blocks < file2blocks) {
			for (unsigned k = 0; k < 16; ++k)
				msg1[k] = xrng128()+xrng128();
			md5compress(ihv1, msg1);
			if (of1) save_block(*of1, msg1);
			++file1blocks;
		}
		// expand length outputfile2 to length outputfile1, by saving and md5'ing
		while (file2blocks < file1blocks) {
			for (unsigned k = 0; k < 16; ++k)
				msg2[k] = xrng128()+xrng128();
			md5compress(ihv2, msg2);
			if (of2) save_block(*of2, msg2);
			++file2blocks;
		}

		if (of1) delete of1;
		if (of2) delete of2;
	} // close if1, if2, of1, of2
	addseed(uint32(time(0)));
	addseed(parameters.modi);

	for (unsigned k = 0; k < 4; ++k)
	{
		parameters.ihv1[k] = ihv1[k];
		parameters.ihv2[k] = ihv2[k];
	}
	for (unsigned k = 0; k < 16; ++k)
	{
		parameters.msg1[k] = msg1[k];
		parameters.msg2[k] = msg2[k];
	}

	cout << "IHV1 = {" << ihv1[0] << "," << ihv1[1] << "," << ihv1[2] << "," << ihv1[3] << "}" << endl;
	cout << "IHV1 = " << hex;
	for (unsigned k = 0; k < 4; ++k)
		for (unsigned c = 0; c < 4; ++c)
		{
			cout.width(2); cout.fill('0');
			cout << ((ihv1[k]>>(c*8))&0xFF);
		}
	cout << dec << endl << endl;

	cout << "IHV2 = {" << ihv2[0] << "," << ihv2[1] << "," << ihv2[2] << "," << ihv2[3] << "}" << endl;
	cout << "IHV2 = " << hex;
	for (unsigned k = 0; k < 4; ++k)
		for (unsigned c = 0; c < 4; ++c)
		{
			cout.width(2); cout.fill('0');
			cout << ((ihv2[k]>>(c*8))&0xFF);
		}
	cout << dec << endl << endl;

	// Start job with given parameters
	birthday(parameters);

	ofstream* of1 = 0;
	ofstream* of2 = 0;
	if (parameters.modi == 0) {
		of1 = new ofstream(parameters.outputfile1.c_str(), ios::binary | ios::app);
		if (!(*of1)) {
			cerr << "Error: cannot open outputfile 1 '" << parameters.outputfile1 << "'!" << endl;
			return 1;
		}
		of2 = new ofstream(parameters.outputfile2.c_str(), ios::binary | ios::app);
		if (!(*of2)) {
			delete of1;
			cerr << "Error: cannot open outputfile 2 '" << parameters.outputfile1 << "'!" << endl;
			return 1;
		}
		msg1[13] = collc1;
		msg1[14] = colla1;
		msg1[15] = collb1;
		msg2[13] = collc2;
		msg2[14] = colla2;
		msg2[15] = collb2;
		if (of1) save_block(*of1, msg1);
		if (of2) save_block(*of2, msg2);
	}
	if (of1) delete of1;
	if (of2) delete of2;

	return 0;
}
