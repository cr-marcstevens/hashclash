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

#include <time.h>
#include <iostream>
#include "rng.hpp"

#if defined(__linux__) || defined (__FreeBSD__)
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>    // open
#include <unistd.h>   // read, close
namespace hashclash {
	void getosrnd(uint32 buf[256])
	{
		int fd;
		if ((fd = open("/dev/urandom", O_RDONLY)) < 0) return;
		read(fd, reinterpret_cast<char*>(buf), 256*4);
		close(fd);
	}
}
#elif defined(WIN32)
#include <windows.h>
#include <wincrypt.h>
namespace hashclash {
	void getosrnd(uint32 buf[256])
	{
		HCRYPTPROV g_hCrypt;
		if(!CryptAcquireContext(&g_hCrypt,NULL,NULL,PROV_RSA_FULL,CRYPT_VERIFYCONTEXT))
			return;
		CryptGenRandom(g_hCrypt,256*4,reinterpret_cast<BYTE*>(buf));
		CryptReleaseContext(g_hCrypt,0);
	}
}
#else
namespace hashclash {
	void getosrnd(uint32 buf[256])
	{
		std::cout << "Warning: no OS randomness!" << std::endl;
		// without OS randomness
		// use fact that buf is uninitialized
	}
}
#endif

namespace hashclash {

	void getosrnd(uint32 buf[256]);

	uint32 seedd;
	uint32 seed32_1;
	uint32 seed32_2;
	uint32 seed32_3;
	uint32 seed32_4;

	void seed(uint32 s)
	{
		seedd = 0;
		seed32_1 = s;
		seed32_2 = 2;
		seed32_3 = 3;
		seed32_4 = 4;
		for (unsigned i = 0; i < 0x1000; ++i)
			xrng128();
	}

	void seed(uint32* sbuf, unsigned len)
	{
		seedd = 0;
		seed32_1 = 1;
		seed32_2 = 2;
		seed32_3 = 3;
		seed32_4 = 4;
		for (unsigned i = 0; i < len; ++i)
		{
			seed32_1 ^= sbuf[i];
			xrng128();
		}
		for (unsigned i = 0; i < 0x1000; ++i)
			xrng128();
	}

	void addseed(uint32 s)
	{
		xrng128();
		seed32_1 ^= s;
		xrng128();
	}

	void addseed(const uint32* sbuf, unsigned len)
	{
		xrng128();
		for (unsigned i = 0; i < len; ++i)
		{
			seed32_1 ^= sbuf[i];
			xrng128();
		}
	}

	struct hashclash_rng__init {
		hashclash_rng__init()
		{
			seed(uint32(time(NULL)));
			uint32 rndbuf[256]; // uninitialized on purpose
			addseed(rndbuf, 256);
			getosrnd(rndbuf);
			addseed(rndbuf, 256);
		}
	};
	hashclash_rng__init hashclash_rng__init__now;

	void hashclash_rng_hpp_init()
	{
		hashclash_rng__init here;
	}

} // namespace hashclash
