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

#ifndef HASHCLASH_SAVELOAD_GZ_HPP
#define HASHCLASH_SAVELOAD_GZ_HPP

#include "saveload.hpp"

#pragma warning(disable: 4996)
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace hashclash {

    template<class value_type>
    void save_gz(const value_type& val, const std::string& filepath, archive_type artype)
    {
		switch (artype) {
			case binary_archive:
				save_gz(val, artype, path(filepath + ".bin.gz"));
				break;
			case xml_archive:
				save_gz(val, artype, path(filepath + ".xml.gz"));
				break;
			case text_archive:
				save_gz(val, artype, path(filepath + ".txt.gz"));
				break;
			default:
				throw std::runtime_error("save_gz(): archive type not supported!");
		}
    }

    template<class value_type>
    void save_gz(const value_type& val, archive_type artype, const path& filepath)
    {
		std::ofstream ofs(filepath.string().c_str(), std::ios::binary);
        if (!ofs) throw std::runtime_error("save_gz(): could not open file!");
        boost::iostreams::filtering_ostream filter;
        filter.push(boost::iostreams::gzip_compressor());
        filter.push(ofs);
        switch (artype) {
			case binary_archive:
				{
					boost::archive::binary_oarchive oa(filter);
					oa << boost::serialization::make_nvp("saveload", val);
				}
				break;

			case xml_archive:
				{
					boost::archive::xml_oarchive oa(filter);
					oa << boost::serialization::make_nvp("saveload", val);
				}
				break;

			case text_archive:
				{
					boost::archive::text_oarchive oa(filter);
					oa << boost::serialization::make_nvp("saveload", val);
				}
				break;

			default:
				throw std::runtime_error("save_gz(): archive type not supported!");
        }
        if (!ofs) throw std::runtime_error("save_gz(): write error!");
    }

    template<class value_type>
    void load_gz(value_type& val, const std::string& filepath, archive_type artype)
    {
		switch (artype) {
			case binary_archive:
				load_gz(val, artype, path(filepath + ".bin.gz"));
				break;
			case xml_archive:
				load_gz(val, artype, path(filepath + ".xml.gz"));
				break;
			case text_archive:
				load_gz(val, artype, path(filepath + ".txt.gz"));
				break;
			default:
				throw std::runtime_error("load_gz(): archive type not supported!");
		}
    }

    template<class value_type>
    void load_gz(value_type& val, archive_type artype, const path& filepath)
    {
		std::ifstream ifs(filepath.string().c_str(), std::ios::binary);
        if (!ifs) throw std::runtime_error("load_gz(): could not open file!");
        boost::iostreams::filtering_istream filter;
        filter.push(boost::iostreams::gzip_decompressor());
        filter.push(ifs);
        switch (artype) {
			case binary_archive:
				{
					boost::archive::binary_iarchive ia(filter);
					ia >> boost::serialization::make_nvp("saveload", val);
				}
				break;

			case xml_archive:
				{
					boost::archive::xml_iarchive ia(filter);
					ia >> boost::serialization::make_nvp("saveload", val);
				}
				break;

			case text_archive:
				{
					boost::archive::text_iarchive ia(filter);
					ia >> boost::serialization::make_nvp("saveload", val);
				}
				break;

			default:
				throw std::runtime_error("load_gz(): archive type not supported!");
        }
        if (!ifs) throw std::runtime_error("load_gz(): read error!");
    }


	/* specialization types to distinguish between differentialpath & vector<differentialpath> */
	template<>
	inline void save_gz(const differentialpath& val, archive_type artype, const path& filepath)
	{
		const saveload_diffpath tmp(val);
		save_gz(tmp, artype, filepath);
	}
	template<>
	inline void load_gz(differentialpath& val, archive_type artype, const path& filepath)
	{
		saveload_diffpath tmp(val);
		load_gz(tmp, artype, filepath);
	}
	template<>
	inline void save_gz(const std::vector<differentialpath>& val, archive_type artype, const path& filepath)
	{
		const saveload_vecdiffpath tmp(val);
		save_gz(tmp, artype, filepath);
	}
	template<>
	inline void load_gz(std::vector<differentialpath>& val, archive_type artype, const path& filepath)
	{
		saveload_vecdiffpath tmp(val);
		load_gz(tmp, artype, filepath);
	}


} // namespace hashclash

#endif // HASHCLASH_SAVELOAD_GZ_HPP
