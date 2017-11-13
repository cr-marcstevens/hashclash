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

#ifndef HASHCLASH_SAVELOAD_HPP
#define HASHCLASH_SAVELOAD_HPP

#include <vector>
#include <string>
#include <stdexcept>

#include "types.hpp"
#include "differentialpath.hpp"

#ifdef NOSERIALIZATION
#error "Required boost::serialization is disabled. Undefine NOSERIALIZATION."
#endif // NOSERIALIZATION

#pragma warning(push)
#pragma warning(disable: 4267)
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/deque.hpp>
#include <boost/serialization/hash_map.hpp>
#include <boost/serialization/hash_set.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/scoped_ptr.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/slist.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#pragma warning(pop)

namespace hashclash {

	namespace fs = boost::filesystem;
	using fs::path;

	enum archive_type
	{
		binary_archive = 0, 
		xml_archive, 
		text_archive
	};

	template<class value_type>
	void save(const value_type& val, const std::string& filepath, archive_type artype)
	{
		switch (artype) {
			case binary_archive:
				save(val, artype, path(filepath + ".bin"));
				break;
			case xml_archive:
				save(val, artype, path(filepath + ".xml"));
				break;
			case text_archive:
				save(val, artype, path(filepath + ".txt"));
				break;
			default:
				throw std::runtime_error("save(): archive type not supported!");
		}
	}

	template<class value_type>
	void save(const value_type& val, archive_type artype, const path& filepath)
	{
		switch (artype) {
			case binary_archive:
				{ 
					std::ofstream ofs(filepath.string().c_str(), std::ios::binary);
					if (!ofs) throw std::runtime_error("save(): could not open file!");
					boost::archive::binary_oarchive oa(ofs);
					oa << boost::serialization::make_nvp("saveload", val);
					if (!ofs) throw std::runtime_error("save(): write error!");
				}
				break;

			case xml_archive:
				{
					std::ofstream ofs(filepath.string().c_str());
					if (!ofs) throw std::runtime_error("save(): could not open file!");
					boost::archive::xml_oarchive oa(ofs);
					oa << boost::serialization::make_nvp("saveload", val);
					if (!ofs) throw std::runtime_error("save(): write error!");
				}
				break;

			case text_archive:
				{
					std::ofstream ofs(filepath.string().c_str());
					if (!ofs) throw std::runtime_error("save(): could not open file!");
					boost::archive::text_oarchive oa(ofs);
					oa << boost::serialization::make_nvp("saveload", val);
					if (!ofs) throw std::runtime_error("save(): write error!");
				}
				break;

			default:
				throw std::runtime_error("save(): archive type not supported!");
		}		
	}

	template<class value_type>
	void load(value_type& val, const std::string& filepath, archive_type artype)
	{
		switch (artype) {
			case binary_archive:
				load(val, artype, path(filepath + ".bin"));
				break;
			case xml_archive:
				load(val, artype, path(filepath + ".xml"));
				break;
			case text_archive:
				load(val, artype, path(filepath + ".txt"));
				break;
			default:
				throw std::runtime_error("load(): archive type not supported!");
		}
	}

	template<class value_type>
	void load(value_type& val, archive_type artype, const path& filepath)
	{
		switch (artype)	{
			case binary_archive:
				{
					std::ifstream ifs(filepath.string().c_str(), std::ios::binary);
					if (!ifs) throw std::runtime_error("load(): could not open file!");
					boost::archive::binary_iarchive ia(ifs);
					ia >> boost::serialization::make_nvp("saveload", val);
					if (!ifs) throw std::runtime_error("load(): read error!");	
				}
				break;

			case xml_archive:
				{
					std::ifstream ifs(filepath.string().c_str());
					if (!ifs) throw std::runtime_error("load(): could not open file!");
					boost::archive::xml_iarchive ia(ifs);
					ia >> boost::serialization::make_nvp("saveload", val);
					if (!ifs) throw std::runtime_error("load(): read error!");	
				}
				break;

			case text_archive:
				{
					std::ifstream ifs(filepath.string().c_str());
					if (!ifs) throw std::runtime_error("load(): could not open file!");
					boost::archive::text_iarchive ia(ifs);
					ia >> boost::serialization::make_nvp("saveload", val);
					if (!ifs) throw std::runtime_error("load(): read error!");	
					break;
				}

			default:
				throw std::runtime_error("load(): archive type not supported!");
		}
	}


	/* specialization types to distinguish between differentialpath & vector<differentialpath> */
	struct saveload_diffpath {
		differentialpath& p;
		uint32 magic;
		saveload_diffpath(const differentialpath& path): p(const_cast<differentialpath&>(path)), magic(0x23580001) {}
		template<class Archive>
		void serialize(Archive& ar, const unsigned int file_version)
		{			
			ar & boost::serialization::make_nvp("magic", magic);
			if (magic != 0x23580001) throw std::runtime_error("serialize(): type is not differentialpath");
			ar & boost::serialization::make_nvp("path", p);
		}
	};
	template<>
	inline void save(const differentialpath& val, archive_type artype, const path& filepath)
	{
		const saveload_diffpath tmp(val);
		save(tmp, artype, filepath);
	}
	template<>
	inline void load(differentialpath& val, archive_type artype, const path& filepath)
	{
		saveload_diffpath tmp(val);
		load(tmp, artype, filepath);
	}

	struct saveload_vecdiffpath {
		std::vector<differentialpath>& p;
		uint32 magic;
		saveload_vecdiffpath(const std::vector<differentialpath>& path): p(const_cast<std::vector<differentialpath>&>(path)), magic(0x23580002) {}
		template<class Archive>
		void serialize(Archive& ar, const unsigned int file_version)
		{			
			ar & boost::serialization::make_nvp("magic", magic);
			if (magic != 0x23580002) throw std::runtime_error("serialize(): type is not vector<differentialpath>");
			ar & boost::serialization::make_nvp("path", p);
		}
	};
	template<>
	inline void save(const std::vector<differentialpath>& val, archive_type artype, const path& filepath)
	{
		const saveload_vecdiffpath tmp(val);
		save(tmp, artype, filepath);
	}
	template<>
	inline void load(std::vector<differentialpath>& val, archive_type artype, const path& filepath)
	{
		saveload_vecdiffpath tmp(val);
		load(tmp, artype, filepath);
	}

} // namespace 

#endif // HASHCLASH_SAVELOAD_HPP
