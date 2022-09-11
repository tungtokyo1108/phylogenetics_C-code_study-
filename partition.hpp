//
//  partition.hpp
//  Phylogenetic_code
//
//  Created by TungDang on 2022/09/02.
//

#ifndef partition_h
#define partition_h

#pragma  once

#include <tuple>
#include <limits>
#include <cmath>
#include <iostream>
#include <regex>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "genetic_code.hpp"
#include "datatype.hpp"
#include "xstrom.hpp"

namespace strom {
    
    class Partition {
        
        public:
            typedef std::match_results<std::string::const_iterator>::const_reference    regex_match_t;
            typedef std::tuple<unsigned, unsigned, unsigned, unsigned> subset_range_t;
            typedef std::vector<subset_range_t> partition_t;
            typedef std::vector<DataType> datatype_vect_t;
            typedef std::vector<unsigned> subset_sizes_vect_t;
            typedef std::vector<std::string> subset_name_vect_t;
            typedef std::shared_ptr<Partition> SharedPtr;
            
            Partition();
            ~Partition();
            
            unsigned getNumSites() const;
            unsigned getNumSubsets() const;
            std::string getSubsetName(unsigned subset) const;
            
            const partition_t & getSubsetRangeVect() const;
            
            unsigned findSubsetByName(const std::string & subset_name) const;
            unsigned findSubsetForSite(unsigned site_index) const;
            bool siteInSubset(unsigned site_index, unsigned subset_index) const;
            DataType getDataTypeForSubset(unsigned subset_index) const;
            const datatype_vect_t & getSubsetDataType() const;
            
            unsigned numSitesInSubset(unsigned subset_index) const;
            subset_sizes_vect_t calcSubsetSize() const;
            
            void defaultPartition(unsigned nsites = std::numeric_limits<unsigned>::max());
            void parseSubsetDefinition(std::string & s);
            void finalize(unsigned nsites);
            
            void clear();
            
        private:
            int extractIntFromRegexMatch(regex_match_t s, unsigned min_value);
            void addSubsetRange(unsigned subset_index, std::string range_definition);
            void addSubset(unsigned subset_index, std::string subset_definition);
            
            unsigned _num_sites;
            unsigned _num_subsets;
            subset_name_vect_t _subset_names;
            partition_t _subset_ranges;
            datatype_vect_t _subset_data_types;
            const unsigned _infinity;
    };
    
    inline Partition::Partition() : _infinity(std::numeric_limits<unsigned>::max()) {
        //std::cout << "Constructing a Partition" << std::endl;
        clear();
    }

    inline Partition::~Partition() {
        //std::cout << "Destroying a Partition" << std::endl;
    }
    
    inline unsigned Partition::getNumSites() const {
        return _num_sites;
    }
    
    inline unsigned Partition::getNumSubsets() const {
        return _num_subsets;
    }
    
    inline std::string Partition::getSubsetName(unsigned subset) const {
        assert(subset < _num_subsets);
        return _subset_names[subset];
    }
    
    inline const Partition::partition_t & Partition::getSubsetRangeVect() const {
        return _subset_ranges;
    }
    
    inline DataType Partition::getDataTypeForSubset(unsigned subset_index) const {
        assert(subset_index < _subset_data_types.size());
        return _subset_data_types[subset_index];
    }

    inline const std::vector<DataType> & Partition::getSubsetDataType() const {
        return _subset_data_types;
    }
    
    inline unsigned Partition::findSubsetByName(const std::string & subset_name) const {
        auto iter = std::find(_subset_names.begin(), _subset_names.end(), subset_name);
        return (unsigned)std::distance(_subset_names.begin(), iter);
    }
    
    inline unsigned Partition::findSubsetForSite(unsigned site_index) const {
        for (auto & t : _subset_ranges) {
            unsigned begin_site = std::get<0>(t);
            unsigned end_site = std::get<1>(t);
            unsigned stride = std::get<2>(t);
            unsigned site_subset = std::get<3>(t);
            bool inside_range = site_index >= begin_site && site_index <= end_site;
            if (inside_range && (site_index - begin_site) % stride == 0)
                return site_subset;
        }
        throw XStrom(boost::format("Site %d not found in any subset of partition") % (site_index + 1));
    }
    
    inline bool Partition::siteInSubset(unsigned site_index, unsigned subset_index) const {
        unsigned which_subset = findSubsetForSite(subset_index);
        return (which_subset == subset_index ? true : false);
    }
    
    inline unsigned Partition::numSitesInSubset(unsigned subset_index) const {
        unsigned nsites = 0;
        for (auto & t : _subset_ranges) {
            unsigned begin_site = std::get<0>(t);
            unsigned end_site = std::get<1>(t);
            unsigned stride = std::get<2>(t);
            unsigned site_subset = std::get<3>(t);
            if (site_subset == subset_index) {
                unsigned n = end_site - begin_site + 1;
                nsites += (unsigned)(floor(n/stride)) + (n % stride == 0 ? 0 : 1);
            }
        }
        return nsites;
    }
    
    inline std::vector<unsigned> Partition::calcSubsetSize() const {
        std::vector<unsigned> nsites_vect(_num_subsets, 0);
        for (auto & t : _subset_ranges) {
            unsigned begin_site = std::get<0>(t);
            unsigned end_site = std::get<1>(t);
            unsigned stride = std::get<2>(t);
            unsigned site_subset = std::get<3>(t);
            unsigned hull = end_site - begin_site + 1;
            unsigned n = (unsigned)(floor(hull/stride)) + (hull % stride == 0 ? 0 : 1);
            nsites_vect[site_subset] += n;
        }
        return nsites_vect;
    }
    
    inline void Partition::clear() {
        _num_sites = 0;
        _num_subsets = 1;
        _subset_data_types.clear();
        _subset_data_types.push_back(DataType());
        _subset_names.clear();
        _subset_names.push_back("default");
        _subset_ranges.clear();
        _subset_ranges.push_back(std::make_tuple(1, _infinity, 1, 0));
    }
    
    inline void Partition::parseSubsetDefinition(std::string & s) {
        std::vector<std::string> v;
        
        // first separate part before colon from part after colon
        boost::split(v, s, boost::is_any_of(":"));
        
        std::string before_colon = v[0];
        std::string subset_definition = v[1];
        
        // now see if before_colon contains a data type specification in square brackets
        const char * pattern_string = R"((.+?)\s*(\[(\S+?)\])*)";
        std::regex re(pattern_string);
        std::smatch match_obj;
        bool matched = std::regex_match(before_colon, match_obj, re);
        
        if (!matched) {
            throw XStrom(boost::format("Could not interpret \"%s\" as a subset label with optional data type in square brackets") % before_colon);
        }
        
        std::string subset_name = match_obj[1].str();
        DataType dt;    // nucleotide by default
        std::string datatype = "nucleotide";
        if (match_obj.size() == 4 && match_obj[3].length() > 0) {
            datatype = match_obj[3].str();
            boost::to_lower(datatype);

            // check for comma plus genetic code in case of codon
            std::regex re(R"(codon\s*,\s*(\S+))");
            std::smatch m;
            if (std::regex_match(datatype, m, re)) {
                dt.setCodon();
                std::string genetic_code_name = m[1].str();
                dt.setGeneticCodeFromName(genetic_code_name);
            }
            else if (datatype == "codon") {
                dt.setCodon();  // assumes standard genetic code
            }
            else if (datatype == "protein") {
                dt.setProtein();
                }
            else if (datatype == "nucleotide") {
                dt.setNucleotide();
                }
            else if (datatype == "standard") {
                dt.setStandard();
                }
            else {
                throw XStrom(boost::format("Datatype \"%s\" specified for subset(s) \"%s\" is invalid: must be either nucleotide, codon, protein, or standard") % datatype % subset_name);
                }
            }
        unsigned end_site = std::get<1>(_subset_ranges[0]);
        if (_num_subsets == 1 && end_site == _infinity) {
            _subset_names.clear();
            _subset_data_types.clear();
            _subset_ranges.clear();
        }
        else if (subset_name == "default") {
            throw XStrom("Cannot specify \"default\" partition subset after already defining other subsets");
        }
        
        _subset_names.push_back(subset_name);
        _subset_data_types.push_back(dt);
        _num_subsets = (unsigned)_subset_names.size();
        addSubset(_num_subsets - 1, subset_definition);
        
        std::cout << boost::str(boost::format("Partition subset %s comprises sites %s and has type %s") % subset_name % subset_definition % datatype) << std::endl;
    }
    
    inline void Partition::addSubset(unsigned subset_index, std::string subset_definition) {
        std::vector<std::string> parts;
        boost::split(parts, subset_definition, boost::is_any_of(","));
        for (auto subset_component : parts) {
            addSubsetRange(subset_index, subset_component);
        }
    }
    
    inline void Partition::addSubsetRange(unsigned subset_index, std::string range_definition) {
        const char * pattern_string = R"((\d+)\s*(-\s*([0-9.]+)(\\\s*(\d+))*)*)";
        std::regex re(pattern_string);
        std::smatch match_obj;
        bool matched = std::regex_match(range_definition, match_obj, re);
        if (!matched) {
            throw XStrom(boost::format("Could not interpret \"%s\" as a range of site indices") % range_definition);
        }
        
        int ibegin = extractIntFromRegexMatch(match_obj[1], 1);
        int iend   = extractIntFromRegexMatch(match_obj[3], ibegin);
        int istep  = extractIntFromRegexMatch(match_obj[5], 1);
        
        // record the triplet
        _subset_ranges.push_back(std::make_tuple(ibegin, iend, istep, subset_index));
        
        // determine last site in subset
        unsigned last_site_in_subset = iend - ((iend - ibegin) % istep);
        if (last_site_in_subset > _num_sites) {
            _num_sites = last_site_in_subset;
        }
    }
    
    inline int Partition::extractIntFromRegexMatch(regex_match_t s, unsigned min_value) {
        int int_value = min_value;
        if (s.length() > 0) {
            std::string str_value = s.str();
            try {
                int_value = std::stoi(str_value);
            }
            catch(std::invalid_argument) {
                throw XStrom(boost::format("Could not interpret \"%s\" as a number in partition subset definition") % s.str());
            }
            
            // sanity check
            if (int_value < (int)min_value) {
                throw XStrom(boost::format("Value specified in partition subset definition (%d) is lower than minimum value (%d)") % int_value % min_value);
            }
        }
        return int_value;
    }
    
    inline void Partition::finalize(unsigned nsites) {
        if (_num_sites == 0) {
            defaultPartition(nsites);
            return;
        }

        // First sanity check:
        //   nsites is the number of sites read in from a data file;
        //   _num_sites is the maximum site index specified in any partition subset.
        //   These two numbers should be the same.
        if (_num_sites != nsites) {
            throw XStrom(boost::format("Number of sites specified by the partition (%d) does not match actual number of sites (%d)") % _num_sites % nsites);
        }
        
        // Second sanity check: ensure that no sites were left out of all partition subsets
        // Third sanity check: ensure that no sites were included in more than one partition subset
        std::vector<int> tmp(nsites, -1);   // begin with -1 for all sites
        for (auto & t : _subset_ranges) {
            unsigned begin_site  = std::get<0>(t);
            unsigned end_site    = std::get<1>(t);
            unsigned stride  = std::get<2>(t);
            unsigned site_subset = std::get<3>(t);
            for (unsigned s = begin_site; s <= end_site; s += stride) {
                if (tmp[s-1] != -1)
                    throw XStrom("Some sites were included in more than one partition subset");
                else
                    tmp[s-1] = site_subset;
            }
        }
        if (std::find(tmp.begin(), tmp.end(), -1) != tmp.end()) {
            throw XStrom("Some sites were not included in any partition subset");
        }
        tmp.clear();
    }
    
    inline void Partition::defaultPartition(unsigned nsites) {
        clear();
        _num_sites = nsites;
        _num_subsets = 1;
        _subset_ranges[0] = std::make_tuple(1, nsites, 1, 0);
    }
    
    
}

#endif /* partition_h */
