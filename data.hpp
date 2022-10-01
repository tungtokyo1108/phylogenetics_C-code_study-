//
//  data.hpp
//  Phylogenetic_code
//
//  Created by TungDang on 2022/09/11.
//

#ifndef data_h
#define data_h

#pragma once

#include <fstream>
#include <regex>
#include <string>
#include <vector>
#include <numeric>
#include <limits>
#include <map>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>
#include "genetic_code.hpp"
#include "datatype.hpp"
#include "partition.hpp"
#include "xstrom.hpp"
#include "nxsmultiformat.h"

namespace strom {
    
    class Data {
      
      public:
        typedef std::vector<std::string> taxon_names_t;
        typedef unsigned long long state_t;
        typedef std::vector<state_t> pattern_vect_t;
        typedef std::vector<state_t> monomorphic_vect_t;
        typedef std::vector<int> partition_key_t;
        typedef std::map<pattern_vect_t,unsigned> pattern_map_t;
        typedef std::vector<pattern_vect_t> data_matrix_t;
        typedef std::vector<pattern_map_t> pattern_map_vect_t;
        typedef std::vector<double> pattern_counts_t;
        typedef std::vector<unsigned> subset_end_t;
        typedef std::vector<unsigned> npattern_vect_t;
        typedef std::pair<unsigned, unsigned> begin_end_pair_t;
        typedef std::shared_ptr<Data> SharedPtr;
        
        Data();
        ~Data();
        
        Partition::SharedPtr getPartition();
        void setPartition(Partition::SharedPtr partition);
        
        void getDataFromFile(const std::string filename);
        
        unsigned getNumSubsets() const;
        std::string getSubsetName(unsigned subset) const;
        
        unsigned getNumTaxa() const;
        const taxon_names_t & getTaxonNames() const;
        
        unsigned getNumPatterns() const;
        npattern_vect_t calcNumPatternsVect() const;
        unsigned getNumPatternsInSubset(unsigned subset) const;
        unsigned getNumStatesForSubset(unsigned subset) const;
        unsigned calcSeqLen() const;
        unsigned calcSeqLenInSubset(unsigned subset) const;
        const data_matrix_t & getDataMatrix() const;
        begin_end_pair_t getSubsetBeginEnd(unsigned subset) const;
        const pattern_counts_t & getPatternCounts() const;
        const monomorphic_vect_t & getMonomorphic() const;
        const partition_key_t & getPartitionKey() const;
        
        void clear();
        
      private:
        
        unsigned storeTaxonNames(NxsTaxaBlock * taxaBlock, unsigned taxa_block_index);
        unsigned storeData(unsigned ntax, unsigned nchar, NxsCharactersBlock * charBlock, NxsCharactersBlock::DataTypesEnum datatype);
        unsigned buildSubsetSpecificMaps(unsigned ntaxa, unsigned seqlen, unsigned nsubsets);
        void updatePatternMap(Data::pattern_vect_t & pattern, unsigned subset);
        void compressPatterns();
        
        Partition::SharedPtr _partition;
        pattern_counts_t _pattern_counts;
        monomorphic_vect_t _monomorphic;
        partition_key_t _partition_key;
        pattern_map_vect_t _pattern_map_vect;
        taxon_names_t _taxon_names;
        data_matrix_t _data_matrix;
        subset_end_t _subset_end;
    
    };
    
    inline Data::Data() {
        clear();
    }
    
    inline Data::~Data() {
    }
    
    inline void Data::setPartition(Partition::SharedPtr partition) {
        _partition = partition;
    }
    
    inline Partition::SharedPtr Data::getPartition() {
        return _partition;
    }
    
    inline unsigned Data::getNumSubsets() const {
        return (_partition ? _partition->getNumSubsets() : 1);
    }
    
    inline std::string Data::getSubsetName(unsigned subset) const {
        return _partition ? _partition->getSubsetName(subset) : std::string("default");
    }
    
    inline const Data::partition_key_t & Data::getPartitionKey() const {
        return _partition_key;
    }
    
    inline const Data::pattern_counts_t & Data::getPatternCounts() const {
        return _pattern_counts;
    }
    
    inline const Data::monomorphic_vect_t & Data::getMonomorphic() const {
        return _monomorphic;
    }

    inline const Data::taxon_names_t & Data::getTaxonNames() const {
        return _taxon_names;
    }

    inline const Data::data_matrix_t & Data::getDataMatrix() const {
        return _data_matrix;
    }
    
    inline Data::begin_end_pair_t Data::getSubsetBeginEnd(unsigned subset) const {
        assert(_subset_end.size() > subset);
        if (subset == 0)
            return std::make_pair(0, _subset_end[0]);
        else
            return std::make_pair(_subset_end[subset-1], _subset_end[subset]);
    }
    
    inline void Data::clear() {
        _partition_key.clear();
        _pattern_counts.clear();
        _monomorphic.clear();
        _pattern_map_vect.clear();
        _taxon_names.clear();
        _data_matrix.clear();
        _subset_end.clear();
    }
    
    inline unsigned Data::getNumPatterns() const {
        if (_data_matrix.size() > 0)
            return (unsigned)_data_matrix[0].size();
        else
            return 0;
    }
    
    inline Data::npattern_vect_t Data::calcNumPatternsVect() const {
        unsigned nsubsets = (unsigned)_subset_end.size();
        std::vector<unsigned> num_patters_vect(nsubsets, 0);
        for (unsigned s = 0; s < nsubsets; s++)
            num_patters_vect[s] = getNumPatternsInSubset(s);
        return num_patters_vect;
    }
    
    inline unsigned Data::getNumStatesForSubset(unsigned subset) const {
        DataType data_type = _partition->getDataTypeForSubset(subset);
        return data_type.getNumStates();
    }
    
    inline unsigned Data::getNumPatternsInSubset(unsigned subset) const {
        assert(_subset_end.size() > subset);
        return (unsigned)_subset_end[subset] - (subset == 0 ? 0 : _subset_end[subset-1]);
    }
    
    inline unsigned Data::getNumTaxa() const {
        return (unsigned)_taxon_names.size();
    }
    
    inline unsigned Data::calcSeqLen() const {
        return std::accumulate(_pattern_counts.begin(), _pattern_counts.end(), 0);
    }
    
    inline unsigned Data::calcSeqLenInSubset(unsigned subset) const {
        begin_end_pair_t s = getSubsetBeginEnd(subset);
        return std::accumulate(_pattern_counts.begin() + s.first, _pattern_counts.begin() + s.second , 0);
    }
    
    inline unsigned Data::buildSubsetSpecificMaps(unsigned ntaxa, unsigned seqlen, unsigned nsubsets) {
        pattern_vect_t pattern(ntaxa);
        
        _pattern_map_vect.clear();
        _pattern_map_vect.resize(nsubsets);
        
        const Partition::partition_t & tuples = _partition->getSubsetRangeVect();
        for (auto & t : tuples) {
            unsigned site_begin = std::get<0>(t);
            unsigned site_end = std::get<1>(t);
            unsigned site_skip = std::get<2>(t);
            unsigned site_subset = std::get<3>(t);
            for (unsigned site = site_begin; site <= site_end; site += site_skip) {
                for (unsigned taxon = 0; taxon < ntaxa; ++taxon) {
                    pattern[taxon] = _data_matrix[taxon][site-1];
                }
                updatePatternMap(pattern, site_subset);
            }
        }
        
        unsigned npatterns = 0;
        for (auto & map : _pattern_map_vect) {
            npatterns += (unsigned)map.size();
        }
        
        return npatterns;
    }
    
    inline void Data::updatePatternMap(Data::pattern_vect_t &pattern, unsigned subset) {
        pattern_map_t::iterator lowb = _pattern_map_vect[subset].lower_bound(pattern);
        if(lowb != _pattern_map_vect[subset].end() && !(_pattern_map_vect[subset].key_comp()(pattern, lowb->first))) {
            lowb->second += 1;
        }
        else {
            _pattern_map_vect[subset].insert(lowb, pattern_map_t::value_type(pattern, 1));
        }
    }
    
    inline void Data::compressPatterns() {
        unsigned ntaxa = (unsigned)_data_matrix.size();
        unsigned seqlen = (unsigned)_data_matrix[0].size();
        
        unsigned nsubsets = getNumSubsets();
        _subset_end.resize(nsubsets);
        _partition->finalize(seqlen);
        
        unsigned npatterns = buildSubsetSpecificMaps(ntaxa, seqlen, nsubsets);
        _pattern_counts.assign(npatterns, 0);
        _monomorphic.assign(npatterns, 0);
        _partition_key.assign(npatterns, -1);
        
        _data_matrix.resize(ntaxa);
        for (auto & row : _data_matrix) {
            row.resize(npatterns);
        }
        
        unsigned p = 0;
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            for (auto & pc : _pattern_map_vect[subset]) {
                _pattern_counts[p] = pc.second;
                _partition_key[p] = subset;
                
                state_t constant_state = pc.first[0];
                unsigned t = 0;
                for (auto sc : pc.first) {
                    constant_state &= sc;
                    _data_matrix[t][p] = sc;
                    ++t;
                }
                _monomorphic[p] = constant_state;
                ++p;
            }
            _subset_end[subset] = p;
            _pattern_map_vect[subset].clear();
        }
    }
    
    inline unsigned Data::storeTaxonNames(NxsTaxaBlock *taxaBlock, unsigned taxa_block_index) {
        unsigned ntax = 0;
        if (taxa_block_index == 0) {
            _taxon_names.clear();
            for (auto s : taxaBlock->GetAllLabels())
                _taxon_names.push_back(s);
            ntax = (unsigned)_taxon_names.size();
            _data_matrix.resize(ntax);
        }
        return ntax;
    }
    
    inline unsigned Data::storeData(unsigned ntax, unsigned nchar_before, NxsCharactersBlock *charBlock, NxsCharactersBlock::DataTypesEnum datatype) {
        unsigned seqlen = 0;
        
        unsigned subset_index = _partition->findSubsetForSite(nchar_before + 1);
        DataType dt = _partition->getDataTypeForSubset(subset_index);
        
        NxsCharactersBlock * block = charBlock;
        if (datatype == NxsCharactersBlock::dna || datatype == NxsCharactersBlock::rna || datatype == NxsCharactersBlock::nucleotide) {
            if (dt.isCodon()) {
                block = NxsCharactersBlock::NewCodonsCharactersBlock(charBlock, true, true, true, NULL, NULL);
            }
        }
        else if (datatype == NxsCharactersBlock::standard) {
            std::string symbols = std::string(charBlock->GetSymbols());
            dt.setStandardNumStates((unsigned)symbols.size());
        }
        else {
            return nchar_before;
        }
        
        unsigned num_states = dt.getNumStates();
        unsigned bits_in_state_t = 8*sizeof(state_t);
        
        for (unsigned t = 0; t < ntax; ++t) {
            const NxsDiscreteStateRow & row = block->GetDiscreteMatrixRow(t);
            if (seqlen == 0)
                seqlen = (unsigned)row.size();
            _data_matrix[t].resize(nchar_before + seqlen);
            
            unsigned k = nchar_before;
            for (int raw_state_code : row) {
                state_t state = std::numeric_limits<state_t>::max();
                bool complete_ambiguity = (!dt.isCodon() && raw_state_code == (int)num_states);
                bool all_missing_or_gaps = (raw_state_code < 0);
                if ((!complete_ambiguity) && (!all_missing_or_gaps)) {
                    int state_code = raw_state_code;
                    if (dt.isCodon())
                        state_code = dt.getGeneticCode()->getStateCode(raw_state_code);
                    if (state_code < (int)num_states) {
                        state = (state_t)1 << state_code;
                    }
                    else {
                        const NxsDiscreteDatatypeMapper * mapper = block->GetDatatypeMapperForChar(k - nchar_before);
                        const std::set<NxsDiscreteStateCell> & state_set = mapper->GetStateSetForCode(raw_state_code);
                        state = 0;
                        for (auto s : state_set) {
                            state |= (state_t)1 << s;
                        }
                    }
                }
                _data_matrix[t][k++] = state;
            }
        }
        return seqlen;
    }
    
    inline void Data::getDataFromFile(const std::string filename) {
        MultiFormatReader nexusReader(-1, NxsReader::WARNINGS_TO_STDERR);
        try {
            nexusReader.ReadFilepath(filename.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
        catch(...){
            nexusReader.DeleteBlocksFromFactories();
            throw;
        }
        
        clear();
        
        assert(_partition);
        
        int numTaxaBlocks = nexusReader.GetNumTaxaBlocks();
        if (numTaxaBlocks == 0)
            throw XStrom("No taxa blocks were found in the data file");
            
        unsigned cum_nchar = 0;
        for (int i = 0; i < numTaxaBlocks; ++i) {
            NxsTaxaBlock * taxaBlock = nexusReader.GetTaxaBlock(i);
            unsigned ntax = storeTaxonNames(taxaBlock, i);
            const unsigned numCharBlocks = nexusReader.GetNumCharactersBlocks(taxaBlock);
            for (unsigned j = 0; j < numCharBlocks; ++j) {
                NxsCharactersBlock * charBlock = nexusReader.GetCharactersBlock(taxaBlock, j);
                NxsCharactersBlock::DataTypesEnum datatype = charBlock->GetOriginalDataType();
                cum_nchar += storeData(ntax, cum_nchar, charBlock, datatype);
            }
        }
        
        nexusReader.DeleteBlocksFromFactories();
        
        if (_data_matrix.empty()) {
            std::cout << "No data were stored from the file \"" << filename << "\"" << std::endl;
            clear();
        }
        else {
            compressPatterns();
        }
        
    }

}

#endif /* data_h */
