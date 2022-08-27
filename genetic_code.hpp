//
//  genetic_code.hpp
//  Phylogenetic_code
//
//  Created by TungDang on 2022/08/27.
//

#ifndef genetic_code_h
#define genetic_code_h

#pragma once
#include "xstrom.hpp"
#include <boost/algorithm/string.hpp>

namespace  strom {
    
    class Data;
    class Model;
    class QMatrix;
    
    class GeneticCode {
        
        friend class Data;
        friend class Model;
        friend class QMatrix;
        
        public:
        
            typedef std::map<int, int> genetic_code_map_t;
            typedef std::map<char, unsigned> amino_acid_map_t;
            typedef std::vector<unsigned> amino_acid_vect_t;
            typedef std::vector<std::string> codon_vect_t;
            typedef std::vector<char> amino_acid_symbol_vect_t;
            typedef std::map<std::string, std::string> genetic_code_definitions_t;
            typedef std::vector<std::string> genetic_code_names_t;
            
            GeneticCode();
            GeneticCode(std::string name);
            ~GeneticCode();
            
            std::string getGeneticCodeName() const;
            void useGeneticCode(std::string name);
            
            unsigned getNumNonStopCodons() const;
            int getStateCode(int triplet_index) const;
            char getAminoAcidAbbrev(unsigned aa_index) const;
            
            void copyCodons(std::vector<std::string> & codon_vect) const;
            void copyAminoAcids(std::vector<unsigned> & aa_vect) const;
            
            static genetic_code_names_t getRecognizedGeneticCodeNames();
            static bool isRecognizedGeneticCodeName(const std::string & name);
            static void ensureGeneticCodeNameIsValid(const std::string & name);
            
        private:
        
            void buildGeneticCodeTranslators();
            
            std::string _genetic_code_name;
            
            genetic_code_map_t _genetic_code_map;
            amino_acid_map_t _amino_acid_map;
            
            amino_acid_vect_t _amino_acids;
            codon_vect_t _codons;
            
            const amino_acid_symbol_vect_t  _all_amino_acids = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
            const std::vector<std::string> _all_codons = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};
            
            static genetic_code_definitions_t _definitions;
        
        public:
            
            typedef std::shared_ptr<GeneticCode> SharePtr;

    };
    
    inline GeneticCode::GeneticCode() {
        useGeneticCode("standard");
    }
    
    inline GeneticCode::GeneticCode(std::string name) {
        useGeneticCode(name);
    }
    
    inline GeneticCode::~GeneticCode(){
        std::cout << "Destroying a GeneticCode" << std::endl;
    }
    
    inline std::string GeneticCode::getGeneticCodeName() const {
        return _genetic_code_name;
    }
    
    inline void GeneticCode::useGeneticCode(std::string name) {
        _genetic_code_name = name;
        buildGeneticCodeTranslators();
    }
    
    inline unsigned GeneticCode::getNumNonStopCodons() const {
        return (unsigned)_codons.size();
    }
    
    inline int GeneticCode::getStateCode(int triplet_index) const {
        return _genetic_code_map.at(triplet_index);
    }
    
    inline char GeneticCode::getAminoAcidAbbrev(unsigned aa_index) const {
        return _all_amino_acids[aa_index];
    }
    
    inline void GeneticCode::copyCodons(std::vector<std::string> & codon_vect) const {
        codon_vect.resize(_codons.size());
        std::copy(_codons.begin(), _codons.end(), codon_vect.begin());
    }
    
    inline void GeneticCode::copyAminoAcids(std::vector<unsigned> & aa_vect) const {
        aa_vect.resize(_amino_acids.size());
        std::copy(_amino_acids.begin(), _amino_acids.end(), aa_vect.begin());
    }
    
    inline void GeneticCode::buildGeneticCodeTranslators() {
        _amino_acid_map.clear();
        for (unsigned i = 0; i < 20; ++i) {
            char aa = _all_amino_acids[i];
            _amino_acid_map[aa] = i;
        }
        
        ensureGeneticCodeNameIsValid(_genetic_code_name);
        std::string gcode_aa = _definitions[_genetic_code_name];
        
        int k = 0;
        int state_code = 0;
        _codons.clear();
        _amino_acids.clear();
        _genetic_code_map.clear();
        for (char ch : gcode_aa) {
            if (ch != '*') {
                _genetic_code_map[k] = state_code++;
                _codons.push_back(_all_codons[k]);
                _amino_acids.push_back(_amino_acid_map[ch]);
            }
            ++k;
        }
    }
    
}

#endif /* genetic_code_h */