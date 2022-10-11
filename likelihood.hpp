//
//  likelihood.hpp
//  Phylogenetic_code
//
//  Created by TungDang on 2022/10/04.
//

#ifndef likelihood_h
#define likelihood_h

#pragma once

#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/range/adaptor/replaced.hpp>
#include "libhmsbeagle/beagle.h"
#include "tree.hpp"
#include "data.hpp"
#include "xstrom.hpp"

namespace  strom {
    
    class Likelihood {
    
        public:
            
            Likelihood();
            ~Likelihood();
            
            void setRooted(bool is_rooted);
            void setPreferGPU(bool prefer_gpu);
            void setAmbiguityEqualsMissing(bool ambig_equals_missing);
            
            bool usingStoredData() const;
            void useStoredData(bool using_data);
            
            std::string beagleLibVersion() const;
            std::string availableResources() const;
            std::string usedResources() const;
            
            void initBeagleLib();
            void finalizeBeagleLib(bool use_exceptions);
            
            double calcLogLikelihood(Tree::SharedPtr t);
            
            Data::SharedPtr getData();
            void setData(Data::SharedPtr d);
            
            void clear();
            
            unsigned calcNumEdgesInFullyResolvedTree() const;
            unsigned calcNumInternalsInFullyResolvedTree() const;
        
        private:
            
            struct InstanceInfo {
                
                int handle;
                int resourcenumbers;
                std::string resourcenames;
                unsigned nstates;
                unsigned nratecateg;
                unsigned npatterns;
                unsigned partial_offset;
                unsigned tmatrix_offset;
                std::vector<unsigned> subsets;
                
                InstanceInfo(): handle(-1), resourcenumbers(-1), resourcenames(""), nstates(0), nratecateg(0), npatterns(0), partial_offset(0), tmatrix_offset(0) {}
            };
            
            typedef std::pair<unsigned, int> instance_pair_t;
            
            unsigned getScalerIndex(Node * nd, InstanceInfo & info) const;
            unsigned getPartialIndex(Node * nd, InstanceInfo & info) const;
            unsigned getTMatrixIndex(Node * nd, InstanceInfo & info, unsigned subset_index) const;
            
            void updateInstanceMap(instance_pair_t & p, unsigned subset);
            void newInstance(unsigned nstates, int nrates, std::vector<unsigned> & subset_indices);
            void setTipStates();
            void setTipPartials();
            void setPatternPartitionAssignments();
            void setPatternWeights();
            void setAmongSiteRateHeterogenetity();
            void setModeRateMatrix();
            void addOperation(InstanceInfo & info, Node * nd, Node * lchild, Node * rchild, unsigned subset_index);
            void queuePartialsRecalculation(Node * nd, Node * lchild, Node * rchild);
            void queueTMatrixRecalculation(Node * nd);
            void defineOperations(Tree::SharedPtr t);
            void updateTransitionMatrices();
            void calculatePartials();
            double calcInstanceLogLikelihood(InstanceInfo & inst, Tree::SharedPtr t);
            
            std::vector<InstanceInfo> _instances;
            std::map<int, std::string> _beagle_error;
            std::map<int, std::vector<int>> _operations;
            std::map<int, std::vector<int>> _pmatrix_index;
            std::map<int, std::vector<double>> _edge_lengths;
            std::map<int, std::vector<int>> _eigen_indices;
            std::map<int, std::vector<int>> _category_rate_indices;
            double _relrate_normalizing_constant;
            
            std::vector<int> _subset_indices;
            std::vector<int> _parent_indices;
            std::vector<int> _child_indices;
            std::vector<int> _tmatrix_indices;
            std::vector<int> _weights_indices;
            std::vector<int> _freqs_indices;
            std::vector<int> _scaling_indices;
            
            Data::SharedPtr _data;
            unsigned _ntaxa;
            bool _rooted;
            bool _prefer_gpu;
            bool _ambiguity_equals_missing;
            bool _using_data;
            
        public:
            
            typedef std::shared_ptr<Likelihood> SharedPtr;
            
    };
    
    inline Likelihood::Likelihood() {
        clear();
    }
    
    inline Likelihood::~Likelihood() {
        finalizeBeagleLib(false);
        clear();
    }
    
    inline unsigned Likelihood::calcNumEdgesInFullyResolvedTree() const {
        assert(_ntaxa > 0);
        return (_rooted ? (2*_ntaxa - 2) : (2*_ntaxa - 3));
    }
    
    inline unsigned Likelihood::calcNumInternalsInFullyResolvedTree() const {
        assert(_ntaxa > 0);
        return (_rooted ? (_ntaxa - 1) : (_ntaxa - 2));
    }
    
    inline void Likelihood::finalizeBeagleLib(bool use_exceptions) {
        // Close down all BeagleLib instances if active
        
        for (auto info : _instances) {
            if (info.handle >= 0) {
                int code = beagleFinalizeInstance(info.handle);
                if (code != 0) {
                    if (use_exceptions)
                        throw XStrom(boost::format("Likelihood failed to finalize BeagleLib instance. BeagleLib error code was %d (%s).") % code % _beagle_error[code]);
                    else
                        std::cerr << boost::format("Likelihood destructor failed to finalize BeagleLib instance. BeagleLib error code was %d (%s).") % code % _beagle_error[code] << std::endl;
                }
            }
        }
        _instances.clear();
    }
    
    inline void Likelihood::clear() {
        finalizeBeagleLib(true);
        
        _ntaxa                      = 0;
        _rooted                     = false;
        _prefer_gpu                 = false;
        _ambiguity_equals_missing   = true;
        _using_data                 = true;
        _data                       = nullptr;
        
        _operations.clear();
        _pmatrix_index.clear();
        _edge_lengths.clear();
        _eigen_indices.clear();
        _category_rate_indices.clear();
        _relrate_normalizing_constant = 1.0;
        _subset_indices.assign(1, 0);
        _parent_indices.assign(1, 0);
        _child_indices.assign(1, 0);
        _tmatrix_indices.assign(1, 0);
        _weights_indices.assign(1, 0);
        _freqs_indices.assign(1, 0);
        _scaling_indices.assign(1, 0);
        
        // Store BeagleLib error codes so that useful
        // error messages may be provided to the user
        _beagle_error.clear();
        _beagle_error[0]  = std::string("success");
        _beagle_error[-1] = std::string("unspecified error");
        _beagle_error[-2] = std::string("not enough memory could be allocated");
        _beagle_error[-3] = std::string("unspecified exception");
        _beagle_error[-4] = std::string("the instance index is out of range, or the instance has not been created");
        _beagle_error[-5] = std::string("one of the indices specified exceeded the range of the array");
        _beagle_error[-6] = std::string("no resource matches requirements");
        _beagle_error[-7] = std::string("no implementation matches requirements");
        _beagle_error[-8] = std::string("floating-point range exceeded");
    }
    
    inline std::string Likelihood::beagleLibVersion() const {
        return std::string(beagleGetVersion());
    }
    
    inline std::string Likelihood::availableResources() const {
        BeagleResourceList * rsrcList = beagleGetResourceList();
        std::string s;
        for (int i = 0; i < rsrcList->length; ++i) {
            std::string desc = rsrcList->list[i].description;
            boost::trim(desc);
            if (desc.size() > 0)
                s += boost::str(boost::format("resource %d: %s (%s)\n") % i % rsrcList->list[i].name % desc);
            else
                s += boost::str(boost::format("resource %d: %s\n") % i % rsrcList->list[i].name);
        }
        boost::trim_right(s);
        return s;
    }
    
    inline std::string Likelihood::usedResources() const {
        std::string s;
        for (unsigned i = 0; i < _instances.size(); i++) {
            s  += boost::str(boost::format("instance %d: %s (resource %d)\n") % _instances[i].handle % _instances[i].resourcenames % _instances[i].resourcenumbers);
        }
        return s;
    }
    
    inline Data::SharedPtr Likelihood::getData() {
        return _data;
    }
    
    inline void Likelihood::setData(Data::SharedPtr data) {
        assert(_instances.size() == 0);
        assert(!data->getDataMatrix().empty());
        _data = data;
    }
    
    inline void Likelihood::setRooted(bool is_rooted) {
        assert(_instances.size() == 0 || _rooted == is_rooted);
        _rooted = is_rooted;
    }
    
    inline void Likelihood::setPreferGPU(bool prefer_gpu) {
        assert(_instances.size() == 0 || _prefer_gpu == prefer_gpu);
        _prefer_gpu = prefer_gpu;
    }
    
    inline void Likelihood::setAmbiguityEqualsMissing(bool ambig_equals_missing) {
        assert(_instances.size() == 0 || _ambiguity_equals_missing == ambig_equals_missing);
        _ambiguity_equals_missing = ambig_equals_missing;
    }
    
    inline void Likelihood::useStoredData(bool using_data) {
        _using_data = using_data;
    }
    
    inline bool Likelihood::usingStoredData() const {
        return _using_data;
    }
    
    inline void Likelihood::initBeagleLib() {
        assert(_data);
        
        finalizeBeagleLib(true);
        
        _ntaxa = _data->getNumTaxa();
        
        unsigned nsubsets = _data->getNumSubsets();
        std::set<instance_pair_t> nstates_ncateg_combinations;
        std::map<instance_pair_t, std::vector<unsigned>> subsets_for_pair;
        
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            // Create a pair comprising number of states and number of rate categories
            unsigned nstates = _data->getNumStatesForSubset(subset);
            int nrates = 1;
            instance_pair_t p = std::make_pair(nstates, nrates);
            
            nstates_ncateg_combinations.insert(p);
            subsets_for_pair[p].push_back(subset);
        }
        
        // Create one instance for each distinct nstates-nrates combination
        _instances.clear();
        for (auto p : nstates_ncateg_combinations) {
            newInstance(p.first, p.second, subsets_for_pair[p]);
            
            InstanceInfo & info = *_instances.rbegin();
            //std::cout << boost::str(boost::format("Created BeagleLib instance %d (%d state, % rate%s, %d subset%s)") % info.handle % info.nstates % info.nratecateg % (info.nratecateg == 1 ? "" : "s") % info.subsets.size() % (info.subsets.size() == 1 ? "" : "s")) << std::endl;
            std::cout << boost::str(boost::format("Created BeagleLib instance %d (%d states, %d rate%s, %d subset%s)") % info.handle % info.nstates % info.nratecateg % (info.nratecateg == 1 ? "" : "s") % info.subsets.size() % (info.subsets.size() == 1 ? "" : "s")) << std::endl;
        }
        
        if (_ambiguity_equals_missing)
            setTipStates();
        else
            setTipPartials();
        
        setPatternWeights();
        setPatternPartitionAssignments();
    }
    
    inline void Likelihood::newInstance(unsigned nstates, int nrates, std::vector<unsigned> &subset_indices) {
        
        unsigned num_subsets = (unsigned)subset_indices.size();
        
        unsigned ngammacat = (unsigned)nrates;
        
        unsigned num_patterns = 0;
        for (auto s : subset_indices) {
            num_patterns += _data->getNumPatternsInSubset(s);
        }
        
        unsigned num_internals = calcNumInternalsInFullyResolvedTree();
        unsigned num_edges = calcNumEdgesInFullyResolvedTree();
        unsigned num_nodes = num_edges + 1;
        unsigned num_transition_probs = num_nodes * num_subsets;
        
        long requirementFlags = 0;
        long preferenceFlags = BEAGLE_FLAG_PRECISION_SINGLE | BEAGLE_FLAG_THREADING_CPP;
        if (_prefer_gpu)
            preferenceFlags |= BEAGLE_FLAG_PROCESSOR_GPU;
        else
            preferenceFlags |= BEAGLE_FLAG_PROCESSOR_CPU;
        
        BeagleInstanceDetails instance_details;
        
        unsigned npartials = num_internals + _ntaxa;
        unsigned nsequences = 0;
        
        if (_ambiguity_equals_missing) {
            npartials -= _ntaxa;
            nsequences += _ntaxa;
        }
        
        int inst = beagleCreateInstance(_ntaxa, npartials, nsequences, nstates, num_patterns, num_subsets, num_subsets*num_transition_probs, ngammacat,
        0, NULL, 0, preferenceFlags, requirementFlags, &instance_details);
        
        if (inst < 0) {
            throw XStrom(boost::str(boost::format("Likelihood init function failed to create BeagleLib instance (BeagleLib error code was %d)") % _beagle_error[inst]));
        }
        
        InstanceInfo info;
        info.handle = inst;
        info.resourcenumbers = instance_details.resourceNumber;
        info.resourcenames = instance_details.resourceName;
        info.nstates = nstates;
        info.nratecateg = ngammacat;
        info.subsets = subset_indices;
        info.npatterns = num_patterns;
        info.partial_offset = num_internals;
        info.tmatrix_offset = num_nodes;
        _instances.push_back(info);
    }
    
    inline void Likelihood::setTipStates() {
        assert(_instances.size() > 0);
        assert(_data);
        Data::state_t one = 1;
        
        for (auto & info : _instances) {
            std::vector<int> states(info.nstates*info.npatterns);
            
            // Loop through all rows of the data matrix, setting the tip states for one taxon each row
            unsigned t = 0;
            for (auto & row : _data->getDataMatrix()) {
                
                // Loop through all subsets assigned to this instance
                unsigned k = 0;
                for (unsigned s : info.subsets) {
                    
                    // Loop through all patterns in this subset
                    auto interval = _data->getSubsetBeginEnd(s);
                    for (unsigned p = interval.first; p < interval.second; p++) {
                        
                        Data::state_t d = row[p];
                        
                        // Handle common nucleotide case separately
                        if (info.nstates == 4) {
                            if (d == 1)
                                states[k++] = 0;
                            else if (d == 2)
                                states[k++] = 1;
                            else if (d == 4)
                                states[k++] = 2;
                            else if (d == 8)
                                states[k++] = 3;
                            else
                                states[k++] = 4;
                        }
                        else {
                            int s = -1;
                            for (unsigned b = 0; b < info.nstates; b++) {
                                if (d == one << b) {
                                    s = b;
                                    break;
                                }
                            }
                            if (s == -1)
                                states[k++] = info.nstates;
                            else
                                states[k++] = s;
                        }
                    }
                }
            
            int code = beagleSetTipStates(info.handle, t, &states[0]);
            
            if (code != 0)
                throw XStrom(boost::format("failed to set tip state for taxon %d (\"%s\"; BeagleLib error code was %d)") % (t + 1) % _data->getTaxonNames()[t] % code % _beagle_error[code]);
            ++t;
            }
        }
    }
    
    inline void Likelihood::setTipPartials() {
        assert(_instances.size() > 0);
        assert(_data);
        
        Data::state_t one = 1;
        
        for (auto & info : _instances) {
            
            std::vector<double> partials(info.nstates*info.npatterns);
            
            // Loop through all rows of data matrix, setting the tip states for one taxon each row
            unsigned t = 0;
            for (auto & row : _data->getDataMatrix()) {
                
                // Loop through all subset assigned to this instance
                unsigned k = 0;
                for (unsigned s : info.subsets) {
                    
                    // Loop through all patterns in this subset
                    auto interval = _data->getSubsetBeginEnd(s);
                    for (unsigned p = interval.first; p < interval.second; p++) {
                        
                        Data::state_t d = row[p];
                        
                        if (info.nstates == 4) {
                            partials[k++] = d & 1 ? 1.0 : 0.0;
                            partials[k++] = d & 2 ? 1.0 : 0.0;
                            partials[k++] = d & 4 ? 1.0 : 0.0;
                            partials[k++] = d & 8 ? 1.0 : 0.0;
                        }
                        else {
                            for (unsigned b = 0; b < info.nstates; b++) {
                                partials[k++] = d & (one << b) ? 1.0 : 0.0;
                            }
                        }
                    }
                }
            int code = beagleSetTipPartials(info.handle, t, &partials[0]);
            
            if (code != 0)
                throw XStrom(boost::format("failed to set tip state for taxon %d (\"%s\"; BeagleLib error code was %d)") % (t + 1) % _data->getTaxonNames()[t] % code % _beagle_error[code]);
            ++t;
            }
        }
    }
    
    inline void Likelihood::setPatternPartitionAssignments() {
        
        assert(_instances.size() > 0);
        assert(_data);
        
        if (_instances.size() == 1 && _instances[0].subsets.size() == 1)
            return;
            
        Data::partition_key_t v;
        
        // Loop through all instances
        for (auto & info : _instances) {
            unsigned nsubsets = (unsigned)info.subsets.size();
            v.resize(info.npatterns);
            unsigned pattern_index = 0;
            
            // Loop through all subsets assigned to this instance
            unsigned instance_specific_subset_index = 0;
            for (unsigned s : info.subsets) {
                // Loop through all patterns in this subset
                auto interval = _data->getSubsetBeginEnd(s);
                for (unsigned p = interval.first; p < interval.second; p++) {
                    v[pattern_index++] = instance_specific_subset_index;
                }
                ++instance_specific_subset_index;
            }
            
            int code = beagleSetPatternPartitions(info.handle, nsubsets, &v[0]);
            
            if (code != 0)
                throw XStrom(boost::format("failed to set pattern partition. BeagleLib error code was %d (%s)") % code % _beagle_error[code]);
        }
    }
    
    inline void Likelihood::setPatternWeights() {
        assert(_instances.size() > 0);
        assert(_data);
        
        Data::pattern_counts_t v;
        auto pattern_counts = _data->getPatternCounts();
        assert(pattern_counts.size() > 0);
        
        // Loop through all instance
        for (auto & info : _instances) {
            v.resize(info.npatterns);
            unsigned pattern_index = 0;
            
            // Loop through all subsets assigned to this instance
            for (unsigned s : info.subsets) {
                
                // Loop through all patterns in this subset
                auto interval = _data->getSubsetBeginEnd(s);
                for (unsigned p = interval.first; p < interval.second; p++) {
                    v[pattern_index++] = pattern_counts[p];
                }
            }
            
            int code = beagleSetPatternWeights(info.handle, &v[0]);
            
            if (code != 0)
                throw XStrom(boost::format("failed to set pattern weights for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]);
            
        }
    }
    
    inline void Likelihood::setAmongSiteRateHeterogenetity() {
        assert(_instances.size() > 0);
        int code = 0;
        
        double rates[1] = {1.0};
        double probs[1] = {1.0};
        
        // Loop throgh all instance
        for (auto & info : _instances) {
            
            // Loop through all subsets assigned to this instance
            for (unsigned instance_specific_subset_index = 0; instance_specific_subset_index < info.subsets.size(); instance_specific_subset_index++) {
                code = beagleSetCategoryRatesWithIndex(info.handle, instance_specific_subset_index, rates);
                
                if (code != 0)
                    throw XStrom(boost::format("failed to set category rate for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]);
                    
                code = beagleSetCategoryWeights(info.handle, instance_specific_subset_index, probs);
                
                if (code != 0)
                    throw XStrom(boost::format("failed to set category probability for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]);
            }
        }
    }
    
    inline void Likelihood::setModeRateMatrix() {
        assert(_instances.size() > 0);
        int code = 0;
        double state_freqs[4] = {0.25, 0.25, 0.25, 0.25};
        
        double eigenvalues[4] = {
            -4.0/3.0,
            -4.0/3.0,
            -4.0/3.0,
            0.0
        };
        
        double eigenvectors[16] = {
            -1, -1, -1, 1,
            0, 0, 1, 1,
            0, 1, 0, 1,
            1, 0, 0, 1
        };
        
        double inverse_eigenvectors[16] = {
            -0.25, -0.25, -0.25, 0.75,
            -0.25, -0.25, 0.75, -0.25,
            -0.25, 0.75, -0.25, -0.25,
            0.25, 0.25, 0.25, 0.25
        };
        
        for (auto & info : _instances) {
            
            // Loop through all subsets assigned to this intanse
            for (unsigned instance_specific_subset_index = 0; instance_specific_subset_index < info.subsets.size(); instance_specific_subset_index++) {
                
                code = beagleSetStateFrequencies(info.handle, instance_specific_subset_index, state_freqs);
                
                if (code != 0)
                    throw XStrom(boost::format("failed to set state frequencies for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]);
                    
                code = beagleSetEigenDecomposition(info.handle, instance_specific_subset_index, (const double *)eigenvectors, (const double *)inverse_eigenvectors, eigenvalues);
                
                if (code != 0)
                    throw XStrom(boost::format("failed to set eigen decomposition for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]);
            }
        }
    }
    
    inline void Likelihood::defineOperations(Tree::SharedPtr t) {
        assert(_instances.size() > 0);
        assert(t);
        assert(t->isRooted() == _rooted);
        
        _relrate_normalizing_constant = 1.0;
        
        for (auto & info : _instances) {
            _operations[info.handle].clear();
            _pmatrix_index[info.handle].clear();
            _edge_lengths[info.handle].clear();
            _eigen_indices[info.handle].clear();
            _category_rate_indices[info.handle].clear();
        }
        
        // Loop through all nodes in reverse level order
        for (auto nd : boost::adaptors::reverse(t->_levelorder)) {
            assert(nd->_number >= 0);
            if (!nd->_left_child) {
                queueTMatrixRecalculation(nd);
            }
            else {
                queueTMatrixRecalculation(nd);
                
                // Internal node have partials to be calculated, so define
                // an operation to compute the partials for this node
                Node * lchild = nd->_left_child;
                assert(lchild);
                Node * rchild = lchild->_right_sib;
                assert(rchild);
                queuePartialsRecalculation(nd, lchild, rchild);
            }
        }
    }
    
    inline void Likelihood::queueTMatrixRecalculation(Node *nd) {
        double subset_relative_rate = 1.0;
        for (auto & info : _instances) {
            unsigned instance_specific_subset_index = 0;
            for (unsigned s : info.subsets) {
                unsigned tindex = getTMatrixIndex(nd, info, instance_specific_subset_index);
                _pmatrix_index[info.handle].push_back(tindex);
                _edge_lengths[info.handle].push_back(nd->_edge_length*subset_relative_rate);
                _eigen_indices[info.handle].push_back(s);
                _category_rate_indices[info.handle].push_back(s);
                
                ++instance_specific_subset_index;
            }
        }
    }
    
    inline void Likelihood::queuePartialsRecalculation(Node *nd, Node *lchild, Node *rchild) {
        for (auto & info : _instances) {
            unsigned instance_specific_subset_index = 0;
            for (unsigned s : info.subsets) {
                addOperation(info, nd, lchild, rchild, instance_specific_subset_index);
                ++instance_specific_subset_index;
            }
        }
    }
    
    inline void Likelihood::addOperation(InstanceInfo &info, Node *nd, Node *lchild, Node *rchild, unsigned subset_index) {
        
        assert(nd);
        assert(lchild);
        assert(rchild);
        
        // 1. destination partial to be calculated
        int partial_dest = getPartialIndex(nd, info);
        _operations[info.handle].push_back(partial_dest);
        
        // 2. destination scaling buffer index to write to
        int scaler_index = getScalerIndex(nd, info);
        _operations[info.handle].push_back(scaler_index);
        
        // 3. destination scaling buffer index to read from
        _operations[info.handle].push_back(BEAGLE_OP_NONE);
        
        // 4. left child partial index
        int partial_lchild = getPartialIndex(lchild, info);
        _operations[info.handle].push_back(partial_lchild);
        
        // 5. left child transition matrix index
        unsigned tindex_lchild = getTMatrixIndex(lchild, info, subset_index);
        _operations[info.handle].push_back(tindex_lchild);
        
        // 6. right child partial index
        int partial_rchild = getPartialIndex(rchild, info);
        _operations[info.handle].push_back(partial_rchild);
        
        // 7. right child transition matrix index
        unsigned tindex_rchild = getTMatrixIndex(rchild, info, subset_index);
        _operations[info.handle].push_back(tindex_rchild);
        
        if (info.subsets.size() > 1) {
            // 8. index of partition subset
            _operations[info.handle].push_back(subset_index);
            
            // 9. cumulative scale index
            _operations[info.handle].push_back(BEAGLE_OP_NONE);
        }
    }
    
    inline unsigned Likelihood::getPartialIndex(Node *nd, InstanceInfo &info) const {
        assert(nd->_number >= 0);
        return nd->_number;
    }
    
    inline unsigned Likelihood::getTMatrixIndex(Node *nd, InstanceInfo &info, unsigned subset_index) const {
        unsigned tindex = subset_index*info.tmatrix_offset + nd->_number;
        return tindex;
    }
    
    inline unsigned Likelihood::getScalerIndex(Node *nd, InstanceInfo &info) const {
        return BEAGLE_OP_NONE;
    }
    
    inline void Likelihood::updateTransitionMatrices() {
        
        assert(_instances.size() > 0);
        
        if (_pmatrix_index.size() == 0)
            return;
            
        for (auto & info : _instances) {
            int code = 0;
            
            unsigned nsubsets = (unsigned)info.subsets.size();
            if (nsubsets > 1) {
                code = beagleUpdateTransitionMatricesWithMultipleModels(info.handle, &_eigen_indices[info.handle][0], &_category_rate_indices[info.handle][0], &_pmatrix_index[info.handle][0],
                NULL, // first derivative matrix to update
                NULL, // second derivative matrix to update
                &_edge_lengths[info.handle][0], (int)_pmatrix_index[info.handle].size());
            } else {
                code = beagleUpdateTransitionMatrices(info.handle, 0, &_pmatrix_index[info.handle][0],
                NULL, // firest derivative matrices to update
                NULL, // second derivative matrices to update
                &_edge_lengths[info.handle][0], (int)_pmatrix_index[info.handle].size());
            }
            if (code != 0)
                    throw XStrom(boost::format("failed to update transition matrices for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]);
        }
    }
    
    inline void Likelihood::calculatePartials() {
        assert(_instances.size() > 0);
        
        if (_operations.size() == 0)
            return;
            
        int code = 0;
        
        for (auto & info : _instances) {
            unsigned nsubsets = (unsigned)info.subsets.size();
            
            if (nsubsets > 1) {
                code = beagleUpdatePartialsByPartition(info.handle,
                (BeagleOperationByPartition *) &_operations[info.handle][0],
                (int)(_operations[info.handle].size()/9));
                if (code != 0)
                    throw XStrom(boost::format("failed to update partials for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]);
            } else {
                code = beagleUpdatePartials(info.handle,
                (BeagleOperation *) &_operations[info.handle][0],
                (int)(_operations[info.handle].size()/7), BEAGLE_OP_NONE);
                if (code != 0)
                    throw XStrom(boost::format("failed to update partials for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]);
            }
        }
    }
    
    inline double Likelihood::calcInstanceLogLikelihood(InstanceInfo & info, Tree::SharedPtr t) {
        int code = 0;
        unsigned nsubsets = (unsigned)info.subsets.size();
        assert(nsubsets > 0);
        
        assert(_pmatrix_index[info.handle].size() == _edge_lengths[info.handle].size());
        
        int state_frequency_index = 0;
        int category_weights_index = 0;
        int cumulative_scale_index = BEAGLE_OP_NONE;
        int child_partials_index = getPartialIndex(t->_root, info);
        int parent_partials_index = getPartialIndex(t->_preorder[0], info);
        int parent_tmatrix_index = getTMatrixIndex(t->_preorder[0], info, 0);
        
        // storage for results of the likelihood calculation
        std::vector<double> subset_log_likelihoods(nsubsets, 0.0);
        double log_likelihood = 0;
        
        if (nsubsets > 1) {
            _parent_indices.assign(nsubsets, parent_tmatrix_index);
            _child_indices.assign(nsubsets, child_partials_index);
            _weights_indices.assign(nsubsets, category_weights_index);
            _scaling_indices.resize(nsubsets);
            _subset_indices.resize(nsubsets);
            _freqs_indices.resize(nsubsets);
            _tmatrix_indices.resize(nsubsets);
            
            for (unsigned s = 0; s < nsubsets; s++) {
                _scaling_indices[s] = BEAGLE_OP_NONE;
                _subset_indices[s] = s;
                _freqs_indices[s] = s;
                _tmatrix_indices[s] = getTMatrixIndex(t->_preorder[0], info, s);
            }
            
            code = beagleCalculateEdgeLogLikelihoodsByPartition(info.handle, &_parent_indices[0], &_child_indices[0], &_tmatrix_indices[0],
            NULL, NULL, &_weights_indices[0], &_freqs_indices[0], &_scaling_indices[0], &_subset_indices[0], nsubsets, 1,
            &subset_log_likelihoods[0], &log_likelihood, NULL, NULL, NULL, NULL);
        }
        else {
            code = beagleCalculateEdgeLogLikelihoods(info.handle, &parent_partials_index, &child_partials_index, &parent_tmatrix_index,
            NULL, NULL, &category_weights_index, &state_frequency_index, &cumulative_scale_index, 1, &log_likelihood, NULL, NULL);
        }
        
        return log_likelihood;
    }
    
    inline double Likelihood::calcLogLikelihood(Tree::SharedPtr t) {
        
        assert(_instances.size() > 0);
        
        if (!_using_data)
            return 0.0;
         
        assert(_data);
        
        if (t->_is_rooted)
            throw XStrom("This version of the program can only compute likelihood of unrooted tree");
            
        setModeRateMatrix();
        setAmongSiteRateHeterogenetity();
        defineOperations(t);
        updateTransitionMatrices();
        calculatePartials();
        
        double log_likelihood = 0.0;
        for (auto & info: _instances) {
            log_likelihood += calcInstanceLogLikelihood(info, t);
        }
        
        return log_likelihood;
    }
}

#endif /* likelihood_h */
