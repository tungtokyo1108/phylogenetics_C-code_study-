//
//  main.cpp
//  Phylogenetic_code
//
//  Created by TungDang on 2022/08/20.
//

#include <iostream>
//#include "node.hpp"
//#include "tree.hpp"
//#include "TreeManip.hpp"
#include "genetic_code.hpp"
#include "datatype.hpp"
#include "partition.hpp"
#include "data.hpp"
#include "tree_summary.hpp"
#include "likelihood.hpp"
//#include "strom.hpp"
//#include "xstrom.hpp"

using namespace strom;

const double Node::_smallest_edge_length = 1.0e-12;

//std::string  Strom::_program_name        = "strom";
//unsigned     Strom::_major_version       = 1;
//unsigned     Strom::_minor_version       = 0;
GeneticCode::genetic_code_definitions_t GeneticCode::_definitions = {
                             // codon order is alphabetical: i.e. AAA, AAC, AAG, AAT, ACA, ..., TTT
    {"standard",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"vertmito",             "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"yeastmito",            "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"moldmito",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"invertmito",           "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"ciliate",              "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF"},
    {"echinomito",           "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"euplotid",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF"},
    {"plantplastid",         "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"altyeast",             "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"},
    {"ascidianmito",         "KNKNTTTTGSGSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"altflatwormmito",      "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVYY*YSSSSWCWCLFLF"},
    {"blepharismamacro",     "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YQYSSSS*CWCLFLF"},
    {"chlorophyceanmito",    "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLYSSSS*CWCLFLF"},
    {"trematodemito",        "NNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"},
    {"scenedesmusmito",      "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*YLY*SSS*CWCLFLF"},
    {"thraustochytriummito", "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWC*FLF"}
};

int main(int argc, const char * argv[]) {

    std::cout << "Starting ... " << std::endl;
    
    /*
    
    // Test tree structure
    
    TreeSummary sumt;
    try {
        sumt.readTreefile("test.tre", 1);
    } catch(NxsException x) {
        std::cerr << "Program aborting due to error encounted reading tree file" << std::endl;
        std::cerr << x.what() << std::endl;
    }
    
    sumt.showSummary();
    */
    
    
    
    // Test input dataset
    
    Partition::SharedPtr _partition;
    Data::SharedPtr _data;
    Likelihood::SharedPtr _likelihood;
    TreeSummary::SharedPtr _tree_summary;
    bool _use_gpu;
    bool _ambig_missing;
    
    std::vector<std::string> partition_subsets;
    partition_subsets.push_back("first:1-20");
    partition_subsets.push_back("second:21-40");
    partition_subsets.push_back("third:41-60");
    //partition_subsets.push_back("rbcL[codon,plantplastid]:1-10");
    //partition_subsets.push_back("rbcL[protein]:11-40");
    //partition_subsets.push_back("morph[standard]:41-45");
    _partition.reset(new Partition());
    for (auto s : partition_subsets) {
        _partition->parseSubsetDefinition(s);
    }
    try {
        std::cout << "\n*** Reading and storing the data in the file" << std::endl;
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        _data->getDataFromFile("rbcL.nex");
        //_data->getDataFromFile("datatest.nex");
        
        // Report information about data partition subsets
        unsigned nsubsets = _data->getNumSubsets();
        std::cout << "\nNumber of taxa: " << _data->getNumTaxa() << std::endl;
        std::cout << "Number of partition subset: " << nsubsets << std::endl;
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            DataType dt = _partition->getDataTypeForSubset(subset);
            std::cout << "Subset " << (subset + 1) << "(" << _data->getSubsetName(subset) << ")" << std::endl;
            std::cout << "data type: " << dt.getDataTypeAsString() << std::endl;
            std::cout << "site: " << _data->calcSeqLenInSubset(subset) << std::endl;
            std::cout << "pattern: " << _data->getNumPatternsInSubset(subset) << std::endl;
        }
        
        std::cout << "\n*** Resource available to BeagleLib" << _likelihood->beagleLibVersion() << ":\n";
        std::cout << _likelihood->availableResources() << std::endl;
        
        std::cout << "\n*** Creating the likelihood calculator" << std::endl;
        _likelihood = Likelihood::SharedPtr(new Likelihood());
        _use_gpu = false;
        _ambig_missing = true;
        _likelihood->setPreferGPU(_use_gpu);
        _likelihood->setAmbiguityEqualsMissing(_ambig_missing);
        _likelihood->setData(_data);
        _likelihood->initBeagleLib();
        
        std::cout << "\n*** Reading and storing the first tree in the file" << std::endl;
        _tree_summary = TreeSummary::SharedPtr(new TreeSummary());
        _tree_summary->readTreefile("rbcLjc.tree", 0);
        Tree::SharedPtr tree = _tree_summary->getTree(0);
        
        if (tree->numLeaves() != _data->getNumTaxa())
            throw XStrom(boost::format("Number of taxa in tree (%d) dost not equal the number of taxa in data matrix (%d)") % tree->numLeaves() % _data->getNumTaxa());
        
        std::cout << "\n*** Calculating the likelihood of the tree" << std::endl;
        double lnL = _likelihood->calcLogLikelihood(tree);
        std::cout << boost::str(boost::format("log likelihood = %0.5f") % lnL) << std::endl;
        //std::cout << lnL << std::endl;
        
    } catch(NxsException x) {
        std::cerr << "Program aborting due to error encounted reading tree file" << std::endl;
        std::cerr << x.what() << std::endl;
    }
    
    
    std::cout << "\n Finished!" << std::endl;
    
    /*
    Strom strom;
    try {
        strom.processCommandLineOptions(argc, argv);
        strom.run();
    }
    catch(std::exception & x) {
        std::cerr << "Exception: " << x.what() << std::endl;
        std::cerr << "Aborted." << std::endl;
    }
    catch(...) {
        std::cerr << "Exception of unknown type!\n";
    }
    */
    
    return 0;
}
