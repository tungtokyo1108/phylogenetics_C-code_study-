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
    TreeSummary sumt;
    try {
        sumt.readTreefile("test.tre", 1);
    } catch(NxsException x) {
        std::cerr << "Program aborting due to error encounted reading tree file" << std::endl;
        std::cerr << x.what() << std::endl;
    }
    
    sumt.showSummary();
    */
    Partition::SharedPtr _partition;
    Data::SharedPtr _data;
    std::vector<std::string> partition_subsets;
    //partition_subsets.push_back("first:1-10");
    //partition_subsets.push_back("second:11-40");
    //partition_subsets.push_back("third:41-60");
    partition_subsets.push_back("rbcL[codon,plantplastid]:1-20");
    partition_subsets.push_back("rbcL[protein]:21-40");
    partition_subsets.push_back("morph[standard]:41-45");
    _partition.reset(new Partition());
    for (auto s : partition_subsets) {
        _partition->parseSubsetDefinition(s);
    }
    try {
        _data = Data::SharedPtr(new Data());
        _data->setPartition(_partition);
        //_data->getDataFromFile("rbcL.nex");
        _data->getDataFromFile("datatest.nex");
        
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
