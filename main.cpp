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
#include "tree_summary.hpp"
//#include "strom.hpp"
//#include "xstrom.hpp"

using namespace strom;

const double Node::_smallest_edge_length = 1.0e-12;

//std::string  Strom::_program_name        = "strom";
//unsigned     Strom::_major_version       = 1;
//unsigned     Strom::_minor_version       = 0;

int main(int argc, const char * argv[]) {

    std::cout << "Starting ... " << std::endl;
    
    TreeSummary sumt;
    try {
        sumt.readTreefile("test.tre", 1);
    } catch(NxsException x) {
        std::cerr << "Program aborting due to error encounted reading tree file" << std::endl;
        std::cerr << x.what() << std::endl;
    }
    
    sumt.showSummary();
    
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
