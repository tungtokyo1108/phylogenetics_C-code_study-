//
//  main.cpp
//  Phylogenetic_code
//
//  Created by TungDang on 2022/08/20.
//

#include <iostream>
#include "node.hpp"
#include "tree.hpp"
#include "TreeManip.hpp"
#include "genetic_code.hpp"

using namespace strom;

const double Node::_smallest_edge_length = 1.0e-12;

int main(int argc, const char * argv[]) {
    // insert code here...
    
    std::cout << "Starting ... " << std::endl;
    TreeManip tm;
    tm.createTestTree();
    std::cout << tm.makeNewick(3) << std::endl;
    std::cout << tm.makeNewick(3, true) << std::endl;
    std::cout << "\n Finished!" << std::endl;
    
    return 0;
}
