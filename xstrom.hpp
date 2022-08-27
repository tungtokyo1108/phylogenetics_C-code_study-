//
//  xstrom.hpp
//  Phylogenetic_code
//
//  Created by TungDang on 2022/08/26.
//

#ifndef xstrom_h
#define xstrom_h

#include <boost/format.hpp>

namespace strom {

    class XStrom : public std::exception {
    
        public:
            
                XStrom() throw() {}
                XStrom(const std::string s) throw() : _msg() {_msg = s;}
                XStrom(const boost::format & f) throw() : _msg() {_msg = boost::str(f);}
            virtual ~XStrom() throw() {}
            const char * what() const throw() {return _msg.c_str();}
        
        private:
            std::string _msg;
    };

}


#endif /* xstrom_h */
