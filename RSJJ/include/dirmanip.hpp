#ifndef __dirmanip_hpp__included__
#define __dirmanip_hpp__included__

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>

class Folder{
	std::string sdir;

	public:
	Folder(){ sdir = "dir"; };
	Folder(std::string str){ sdir = str; };
	~Folder(){};
	void open(std::string str){ sdir = str; };
	void mkdir(){
		std::stringstream ss;
		ss.str(""); ss<<"test -d "<<sdir<<"-"<<0;

		int iterdir(0);
		while( system( ss.str().c_str() )==0 ){
			iterdir++;
			ss.str(""); ss<<"test -d "<<sdir<<"-"<<iterdir;
		}
		ss.str(""); ss<<"mkdir "<<sdir<<"-"<<iterdir;
		system( ss.str().c_str() );

		ss.str(""); ss<<sdir<<"-"<<iterdir;
		sdir = ss.str();
	};
	std::string getdir(){ return sdir; };
};

#endif
