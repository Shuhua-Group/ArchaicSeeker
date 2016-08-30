/*
 * ArchTools.h
 *
 *  Created on: 2015.11.30
 *      Author: Yuan Kai
 */

#ifndef ARCHTOOLS_H_
#define ARCHTOOLS_H_

# include <iostream>
# include <fstream>
# include <cstdlib>

inline void err_print(const std::string & err_info)
{
	std::cerr << err_info << std::endl;
	exit(1);
}

inline bool FileCheck(const std::string & path)
{
	std::ifstream fp(path.c_str());
	if(fp)
		return true;
	else
		return false;
}

#endif /* ARCHTOOLS_H_ */
