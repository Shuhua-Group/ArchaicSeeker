/*
 * ArchPar.h
 *
 *  Created on: 2015.11.30
 *      Author: Yuan Kai
 */

#ifndef ARCHPAR_H_
#define ARCHPAR_H_

# include <iostream>
# include <vector>

class ArchPar
{
public:
	ArchPar(int & argc, char ** & argv);
	ArchPar() : gzInput(false), AfrRank(0), MinBinSites(10), BinSize(10000),
			PathAfrInput(""), PathArchInput(""), PathTestInput(""), PathSumOut(""), PathSegOut("") {};

	bool gzInput;
	int AfrRank, MinBinSites;
	long BinSize;
	std::string PathAfrInput, PathArchInput, PathTestInput, PathSumOut, PathSegOut;
};



#endif /* ARCHPAR_H_ */
