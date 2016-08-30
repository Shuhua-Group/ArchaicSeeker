/*
 * ArchData.h
 *
 *  Created on: 2015.11.30
 *      Author: Yuan Kai
 */

#ifndef ARCHDATA_H_
#define ARCHDATA_H_

# include "ArchPar.h"
# include <vector>
# include <algorithm>
# include <cstring>

void gzSeeker (const ArchPar & par);

void Seeker (const ArchPar & par);

void OutPut (const ArchPar & par, const std::vector<std::vector<std::vector<bool> > > & similarity,
		const std::vector<std::vector<long> > & eallele, const std::vector<std::string>  & ids, const int &valid);

const int BUFFSIZE = 4096;

const int NCHR = 22;

#endif /* ARCHDATA_H_ */
