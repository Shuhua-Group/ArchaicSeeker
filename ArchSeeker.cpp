/*
 * ArchSeeker.cpp
 *
 *  Created on: 2015.11.30
 *      Author: Yuan Kai
 */

# include "ArchPar.h"
# include "ArchData.h"

using namespace std;

int main(int argc,char **argv)
{
	ArchPar par(argc, argv);
	if(par.gzInput)
		gzSeeker(par);
	else
		Seeker(par);
}


