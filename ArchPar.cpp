/*
 * ArchPar.cpp
 *
 *  Created on: 2015.11.30
 *      Author: Yuan Kai
 */

# include "ArchPar.h"
# include "ArchTools.h"
# include <cstring>
# include <cstdlib>

ArchPar::ArchPar(int & argc, char ** &argv)
{
	bool IsGzInput(false), IsAfrRank(false), IsMinBinSites(false), IsBinSize(false),
			IsPathAfrInput(false), IsPathArchInput(false), IsPathTestInput(false), IsPathSumOut(false), IsPathSegOut(false);
	for(int i = 1; i < argc; ++i)
	{
		if(strcmp(argv[i],"-gzInput") == 0)
		{
			if(IsGzInput)
				err_print("Redundant input of \"gzInput\" parameter!");
			else
			{
				IsGzInput = true;
				gzInput = true;
			}
		}
		else if(strcmp(argv[i],"-AfrRank") == 0)
		{
			if(IsAfrRank)
				err_print("Redundant input of \"AfrRank\" parameter!");
			else
			{
				IsAfrRank = true;
				AfrRank = atoi(argv[++i]);
			}
		}
		else if(strcmp(argv[i],"-MinBinSites") == 0)
		{
			if(IsMinBinSites)
				err_print("Redundant input of \"MinBinSites\" parameter!");
			else
			{
				IsMinBinSites = true;
				MinBinSites = atoi(argv[++i]);
			}
		}
		else if(strcmp(argv[i],"-BinSize") == 0)
		{
			if(IsBinSize)
				err_print("Redundant input of \"BinSize\" parameter!");
			else
			{
				IsBinSize = true;
				BinSize = atol(argv[++i]);
			}
		}
		else if(strcmp(argv[i],"-Afr") == 0)
		{
			if(IsPathAfrInput)
				err_print("Redundant input of \"Afr\" parameter!");
			else
			{
				IsPathAfrInput = true;
				PathAfrInput = argv[++i];
			}
		}
		else if(strcmp(argv[i],"-Arch") == 0)
		{
			if(IsPathArchInput)
				err_print("Redundant input of \"Arch\" parameter!");
			else
			{
				IsPathArchInput = true;
				PathArchInput = argv[++i];
			}
		}
		else if(strcmp(argv[i],"-Test") == 0)
		{
			if(IsPathTestInput)
				err_print("Redundant input of \"Test\" parameter!");
			else
			{
				IsPathTestInput = true;
				PathTestInput = argv[++i];
			}
		}
		else if(strcmp(argv[i],"-Summary") == 0)
		{
			if(IsPathSumOut)
				err_print("Redundant input of \"Summary\" parameter!");
			else
			{
				IsPathSumOut = true;
				PathSumOut = argv[++i];
			}
		}
		else if(strcmp(argv[i],"-Seg") == 0)
		{
			if(IsPathSegOut)
				err_print("Redundant input of \"Seg\" parameter!");
			else
			{
				IsPathSegOut = true;
				PathSegOut = argv[++i];
			}
		}
		else
			err_print((std::string)("Cannot identify ")+(std::string)argv[i]);
	}
	if(!IsGzInput)
		gzInput = false;
	if(!IsAfrRank)
		AfrRank = 0;
	if(!IsMinBinSites)
		MinBinSites = 10;
	if(!IsBinSize)
		BinSize = 10000;
	if(!IsPathAfrInput || !FileCheck(PathAfrInput))
		err_print("Invalid path of African data!");
	if(!IsPathArchInput || !FileCheck(PathArchInput))
		err_print("Invalid path of Archaic data!");
	if(!IsPathTestInput || !FileCheck(PathTestInput))
		err_print("Invalid path of Test data!");
	if(!IsPathSumOut)
		err_print("Invalid path of output summary file!");
	if(!IsPathSegOut)
		err_print("Invalid path of segment file!");
}


