/*
 * ArchData.cpp
 *
 *  Created on: 2015.11.30
 *      Author: Yuan Kai
 */

# include "ArchData.h"
# include "ArchTools.h"
# include "ArchPar.h"
# include <zlib.h>
# include <vector>
# include <fstream>
# include <sstream>

using namespace std;

void gzSeeker(const ArchPar & par)
{
	gzFile FpAfr, FpArch, FpTest;
	FpAfr = gzopen(par.PathAfrInput.c_str(), "r");
	FpArch = gzopen(par.PathArchInput.c_str(), "r");
	FpTest = gzopen(par.PathTestInput.c_str(), "r");
	ofstream FpSum(par.PathSumOut.c_str()), FpSeg(par.PathSegOut.c_str());

	char buff[BUFFSIZE];
	bool EndLine(false);
	string tmp;
	string head;
	while(gzgets(FpAfr,buff,BUFFSIZE))
	{
		EndLine = false;
		if(buff[0] == '#' && buff[1] != '#')
		{
			EndLine = true;
			head = buff;
		}
		while(strlen(buff) == BUFFSIZE-1 && buff[BUFFSIZE - 2] != '\n')
		{
			gzgets(FpAfr,buff,BUFFSIZE);
			if(EndLine)
				head += (string) buff;
		}
		if(EndLine)
			break;
	}
	istringstream HeadStream (head);
	for(int i = 0;i < 9 ;++i)
		HeadStream >> tmp;
	int nafr(0);
	while(HeadStream >> tmp)
		nafr += 2;
	while(gzgets(FpArch,buff,BUFFSIZE))
	{
		EndLine = false;
		if(buff[0] == '#' && buff[1] != '#')
		{
			EndLine = true;
			head = buff;
		}
		while(strlen(buff) == BUFFSIZE-1 && buff[BUFFSIZE - 2] != '\n')
		{
			gzgets(FpArch,buff,BUFFSIZE);
			if(EndLine)
				head += (string) buff;
		}
		if(EndLine)
			break;
	}
	HeadStream.clear();
	HeadStream.str(buff);
	for(int i = 0; i < 9; ++i)
		HeadStream >> tmp;
	int narch(0);
	while(HeadStream >> tmp)
		narch += 2;
	while(gzgets(FpTest,buff,BUFFSIZE))
	{
		EndLine = false;
		if(buff[0] == '#' && buff[1] != '#')
		{
			EndLine = true;
			head = buff;
		}
		while(strlen(buff) == BUFFSIZE-1 && buff[BUFFSIZE -2] != '\n')
		{
			gzgets(FpTest,buff,BUFFSIZE);
			if(EndLine)
				head += (string) buff;
		}
		if(EndLine)
			break;
	}
	HeadStream.clear();
	HeadStream.str(buff);
	for(int i=0; i < 9; ++i)
		HeadStream >> tmp;
	string id;
	vector<string> ids;
	while(HeadStream >> id)
	{
		ids.push_back(id+"_1");
		ids.push_back(id+"_2");
	}
	int nad = ids.size();

	vector < vector < vector < bool > > > similarity;
	int valid(0);
	vector<vector<long> >eallele;
	similarity.resize(NCHR);
	for(int i = 0 ; i < NCHR ; ++i)
		similarity[i].resize(nad);

	eallele.resize(NCHR);

	vector<vector<int> >DifAfr, DifArch;
	DifAfr.resize(nad);
	DifArch.resize(nad);
	for(int i = 0;i < nad; ++i)
	{
		DifAfr[i].resize(nafr);
		DifArch[i].resize(narch);
		for(int j = 0;j < nafr; ++j)
		{
			DifAfr[i][j] = 0;
		}
		for(int j = 0;j < narch; ++j)
		{
			DifArch[i][j] = 0;
		}
	}

	vector<char>HapAfr, HapArch, HapTest;
	HapAfr.resize(nafr);
	HapArch.resize(narch);
	HapTest.resize(nad);

	int nchr(1), CountSite(0);
	long from(0), to(par.BinSize);
	while(gzgets(FpAfr, buff, BUFFSIZE))
	{
		int NumAfr(0), NumArch(0);	//count for e allele
		string line;
		bool ExLen(false);
		istringstream stream;
		int TChrAfr, TChrArch, TChrTest;
		long TPosAfr, TPosArch, TPosTest;
		bool OKAfr(false), OKArch(false), OKTest(false);
		string TA1Afr, TA2Afr, TA1Arch, TA2Arch, TA1Test, TA2Test;
		if(strlen(buff) == BUFFSIZE - 1 && buff[BUFFSIZE - 2] != '\n')
		{
			ExLen = true;
			line = buff;
			while(strlen(buff) == BUFFSIZE - 1 && buff[BUFFSIZE - 2] != '\n')
			{
				gzgets(FpAfr, buff, BUFFSIZE);
				line += (string) buff;
			}
		}
		if(ExLen)
			stream.str(line);
		else
			stream.str(buff);
		stream >> TChrAfr >> TPosAfr >> tmp >> TA1Afr >> TA2Afr >> tmp >> tmp >> tmp >> tmp;
		for(int i = 0;i < nafr/2 ;++i)
		{
			stream >> tmp;
			HapAfr[i * 2] = tmp[0];
			HapAfr[i * 2 + 1] = tmp[2];
			if(tmp[0] != '.' || tmp[2] != '.')
				OKAfr = true;
			if(tmp[1] != '|')
			{
				cerr << "Afr: " << TChrAfr << '\t' << TPosAfr << endl;
				cerr << tmp <<endl;
				err_print("Data should be phased!");
			}
			if(tmp[0] == '1')
				++NumAfr;
			if(tmp[2] == '1')
				++NumAfr;
		}
		if(!stream)
			err_print("Insufficient number of haplotypes in African data!");
		stream >> tmp;
		if(stream)
			err_print("Redundant information in African data!");

		gzgets(FpArch, buff, BUFFSIZE);
		ExLen = false;
		stream.clear();
		if(strlen(buff) == BUFFSIZE - 1 && buff[BUFFSIZE - 2] != '\n')
		{
			ExLen = true;
			line = buff;
			while(strlen(buff) == BUFFSIZE - 1 && buff[BUFFSIZE - 2] != '\n')
			{
				gzgets(FpArch, buff, BUFFSIZE);
				line += (string) buff;
			}
		}
		if(ExLen)
			stream.str(line);
		else
			stream.str(buff);
		stream >> TChrArch >> TPosArch >> tmp >> TA1Arch >> TA2Arch >> tmp >> tmp >> tmp >> tmp;
		if(TChrAfr != TChrArch || TPosAfr != TPosArch || TA1Afr != TA1Arch || TA2Afr != TA2Arch)
		{
			cerr << "Afr: \t" << TChrAfr << '\t' << TPosAfr << '\t' << TA1Afr << '\t' << TA2Afr <<endl;
			cerr << "Arch: \t" << TChrArch << '\t' << TPosArch << '\t' << TA1Arch << '\t' << TA2Arch <<endl;
			err_print("Inconsistent SNP information between African and Archaic data!");
		}
		for(int i = 0;i < narch/2 ;++i)
		{
			stream >> tmp;
			HapArch[i * 2] = tmp[0];
			HapArch[i * 2 + 1] = tmp[2];
			if(tmp[0] != '.' || tmp[2] != '.')
				OKArch = true;
			if(tmp[1] != '|')
			{
				cerr << "Arch: " << TChrAfr << '\t' << TPosAfr << endl;
				cerr << tmp <<endl;
				err_print("Data should be phased!");
			}
			if(tmp[0] == '1')
				++NumArch;
			if(tmp[2] == '1')
				++NumArch;
		}
		if(!stream)
			err_print("Insufficient number of haplotypes in Archaic data!");
		stream >> tmp;
		if(stream)
			err_print("Redundant information in Archaic data!");

		gzgets(FpTest, buff, BUFFSIZE);
		ExLen = false;
		stream.clear();
		if(strlen(buff) == BUFFSIZE - 1 && buff[BUFFSIZE - 2] != '\n')
		{
			ExLen = true;
			line = buff;
			while(strlen(buff) == BUFFSIZE - 1 && buff[BUFFSIZE - 2] != '\n')
			{
				gzgets(FpTest, buff, BUFFSIZE);
				line += (string) buff;
			}
		}
		if(ExLen)
			stream.str(line);
		else
			stream.str(buff);
		stream >> TChrTest >> TPosTest >> tmp >> TA1Test >> TA2Test >> tmp >> tmp >> tmp >> tmp;
		if(TChrAfr != TChrTest || TPosAfr != TPosTest || TA1Afr != TA1Test || TA2Afr != TA2Test)
		{
			cerr << "Afr: \t" << TChrAfr << '\t' << TPosAfr << '\t' << TA1Afr << '\t' << TA2Afr <<endl;
			cerr << "Test: \t" << TChrTest << '\t' << TPosTest << '\t' << TA1Test << '\t' << TA2Test <<endl;
			err_print("Inconsistent SNP information between African and Test data!");
		}
		for(int i = 0;i < nad/2 ;++i)
		{
			stream >> tmp;
			HapTest[i * 2] = tmp[0];
			HapTest[i * 2 + 1] = tmp[2];
			if(tmp[0] != '.' || tmp[2] != '.')
				OKTest = true;
			if(tmp[1] != '|')
			{
				cerr << "Test: " << TChrAfr << '\t' << TPosAfr << endl;
				cerr << tmp <<endl;
				err_print("Data should be phased!");
			}
		}
		if(!stream)
			err_print("Insufficient number of haplotypes in Test data!");
		stream >> tmp;
		if(stream)
			err_print("Redundant information in Test data!");

		if(TChrAfr != nchr || (TChrAfr == nchr && TPosAfr > to))
		{
			if(CountSite < par.MinBinSites)
			{
				for(int i = 0; i < nad ; ++i)
					similarity[nchr-1][i].push_back(false);
				++valid;
			}
			else
			{
				vector<double> tafr, tarch;
				tafr.resize(nafr);
				int rank = par.AfrRank;
				tarch.resize(narch);
				for(int i = 0; i < nad; ++i)
				{
					for(int j = 0; j < nafr; ++j)
						tafr[j] = (double) DifAfr[i][j] / CountSite;
					sort(tafr.begin(), tafr.end());
					for(int j = 0; j < narch; ++j)
						tarch[j] = (double) DifArch[i][j] / CountSite;
					sort(tarch.begin(), tarch.end());
					if(tafr[rank] - tarch[0] > 0)
						similarity[nchr-1][i].push_back(true);
					else
						similarity[nchr-1][i].push_back(false);
				}
				++valid;
			}
			if(nchr != TChrAfr)
			{
				nchr = TChrAfr;
				from = 0;
				to = par.BinSize;
			}
			else
			{
				to += par.BinSize;
				from += par.BinSize;
			}		
			while(TPosAfr > to)
			{
				for(int i = 0; i < nad ; ++i)
					similarity[nchr-1][i].push_back(false);
				to += par.BinSize;
				from += par.BinSize;
				++valid;
			}
			CountSite = 0;
			for(int i = 0; i < nad; ++i)
			{
				for(int j = 0; j < nafr; ++j)
					DifAfr[i][j] = 0;
				for(int j = 0; j < narch; ++j)
					DifArch[i][j] = 0;
			}
		}

		if(!OKAfr || !OKArch || !OKTest)
			continue;

		for(int i = 0; i < nad; ++i)
		{
			for(int j = 0; j < nafr; ++j)
			{
				if(HapTest[i] != HapAfr[j] && HapTest[i] != '.' && HapAfr[j] != '.')
					++DifAfr[i][j];
			}
			for(int j = 0; j < narch; ++j)
			{
				if(HapTest[i] != HapArch[j] && HapTest[i] != '.' && HapArch[j] != '.')
					++DifArch[i][j];
			}
		}
		++CountSite;
		if(NumAfr < 2 && NumArch >=1)
			eallele[TChrAfr - 1].push_back(TPosAfr);
	}
	if(CountSite < par.MinBinSites)
	{
		for(int i = 0; i < nad ; ++i)
			similarity[nchr-1][i].push_back(false);
	}
	else
	{
		vector<double> tafr, tarch;
		tafr.resize(nafr);
		int rank = par.AfrRank;
		tarch.resize(narch);
		for(int i = 0; i < nad; ++i)
		{
			for(int j = 0; j < nafr; ++j)
				tafr[j] = (double) DifAfr[i][j] / CountSite;
			sort(tafr.begin(), tafr.end());
			for(int j = 0; j < narch; ++j)
				tarch[j] = (double) DifArch[i][j] / CountSite;
			sort(tarch.begin(), tarch.end());
			if(tafr[rank] - tarch[0] > 0)
				similarity[nchr-1][i].push_back(true);
			else
				similarity[nchr-1][i].push_back(false);
		}
		++valid;
	}
	OutPut(par, similarity, eallele,ids, valid);
}

void Seeker(const ArchPar & par)
{
	ifstream FpAfr(par.PathAfrInput.c_str()), FpArch(par.PathArchInput.c_str()), FpTest(par.PathTestInput.c_str());
	ofstream FpSum(par.PathSumOut.c_str()), FpSeg(par.PathSegOut.c_str());

	string buff;
	string tmp;
	while(getline(FpAfr,buff))
	{
		if(buff[0] == '#' && buff[1] != '#')
			break;
	}
	istringstream HeadStream (buff);
	for(int i = 0;i < 9 ;++i)
		HeadStream >> tmp;
	int nafr(0);
	while(HeadStream >> tmp)
		nafr += 2;
	while(getline(FpArch,buff))
	{
		if(buff[0] == '#' && buff[1] != '#')
			break;
	}
	HeadStream.clear();
	HeadStream.str(buff);
	for(int i = 0; i < 9; ++i)
		HeadStream >> tmp;
	int narch(0);
	while(HeadStream >> tmp)
		narch += 2;
	while(getline(FpTest,buff))
	{
		if(buff[0] == '#' && buff[1] != '#')
			break;
	}
	HeadStream.clear();
	HeadStream.str(buff);
	for(int i=0; i < 9; ++i)
		HeadStream >> tmp;
	string id;
	vector<string> ids;
	while(HeadStream >> id)
	{
		ids.push_back(id+"_1");
		ids.push_back(id+"_2");
	}
	int nad = ids.size();

	vector < vector < vector < bool > > > similarity;
	int valid(0);
	vector<vector<long> >eallele;
	similarity.resize(NCHR);
	for(int i = 0 ; i < NCHR ; ++i)
		similarity[i].resize(nad);

	eallele.resize(NCHR);

	vector<vector<int> >DifAfr, DifArch;
	DifAfr.resize(nad);
	DifArch.resize(nad);
	for(int i = 0;i < nad; ++i)
	{
		DifAfr[i].resize(nafr);
		DifArch[i].resize(narch);
		for(int j = 0;j < nafr; ++j)
		{
			DifAfr[i][j] = 0;
		}
		for(int j = 0;j < narch; ++j)
		{
			DifArch[i][j] = 0;
		}
	}

	vector<char>HapAfr, HapArch, HapTest;
	HapAfr.resize(nafr);
	HapArch.resize(narch);
	HapTest.resize(nad);

	int nchr(1), CountSite(0);
	long from(0), to(par.BinSize);
	while(getline(FpAfr, buff))
	{
		int NumAfr(0), NumArch(0);	//count for e allele
		string line;
		istringstream stream;
		int TChrAfr, TChrArch, TChrTest;
		long TPosAfr, TPosArch, TPosTest;
		bool OKAfr(false), OKArch(false), OKTest(false);
		string TA1Afr, TA2Afr, TA1Arch, TA2Arch, TA1Test, TA2Test;
		stream.str(buff);
		stream >> TChrAfr >> TPosAfr >> tmp >> TA1Afr >> TA2Afr >> tmp >> tmp >> tmp >> tmp;
		for(int i = 0;i < nafr/2 ;++i)
		{
			stream >> tmp;
			HapAfr[i * 2] = tmp[0];
			HapAfr[i * 2 + 1] = tmp[2];
			if(tmp[0] != '.' || tmp[2] != '.')
				OKAfr = true;
			if(tmp[1] != '|')
			{
				cerr << "Afr: " << TChrAfr << '\t' << TPosAfr << endl;
				cerr << tmp <<endl;
				err_print("Data should be phased!");
			}
			if(tmp[0] == '1')
				++NumAfr;
			if(tmp[2] == '1')
				++NumAfr;
		}
		if(!stream)
			err_print("Insufficient number of haplotypes in African data!");
		stream >> tmp;
		if(stream)
			err_print("Redundant information in African data!");

		getline(FpArch, buff);
		stream.clear();
		stream.str(buff);
		stream >> TChrArch >> TPosArch >> tmp >> TA1Arch >> TA2Arch >> tmp >> tmp >> tmp >> tmp;
		if(TChrAfr != TChrArch || TPosAfr != TPosArch || TA1Afr != TA1Arch || TA2Afr != TA2Arch)
		{
			cerr << "Afr: \t" << TChrAfr << '\t' << TPosAfr << '\t' << TA1Afr << '\t' << TA2Afr <<endl;
			cerr << "Arch: \t" << TChrArch << '\t' << TPosArch << '\t' << TA1Arch << '\t' << TA2Arch <<endl;
			err_print("Inconsistent SNP information between African and Archaic data!");
		}
		for(int i = 0;i < narch/2 ;++i)
		{
			stream >> tmp;
			HapArch[i * 2] = tmp[0];
			HapArch[i * 2 + 1] = tmp[2];
			if(tmp[0] != '.' || tmp[2] != '.')
				OKArch = true;
			if(tmp[1] != '|')
			{
				cerr << "Arch: " << TChrAfr << '\t' << TPosAfr << endl;
				cerr << tmp <<endl;
				err_print("Data should be phased!");
			}
			if(tmp[0] == '1')
				++NumArch;
			if(tmp[2] == '1')
				++NumArch;
		}
		if(!stream)
			err_print("Insufficient number of haplotypes in Archaic data!");
		stream >> tmp;
		if(stream)
			err_print("Redundant information in Archaic data!");

		getline(FpTest, buff);
		stream.clear();
		stream.str(buff);
		stream >> TChrTest >> TPosTest >> tmp >> TA1Test >> TA2Test >> tmp >> tmp >> tmp >> tmp;
		if(TChrAfr != TChrTest || TPosAfr != TPosTest || TA1Afr != TA1Test || TA2Afr != TA2Test)
		{
			cerr << "Afr: \t" << TChrAfr << '\t' << TPosAfr << '\t' << TA1Afr << '\t' << TA2Afr <<endl;
			cerr << "Test: \t" << TChrTest << '\t' << TPosTest << '\t' << TA1Test << '\t' << TA2Test <<endl;
			err_print("Inconsistent SNP information between African and Test data!");
		}
		for(int i = 0;i < nad/2 ;++i)
		{
			stream >> tmp;
			HapTest[i * 2] = tmp[0];
			HapTest[i * 2 + 1] = tmp[2];
			if(tmp[0] != '.' || tmp[2] != '.')
				OKTest = true;
			if(tmp[1] != '|')
			{
				cerr << "Test: " << TChrAfr << '\t' << TPosAfr << endl;
				cerr << tmp <<endl;
				err_print("Data should be phased!");
			}
		}
		if(!stream)
			err_print("Insufficient number of haplotypes in Test data!");
		stream >> tmp;
		if(stream)
			err_print("Redundant information in Test data!");

		if(TChrAfr != nchr || (TChrAfr == nchr && TPosAfr > to))
		{
			if(CountSite < par.MinBinSites)
			{
				for(int i = 0; i < nad ; ++i)
					similarity[nchr-1][i].push_back(false);
				++valid;
			}
			else
			{
				vector<double> tafr, tarch;
				tafr.resize(nafr);
				int rank = par.AfrRank;
				tarch.resize(narch);
				for(int i = 0; i < nad; ++i)
				{
					for(int j = 0; j < nafr; ++j)
						tafr[j] = (double) DifAfr[i][j] / CountSite;
					sort(tafr.begin(), tafr.end());
					for(int j = 0; j < narch; ++j)
						tarch[j] = (double) DifArch[i][j] / CountSite;
					sort(tarch.begin(), tarch.end());
					if(tafr[rank] - tarch[0] > 0)
						similarity[nchr-1][i].push_back(true);
					else
						similarity[nchr-1][i].push_back(false);
				}
				++valid;
			}
			if(nchr != TChrAfr)
			{
				nchr = TChrAfr;
				from = 0;
				to = par.BinSize;
			}
			else
			{
				to += par.BinSize;
				from += par.BinSize;
			}
			while(TPosAfr > to)
			{
				for(int i = 0; i < nad ; ++i)
					similarity[nchr-1][i].push_back(false);
				to += par.BinSize;
				from += par.BinSize;
				++valid;
			}
			CountSite = 0;
			for(int i = 0; i < nad; ++i)
			{
				for(int j = 0; j < nafr; ++j)
					DifAfr[i][j] = 0;
				for(int j = 0; j < narch; ++j)
					DifArch[i][j] = 0;
			}
		}

		if(!OKAfr || !OKArch || !OKTest)
			continue;

		for(int i = 0; i < nad; ++i)
		{
			for(int j = 0; j < nafr; ++j)
			{
				if(HapTest[i] != HapAfr[j] && HapTest[i] != '.' && HapAfr[j] != '.')
					++DifAfr[i][j];
			}
			for(int j = 0; j < narch; ++j)
			{
				if(HapTest[i] != HapArch[j] && HapTest[i] != '.' && HapArch[j] != '.')
					++DifArch[i][j];
			}
		}
		++CountSite;
		if(NumAfr < 2 && NumArch >=1)
			eallele[TChrAfr - 1].push_back(TPosAfr);
	}
	if(CountSite < par.MinBinSites)
	{
		for(int i = 0; i < nad ; ++i)
			similarity[nchr-1][i].push_back(false);
	}
	else
	{
		vector<double> tafr, tarch;
		tafr.resize(nafr);
		int rank = par.AfrRank;
		tarch.resize(narch);
		for(int i = 0; i < nad; ++i)
		{
			for(int j = 0; j < nafr; ++j)
				tafr[j] = (double) DifAfr[i][j] / CountSite;
			sort(tafr.begin(), tafr.end());
			for(int j = 0; j < narch; ++j)
				tarch[j] = (double) DifArch[i][j] / CountSite;
			sort(tarch.begin(), tarch.end());
			if(tafr[rank] - tarch[0] > 0)
				similarity[nchr-1][i].push_back(true);
			else
				similarity[nchr-1][i].push_back(false);
		}
		++valid;
	}
	OutPut(par, similarity, eallele,ids, valid);
}

void OutPut (const ArchPar & par, const std::vector<std::vector<std::vector<bool> > > & similarity,
		const std::vector<std::vector<long> > & eallele, const std::vector<std::string>  & ids, const int &valid)
{

	ofstream fpseg(par.PathSegOut.c_str()), fpsum(par.PathSumOut.c_str());
	for(int i = 0 ; i < ids.size() ; ++i)
	{
		long len(0);
		for(int j = 0; j < similarity.size() ;++j )
		{
			if(eallele[j].size() == 0)
				continue;
			for(int k = 0; k < eallele[j].size() - 1 ; ++k)
			{
				int p,q,kstart(k);
				p = eallele[j][k] / par.BinSize;
				if(!similarity[j][i][p])
					continue;
				q = eallele[j][k + 1] / par.BinSize;
				if(!similarity[j][i][q])
				{
					++k;
					continue;
				}
				if(q - p >= 2)
				{
					int tmp = p + 1;
					bool ok(true);
					while(tmp != q)
					{
						if(!similarity[j][i][tmp])
						{
							ok = false;
							break;
						}
						++tmp;
					}
					if(!ok)
						continue;
				}
				k += 2;
				while(k < eallele[j].size())
				{
					int tmpq = eallele[j][k] / par.BinSize;
					if(!similarity[j][i][tmpq])
					{
						break;
					}
					if(tmpq - q >= 2 )
					{
						int tmp = q + 1;
						bool ok(true);
						while(tmp != tmpq)
						{
							if(!similarity[j][i][tmp])
							{
								ok = false;
								break;
							}
							++tmp;
						}
						if(!ok)
							break;
					}
					q = tmpq;
					++k;
				}
				--k;
				fpseg << j + 1 << '\t' << eallele[j][kstart] << '\t' << eallele[j][k] << '\t' << ids[i] << endl;
				len += eallele[j][k] - eallele[j][kstart];
			}
		}
		fpsum << ids[i] << '\t' << (double) len / valid / par.BinSize << endl;
	}
}
