/*
 * ancSegs.cpp
 *
 *  Created on: June 8, 2021
 *      Author: yuankai & zhangrui
 */

# include "ancSegs.h"
# include "exception.h"

# include <fstream>
# include <cstdlib>
# include <map>
# include <set>

ancSegs::ancSegs(const std::string & path, double cutoff) : totNumSegs(0), minLen(cutoff), \
		totLen(0)
{
#ifdef DEBUG
	std::cout << "Start loading data" << std::endl;
	std::cout << "\tcutoff: " << cutoff << std::endl;
#endif //DEBUG
	std::ifstream fps(path.c_str());
	if(!fps)
	{
		except e;
		e << "Error: Cannot open the ancestral segments file!\n";
		e << "    The input path is \"" << path << '\"';
		throw e;
	}

	double start, end;

	std::string lab;
	std::map<std::string, int> popCount;
	std::map<std::string, double> ancLen, tancLen;
	std::set<std::string> popcheck;

	while(fps >> start >> end >> lab)
	{
		double len(end - start);
		if(len > cutoff)
			popcheck.insert(lab);
	}
	int tpop(0);
	for(std::set<std::string>::iterator it = popcheck.begin(); it != popcheck.end(); ++it)
	{
		popCount[*it] = tpop;
		popLabels.push_back(*it);
		++tpop;
	}
	segs.resize(tpop);
	fps.clear();
	fps.close();
	fps.open(path.c_str());

	while(fps >> start >> end >> lab)
	{
		double len(end - start);
		if(len > cutoff)
		{
			segs[popCount[lab]].push_back(len);
		}
		segs4bootstrap.push_back(len);
		segMap.push_back(segs[popCount[lab]].size()-1);
		popLabels4bootstrap.push_back(popCount[lab]);
		++totNumSegs;
		ancLen[lab] += len;
		totLen += len;
	}
	ancProp.resize(popLabels.size(), 0);
	for(unsigned int i = 0 ; i < popLabels.size(); ++i)
		ancProp[i] = ancLen[popLabels[i]] / totLen;
	
#ifdef DEBUG
	std::cout << "TotalNumSegs:\t" << totNumSegs << std::endl;
	std::cout << "TotalLength:\t" << totLen << std::endl;
	for(unsigned int i = 0 ; i < popLabels.size(); ++i)
		std::cout << ancProp[i] << std::endl;
#endif //DEBUG
}


