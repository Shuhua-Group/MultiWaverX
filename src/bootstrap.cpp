/*
 * bootstrap.cpp
 *
 *  Created on: June 8, 2021
 *      Author: yuankai & zhangrui
 */

# include "ancSegs.h"

# include <cstdlib>

void ancBootstrap::bootstrap(const ancSegs & dat)
{
	totNumSegs = dat.totNumSegs;
	minLen = dat.minLen;
	totLen = 0;
	segs.resize(dat.segs.size());
	for(std::vector<std::vector<double> >::iterator it = segs.begin(); \
			it != segs.end(); ++it)
		it->reserve(dat.totNumSegs);
	ancProp.resize(dat.ancProp.size(), 0);
	popLabels.insert(popLabels.begin(), dat.popLabels.begin(), dat.popLabels.end());

	double ttotLen(0);
	for(int i = 0 ; i < dat.totNumSegs; ++i)
	{
		int index = (double) rand() / RAND_MAX * dat.totNumSegs;
		double tseg = dat.segs4bootstrap[index];
		int lab = dat.popLabels4bootstrap[index];
		if(tseg > minLen)
		{
			segs[lab].push_back(tseg);
		}	
		ancProp[lab] += tseg;
		totLen += tseg;
	}

	for(unsigned int i = 0 ; i < ancProp.size(); ++i)
	{
		ancProp[i] /= totLen;
	}
}

void ancBootstrap::bootstrapMultiX(const ancSegs &dat, const std::vector<int> &Nwave)
{
	totNumSegs = dat.totNumSegs;
	minLen = dat.minLen;
	totLen = 0;
	segs.resize(dat.segs.size());
	segProp.resize(dat.segs.size());
	for (int i = 0; i < Nwave.size(); i++)
		segProp[i].resize(Nwave[i]);
	for(std::vector<std::vector<double> >::iterator it = segs.begin(); \
			it != segs.end(); ++it)
		it->reserve(dat.totNumSegs);
	ancProp.resize(dat.ancProp.size(), 0);
	popLabels.insert(popLabels.begin(), dat.popLabels.begin(), dat.popLabels.end());

	double ttotLen(0);
	for(int i = 0 ; i < dat.totNumSegs; ++i)
	{
		int index = (double) rand() / RAND_MAX * dat.totNumSegs;
		double tseg = dat.segs4bootstrap[index];
		int lab = dat.popLabels4bootstrap[index];
		
		if(tseg > minLen)
		{
			segs[lab].push_back(tseg);
		/*	double tempRand = (double) rand()/RAND_MAX;
			double tempSum = 0;
			for (int j = 0; j < dat.EMprop[lab][dat.segMap[index]].size(); j++)
			{
				tempSum += dat.EMprop[lab][dat.segMap[index]][j];
				if (tempRand < tempSum)
				{
					segProp[lab][j] += tseg;
					break;
				}
			} */

			for (int j = 0; j < dat.EMprop[lab][dat.segMap[index]].size(); j++)
				segProp[lab][j] += dat.EMprop[lab][dat.segMap[index]][j]*tseg;
		}
		ancProp[lab] += tseg;
		totLen += tseg; 	
	}

	for(unsigned int i = 0 ; i < ancProp.size(); ++i)
		ancProp[i] /= totLen;
}

void ancBootstrap::clean()
{
	segs.clear();
	std::vector<std::vector<double> > ssd;
	segs = ssd;
}

