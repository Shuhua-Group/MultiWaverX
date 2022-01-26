/*
 * ancSegs.h
 *
 *  Created on: June 8, 2021
 *      Author: yuankai & zhangrui
 */

#ifndef ANCSEGS_H_
#define ANCSEGS_H_

# include <iostream>
# include <vector>

# include "bootstrap.h"

class ancSegs
{
public:
	ancSegs(const std::string & path, double cutoff);

	friend void ancBootstrap::bootstrap(const ancSegs & dat);
	friend void ancBootstrap::bootstrapMultiX(const ancSegs &dat, const std::vector<int> &Nwave);

	std::vector<std::vector<double> > segs;
	std::vector<std::vector<std::vector<double>>> EMprop;
	std::vector<int> segMap;
	std::vector<double> ancProp;
	std::vector<std::string> popLabels;

	int totNumSegs;
	double minLen, totLen;

private:
	std::vector<double> segs4bootstrap;
	std::vector<int> popLabels4bootstrap;
};

#endif /* ANCSEGS_H_ */
