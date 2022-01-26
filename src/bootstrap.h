/*
 * bootstrap.h
 *
 *  Created on: June 8, 2021
 *      Author: yuankai & zhangrui
 */

#ifndef BOOTSTRAP_H_
#define BOOTSTRAP_H_

class ancSegs;

class ancBootstrap
{
public:
	void bootstrap(const ancSegs & dat);
	void bootstrapMultiX(const ancSegs &dat, const std::vector<int> &Nwave); 
	void clean();

	std::vector<std::vector<double> > segs;
	std::vector<double> ancProp;
	std::vector<std::string> popLabels;
	std::vector<std::vector<double>> segProp;

	int totNumSegs;
	double minLen, totLen;
};



#endif /* BOOTSTRAP_H_ */
