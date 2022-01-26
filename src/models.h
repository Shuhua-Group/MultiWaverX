/*
 * models.h
 *
 *  Created on: June 8, 2021
 *      Author: yuankai & zhangrui
 */

#ifndef MODELS_H_
#define MODELS_H_

# include "ancSegs.h"

# include <vector>
# include <iostream>
# include <algorithm>

class models
{
public:
	models(ancSegs & dat) : segs(dat.segs), prop(dat.ancProp), \
			popLabels(dat.popLabels), hiT(0), cgfdT(0), cgfrT(0), gaT(0), \
			hiLlk(0), cgfdLlk(0), cgfrLlk(0), gaLlk(0), multiLlk(0), intMultiLlk(0), \
			minLen(dat.minLen), multiScenarioCount(0){};
	models(ancBootstrap & dat) : segs(dat.segs), prop(dat.ancProp), segProp(dat.segProp), \
			popLabels(dat.popLabels), hiT(0), cgfdT(0), cgfrT(0), gaT(0), \
			hiLlk(0), cgfdLlk(0), cgfrLlk(0), gaLlk(0), multiLlk(0), intMultiLlk(0), \
			minLen(dat.minLen), multiScenarioCount(0){};

	void solve(const std::string &chr, std::vector<int> &Aorder, const std::string &Amodel);
	void solveMultiX(const std::vector<int> &Aorder, const std::vector<double> &AT, const std::vector<int> &APop);

	const std::vector<std::vector<double> > & segs;
	const std::vector<double> & prop;
	const std::vector<std::string> & popLabels;
	std::vector<std::vector<std::vector<double>>> EMprop;
	std::vector<std::vector<double>> segProp;

	int hiT, cgfdT, cgfrT, gaT;
	double hiLlk, cgfdLlk, cgfrLlk, gaLlk, multiLlk, intMultiLlk;
	double minLen;
	int multiScenarioCount;
	std::vector<double> multiT, multiProp;
	std::vector<int> multiPop;

	static bool isHI, isCGF, isGA, isMulti;
//	static int MultiMaxWave, MultiMaxPop, MaxInter;
	static int MultiMaxWave, MaxInter;
	static double epsilon, critical, minProp;
private:
	inline double HI(int t) const;
	inline double CGFR(int t) const;
	inline double CGFD(int t) const;
	inline double CGF(const std::vector<double> & donor, const std::vector<double> & recipient,\
			double m, int t) const;
	inline double GA(int t) const;
	inline double Multi(const std::vector<double> & psegs, double* lambda, \
			double* prop, int nwave) const;
	double MultiEM(const std::vector<double> & psegs, const double ancProp, \
			double* & curLambda, double* & curProp, int & curNwave, \
			const std::string &chr, const int &popIndex, const int &npop);

	void bisectionSearch(int & t, double & llk, double (models::*getLlk)(int) const);
	void multiWaveInfer(const std::string &chr, std::vector<int> &Aorder);

public:
	std::string info(const std::string &chr) const;
	std::string summary(std::vector<int> & pop, std::vector<int> & t, \
			std::vector<double> & pro, const std::string &chr) const;
};

template<class T>
std::vector<std::vector<T> > perm(std::vector<T> &seq)
{
	std::vector<std::vector<T> > result;
	int size = seq.size();
	std::sort(seq.begin(), seq.end());
	do
	{
		if (seq.at(size - 1) == seq.at(size - 2))
			continue;
		int rsize = result.size();
		if (rsize > 0 && result.at(rsize - 1).at(size - 2) == seq.at(size - 1))
			continue;
		result.push_back(seq);
	} while (std::next_permutation(seq.begin(), seq.end()));
	return result;
}

void Summary(const models& md, const std::vector<models *>& bootstrap, \
		const std::vector<std::string>&  popLabels, double ci, \
		std::ofstream & fpo, std::ofstream & fplog, \
		std::vector<double>& Aprop, std::string &Amodel, \
		const std::string &chr, std::vector<int> &AT);

class SegData
{
public:
	SegData(const std::vector<double> &_psegs, int &_k): psegs(_psegs), k(_k){};
	const std::vector<double> &psegs;
	int k;
};

class event
{
public:
	int pop;
	double t, alpha;
};
bool cmp(const event & a, const event &b );
double direcRatio(const double &mA, const std::vector<double> &mX, const std::string &sex, const int &T, const std::string &model);

#endif /* MODELS_H_ */
