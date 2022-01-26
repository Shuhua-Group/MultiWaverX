/*
 * models.cpp
 *
 *  Created on: June 8, 2021
 *      Author: yuankai & zhangrui
 */


# include "models.h"

# include <cassert>
# include <cmath>
# include <cstdlib>
# include <map>
# include <sstream>
# include <fstream>
# include <iostream> 	//
# include <set>
# include "boost/algorithm/string.hpp"

bool models::isHI = true;
bool models::isCGF = true;
bool models::isGA = true;
bool models::isMulti = true;

int models::MultiMaxWave = 2;
//int models::MultiMaxPop = 2;
//int models::MultiMaxPop = 5;
int models::MaxInter = 10000;
double models::epsilon = 1e-6;
double models::critical = 0.001;
double models::minProp = 0.05;

void models::solve(const std::string &chr, std::vector<int> &Aorder, const std::string &Amodel)
{
	if (chr == "X")
	{
		if (Amodel == "No solution")
		{
			isHI = isCGF = isGA = isMulti = false;
		}
		else if (Amodel == "HI")
		{
			isHI = true;
			isCGF = isGA = isMulti = false;
		}
		else if (Amodel == "GA")
		{
			isGA = true;
			isHI = isCGF = isMulti = false;
		}
		else if (Amodel.substr(0, 3) == "CGF")
		{
			isCGF = true;
			isHI = isGA = isMulti = false;
		}
		else
		{
			isMulti = true;
			isHI = isGA = isCGF = false;
		}
	}

	if(chr == "A" && isHI)
	{
#ifdef DEBUG
		std::cout << "HI start!" << std::endl;
#endif //DEBUG
		bisectionSearch(hiT, hiLlk, &models::HI);
#ifdef DEBUG
		std::cout << "HI Model:" << std::endl;
		std::cout << "\tT: " << hiT << std::endl;
		std::cout << "\tLlk: " << hiLlk << std::endl;
#endif //DEBUG
	}
	if(chr == "A" && isCGF)
	{
#ifdef DEBUG
		std::cout << "CGFD start!" << std::endl;
#endif //DEBUG
		bisectionSearch(cgfdT, cgfdLlk, &models::CGFD);
#ifdef DEBUG
		std::cout << "CGFD Model:" << std::endl;
		std::cout << "\tT: " << cgfdT << std::endl;
		std::cout << "\tLlk: " << cgfdLlk << std::endl;
		std::cout << "CGFR start!" << std::endl;
#endif //DEBUG
		bisectionSearch(cgfrT, cgfrLlk, &models::CGFR);
#ifdef DEBUG
		std::cout << "CGFR Model:" << std::endl;
		std::cout << "\tT: " << cgfrT << std::endl;
		std::cout << "\tLlk: " << cgfrLlk << std::endl;
#endif //DEBUG
	}
	if(chr == "A" && isGA)
	{
#ifdef DEBUG
		std::cout << "GA start!" << std::endl;
#endif //DEBUG
		bisectionSearch(gaT, gaLlk, &models::GA);
#ifdef DEBUG
		std::cout << "GA Model:" << std::endl;
		std::cout << "\tT: " << gaT << std::endl;
		std::cout << "\tLlk: " << gaLlk << std::endl;
#endif //DEBUG
	}
	if(isMulti)
	{
#ifdef DEBUG
		std::cout << "MultiWaveInfer start!" << std::endl;
#endif //DEBUG
		multiWaveInfer(chr, Aorder);
#ifdef DEBUG
		std::cout << "multiScenarioCount: " << multiScenarioCount << std::endl;
		for(int i = 0 ; i < multiScenarioCount; ++i)
		{
			std::cout << "Scenario " << i + 1 << std::endl;
			int nwave = multiPop.size();
			for(int j = 0 ; j < nwave; ++j)
			{
				std::cout << "\tPop " << multiPop[j] + 1 << ":\t" << multiT[j] << \
						"\t" << multiProp[j] << std::endl;
			}
		}
		std::cout << "\tLlk: " << multiLlk << std::endl;
#endif //DEBUG
	}
}

inline double models::HI(int t) const
{
	double llk(0);
	unsigned int npop = segs.size();
	assert(npop == 2);

	for(unsigned int ipop = 0 ; ipop < npop; ++ipop)
	{
		const std::vector<double> & psegs = segs[ipop];
		const double lambda = (1 - prop[ipop]) * t;
		const double logLambda = log(lambda);
		int nsegs = psegs.size();
		double pllk(0);
		for(int i = 0 ; i < nsegs; ++i)
			pllk += psegs[i];
		pllk -= minLen * nsegs;
		pllk *= - lambda;
		pllk += logLambda * nsegs;
		llk += pllk;
	}

	return llk;
}

inline double models::CGFR(int t) const
{
	unsigned int npop = segs.size();
	assert(npop == 2);

	const std::vector<double> & donor = segs[1];
	const std::vector<double> & recipient = segs[0];
	const double m = prop[0];

	return CGF(donor, recipient, m, t);
}

inline double models::CGFD(int t) const
{
	unsigned int npop = segs.size();
	assert(npop == 2);

	const std::vector<double> & donor = segs[0];
	const std::vector<double> & recipient = segs[1];
	const double m = prop[1];

	return CGF(donor, recipient, m, t);
}

inline double models::CGF(const std::vector<double> & donor, \
		const std::vector<double> & recipient, double m, int t) const
{
	double llk(0);
	const double mpow_1_t = pow(m, 1.0 / t);
	const double lambda = (double)t - ( 1.0 - m ) * mpow_1_t/ ( 1 - mpow_1_t );
	const double logLambda = log(lambda);
	int nrsegs = recipient.size();
	for(int i = 0 ; i < nrsegs; ++i)
		llk += recipient[i];
	llk -= minLen * nrsegs;
	llk *= - lambda;
	llk += logLambda * nrsegs;

	const double mpow_tp1_t = pow(m, (1.0 + t) / t);
	double *tfactor, *expfactor, *tdefactor;

	tfactor = new double[t + 1];
	expfactor = new double[t + 1];
	tdefactor = new double[t + 1];

	for(int i = 1 ; i <= t ; ++i)
	{
		const double mpow_i_t = pow(m, (double) i / t);
		const double tmp = mpow_i_t - mpow_tp1_t;
		expfactor[i] = - tmp / (1.0 - mpow_1_t);
		tdefactor[i] = tmp / mpow_i_t;
		tfactor[i] = tdefactor[i] * tmp;
	}

	double factor(0);
	for(int i = 1 ; i <= t ; ++i)
		factor += tdefactor[i] * exp(expfactor[i] * minLen);
	factor *= 1.0 - mpow_1_t;
	factor = 1.0 / factor;
	const double logfactor = log(factor);

	int ndsegs = donor.size();
	double dllk(0);
	for(int i = 0 ; i < ndsegs; ++i)
	{
		double tlk(0);
		for(int j = 1; j <= t ; ++j)
			tlk += tfactor[j] * exp(expfactor[j] * donor[i]);
		dllk += log(tlk);
	}
	dllk += logfactor * ndsegs;
	llk += dllk;

	delete []tfactor;
	delete []expfactor;
	delete []tdefactor;

	return llk;
}

inline double models::GA(int t) const
{
	double llk(0);
	unsigned int npop = segs.size();
	assert(npop == 2);

	double *tfactor, *expfactor, *tdefactor;

	tfactor = new double[t + 1];
	expfactor = new double[t + 1];
	tdefactor = new double[t + 1];
	for(int i = 2 ; i <= t ; ++i)
	{
		double ret = t - i + 1;
		double tpowt = pow(1.0 - 1.0 / t, 1.0 - i);
		tdefactor[i] = ret / t * tpowt;
		tfactor[i] = tdefactor[i] * ret;
		expfactor[i] = -ret;
	}

	for(unsigned int ipop = 0 ; ipop < npop; ++ipop)
	{
		const std::vector<double> & psegs = segs[ipop];
		double m = 1 - prop[ipop];
		double factor(0);
		for(int i = 2 ; i <= t ; ++i)
			factor += tdefactor[i] * exp (expfactor[i] * m * minLen);
		factor += t * exp( - t * m * minLen);
		factor = m / factor;
		const double logfactor = log(factor);
		double pllk(0);
		int nsegs = psegs.size();
		for(int i = 0 ; i < nsegs; ++i)
		{
			double plk(0);
			for(int j = 2 ; j <= t ; ++j)
				plk += tfactor[j] * exp(expfactor[j] * m * psegs[i]);
			plk += t * t * exp(-t * m * psegs[i]);
			pllk += log(plk);
		}
		pllk += nsegs * logfactor;
		llk += pllk;
	}

	delete []tfactor;
	delete []expfactor;
	delete []tdefactor;

	return llk;
}  

inline double models::Multi(const std::vector<double> & psegs, double* lambda, \
		double* prop, int nwave) const
{
	double llk(0);
	int nseg = psegs.size();
	for(int i = 0 ; i < nseg; ++i)
	{
		double tlk(0);
		for(int j = 0 ; j < nwave ; ++j)
			tlk += prop[j] * lambda[j] * exp(- lambda[j] * psegs[i]);
		llk += log(tlk);
	}
	return llk;
}


double models::MultiEM(const std::vector<double> & psegs, const double ancProp, \
		double* & curLambda, double* & curProp, int & curNwave, \
		const std::string &chr, const int &popIndex, const int &npop)
{
	int k(1);
	if (chr == "X")
	{
		k = curNwave;
	}
	int nsegs = psegs.size();
	double llk(0);
	while(k <= MultiMaxWave)
	{
		double *iterLambda, *iterProp;
		iterLambda = new double[k];
		iterProp = new double[k];

		for(int i = 0 ; i < k ; ++i)
		{
			iterLambda[i] = (double) rand() / RAND_MAX;
			iterProp[i] = 1.0 / k;
		}

		int nit(0);
		double *iterNLambda, *iterNProp;
		iterNLambda = new double[k];
		iterNProp = new double[k];

		double **pval, **pxval;
		pval = new double *[k];
		pxval = new double *[k];
		for(int i = 0 ; i < k; ++i)
		{
			pval[i] = new double[nsegs];
			pxval[i] = new double[nsegs];
		}
		while(++nit < MaxInter)
		{
			//E-step
			for(int i = 0 ; i < nsegs; ++i)
			{
				double denon(0);
				for(int j = 0 ; j < k ; ++j)
				{
					pval[j][i] = iterProp[j] * iterLambda[j] * exp(- iterLambda[j] * psegs[i]);
					denon += pval[j][i];
				}
				for(int j = 0 ; j < k ; ++j)
				{
					pval[j][i] /= denon;
					pxval[j][i] = pval[j][i] * psegs[i];
				}
			}
			//M-step
			for(int i = 0 ; i < k ; ++i)
			{
				double sump(0), sumpx(0);
				for(int j = 0 ; j < nsegs; ++j)
				{
					sump += pval[i][j];
					sumpx += pxval[i][j];
				}
				iterNProp[i] = sump / nsegs;		
				iterNLambda[i] = sump / sumpx;
			}

			//check converge
			bool converge(true);
			for(int i = 0 ; i < k ; ++i)
			{
				if( (fabs(iterLambda[i] - iterNLambda[i]) > epsilon) || \
						(fabs(iterProp[i] - iterNProp[i]) > epsilon) )
					converge = false;
				iterLambda[i] = iterNLambda[i];
				iterProp[i] = iterNProp[i];
			}
			if(converge)
				break;
		}

		double tllk = Multi(psegs, iterLambda, iterProp, k);

		if (chr == "X")
		{
			llk = tllk;
			curLambda = iterLambda;
			curProp = iterProp;
			k = MultiMaxWave+1;
		}
		else
		{
			if(k == 1)
			{
				llk = tllk;
				curLambda = iterLambda;
				curProp = iterProp;
				curNwave = k;
			}
			else
			{
				double findopt = false;
				if( 2 * (tllk - llk) < critical)
					findopt =true;
				else
				{
					double tempSum(0);
					std::vector<double> temp;
					temp.resize(k, 0);
					for(int j = 0 ; j < k ; ++j)
					{
						temp[j] = iterProp[j] / iterLambda[j];
						tempSum += temp[j] / temp[0];
					}
					tempSum = ancProp / tempSum;
					for(int j = 0 ; j < k ; ++j)
					{
						if(tempSum * temp[j] / temp[0] < minProp)
						{
							findopt = true;
							break;
						}
					}
				}
				if(!findopt)
				{
					delete[] curLambda;
					delete[] curProp;
					llk = tllk;
					curLambda = iterLambda;
					curProp = iterProp;
					curNwave = k;
				}
				else
				{
					delete[] iterLambda;
					delete[] iterProp;
				}
			}
			++k;
		}
		delete[] iterNLambda;
		delete[] iterNProp;
		for(int i = 0 ; i < curNwave; ++i)
		{
			delete[] pval[i];
			delete[] pxval[i];
		}
		delete[] pval;
		delete[] pxval;
	}

	std::map<double, double> SortPar;
	for(int j = 0 ; j < curNwave; ++j)
		SortPar[curLambda[j]] = curProp[j];
	int p (0);
	for(std::map<double, double>::iterator it = SortPar.begin(); it != SortPar.end(); ++it)
	{
		curLambda[p] = it->first;
		curProp[p] = it->second;
		++p;
	}

	if (chr == "X")
	{
		EMprop.resize(npop);
		EMprop[popIndex].resize(nsegs);
		for (int i = 0; i < nsegs; i++)
			EMprop[popIndex][i].resize(curNwave);
		for (int i = 0; i < nsegs; i++)
		{
			double sumE = 0;
			for (int j = 0; j < curNwave; j++)
				sumE += curProp[j]*curLambda[j]*exp(-curLambda[j]*psegs[i]);

			for (int j = 0; j < curNwave; j++)
				EMprop[popIndex][i][j] =  curProp[j]*curLambda[j]*exp(-curLambda[j]*psegs[i])/sumE;
		}
	}

	std::vector<double> temp;
	temp.resize(curNwave, 0);
	temp[0] = 1.0;
	double tempSum = temp[0];
	for(int j = 1; j < curNwave ; ++j)
	{
		temp[j] = curProp[j] * exp(( curLambda[j] - curLambda[0]) * minLen) / curProp[0];
		tempSum += temp[j];
	}
	curProp[0] = 1.0 / tempSum;
	for(int j = 1 ; j < curNwave; ++j)
		curProp[j] = curProp[0] * temp[j];

	tempSum = 0;
	for(int j = 0 ; j < curNwave; ++j)
	{
		temp[j] = curProp[j] / curLambda[j];
		tempSum += temp[j] / temp[0];
	}
	tempSum = ancProp / tempSum;
	for(int j = 0 ; j < curNwave; ++j)
		curProp[j] = tempSum * temp[j] / temp[0];

	return llk;
}

void models::bisectionSearch(int & t, double & llk, double (models::*getLlk)(int) const)
{
	int leftBound = 1;
	int rightBound = 2;
	double llk1 = 0;
	double llk2 = 0;
	bool findBound = false;
	while (!findBound)
	{
		llk1 = (this->*getLlk)(rightBound);
		llk2 = (this->*getLlk)(rightBound - 1);
		if (llk1 > llk2)
		{
			leftBound = rightBound;
			rightBound *= 2;
		}
		else
		{
			rightBound = rightBound - 1;
			findBound = true;
		}
	}

	int left = leftBound;
	int right = rightBound;
	bool isOpt = false;
	while (!isOpt)
	{
		//if only two points lefted, just compare it
		if (right - left < 2)

		{
			llk1 = (this->*getLlk)(right);
			llk2 = (this->*getLlk)(left);
			if (llk1 > llk2)
			{
				llk = llk1;
				t = right;
			}
			else
			{
				llk = llk2;
				t = left;
			}
			isOpt = true;
		}
		else
		{
			int middle = (right + left) / 2;
			llk1 = (this->*getLlk)(middle);
			llk2 = (this->*getLlk)(middle - 1);
			if (llk1 < llk2)
			{
				right = middle - 1;
			}
			else
			{
				left = middle;
			}
		}
	}
}

void models::multiWaveInfer(const std::string &chr, std::vector<int> &Aorder)
{
	int npop = segs.size();
//	assert(npop <= MultiMaxPop && npop >= 2);
	assert(npop >= 2);

	double **multilambda, **multiprop;
	multilambda = new double* [npop];
	multiprop = new double* [npop];

	int *nwave;
	nwave = new int[npop];

	if (chr == "X")
	{
		for (int i = 0; i < npop; i++)
		{
			nwave[i] = count(Aorder.begin(), Aorder.end(), i);
		}
	}

	multiLlk = 0;
	int totalWaves(0);
	std::vector<int> popOrder;
	for(int i = 0; i < npop; ++i)
	{
		std::vector<double> csegs;
		int pnsegs = segs[i].size();
		const std::vector<double> & curSegs = segs[i];
		csegs.resize(pnsegs);
		for(int j = 0 ; j < pnsegs; ++j)
			csegs[j] = curSegs[j] - minLen;
		multiLlk += MultiEM(csegs, prop[i], multilambda[i], multiprop[i], nwave[i], chr, i, npop);
#ifdef DEBUG
		std::cout << "Pop" << i + 1 << ": waves:\t" << nwave[i] << std::endl;
		for(int j = 0 ; j < nwave[i]; ++j)
			std::cout << "\twave " << j + 1 << '\t' << multilambda[i][j] << '\t' \
			<< multiprop[i][j] << std::endl;
#endif //DEBUG
		totalWaves += nwave[i];
		for(int j = 0 ; j < nwave[i]; ++j)
			popOrder.push_back(i);
	}

	std::vector<std::vector<int> > allOrder;
	if (chr == "A")
	{
		allOrder = perm(popOrder);
	}
	else
	{
		allOrder.resize(1);
		allOrder[0] = Aorder;
	}
	for (std::vector< std::vector<int> >::iterator iter = allOrder.begin(); \
			iter != allOrder.end(); ++iter)
	{
		/*for a given order, calculate alpha, H and time T */
		std::vector<double> mInOrder;
		mInOrder.resize(totalWaves);
		for (int i = 0; i < npop; ++i)
		{
			int index = 0;
			for (int j = 0; j < totalWaves; ++j)
			{
				if (iter->at(j) == i)
				{
					mInOrder[j] = multiprop[i][index];
					index++;
				}
			}
		}

		/* calculate alpha */
		std::vector<double> alphaInOrder;
		alphaInOrder.resize(totalWaves);
		double multiplier = 1.0;
		alphaInOrder[0] = mInOrder[0];
		for (int i = 1; i < totalWaves - 1; ++i)
		{
			multiplier *= (1 - alphaInOrder[i - 1]);
			alphaInOrder[i] = mInOrder[i] / multiplier;
		}
		//here make sure the proportions of first two populations sum to 1
		alphaInOrder[totalWaves - 1] = 1.0 - alphaInOrder[totalWaves - 2];

		/* calculate total ancestry proportion of kth ancestral population at t generation H_k(t)*/
		std::vector< std::vector<double> > hInOrder;
		hInOrder.resize(npop);
		for(int i = 0 ; i < npop; ++i)
			hInOrder[i].resize(totalWaves);
		/*
		 * Initial last element of H for each population.
		 * look at the last two positions of population order: if the second last
		 * position is from population k, then set the last value of H for population
		 * k as the alpha in order of the second last value; if the last position is
		 * from population k, then set the last value of H for population k as the alpha
		 * in order of the last value; otherwise, set the last value of H for population
		 * k as 0.0
		 */
		for (int i = 0; i < npop; ++i)
		{
			int lastIndexOfH = totalWaves - 2;
			for (int j = 0; j < npop; ++j)
			{
				if (iter->at(lastIndexOfH) == j)
					hInOrder[j][lastIndexOfH] = alphaInOrder[lastIndexOfH];
				else if (iter->at(lastIndexOfH + 1) == j)
					hInOrder[j][lastIndexOfH] = alphaInOrder[lastIndexOfH + 1];
				else
					hInOrder[j][lastIndexOfH] = 0;
			}
			/*
			 * Recursively update H
			 * if alpha in order from population k, then update H as H*(1-alpha)+alpha;
			 * else update H as H*(1-alpha)
			 */
			for (int j = totalWaves - 3; j >= 0; --j)
			{
				if (iter->at(j) == i)
					hInOrder[i][j] = hInOrder[i][j + 1] * \
							(1 - alphaInOrder[j]) + alphaInOrder[j];
				else
					hInOrder[i][j] = hInOrder[i][j + 1] * (1 - alphaInOrder[j]);
			}
		}

        for (int i = 0; i < npop; ++i)
            hInOrder[i][totalWaves - 1] = hInOrder[i][totalWaves - 2];

		std::vector<double> admixTime; //here is delta T between waves, in reverse order
		admixTime.resize(totalWaves);
		std::vector<int> indexes;
		indexes.resize(npop, 0);

		for (int i = 0; i < totalWaves; ++i)
		{
			for (int j = 0; j < npop; ++j)
			{
				if (iter->at(i) == j)
				{
					double tempSum = 0;
					for (int k = 0; k < i; ++k)
						tempSum += (1 - hInOrder[j][k]) * admixTime[k];

					if (chr == "X")
						tempSum *= 2.0/3;
					
					double rate = multilambda[j][indexes[j]] - tempSum;

					admixTime[i] = rate / (1 - hInOrder[j][i]);
					if (chr == "X")
						admixTime[i] *= 1.5;
					++indexes[j];
					break;
				}
			}
		}

		/*
         * Check whether the results is reasonable or not
         * Assume total N waves of admixture events, let T[i] denotes the admixture time of i+1 wave
         *  t[i] denotes the time difference between i and i+1 wave
         * then T[i] = sum_{j=0}^i (t[j]}
         * if the order is reasonable, we must make sure T[0] <= T[1] <= ... <= T[N-2] and T[N-3] <= T[N-1]
         * T[0] <= T[1] <= ... <= T[N-2] => t[i] >= 0 for i = 0, 1, 2, ..., N-2
         * and T[N-3] <= T[N-1] => T[N-3] <= T[N-3] + t[N-2] + t[N-1] => t[N-2] + t[N-1] >= 0
         */

		bool isReasonable = true;

		for (int i = 0; i < totalWaves - 1; ++i)
		{
			if (admixTime[i] < 0)
			{
				isReasonable = false;
				break;
			}
		}

        if (admixTime[totalWaves - 2] + admixTime[totalWaves - 1] < 0)
   	        isReasonable = false;

		for (int i = 1; i < totalWaves; ++i)
			admixTime[i] += admixTime[i - 1];

        /*
         * This part is to check whether the times of first admixture events are closer enough
         * In theory, these two times should be the same, in order to keep at least one possible result,
         * The filtering is toggered off now
         */


		if(isReasonable)
		{
			double meanTime = (admixTime[totalWaves - 1] + admixTime[totalWaves - 2]) / 2.0;
			admixTime[totalWaves - 1] = meanTime;
			admixTime[totalWaves - 2] = meanTime;

			//inter update llk

			std::vector<double> ht, uk, s;
			ht.resize(npop, 0);
			uk.resize(totalWaves, 0);
			s.resize(totalWaves, 0);
			ht[iter->at(totalWaves - 1)] = 1;
			s.back() = 1;
			for(int i = totalWaves - 2; i >= 0; --i)
			{
				int pop = iter->at(i);
				for(int j = 0 ; j < npop; ++j)
				{
					ht[j] *= 1 - alphaInOrder[i];
					if(pop == j)
						ht[j] += alphaInOrder[i];
				}
				for(int j = i ; j < totalWaves; ++j)
				{
					if(i == 0)
						uk[j] += (1 - ht[iter->at(j)]) * ((int)(admixTime[i] + 0.5));
					else
						uk[j] += (1 - ht[iter->at(j)]) * ((int)(admixTime[i] + 0.5) \
								 - (int)(admixTime[i - 1] + 0.5));
				}
				if (chr == "X")
				{
					for (int j = i; j < totalWaves; ++j)
						uk[j] *= 2.0/3;
				}

				s[i] = alphaInOrder[i];
				for(int j = i + 1 ; j < totalWaves; ++j)
					s[j] *= 1 - alphaInOrder[i];
			}
//			for(int i = 0 ; i < totalWaves; ++i)
//				std::cout << uk[i] << '\t' << s[i] << '\t' << iter->at(i) << std::endl;
			std::vector<double> w, sumw;
			w.resize(totalWaves, 0);
			sumw.resize(npop, 0);

			for(int i = 0 ; i < totalWaves; ++i)
			{
				w[i] = s[i] * uk[i];
				sumw[iter->at(i)] += w[i];
			}
			for(int i = 0 ; i < totalWaves; ++i)
				w[i] /= sumw[iter->at(i)];

			std::vector<double> factor, kfactor;
			factor.resize(totalWaves, 0);
			kfactor.resize(npop, 0);
			for(int i = 0 ; i < totalWaves; ++i)
			{
				double t = w[i] * exp( - uk[i] * minLen);
				factor[i] = t;
				kfactor[iter->at(i)] += t;
			}
			for(int i = 0 ; i < totalWaves; ++i)
				factor[i] /= kfactor[iter->at(i)];

			double tllk(0);
			for(int i = 0 ; i < npop; ++i)
			{
				const std::vector<double> & curSegs = segs[i];
				int curNseg = curSegs.size();
				for(int j = 0 ; j < curNseg; ++j)
				{
					double tlk(0);
					for(int k = 0 ; k < totalWaves; ++k)
					{
						if(iter->at(k) == i)
						{
							double curUk = uk[k];
							tlk += exp( - curUk * (curSegs[j] - minLen)) * \
									factor[k] * curUk;
						}
					}
					tllk += log(tlk);
				}
			}

			if(multiScenarioCount)
			{
				if(tllk < intMultiLlk)
					continue;
			}	

			intMultiLlk = tllk;
			if (chr == "A")
			{
				Aorder = *iter;	
			}

			std::vector<event> event4sort;
			for(int i = totalWaves - 1; i >= 0; --i)
			{
				event t;
				t.pop = iter->at(i);
				t.t = admixTime[i];
				t.alpha = alphaInOrder[i];
				event4sort.push_back(t);
			}
			sort(event4sort.begin(), event4sort.end(), cmp);

			for(int i = 0; i < totalWaves; ++i)
			{
				multiT.push_back(event4sort[i].t + 0.5);
				multiProp.push_back(event4sort[i].alpha);
				multiPop.push_back(event4sort[i].pop);
			}
			++multiScenarioCount;
		}
	}

	for(int i = 0 ; i < npop; ++i)
	{
		delete[] multilambda[i];
		delete[] multiprop[i];
	}
	delete[] multilambda;
	delete[] multiprop;
	delete[] nwave;
}

void models::solveMultiX(const std::vector<int> &Aorder, const std::vector<double> &AT, const std::vector<int> &APop)
{
	double sum = 0;
	for (int i = 0; i < segProp.size(); i++)
	{
		for (int j = 0; j < segProp[i].size(); j++)
			sum += segProp[i][j];		
	}

	for (int i = 0; i < segProp.size(); i++)
	{
		for (int j = 0; j < segProp[i].size(); j++)
			segProp[i][j] /= sum;
	}

	std::vector<int> tempOrder;
	tempOrder.resize(Aorder.size());
	for (int i = 0; i < Aorder.size(); i++)
	{
		int n = count(Aorder.begin(), Aorder.begin()+i+1, Aorder[i]);
		tempOrder[i] = n-1;			
	}	

	int npop = segs.size();
	std::vector<std::vector<double>> tempMultiProp;
	tempMultiProp.resize(npop);

	double alpha = segProp[Aorder[0]][0];
	double multiplier = 1-alpha;
	tempMultiProp[Aorder[0]].push_back(alpha);

	for (int i = 1; i < Aorder.size()-1; i++)
	{		
		alpha = segProp[Aorder[i]][tempOrder[i]]/multiplier;
		tempMultiProp[Aorder[i]].push_back(alpha);
		multiplier *= (1-alpha);
	}	

	tempMultiProp[Aorder.back()].push_back(1-alpha);

	multiProp.resize(Aorder.size());
	int k = 0;
	for (int i = 0; i < tempMultiProp.size(); i++)
	{
		for (int j = 0; j < tempMultiProp[i].size(); j++)
		{
			multiProp[k] = tempMultiProp[i][j];
			k++;   
		}	
	}

	multiPop = APop;
	multiT = AT;
	multiScenarioCount = 1;
}

std::string models::info(const std::string &chr) const
{
	std::ostringstream str;
	if (chr == "A" or isMulti)
	{
		if(isHI)
			str << "HI, generation:\t" << hiT << "\t, log(likelihood):\t" \
				<< hiLlk << "\t; ";
		if(isCGF)
		{
			str << "CGFD, generation:\t" << cgfdT << "\t, log(likelihood):\t" \
				<< cgfdLlk << "\t; ";
			str << "CGFR, generation:\t" << cgfrT << "\t, log(likelihood):\t" \
					<< cgfrLlk << "\t; ";
		}
		if(isGA)
			str << "GA, generation:\t" << gaT << "\t, log(likelihood):\t" << gaLlk \
				<< "\t; ";
		if(isMulti)
		{
			if(multiScenarioCount)
			{
				int numWave = (multiPop.size() + 1) / multiScenarioCount - 1;
				str << "MultiWave, NumberWaves:\t" << numWave << "\t, log(likelihood):\t"<< \
					multiLlk << '\t' << intMultiLlk << "\t,\t";
				for(int n = 0 ; n < multiScenarioCount; ++n)
				{
					for(int i = (numWave + 1) * n ; i < (numWave + 1) * (n + 1) - 1; ++i)
					{
						str << "(" << popLabels[multiPop[i]] << ", generation:\t" << multiT[i] << \
							"\t" << multiProp[i] << "), ";
					}
					str << ";";
				}
			}
			else
				str << "MultiWave, No solution;";
		}
	}
	return str.str();
}

std::string models::summary(std::vector<int> & pop, std::vector<int> & t,\
	std::vector<double> & pro, const std::string &chr) const
{
#ifdef DEBUG
	std::cout << "Start summary" << std::endl;
	std::cout << prop.size() << std::endl;
	std::cout << "prop: " << prop[0] << std::endl;
#endif //DEBUG
	std::string modeldes("No solution");
	double maxllk;

	if (chr == "X" && (isHI || isGA || isCGF))
	{
		pop.resize(1);
		pop[0] = 0;
		pro.resize(1);
		pro[0] = prop[0];
		if (isHI)
			modeldes = "HI";
		if (isGA)
			modeldes = "GA";
		if (isCGF)
			modeldes = "CGF";
	}
	else
	{
		if(isHI)
		{
			maxllk = hiLlk;
			modeldes = "HI";
			pop.resize(1);
			pop[0] = 0;
			t.resize(1);
			t[0] = hiT;
			pro.resize(1);
			pro[0] = prop[0];

			if(maxllk < gaLlk)
			{
				maxllk = gaLlk;
				modeldes = "GA";
				t[0] = gaT;
			}
			if(maxllk < cgfdLlk)
			{
				maxllk = cgfdLlk;
				modeldes = "CGF (\"" + popLabels[0] + "\" as donor)";
				t[0] = cgfdT;
			}
			if(maxllk < cgfrLlk)
			{
				maxllk = cgfrLlk;
				modeldes = "CGF (\"" + popLabels[0] + "\" as recipient)";
				t[0] = cgfrT;
			}
		}
	}

//	std::cout << "HI: " << hiLlk << std::endl;	//
//	std::cout << "GA: " << gaLlk << std::endl;	//
//	std::cout << "CGFD: " << cgfdLlk << std::endl;	//	
//	std::cout << "CGFR: " << cgfrLlk << std::endl;	//
//	std::cout << "Multi: " << multiLlk << std::endl;	//

#ifdef DEBUG
	std::cout << "AI best\t" << modeldes << std::endl;
	std::cout << (multiPop.size() + 1) / multiScenarioCount - 1 << std::endl;
	std::cout << "multiT " << multiT.size() << std::endl;
	std::cout << "multiProp " << multiProp.size() << std::endl;
	std::cout << "multiPop " << multiPop.size() << std::endl;
#endif //DEBUG
	if(isMulti && multiScenarioCount)
	{
		int numWave = (multiPop.size() + 1) / multiScenarioCount - 1;
//		double aic = (numWave - 1.0) * 2 - multiLlk;
//		double aic = (numWave - 1.0) * 2 - intMultiLlk;
//		double paic = 2.0 - maxllk;
		int totalSegs(0);
		for (unsigned int i = 0 ; i < segs.size(); ++i)
			totalSegs += segs[i].size();
		double bic = log(totalSegs) * (numWave - 1.0) * 2 - 2 * intMultiLlk;
		double pbic = log(totalSegs) * 2 - 2 * maxllk;

		if(!isHI || bic < pbic)
		{
			maxllk = multiLlk;
			int npop = segs.size();
			if(npop == 2 && numWave == 2)
			{
				modeldes = "HI";
				if(!isHI)
				{
					t.resize(1);
					pop.resize(1);
					pro.resize(1);
				}
				t[0] = multiT[0];
				pop[0] = 0;
				pro[0] = prop[0];
			}
			else
			{
				std::ostringstream mddes;
				std::vector<int> waveCount;
				waveCount.resize(npop, 0);
				for(unsigned int i = 0 ; i < multiPop.size(); ++i)
					if(multiPop[i] >= 0)
						++ waveCount[multiPop[i]];
				mddes << waveCount[0];
				for(int i = 1 ; i < npop; ++i)
					mddes << "-" << waveCount[i];
				modeldes = "Multi " + mddes.str();
				t.clear();
				pop.clear();
				pro.clear();
				t.insert(t.begin(), multiT.begin(), multiT.end());
				pop.insert(pop.begin(), multiPop.begin(), multiPop.end());
				pro.insert(pro.begin(), multiProp.begin(), multiProp.end());
			}
		}
	}

	return modeldes;
}

void Summary(const models& md, const std::vector<models *>& bootstrap, \
		const std::vector<std::string>& popLabels, double ci, \
		std::ofstream & fpo, std::ofstream & fplog, \
		std::vector<double>& Aprop, std::string &Amodel, \
		const std::string &chr, std::vector<int> &AT)
{
	std::vector<int> mdPop, mdT;
	std::vector<double> mdProp;
	std::string mdMod;
	mdMod = md.summary(mdPop, mdT, mdProp, chr);
	
	std::vector<std::string> pDirec; 

	if (chr == "A")
	{
		Amodel = mdMod;
		AT = mdT;
	}

	if (chr == "X")
	{
		if (mdMod != "No solution")
		{
			mdMod = Amodel;
			mdT = AT;		//set X chromosomal Time as that of autosomes
		}
	}		

	std::vector<std::vector<int> > bootstrapPop, bootstrapT;
	std::vector<std::vector<double> > bootstrapProp;
	std::vector<std::string> bootstrapMd;
	int nbootstrap = bootstrap.size();

	if (nbootstrap)
	{
		bootstrapPop.resize(nbootstrap);
		bootstrapT.resize(nbootstrap);
		bootstrapProp.resize(nbootstrap);
		bootstrapMd.resize(nbootstrap);
	}
	for (int i = 0 ; i < nbootstrap; ++i)
	{
		bootstrapMd[i] = bootstrap[i]->summary(bootstrapPop[i], bootstrapT[i],\
				bootstrapProp[i], chr);	
	}

	if (chr == "X")
	{
		for (int i = 0 ; i < nbootstrap; ++i)
		{
			if (bootstrapMd[i] != "No solution")
			{
				bootstrapMd[i] = Amodel;
				bootstrapT[i] = AT;
			}
		}
	}

	std::map<std::string, std::vector<std::vector<double> > > sumProp;
	std::map<std::string, std::vector<std::vector<int> > > sumT, sumPop;
	std::map<std::string, int> addtion;


	for(int i = 0 ; i < nbootstrap; ++i)
	{
		std::vector<std::vector<double> > & curProp = sumProp[bootstrapMd[i]];
		std::vector<std::vector<int> > & curT = sumT[bootstrapMd[i]];
		std::vector<std::vector<int> > & curPop = sumPop[bootstrapMd[i]];
		
		if(bootstrapPop[i].size() == 1)
		{
			if(curProp.size() == 0)
			{
				curProp.resize(1);
				curT.resize(1);
				curPop.resize(1);
			}
			curProp[0].push_back(bootstrapProp[i][0]);
			curT[0].push_back(bootstrapT[i][0]);
			curPop[0].push_back(0);
		}
		else if(bootstrapPop[i].size())
		{
			int nwave;
			if(curProp.size() == 0)
			{
				nwave = 0;
				for(unsigned int j = 0 ; j < bootstrapPop[i].size(); ++j)
					if(bootstrapPop[i][j] != -1)
						++nwave;
				curProp.resize(nwave);
				curT.resize(nwave);
				curPop.resize(nwave);
			}
			else
				nwave = curProp.size();
			int nscenario = ( bootstrapPop[i].size() + 1 ) / ( nwave + 1 );
			if(nscenario > 1)
				++addtion[bootstrapMd[i]];
			for(int j = 0 ; j < nscenario; ++j)
			{
				for(int k = j * ( nwave + 1) ; k < ( nwave + 1)  * ( j + 1 ) - 1; ++k)
				{
					int index = k - j * ( nwave + 1 );
					curProp[index].push_back(bootstrapProp[i][k]);
					curT[index].push_back(bootstrapT[i][k]);
					curPop[index].push_back(bootstrapPop[i][k]);
				}
			}
		}
	}

	if (chr == "X")
	{
		fpo << "Best Model for X:\t";
	}
	else
	{
		fpo << "Best Model:\t";
	}
	fpo << mdMod << std::endl;
	if (nbootstrap && mdMod != "No solution")
	{
		if (sumPop.find(mdMod) != sumPop.end())
			fpo << "Bootstrapping support ratio: " << (double) sumPop[mdMod].front().size() \
					/ nbootstrap * 100 << "% (" << sumPop[mdMod].front().size() << "/" \
					<< nbootstrap << ")" << std::endl;
		else
			fpo << "Bootstrapping support ratio: " << "0% (0/" << nbootstrap << ")" << std::endl;
	}
	if (mdPop.size() == 1)
	{
		if (chr == "A")
		{
			Aprop.push_back(mdProp[0]);
			Aprop.push_back(1-mdProp[0]);
		}

		fpo << "\t" << popLabels[0] << ": " << mdT[0] << "(G) " << mdProp[0];
		if (chr == "X")
		{
			if (mdMod.find("CGF") != std::string::npos)		//CGF
			{
				if (mdMod.find("recipient") != std::string::npos)  //pop1 as recipient
				{
					double mA = Aprop.at(0);
					double mX = mdProp[0];
					double p = 2-1.5*(pow(mX,1.0/mdT[0])*1.0 / pow(mA,1.0/mdT[0]));	
					if (p < 0)
						p = 0;
					else if (p > 1)
						p = 1;
					fpo << "\tp = " << p;
					if (p < 0.5)
						pDirec.push_back("F");
					else
						pDirec.push_back("M");

				}
				else	//pop1 as donor
				{
					double mA = 1-Aprop.at(0);
					double mX = 1-mdProp[0];
					double p = 2-1.5*((1-pow(mX,1.0/mdT[0]))*1.0/(1-pow(mA,1.0/mdT[0])));
					if (p < 0)
						p = 0;
					else if (p > 1)
						p = 1;	
					fpo << "\tp = " << p;
					if (p < 0.5)
						pDirec.push_back("F");
					else
						pDirec.push_back("M");
				}
			}
			else
			{
				double alphaA = Aprop.at(0);
				double alphaX = mdProp[0];
				double p = 2-1.5*(alphaX/alphaA);
				if (p <= 0)
					p = 0;
				else if (p > 1)
					p = 1;
				fpo << "\tp = " << p;
				if (p < 0.5)
					pDirec.push_back("F");
				else
					pDirec.push_back("M");
			}	
		} 
		fpo << std::endl;
				
		fpo << "\t" << popLabels[1] << ": " << mdT[0] << "(G) " << 1 - mdProp[0];
		if (chr == "X")
		{
			if (mdMod.find("CGF") != std::string::npos)		//CGF
			{
				if (mdMod.find("recipient") != std::string::npos)  // pop1 as recipient
				{
					double mA = Aprop.at(0);
					double mX = mdProp[0];
					double p = 2-1.5*((1-pow(mX,1.0/mdT[0]))*1.0/(1-pow(mA,1.0/mdT[0])));
					if (p < 0)
						p = 0;
					else if (p > 1)
						p = 1;
					fpo << "\tp = " << p;
					if (p < 0.5)
						pDirec.push_back("F");
					else
						pDirec.push_back("M");	
				}
				else		// pop1 as donor
				{
					double mA = 1-Aprop.at(0);
					double mX = 1-mdProp[0];
					double p = 2-1.5*(pow(mX,1.0/mdT[0])*1.0 / pow(mA,1.0/mdT[0]));
					if (p < 0)
						p = 0;
					else if (p > 1)
						p = 1;
					fpo << "\tp = " << p;
					if (p < 0.5)
						pDirec.push_back("F");
					else
						pDirec.push_back("M");
				}
			}
			else
			{
				double alphaA = Aprop.at(1);
				double alphaX = 1-mdProp[0];
				double p = 2-1.5*(alphaX/alphaA);
				if (p < 0)
					p = 0;
				else if (p > 1)
					p = 1;
				fpo << "\tp = " << p;
				if (p < 0.5)
					pDirec.push_back("F");
				else
					pDirec.push_back("M");
			}
		}
 		fpo	<< std::endl;
	}
	else if(md.multiScenarioCount)
	{
	//	std::cout << "mdPop.size() = " << mdPop.size() << std::endl;	//
		int nscenario = md.multiScenarioCount;
	//	std::cout << "nscenario = " << nscenario << std::endl;	//
		int nwave = (mdPop.size() + 1) / nscenario - 1;
		for(int j = 0 ; j < nscenario; ++j)
		{
			if(nscenario > 1)
				fpo << "Possible Scenario " << j + 1 << std::endl;
			for(int k = j * ( nwave + 1) ; k < ( nwave + 1)  * ( j + 1 ) - 1; ++k)
			{
				fpo << "\t" << popLabels[mdPop[k]] << ": " << mdT[k] << "(G) " << mdProp[k];
				if (chr == "A")
				{
					Aprop.push_back(mdProp[k]);
				}
				else
				{
					double alphaA = Aprop.at(k);
					double alphaX = mdProp[k];
					double p = 2-1.5*(alphaX/alphaA);
					if (p < 0)
						p = 0;
					else if (p > 1)
						p = 1;
					fpo << "\tp = " << p;
					if (p < 0.5)
						pDirec.push_back("F");
					else
						pDirec.push_back("M");
				}
				fpo << std::endl;
			}
		}
	}

	fpo << "-----------------------------------------------------------" << std::endl;
	if (mdMod != "No solution")
		fpo << "Bootstrapping details" << std::endl << std::endl;
	for(std::map<std::string, std::vector<std::vector<int> > >::iterator \
			it = sumPop.begin(); it != sumPop.end(); ++it)
	{
		std::string m = it->first;
		std::vector<std::vector<int> > & curT = sumT[m];
		std::vector<std::vector<int> > & curPop = sumPop[m];
		std::vector<std::vector<double> > & curProp = sumProp[m];

		int nwave = curT.size();
		if (m == "No solution")
			continue;
		else if(nwave == 1)
		{
			sort(curT[0].begin(), curT[0].end());
			sort(curProp[0].begin(), curProp[0].end());
			int n = curT[0].size();
			int p = n * (1.0 - ci) / 2;
			int q = n - p;
			if(p < 0)
				p = 0;
			if(q < 0)
				q = 0;
			if(p >= n)
				p = n - 1;
			if(q >= n)
				q = n - 1;
			fpo << m << std::endl;
			fpo << "Bootstrapping supporting ratio: " << (double) n / nbootstrap * 100 \
					<< "% (" << n << "/" << nbootstrap << ")" << std::endl;
			fpo << "\t" << popLabels[0] << ": ";
			if (chr == "A")
				fpo << curT[0][p] << "~" << curT[0][q];
			else
				fpo << AT[0] << "~" << AT[0];
			fpo << "(G) " << curProp[0][p] << "~" << curProp[0][q];

			if (chr == "X")
			{
				if (m.find("CGF") != std::string::npos)		//CGF
				{
					if (m.find("recipient") != std::string::npos)	// pop1 as recipient
					{
						double mA = Aprop.at(0);
						double mXl = curProp[0][p];
						double mXr = curProp[0][q];
						double pl = 2-1.5*(pow(mXr,1.0/AT[0])*1.0 / pow(mA,1.0/AT[0]));
						double pr = 2-1.5*(pow(mXl,1.0/AT[0])*1.0 / pow(mA,1.0/AT[0]));
						if (pl < 0)
							pl = 0;
						else if (pl > 1)
							pl = 1;
						if (pr < 0)
							pr = 0;
						else if (pr > 1)
							pr = 1;
						fpo << "\tp: " << pl << "~" << pr;
						if (m == mdMod)
						{
							double N = direcRatio(mA, curProp[0], pDirec[0], AT[0], "CGF1");
							fpo << "\t" << N/n*100 << "% (" << N << "/" << n << ")";
						}
					}
					else		// pop1 as donor	
					{
						double mA = Aprop.at(1);
						double mXl = 1-curProp[0][q];
						double mXr = 1-curProp[0][p];
						double pl = 2-1.5*((1-pow(mXl,1.0/AT[0]))*1.0 / (1-pow(mA,1.0/AT[0])));
						double pr = 2-1.5*((1-pow(mXr,1.0/AT[0]))*1.0 / (1-pow(mA,1.0/AT[0])));
						if (pl < 0)
							pl = 0;
						else if (pl > 1)
							pl = 1;
						if (pr < 0)
							pr = 0;
						else if (pr > 1)
							pr = 1;
						fpo << "\tp: " << pl << "~" << pr;
						if (m == mdMod)
						{
							double N = direcRatio(mA, curProp[0], pDirec[0], AT[0], "CGF2");
							fpo << "\t" << N/n*100 << "% (" << N << "/" << n << ")";
						}
					}
				}
				else
				{
					double alphaA = Aprop.at(0);
					double alphaXl = curProp[0][p];
					double alphaXr = curProp[0][q];
					double pl = 2-1.5*alphaXr/alphaA;
					double pr = 2-1.5*alphaXl/alphaA;
					if (pl < 0)
						pl = 0;
					else if (pl > 1)
						pl = 1;
					if (pr < 0)
						pr = 0;
					else if (pr > 1)
						pr = 1;
					fpo << "\tp: " << pl << "~" << pr;
					if (m == mdMod)
					{
						double N = direcRatio(alphaA, curProp[0], pDirec[0], 100, "HIGA1");
						fpo << "\t" << N/n*100 << "% (" << N << "/" << n << ")";
					}
				}
			}
			fpo << std::endl;
			fpo << "\t" << popLabels[1] << ": ";
			if (chr == "A")
 				fpo << curT[0][p] << "~" << curT[0][q];
			else
				fpo << AT[0] << "~" << AT[0];
			fpo << "(G) " << 1 - curProp[0][q] << "~" << 1 - curProp[0][p];
			if (chr == "X")
			{
				if (m.find("CGF") != std::string::npos)		//CGF
				{
					if (m.find("recipient") != std::string::npos)   // pop1 as recipient
					{
						double mA = Aprop.at(0);
						double mXl = curProp[0][p];
						double mXr = curProp[0][q];
						double pl = 2-1.5*((1-pow(mXl,1.0/AT[0]))*1.0 / (1-pow(mA,1.0/AT[0])));
						double pr = 2-1.5*((1-pow(mXr,1.0/AT[0]))*1.0 / (1-pow(mA,1.0/AT[0])));
						if (pl < 0)
							pl = 0;
						else if (pl > 1)
							pl = 1;
						if (pr < 0)
							pr = 0;
						else if (pr > 1)
							pr = 1;
						fpo << "\tp: " << pl << "~" << pr;
						if (m == mdMod)
						{
							double N = 0;
							if (pDirec[0] == "M")
								N = direcRatio(mA, curProp[0], "F", AT[0], "CGF3");
							else
								N = direcRatio(mA, curProp[0], "M", AT[0], "CGF3");
							fpo << "\t" << N/n*100 << "% (" << N << "/" << n << ")";	
						}
					}
					else
					{
						double mA = Aprop.at(1);
						double mXl = 1-curProp[0][q];
						double mXr = 1-curProp[0][p];
						double pl = 2-1.5*(pow(mXr,1.0/AT[0])*1.0/pow(mA,1.0/AT[0]));
						double pr = 2-1.5*(pow(mXl,1.0/AT[0])*1.0/pow(mA,1.0/AT[0]));
						if (pl < 0)
							pl = 0;
						else if (pl > 1)
							pl = 1;
						if (pr < 0)
							pr = 0;
						else if (pr > 1)
							pr = 1;
						fpo << "\tp: " << pl << "~" << pr;
						if (m == mdMod)
						{	
							double N = 0;
							if (pDirec[0] == "M")
								N = direcRatio(mA, curProp[0], "F", AT[0], "CGF4");
							else
								N = direcRatio(mA, curProp[0], "M", AT[0], "CGF4");
							fpo << "\t" << N/n*100 << "% (" << N << "/" << n << ")";	
						}
					}
				}	
				else
				{
					double alphaA = Aprop.at(1);
					double alphaXl = 1-curProp[0][q];
					double alphaXr = 1-curProp[0][p];
					double pl = 2-1.5*alphaXr/alphaA;
					double pr = 2-1.5*alphaXl/alphaA;
					if (pl < 0)
						pl = 0;
					else if (pl > 1)
						pl = 1;
					if (pr < 0)
						pr = 0;
					else if (pr > 1)
						pr = 1;
					fpo << "\tp: " << pl << "~" << pr;
					if (m == mdMod)
					{
						double N = 0;
						if (pDirec[0] == "M")
							N = direcRatio(alphaA, curProp[0], "F", 100, "HIGA2");
						else
							N = direcRatio(alphaA, curProp[0], "M", 100, "HIGA2");
						fpo << "\t" << N/n*100 << "% (" << N << "/" << n << ")";
					}
				}
			}
			fpo << std::endl;
		}
		else
		{
			for(int i = 0 ; i < nwave; ++i)
			{
				sort(curT[i].begin(), curT[i].end());
				sort(curProp[i].begin(), curProp[i].end());
			}
			int n = curT[0].size();
			int p = n * (1.0 - ci) / 2;
			int q = n - p;
			if(p < 0)
				p = 0;
			if(q < 0)
				q = 0;
			if(p >= n)
				p = n - 1;
			if(q >= n)
				q = n - 1;
			fpo << m << std::endl;
			fpo << "Bootstrapping supporting ratio: " << (double) n / nbootstrap * 100 \
					<< "% (" << n << "/" << nbootstrap << ")" << std::endl;
			for(int i = 0 ; i < nwave; ++i)
			{
				fpo << "\t" << popLabels[curPop[i].front()] << ": ";
			 	if (chr == "A")
					fpo << curT[i][p] << "~" << curT[i][q];
				else
					fpo << AT[i] << "~" << AT[i];
				fpo << "(G) " << curProp[i][p] << "~" << curProp[i][q];
				if (chr == "X")
				{
					double alphaA = Aprop.at(i);
					double alphaXl = curProp[i][p];
					double alphaXr = curProp[i][q];
					double pl = 2-1.5*alphaXr/alphaA;
					double pr = 2-1.5*alphaXl/alphaA;
					if (pl < 0)
						pl = 0;
					else if (pl > 1)
						pl = 1;
					if (pr < 0)
						pr = 0;
					else if (pr > 1)
						pr = 1;
					fpo << "\tp: " << pl << "~" << pr;
					if (m == mdMod)
					{
						double N = direcRatio(alphaA, curProp[i], pDirec[i], 100, "Multi");
						fpo << "\t" << N/n*100 << "% (" << N << "/" << n << ")";	
					}							
				}
				fpo << std::endl;
			}
		}
		fpo << "-----------------------------------------------------------" << std::endl;
	}
}

bool cmp(const event & a, const event & b)
{
	if(a.pop != b.pop)
		return a.pop < b.pop;
	else if(a.t != b.t)
		return a.t < b.t;
	else
		return a.alpha < b.alpha;
}

double direcRatio(const double &mA, const std::vector<double> &mX, const std::string &sex, const int &T, const std::string &model)
{
	double N = 0;
	for (int i = 0; i < mX.size(); i++)
	{
		double p = 0;
		if (model == "CGF1")
			p = 2-1.5*(pow(mX[i], 1.0/T)*1.0 / pow(mA, 1.0/T));
		else if (model == "CGF2")
		{
			double Xprop = 1-mX[i];
			p = 2-1.5*((1-pow(Xprop, 1.0/T))*1.0 / (1-pow(mA, 1.0/T)));
		}
		else if (model == "CGF3")
			p = 2-1.5*((1-pow(mX[i], 1.0/T))*1.0 / (1-pow(mA, 1.0/T)));
		else if (model == "CGF4")
		{
			double Xprop = 1-mX[i];
			p = 2-1.5*(pow(Xprop, 1.0/T)*1.0/pow(mA, 1.0/T));
		}
		else if (model == "HIGA1" || model == "Multi")
			p = 2-1.5*mX[i]/mA;
		else if (model == "HIGA2")
			p = 2-1.5*(1-mX[i])/mA;

		if (p < 0.5 && sex == "F")
			N += 1;
		else if (p > 0.5 && sex == "M")
			N += 1;
	}
	return N;
}

