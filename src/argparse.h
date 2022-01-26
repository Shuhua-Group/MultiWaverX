/*
 * argparse.h
 *
 *  Created on: June 8, 2021
 *      Author: yuankai & zhangrui
 */

#ifndef ARGPARSE_H_
#define ARGPARSE_H_

# include <iostream>

void print_help();

void argparse(int argc, char **argv);

class par
{
public:
	par() : inpath(""), outpath(""), logpath(""), \
			minLen(0), ci(0.95), nbootstrap(0), nthread(1), chr("A"), inpathX(""){};
	std::string inpath, outpath, logpath;
	double minLen, ci;
	int nbootstrap, nthread;
	std::string chr, inpathX;
};

extern par SoftPar;

#endif /* ARGPARSE_H_ */
