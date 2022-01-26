/*
 * MultiWaveInferX.cpp
 *
 *  Created on: June 8, 2021
 *      Author: yuankai & zhangrui
 */

# include "models.h"
# include "argparse.h"
# include "exception.h"
# include "ancSegs.h"

# include <fstream>

# ifndef NOOMP
# include <omp.h>
# endif//NOOMP

using namespace std;

int main(int argc, char **argv)
{
	try
	{
		argparse(argc, argv);
		ofstream fpo(SoftPar.outpath.c_str());
		ofstream fplog;
		if(!fpo)
		{
			except e;
			e << "Error: Cannot open the output file!\n";
			e << "    The input path is \"" << SoftPar.outpath << '\"';
			throw e;
		}
		if(SoftPar.logpath != "")
		{
			fplog.open(SoftPar.logpath.c_str());
			if(!fplog)
			{
				except e;
				e << "Error: Cannot open the log file!\n";
				e << "    The input path is \"" << SoftPar.logpath << '\"';
				throw e;
			}
		}
		ancSegs dat(SoftPar.inpath, SoftPar.minLen);
		models m(dat);
		std::vector<double> Aprop;
		std::string Amodel;
		std::string Achr = "A";
		std::vector<int> Aorder;
		std::vector<int> AT;
		m.solve(Achr, Aorder, Amodel);
		std::cout << "Analysis on autosomes finished!\n" << std::endl;
		std::vector<models*> bootstrap;
		std::vector<ancBootstrap> bot;
		bootstrap.resize(SoftPar.nbootstrap, NULL);
		bot.resize(SoftPar.nbootstrap);
		srand(time(0));
		for(int i = 0 ; i < SoftPar.nbootstrap; ++i)
		{
			bot[i].bootstrap(dat);
			bootstrap[i] = new models(bot[i]);
		}
#ifndef NOOMP
		omp_set_num_threads(SoftPar.nthread);

#pragma omp parallel for
#endif//NOOMP
		for(int i = 0 ; i < SoftPar.nbootstrap; ++i)
		{
			std::vector<int> bootstrapOrder;
			bootstrap[i]->solve(Achr, bootstrapOrder, Amodel);
		}
		if(fplog)
		{
			fplog << "Results:\t" << m.info(Achr) << std::endl;
			for(int i = 0 ; i < SoftPar.nbootstrap; ++i)
				fplog << "Bootstrapping " << i + 1 << "\t" << \
						bootstrap[i]->info(Achr) << std::endl;
		}
		Summary(m, bootstrap, dat.popLabels, SoftPar.ci, fpo, fplog, Aprop, Amodel, Achr, AT);
		std::cout << "Bootstrapping on autosomal data finished!\n" << std::endl;

		if (SoftPar.chr == "X")
		{
			fplog << std::endl;
			fpo << std::endl;
			ancSegs datX(SoftPar.inpathX, SoftPar.minLen);
			models mX(datX);
			std::string Xchr = "X";
			mX.solve(Xchr, Aorder, Amodel);
			std::cout << "Analysis on X chromosomes finished!\n" << std::endl;
			datX.EMprop = mX.EMprop;
			std::vector<models*> bootstrapX;
			std::vector<ancBootstrap> botX;
	
			if (mX.info(Xchr) == "MultiWave, No solution;")
				SoftPar.nbootstrap = 0;

			bootstrapX.resize(SoftPar.nbootstrap, NULL);
			botX.resize(SoftPar.nbootstrap);
			srand(time(0));

			for (int i = 0; i < SoftPar.nbootstrap; ++i)
			{
				if (Amodel.find("Multi") != std::string::npos)
				{
					std::vector<int> Nwave;
					for (int j = 0; j < Amodel.size(); j++)
					{
						if (Amodel[j] == '-')
							Nwave.push_back(Amodel[j-1]-'0');
					}			
					Nwave.push_back(Amodel.back()-'0');
					botX[i].bootstrapMultiX(datX, Nwave);
				}
				else
					botX[i].bootstrap(datX);
				bootstrapX[i] = new models(botX[i]);
			}

			for (int i = 0; i < SoftPar.nbootstrap; ++i)
			{
				if (Amodel.find("Multi") != std::string::npos)
					bootstrapX[i] -> solveMultiX(Aorder, mX.multiT, mX.multiPop);
				else
					bootstrapX[i] -> solve(Xchr, Aorder, Amodel);
			}
			if (fplog)
			{
				fplog << "Bootstrap results of X chromosome" << std::endl;
				if (mX.info(Xchr) == "")
					fplog << "Under HI/GA/CGF model, the estimation time of X chromosome is identical to that of Autosomes" << std::endl;
				else
				{
					fplog << "Results:\t" << mX.info(Xchr) << std::endl;
					for (int i = 0; i < SoftPar.nbootstrap; ++i)
						fplog << "Bootstrapping " << i + 1 << "\t" << \
								bootstrapX[i]->info(Xchr) << std::endl;
				}
			}
			Summary(mX, bootstrapX, datX.popLabels, SoftPar.ci, fpo, fplog, Aprop, Amodel, Xchr, AT);	
			std::cout << "Bootstrapping on X chromosomal data finished!" << std::endl;
		}
	}
	catch(except & e)
	{
		e.err();
	}
}


