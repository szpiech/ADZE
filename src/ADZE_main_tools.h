#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <gsl/gsl_combination.h>
#include <chrono>
#include <iomanip>
#include "ADZE_binom.h"
#include "ADZE_stats.h"
#include "ADZE_pfile.h"
#include "ADZE_pop.h"
#include "ADZE_pbar.h"


typedef map<int,string> IntStrMap;
const int BARLEN[] = {50,100,500,500};

int min(int,int);



void buildQTable(Population pop[],int numDivs,int locus,int numAlleles,
		 int gMax,int gStride,vector<double>& q);
void calcAllAgs(Population[],int,const ParamSet&,bool,string);
void calcAllPgs(Population[],int,const ParamSet&,bool,string);



string nameCreate(string,string);
double displayTime(ostream& out);
void calcAllPgComb(Population pop[], int numDivs, int k, const ParamSet& p,
		   bool full_comb, string comb_out);
string combineNames(string names[],int k);
bool validK(int numDivs,list<int> k);
void filterLoci(Population pop[],int numDivs,double tol,string private_out,
		bool pp);
Population* readDataset(ParamSet& p,vector<string>& groupNames,int& numDivs);
list<int> parseKVals(string str);




