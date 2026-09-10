#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <chrono>
#include <iomanip>
#include "ADZE_binom.h"
#include "ADZE_stats.h"
#include "ADZE_pfile.h"
#ifdef ADZE_HAVE_ZLIB
#include <zlib.h>
#endif

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
void buildKTuples(int numDivs,int k,vector< vector<int> >& out);
bool readTupleFile(const string& file,Population pop[],int numDivs,
		   vector< vector<int> >& out);
void calcPgTuples(Population pop[],int numDivs,
		  const vector< vector<int> >& tuples,const ParamSet& p,
		  bool full_comb,string comb_out,bool namedTuples);
extern bool ADZE_QUIET;
ostream& adzelog();
bool stderrIsTerminal();
string combineNames(string names[],int k,char sep = ' ');
bool validK(int numDivs,list<int> k);
void filterLoci(Population pop[],int numDivs,double tol,string private_out,
		bool pp);
Population* readDataset(ParamSet& p,vector<string>& groupNames,int& numDivs);
bool wantsVCF(const string& format,const string& path);
list<int> parseKVals(string str);




