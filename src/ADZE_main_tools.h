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


/*
 * Where each locus sits in the genome.  One table for the whole run rather
 * than one per grouping: the missing-data filter condemns a locus for every
 * grouping at once, so the surviving locus set -- and its order -- is the same
 * everywhere, and a single table stays aligned with all of them.
 *
 * Coordinates come from CHROM and POS for VCF input and from --loci-map for
 * the STRUCTURE layout, which carries no positions of its own.  Empty when
 * neither supplied them; only the windowed statistics need them.
 */
struct LocusMap
{
  vector<string> chromName;   //chromosomes, in order of first appearance
  vector<int> chrom;          //per locus: index into chromName
  vector<long long> pos;      //per locus: position

  bool empty() const { return chrom.empty(); }
  size_t size() const { return chrom.size(); }

  int chromIndex(const string& name);          //find or append
  void compact(const vector<char>& del);       //apply a keep-mask
  bool ascending(int& badLocus) const;         //positions increase per chrom
};

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
		bool pp,LocusMap& lmap);
Population* readDataset(ParamSet& p,vector<string>& groupNames,int& numDivs,
			LocusMap& lmap);
bool wantsVCF(const string& format,const string& path);
list<int> parseKVals(string str);




