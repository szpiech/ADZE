#include <fstream>
#include <map>
#include <unordered_map>
#include <gsl/gsl_combination.h>
#include <time.h>
#include "ADZE_binom.h"
#include "ADZE_stats.h"
#include "ADZE_pfile.h"
#include "ADZE_pop.h"
#include "ADZE_pbar.h"


typedef map<int,string> IntStrMap;
const int BARLEN[] = {50,100,500,500};

int min(int,int);
bool sameStr(string, string);
int seenLabel(vector<string>&, string);
int seenBefore(vector<int>&,int);

//bool openFile(ifstream&, string);
void getDivLines(ifstream&,vector<string>&,vector<int>&,const ParamSet&);
void getLociNames(ifstream&, string[], const ParamSet&);
void readData(/*ifstream&,*/ Population[],vector<string>&,int,const ParamSet&);


int calcG(vector<int>[],int);
void calcNji(Population[],int,const string& missing);
void calcNj(Population[],int);
//void storeG(Population[],int,int,int);
void calcAg(Population[],int,int);
void calcPjig(Population pop[],vector<int> Nji[],
	      int Nj[],int locus,int numDivs);
double calcPg(Population[],int,int,int,int);
void calcAllAgs(Population[],int,const ParamSet&,bool,string);
void calcAllPgs(Population[],int,const ParamSet&,bool,string);
int countCols(string);
void checkDatafile(const ParamSet&);


bool isIn(int, gsl_combination*);

string nameCreate(string,string);
double displayTime(ostream& out);
void calcAllPgComb(Population pop[], int numDivs, int k, const ParamSet& p,
		   bool full_comb, string comb_out);
string combineNames(string names[],int k);
bool validK(int numDivs,list<int> k);
void filterLoci(Population pop[],int numDivs,double tol,string private_out,
		string missing, bool pp);
list<int> parseKVals(string str);
void warning(int loci,int numDivs,list<int> k,Population pop[]);



void estimate(ofstream& out, int J, int k, Population pop[],const ParamSet &p, 
	      double time);

