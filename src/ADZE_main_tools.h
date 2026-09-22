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

  /*
   * Windows are runs of consecutive loci, so the loci must arrive grouped by
   * chromosome and in increasing position.  Returns an empty string when they
   * do and a description of the first violation when they do not.
   */
  string checkOrder(const vector<string>& locusName) const;
};

/*
 * One window: a closed basepair interval on one chromosome, and the half-open
 * range [first,last) of surviving loci inside it.
 */
struct Window
{
  int chrom;
  long long start;
  long long end;
  int first;
  int last;

  int numLoci() const { return last - first; }
};

/*
 * Lay out the windows over the surviving loci.  Returns the number of windows
 * that held at least one locus but fewer than --min-window-loci, which the dry
 * run reports.
 */
long long buildWindows(const LocusMap& lmap, const ParamSet& p,
		       vector<Window>& out);

typedef map<int,string> IntStrMap;
const int BARLEN[] = {50,100,500,500};

int min(int,int);



void buildQTable(Population pop[],int numDivs,int locus,int numAlleles,
		 int gMax,int gStride,vector<double>& q);
/*
 * What one pass over the input yields without keeping any allele.
 *
 * The streaming sweep needs three whole-dataset facts before it can compute
 * anything: which loci survive --tolerance, what the default MAX_G resolves
 * to, and where the windows fall. All three follow from how many gene copies
 * each grouping scored at each locus, which is countable without interning a
 * single allele label -- so the scan parses genotypes only far enough to tell
 * a call from a missing one.
 *
 * observed[g][l] is that count. The grouping's total gene copies minus it is
 * the missing count the filter compares against, and it is itself the Nj the
 * ceilings are taken from.
 */
struct ScanResult
{
  vector<string> groupName;           //groupings, in output order
  vector<long long> groupRows;        //gene copies per grouping
  vector< vector<int> > observed;     //[grouping][locus]: non-missing copies
  vector<string> locusName;           //only when something needs names
  LocusMap lmap;                      //coordinates, when the input has them
  long long numLoci;
  long long geneCopies;

  /*
   * What resolve() derives, once every gene copy has been counted. The
   * filter's denominator is a grouping's total, which is only final at the
   * end of the pass, so the survival decision cannot be taken locus by locus
   * while reading -- that is the one thing this pass has to hold.
   */
  vector<char> dropped;               //per locus: filtered out by --tolerance
  long long survivors;
  vector<int> ceiling;                //per grouping: fewest copies at a survivor
  vector<int> binding;                //per grouping: loci sitting at that floor
  vector<int> emptyAt;                //per grouping: survivors it did not score
  int feasibleG;                      //the ceiling over all groupings
  int feasibleGAll;                   //the same before any locus is dropped

  ScanResult() : numLoci(0), geneCopies(0), survivors(0),
		 feasibleG(0), feasibleGAll(0) {}

  int nj(int g,long long l) const {return observed[g][size_t(l)];};
  long long missing(int g,long long l) const
  {return groupRows[g] - observed[g][size_t(l)];};

  void resolve(double tol);
};

void scanDataset(ParamSet& p,ScanResult& out);

void dumpCounts(const char* path,Population pop[],int numDivs,int numLoci,
		const LocusMap& lmap);
void writeFullDataHeader(ostream& out,const string& labelCols,int gFrom,int gTo);
void writeFullDataRows(ostream& out,const string& label,const vector<double>& vals,
		       int lo,int hi,int gFrom,int gTo,int gStride,
		       const Population& pop);
void writeWindowHeader(ostream& out,const char* groupColumn);
void writeWindowStats(ostream& out,vector<double>& perLocus,int numLoci,
		      const vector<Window>& windows,const LocusMap& lmap,
		      const string& label,int gFirst,int gLast,int gBase);
void calcAllAgs(Population[],int,const ParamSet&,bool,string,
		const vector<Window>&,const LocusMap&);
void calcAllPgs(Population[],int,const ParamSet&,bool,string,
		const vector<Window>&,const LocusMap&);



string nameCreate(string,string);
double displayTime(ostream& out);
void buildKTuples(int numDivs,int k,vector< vector<int> >& out);
bool readTupleFile(const string& file,Population pop[],int numDivs,
		   vector< vector<int> >& out);
void calcPgTuples(Population pop[],int numDivs,
		  const vector< vector<int> >& tuples,const ParamSet& p,
		  bool full_comb,string comb_out,bool namedTuples,
		  const vector<Window>& windows,const LocusMap& lmap);
extern bool ADZE_QUIET;
ostream& adzelog();
bool stderrIsTerminal();
string combineNames(string names[],int k,char sep = ' ');
bool validK(int numDivs,list<int> k);
void filterLoci(Population pop[],int numDivs,double tol,string private_out,
		bool pp,LocusMap& lmap,LocusTable& loci);
Population* readDataset(ParamSet& p,vector<string>& groupNames,int& numDivs,
			LocusTable& loci,
			LocusMap& lmap);
bool wantsVCF(const string& format,const string& path);
list<int> parseKVals(string str);




