#include "ADZE_pop.h"
#include <algorithm>

using namespace std;

/*
 * Drop the flagged loci from the name table, recording them from the highest
 * index down -- the order 1.0's reverse-iterating delete loop wrote to the
 * _deletedloci file. Called once for the run, not once per grouping, because
 * the filter condemns a locus for every grouping at once.
 */
void LocusTable::compact(const vector<char>& del)
{
  for(int l = int(name.size()) - 1; l >= 0; l--)
    {
      if(del[l]) deleted.push_back(name[l]);
    }

  size_t keep = 0;
  for(size_t l = 0; l < name.size(); l++)
    {
      if(del[l]) continue;
      if(keep != l) name[keep].swap(name[l]);
      offset[keep+1] = offset[keep] + (offset[l+1] - offset[l]);
      keep++;
    }
  name.resize(keep);
  offset.resize(keep+1);

  return;
}

/*
 * One-line account of what the missing-data filter removed.  1.0 built this
 * inside printDeleted and wrote it to cout as well as to the file, with
 * different stream flags on each, so the console and the _deletedloci file
 * disagreed on how the percentage was formatted.
 */
string Population::deletedSummary() const
{
  ostringstream out;

  const size_t dropped = loci ? loci->deleted.size() : 0;

  out << dropped;
  out << ((dropped == 1) ? " locus has " : " loci have ");
  out << "at least one grouping with";
  if(tol > 0) out << " more than " << 100*tol << "%";
  out << " missing data.\n";

  return out.str();
}

void Population::printDeleted(ostream& out)
{
  if(!loci) return;

  for(vector<string>::const_iterator i = loci->deleted.begin();
      i != loci->deleted.end(); i++)
    {
      out << *i << endl;
    }

  return;
}


/*
 * Drop every locus flagged in del[] in a single pass.
 *
 * ADZE 1.0 deleted one locus at a time, shifting every higher-indexed locus
 * down one slot and copying all of its genotype strings, so filtering D of L
 * loci cost O(D*L*rows) std::string assignments -- quadratic in the locus
 * count.  Compacting once is O(L + slots).  The counts being one block, the
 * kept loci move down within it and the block then shrinks, which returns the
 * dropped loci's memory here rather than at exit.
 *
 * The names and the new offsets are compacted once for the run, in
 * filterLoci; oldOffset is where each locus used to start.
 */
void Population::deleteLoci(const vector<char>& del, const vector<int>& oldOffset)
{
  int keep = 0;
  int at = 0; //write position in the compacted count block

  for(int l = 0; l < numLoci; l++)
    {
      if(del[l]) continue;

      const int from = oldOffset[l];
      const int slots = oldOffset[l+1] - oldOffset[l];

      //The kept prefix is never longer than what it came from, so the write
      //position trails the read position and a forward copy is safe in place.
      if(at != from) copy(nji.begin()+from, nji.begin()+from+slots, nji.begin()+at);
      at += slots;

      if(keep != l)
	{
	  Nj[keep] = Nj[l];
	  missing[keep] = missing[l];
	}
      keep++;
    }

  nji.resize(at);
  numLoci = keep;
  return;
}

/*
 * Flag loci whose missing-data fraction exceeds the tolerance in this
 * grouping; del[] accumulates the union across groupings.  The missing counts
 * were tallied while the data file was read, so no genotype is revisited --
 * or even retained.
 */
void Population::recLociDelete(double tolerance, vector<char>& del)
{
  tol = tolerance;

  if(rows == 0) return; //no gene copies in this grouping

  for(int l = 0; l < numLoci; l++)
    {
      if(del[l]) continue; //already condemned by another grouping

      if(double(missing[l])/double(rows) > tolerance) del[l] = 1;
    }

  return;
}

//Returns the column length of the Nji matrix at the specified locus
int Population::getNjiColLength(int locus)
{
  if(locus > numLoci-1 || locus < 0)
    {
      //bad bounds
      return -9;
    }
  else
    {
      return (*loci).offset[locus+1] - (*loci).offset[locus];
    }
}

//Initializes Nj with specified number
void Population::fillNj(int n)
{
  for(int locus = 0; locus < numLoci; locus++)
    {
      Nj[locus] = n;
    }
  return;
}

//Returns the number of gene copies scored in this grouping at a given locus
int Population::getNj(int locus)
{
  if(locus > numLoci-1 || locus < 0)
    {
      //bad bounds
      return -9;
    }
  else
    {
      return Nj[locus];
    }
}

//Stores a number in Nj
bool Population::putNj(int num, int locus)
{
  if(locus > numLoci-1 || locus < 0)
    {
      //bad bounds
      return 0;
    }
  else
    {
      if(num < minG) minG = num;
      Nj[locus] = num;
      return 1;
    }
}

/* CALCULATE Nj
 *      m
 *      __
 * Nj = \  Nji
 *      /_
 *     i = 0
 */
void Population::sumNj()
{
  for(int locus = 0; locus < numLoci; locus++)
    {
      int total = 0;
      const int base = (*loci).offset[locus];
      const int slots = (*loci).offset[locus+1] - base;
      for(int i = 0; i < slots; i++)
	{
	  total += nji[base+i];
	}
      putNj(total,locus);
    }
  return;
}

//Returns the count of allele i in this grouping at a given locus
int Population::getNji(int row, int locus)
{
   if(locus > numLoci-1 || locus < 0)
    {
      //bad bounds
      return -9;
    }
  else if(row > getNjiColLength(locus)-1 || row < 0)
    {
      //bad bounds
      return -9;
    }
  else
    {
      return nji[(*loci).offset[locus] + row];
    }

}

//Stores a number in to the Nji matrix
bool Population::putNji(int num, int row, int locus)
{
  if(locus > numLoci-1 || locus < 0)
    {
      //bad bounds
      return 0;
    }
  else if(row > getNjiColLength(locus)-1 || row < 0)
    {
      //bad bounds
      return 0;
    }
  else
    {
      nji[(*loci).offset[locus] + row] = num;
      return 1;
    }
}

/*
 * Size the count block from the shared offsets. One allocation per grouping,
 * where a row per locus meant one per locus per grouping.
 */
void Population::allocNji()
{
  nji.assign(size_t(loci ? loci->offset.back() : 0),0);
  return;
}

//Constructor: nothing is sized until setLoci() is called
Population::Population()
{
  loci = NULL;
  Nj = NULL;
  missing = NULL;
  numLoci = 0;
  rows = 0;
  name = "UNDEF";
  minG = INT_MAX;
  tol = 1;
}

//Memclean
Population::~Population()
{
  if(Nj) delete [] Nj;
  if(missing) delete [] missing;
}

/*
 * Size the per-locus arrays.
 *
 * ADZE 1.0's setRowsLoci also allocated a numLoci x rows matrix holding every
 * genotype as a std::string -- 34 bytes per gene copy, 8.5x the size of the
 * input file.  Nothing needs it: the allele counts and missing-data tallies
 * the estimators use are accumulated while the file is read.
 */
void Population::setLoci(int l)
{
  numLoci = l;

  Nj = new int[numLoci];
  fillNj(0);
  missing = new int[numLoci];
  for(int i = 0; i < numLoci;i++)
    {
      missing[i] = 0;
    }
  return;
}

void Population::putMissing(int count, int locus)
{
  if(locus >= 0 && locus < numLoci) missing[locus] = count;
  return;
}

const string& Population::getLocusName(int pos) const
{
  static const string bad = "BAD REF";

  if(!loci || pos < 0 || pos > numLoci-1 || size_t(pos) >= loci->name.size())
    {
      return bad;
    }

  return loci->name[pos];
}
