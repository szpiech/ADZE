#include "ADZE_pop.h"

using namespace std;

void Population::printDeleted(ostream& out)
{
  cout << deletedLocus.size();
  out << deletedLocus.size();
  if(deletedLocus.size() == 1)
    {
      cout << " locus has ";
      out << " locus has ";
    }
  else
    {
      cout << " loci have ";
      out << " loci have ";
    }
  cout << "at least one grouping with";
  cout.setf(ios::fixed,ios::floatfield);
  if(tol > 0) cout << " more than " << 100*tol << "%";
  cout << " missing data.\n";

  out << "at least one grouping with";
  if(tol > 0) out << " more than " << 100*tol << "%";
  out << " missing data.\n";
 
  for(vector<string>::iterator i = deletedLocus.begin();
      i != deletedLocus.end(); i++)
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
 * count.  Compacting once is O(L), and the per-locus allele-count blocks move
 * by swapping pointers.  Swapping rather than assigning keeps every allocated
 * block reachable from Nji[0..numLociOrig), so the destructor still frees each
 * one exactly once.
 *
 * Deleted names are recorded from the highest index down, matching the order
 * that 1.0's reverse-iterating delete loop wrote to the _deletedloci file.
 */
void Population::deleteLoci(const vector<char>& del)
{
  for(int l = numLoci-1; l >= 0; l--)
    {
      if(del[l]) deletedLocus.push_back(locusName[l]);
    }

  int keep = 0;
  for(int l = 0; l < numLoci; l++)
    {
      if(del[l]) continue;
      if(keep != l)
	{
	  locusName[keep].swap(locusName[l]);
	  int* tmp = Nji[keep];
	  Nji[keep] = Nji[l];
	  Nji[l] = tmp;
	  NjiColLength[keep] = NjiColLength[l];
	  Nj[keep] = Nj[l];
	  missing[keep] = missing[l];
	}
      keep++;
    }

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
      return NjiColLength[locus];
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
      for(int i = 0; i < NjiColLength[locus]; i++)
	{
	  total += Nji[locus][i];
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
  else if(row > NjiColLength[locus]-1 || row < 0)
    {
      //bad bounds
      return -9;
    }
  else
    {
      return Nji[locus][row];
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
  else if(row > NjiColLength[locus]-1 || row < 0)
    {
      //bad bounds
      return 0;
    }
  else
    {
      Nji[locus][row] = num;
      return 1;
    }
}

//allocate size of Nji columns
bool Population::setNjiColLength(int size,int locus)
{
  if(locus > numLoci-1 || locus < 0)
    {
      //bad bounds
      return 0;
    }
  else
    {
      NjiColLength[locus] = size;
      if(Nji[locus]) delete [] Nji[locus];
      Nji[locus] = new int[size];
      for(int i = 0; i < size;i++)
	{
	  Nji[locus][i] = 0;
	}
      return 1;
    }
}

//Constructor: nothing is sized until setLoci() is called
Population::Population()
{
  locusName = NULL;
  Nji = NULL;
  NjiColLength = NULL;
  Nj = NULL;
  missing = NULL;
  numLoci = 0;
  numLociOrig = 0;
  rows = 0;
  name = "UNDEF";
  minG = INT_MAX;
  tol = 1;
}

//Memclean
Population::~Population()
{
  if(Nji)
    {
      for (int i = 0; i < numLociOrig; i++)
	{
	  if(Nji[i]) delete [] Nji[i];
	}
    }

  if(Nj) delete [] Nj;
  if(NjiColLength) delete [] NjiColLength;
  if(Nji) delete [] Nji;
  if(missing) delete [] missing;
  if(locusName) delete [] locusName;
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
  numLociOrig = l;

  Nj = new int[numLoci];
  fillNj(0);
  NjiColLength = new int[numLoci];
  missing = new int[numLoci];
  Nji = new int*[numLoci];
  for(int i = 0; i < numLoci;i++)
    {
      Nji[i] = NULL;
      NjiColLength[i] = 0;
      missing[i] = 0;
    }
  locusName = new string[numLoci];
  return;
}

void Population::putMissing(int count, int locus)
{
  if(locus >= 0 && locus < numLoci) missing[locus] = count;
  return;
}

bool Population::setLocusName(string str, int pos)
{
  if (pos > numLoci-1 || pos < 0)
    {
      //error
      return 0;
    }
  else
    {
      locusName[pos] = str;
      return 1;
    }
}

string Population::getLocusName(int pos)
{
  if (pos > numLoci-1 || pos < 0)
    {
      return "BAD REF";
    }
  else
    {
      return locusName[pos];
    }
}
