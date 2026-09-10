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
 * ADZE 1.0 deleted one locus at a time, shifting all higher-indexed loci down
 * by one and copying every genotype string in them, so filtering D of L loci
 * cost O(D*L*rows) string assignments -- quadratic in the locus count.
 * Compacting once is O(L): each locus's row block is moved by swapping the
 * column pointer, so no genotype is copied at all.  Swapping rather than
 * assigning keeps every allocated block reachable from data[0..numLociOrig),
 * so the destructor still frees each one exactly once.
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
	  string* tmp = data[keep];
	  data[keep] = data[l];
	  data[l] = tmp;
	}
      keep++;
    }

  numLoci = keep;
  return;
}

/*
 * Flag loci whose missing-data fraction exceeds the tolerance in this grouping.
 * Called once per grouping; del[] accumulates the union across groupings.
 */
void Population::recLociDelete(double tolerance, const string& miss,
			       vector<char>& del)
{
  tol = tolerance;

  for(int l = 0; l < numLoci; l++)
    {
      if(del[l]) continue; //already condemned by another grouping

      int missing = 0;
      for(int r = 0; r < rows; r++)
	{
	  if(data[l][r].compare(miss) == 0) missing++;
	}

      if(double(missing)/double(rows) > tolerance) del[l] = 1;
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

//Calculated the probability of not finding allele i 
//in a sample of size g from population j
double Population::calcQjig(int i, int g, int locus)
{
  double Q = 1;
  
  
  /* Calculate Qjig
   *          / Nj - Nji \       g-1
   *          \    g     /      ----    Nj - Nji - u
   *  Qjig = --------------  =  |  | -----------------
   *            /  Nj \         |  |      Nj - u
   *            \  g  /         u = 0
   *
   */
  for (int u = 0; u < g; u++)
    {
      Q *= double(Nj[locus]-Nji[locus][i]-u)/double(Nj[locus]-u);
    }
  
  //Slower calculation...
  //Q = nCk(Nj[locus]-Nji[locus][i],g)/nCk(Nj[locus],g);


  return Q;
}

//Implements calcQjig to calculate the allelic richness of population j
//at a sample size of g
double Population::calcAg(int g, int locus)
{
  double ag = 0, P, Q = 1;
  bool good;
  
  if(locus > numLoci-1 || locus < 0)
    {
      //bad bounds
      return -9;
    }
  else if(g > Nj[locus] || g < 2)
    {
      //Can't have a g greater than your actual sample size or less than 2
      return -9;
    }
  else
    {
           
      /* Calculate Qjig
       *          / Nj - Nji \       g-1
       *          \    g     /      ----    Nj - Nji - u
       *  Qjig = --------------  =  |  | -----------------
       *            /  Nj \         |  |      Nj - u
       *            \  g  /         u = 0
       *
       */
      for (int i = 0; i < NjiColLength[locus]; i++)
	{
	  for (int u = 0; u < g; u++)
	    {
	      Q *= double(Nj[locus]-Nji[locus][i]-u)/double(Nj[locus]-u);
	    }
	  
	  //Calculate Pjig
	  P = 1 - Q;
	  //Calc ag
	  ag += P;
	  
	  /*TESTING OUTPUT
	    cout << "Qjig[" << j << "][" << i << "] = " 
	    << Q << " ";
	  
	  cout << "Pj" << i << g << " = " << P << "\n";
	  TESTING*/
	  Q = 1;
	}
      
      return ag;
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

//Returns the number of alleles in population j at given locus
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
    }
}

//Returns the number of i alleles in jth population at given locus
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
      Nji[locus] = new int[size];
      for(int i = 0; i < size;i++)
	{
	  Nji[locus][i] = 0;
	}
      return 1;
    }
}

//Constructor with no args: initializes numLoci, rows, and name
Population::Population()
{
  data = NULL;
  locusName = NULL;
  Nji = NULL;
  NjiColLength = NULL;
  Nj = NULL;
  numLoci = 0;
  numLociOrig = 0;
  rows = 0;
  name = "UNDEF";
  minG = INT_MAX;
}

//Constructor with three args: initialzes
//name, numLoci, rows, data, Nj,NjiColLength, Nji, locusName
Population::Population(string n, int numL, int r)
{
  name = n;
  numLoci = numL;
  numLociOrig = numL;
  rows = r;
  minG = INT_MAX;

  //Allocate the data matrix with dimentions, numLoci, rows
  data = new string*[numLoci];
  for(int i = 0; i < numLoci; i++)
    {
      data[i] = new string[rows];
    }
  
  fillData("0");//Initialize with 0's


  Nj = new int[numLoci];//Allocate Nj vector of length numLoci
  fillNj(0);//Initialize
  NjiColLength = new int[numLoci];//Allocate NjiColLength vector
  Nji = new int*[numLoci];
  for(int i = 0; i < numLoci;i++)
    {
      Nji[i] = NULL;
    }

  locusName = new string[numLoci];

}

//Memclean
Population::~Population()
{
  for (int i = 0; i < numLociOrig; i++)
    {
      if(data[i]) delete [] data[i];
      if(Nji[i]) delete [] Nji[i];
    }
  
  if(Nj) delete [] Nj;
  if(NjiColLength) delete [] NjiColLength;
  if(Nji) delete [] Nji;
  if(data) delete [] data;
  if(locusName) delete [] locusName;
}


void Population::setRowsLoci(int r, int l)
{
  numLoci = l;
  numLociOrig = l;
  rows = r;

  //cout << rows << " " << numLoci << endl;

  data = new string*[numLoci];
  for(int i = 0; i < numLoci; i++)
  {
	data[i] = new string[rows];
  }

  fillData("0");

  Nj = new int[numLoci];
  fillNj(0);
  NjiColLength = new int[numLoci];
  Nji = new int*[numLoci];
  for(int i = 0; i < numLoci;i++)
    {
      Nji[i] = NULL;
    }
  locusName = new string[numLoci];
  return;
}

/*
 * Returns a reference, not a copy: this is called once per genotype per locus
 * by the binning pass, and returning by value copied a std::string on every
 * access.  Out-of-range indices used to `return 0`, i.e. construct a
 * std::string from a null pointer; they now yield a reference to an empty
 * string, which callers already treat as "no allele here".
 */
const string& Population::getDataElement(int line, int locus) const
{
  static const string outOfRange;

  if (locus > numLoci-1 || locus < 0) return outOfRange;
  if (line > rows-1 || line < 0) return outOfRange;

  return data[locus][line];
}

bool Population::putDataElement(string dataElem, int line, int locus)
{
  if (locus > numLoci-1 || locus < 0)
    {
      //error
      return 0;
    }
  else if (line > rows-1 || line < 0)
    {
      //error
      return 0;
    }
  else 
    {
      data[locus][line] = dataElem;
      return 1;
    }
}

void Population::fillData(string x)
{
  for (int i = 0; i < numLoci; i++)
    {
      for (int j = 0; j < rows; j++)
	{
	  data[i][j] = x;
	}
    }

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
