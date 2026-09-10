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


void Population::deleteLocus(int l)
{
  deletedLocus.push_back(locusName[l]);

  for(l; l < numLoci-1; l++)
    {
      locusName[l] = locusName[l+1];
      for(int r = 0; r < rows; r++)
	{
	  data[l][r] = data[l+1][r];
	}
    }

  numLoci--;
  return;
}

list<int> Population::recLociDelete(double tolerance, string miss)
{
  tol = tolerance;
  double* missing = new double[numLoci];
  list<int> rec;
  
  //tally the number of missing at each locus
  for(int l = 0; l < numLoci; l++)
    {
      missing[l] = 0;
      for(int r = 0; r < rows; r++)
	{
	  if(data[l][r].compare(miss) == 0) missing[l]++;
	}
    }

  //change to a proportion of data missing and push to vector if > tolerance
  for(int l = 0; l < numLoci; l++)
    {
      missing[l] /= double(rows);
      if(missing[l] > tolerance) rec.push_back(l);
    }
  
  delete [] missing;

  return rec;  
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

string Population::getDataElement(int line, int locus)
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
      //good vals
      return data[locus][line];
    }      
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
