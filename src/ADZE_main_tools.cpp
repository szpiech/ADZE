#include "ADZE_main_tools.h"

using namespace std;

list<int> parseKVals(string str)
{
  int size = 10;
  char val[size];
  list<int> nums;
  string err = "ERROR: \"";
  err += str;
  err += "\" not a valid K_RANGE definition.\n";

  //cout << str.size() << endl;

  for(int i = 0; i < str.size(); i++)
    {
      //cout << i << " ";
      if(i >= str.size()) break;

      for(int n = 0; n < size-1; n++) val[n] = '\0';
 
      if(str.at(i) == '#') break;

      if(str.at(i) == ' ' ||
	 //	 str.at(i) == '0' ||
	 isalpha(str.at(i))) i++;
      if(i >= str.size()) break;

      if(isdigit(str.at(i)))
	{
	  int a = 0;
	  do
	    {
	      val[a] = str.at(i);
	      i++;
	      if(i >= str.size()) break;
	      a++;
	      if(a >= size) throw err;
	    }while(isdigit(str.at(i)) && i < str.size());

	  nums.push_back(atoi(val));
	}

      if(i >= str.size()) break;
      if(str.at(i) == '-')
	{
	  if(i+1 >= str.size()) break;
	  if(!isdigit(str.at(i+1))) throw err;

	  int a = 0;
	  do
	    {
	      val[a] = str.at(i);
	      i++;
	      if(i >= str.size()) break;
	      a++;
	      if(a >= size) throw err;
	    }while(isdigit(str.at(i)));
	  nums.push_back(atoi(val));
	}
      
    }

  int prev;

  for(list<int>::iterator i = nums.begin(); i != nums.end(); i++)
    {
      //if(*i == 0) throw err;
      if(*i < 0)
	{
	  if((0-(*i)) < prev) throw err;
	  for(int n = prev; n < (0-(*i)); n++)
	    {
	      nums.push_back(n);
	    }
	  //nums.erase(i);
	  *i = (0-(*i));
	}
      prev = *i;
    }

  if(nums.empty()) throw err;

  nums.sort();
  nums.unique();

  return nums;
  
}



void filterLoci(Population pop[],int numDivs, double tol, string file,
		string missing, bool pp)
{
  vector<char> toDelete(pop[0].getNumLoci(),0);

  for(int n = 0; n < numDivs; n++)
    {
      pop[n].recLociDelete(tol,missing,toDelete);
    }

  int size1 = 0;
  for(size_t l = 0; l < toDelete.size(); l++)
    {
      if(toDelete[l]) size1++;
    }

  int size = (size1 == 0) ? 1 : size1;

  ProgressBar bar(&cout,size*numDivs,BARLEN[0]);
  if(pp)
    {
      bar.init();
    }
  
  for(int n = 0; n < numDivs; n++)
    {
      pop[n].deleteLoci(toDelete);
      if(pp) bar.adv(size);
    }
  if(pp) bar.done();
  cout << endl;
  
  file = nameCreate(file,"_deletedloci");

  ofstream lout;
  lout.open(file.c_str());

  pop[0].printDeleted(lout);

  return;
}




double displayTime(ostream& out)
{
  clock_t tempClock = clock();
  double CPU_sec = tempClock/CLOCKS_PER_SEC;
  double sec = CPU_sec;
  
  out.unsetf(ios::floatfield);
  out.unsetf(ios::showpoint);

  out << floor(CPU_sec/86400.0) << ":";
  CPU_sec = CPU_sec - floor(CPU_sec/86400.0)*86400.0;
  out << floor(CPU_sec/3600.0) << ":";
  CPU_sec = CPU_sec - floor(CPU_sec/3600.0)*3600.0;
  out << floor(CPU_sec/60.0) << ":";
  CPU_sec = CPU_sec - floor(CPU_sec/60.0)*60.0;
  out << CPU_sec;
  
  return (sec);
}

bool validK(int n, list<int> k)
{
  int kback = k.back();
  int kfront = k.front();
  
  if(kback <= 0 || kback > n || kfront <= 0 || kfront > n)
    {
      cout << "ERROR: For these parameters, K_RANGE values must be between "
	   << "1 and " << n << " inclusive.\n";
      return 0;
    }
  else
    {
      return 1;
    }
}


string combineNames(string names[],int k)
{
  string tmp = "";

  for(int i = 0; i < k; i++)
    {
      tmp += names[i];
      if(i != k-1) tmp += ' ';
    }

  return tmp;
}

int min(int a, int b)
{
  if(a<b) return a;
  else return b;
}


void calcAllPgComb(Population pop[], int numDivs, int k, const ParamSet &param,
		   bool full_comb, string comb_out)
{
  int maxG = param.g.val+1;
  int numLoci = param.loci.val;
  int tot_m = int(nCk(numDivs,k));
  double* pgcomb;
  string* names;
  //bool gotNames;

  ofstream full_out, reg_out;

  if(full_comb)
    {
      string name;
      name = nameCreate(comb_out,"_fulldata");
      full_out.open(name.c_str());
      
      for(int i = 1; i <= k; i++)
	{
	  full_out << "POP_GROUPING" << i << " ";
	}
      full_out << "G NUM_LOCI ";
      for(int l = 0; l < numLoci; l++)
	{
	  full_out << pop[0].getLocusName(l) << " ";
	}
      
      full_out << "MEAN VAR STD_ERR\n";
    }

  reg_out.open(comb_out.c_str());


  gsl_combination* c = gsl_combination_calloc(numDivs,k);
  gsl_combination_init_first(c);

  //for every instance of nCk
  //for(int m = 0; m < tot_m; m++)

  ProgressBar bar(&cout,tot_m*(maxG-2)*numLoci,BARLEN[min(k-1,3)]);
  if(param.pp.val)
    {
      bar.init();
    }
  
  int m = 0;
  do
    {
  
      //array for the names
      names = new string[k];
      //gotNames = 0;

      //for all g
      for(int g = 2; g < maxG; g++)
	{
	  bool breakG = 0;
	  //allocate 1-dim array
	  pgcomb = new double[numLoci];
	  
	  for(int l = 0; l < numLoci; l++)
	    {
	      double pg = 0;
	      int maxI = pop[0].getNjiColLength(l);

	      //Sum over all alleles
	      for(int i = 0; i < maxI; i++)
		{
		  double P = 1, Q = 1;
		  //Calc P's
		  for(size_t j = 0; j < gsl_combination_k(c); j++)
		    {
		      int c_kmj = int(gsl_combination_get(c,j));
		      
		      P *= (1-pop[c_kmj].calcQjig(i,g,l));
		      
		      names[int(j)] = pop[c_kmj].getName();
		    }
		  
		  //Calc Q's
		  for(size_t j_p = 0; j_p < size_t(numDivs); j_p++)
		    {
		      if(!isIn(j_p,c))
			{
			  Q *= pop[j_p].calcQjig(i,g,l);
			}
		    }
		  
		  //Multiply together and add to total
		  pg += (P*Q);
		}
	      //cout << "\npgtotal = " << pg << endl;
	      pgcomb[l] = pg;
	      
	      //Good to go for next g?
	      for(int p = 0; p < numDivs; p++)
		{
		  int Nj = pop[p].getNj(l);
		  if(Nj < g+1) breakG = 1; 
		} 

	      if(param.pp.val) ++bar;

	    }
	  
	  Stats comb_stats;
	  comb_stats.putData(pgcomb,numLoci);

	  //Calc Avg
	  comb_stats.calcAvg();
	  
	  //Calc Var
	  comb_stats.calcVar();
	  
	  //Calc Std_err
	  comb_stats.calcStdErr();
	  
	  //output
	  string all_names = combineNames(names,k);
	  comb_stats.printStats(reg_out,all_names,g);

	  if(full_comb)
	    {
	      comb_stats.printData(full_out,all_names,g);
	    }
	  
	  //de-allocate array
	  delete [] pgcomb;
	  
	  if(breakG)
	    {
	      if(param.pp.val) bar.adv((maxG-g-1)*numLoci);
	      break;
	    }
	}

      delete [] names;

      reg_out << endl;
      full_out << endl;
      
      m++;
    }while(gsl_combination_next(c) == GSL_SUCCESS);
  
  if(param.pp.val) bar.done();

  gsl_combination_free(c);

  return;
}

bool isIn(int j, gsl_combination* c)
{
  size_t size = gsl_combination_k(c);

  for(size_t i = 0; i < size; i++)
    {
      if(size_t(j) == gsl_combination_get(c,i))
	{
	  return 1;
	}
    }
  
  return 0;
}


void calcAllPgs(Population pop[],int numDivs,const ParamSet &param,
		bool full_priv,string private_out)
{

  int numLoci = param.loci.val;
  int maxG = param.g.val+1;
  double* pg; //private allelic richness
  ofstream pg_full_out,pg_out;

  if(full_priv)
    {
      string name;
      name = nameCreate(private_out,"_fulldata");
      pg_full_out.open(name.c_str());
      
      pg_full_out << "POP_GROUPING G NUM_LOCI ";
      for(int l = 0; l < numLoci; l++)
	{
	  pg_full_out << pop[0].getLocusName(l) << " ";
	}
      
      pg_full_out << "MEAN VAR STD_ERR\n";
    }

  pg_out.open(private_out.c_str());


  ProgressBar bar(&cout,numDivs*(maxG-2)*numLoci,BARLEN[0]);
  if(param.pp.val)
    {
      bar.init();
    }

  //CALCULATE pg's
  for(int j = 0; j < numDivs; j++)
    {
      for(int g = 2; g < maxG; g++)
	{
	  bool breakG = 0;
	  pg = new double[numLoci];

	  for(int locus = 0;locus < numLoci; locus++)
	    {
	      pg[locus] = calcPg(pop,j,locus,g,numDivs);

	      //Good to go for next g?
	      for(int p = 0; p < numDivs; p++)
		{
		  int Nj = pop[p].getNj(locus);
		  if(Nj < g+1) breakG = 1; 
		}

	      if(param.pp.val) ++bar;
	    }
	  
	  Stats pg_stats;
	  pg_stats.putData(pg,numLoci);

	  //Calc Avg
	  pg_stats.calcAvg();
	  
	  //Calc Var
	  pg_stats.calcVar();
	  
	  //Calc Std_err
	  pg_stats.calcStdErr();

	  //Ouput
	  pg_stats.printStats(pg_out,pop[j].getName(),g);

	  if(full_priv)
	    {
	      pg_stats.printData(pg_full_out,pop[j].getName(),g);
	    }

	  //Delete
	  delete [] pg;
	  
	  if(breakG)
	    {
	      if(param.pp.val) bar.adv((maxG-g-1)*numLoci);
	      break;
	    }
	}
    
      pg_out << endl;
      pg_full_out << endl;
    }

  if(param.pp.val) bar.done();


  if(pg_full_out.is_open())
    {
      pg_full_out.close();
    }

  pg_out.close();
  

  return;
}


void calcAllAgs(Population pop[],int numDivs,const ParamSet &param,
		bool full_rich,string richness_out)
{

  int numLoci = param.loci.val;
  int maxG = param.g.val+1;
  double* ag; //allelic richness
  ofstream ag_full_out,ag_out;

  if(full_rich)
    {
      string name;
      name = nameCreate(richness_out,"_fulldata");
      ag_full_out.open(name.c_str());

      ag_full_out << "POP_GROUPING G NUM_LOCI ";
      for(int l = 0; l < numLoci; l++)
	{
	  ag_full_out << pop[0].getLocusName(l) << " ";
	}
      
      ag_full_out << "MEAN VAR STD_ERR\n";
    }

  ag_out.open(richness_out.c_str());

  ProgressBar bar(&cout,numDivs*(maxG-2)*numLoci,BARLEN[0]);
  if(param.pp.val)
    {
      bar.init();
    }
  
  //CALCULATE ag's
  for(int j = 0; j < numDivs; j++)
    {
      for(int g = 2; g < maxG; g++)
	{
	  ag = new double[numLoci];

	  for(int locus = 0;locus < numLoci; locus++)
	    {
	      
	      ag[locus] = pop[j].calcAg(g,locus);
	      if(param.pp.val) ++bar;
	    }
	  
	  Stats ag_stats;
	  ag_stats.putData(ag,numLoci);

	  //Calc Avg
	  ag_stats.calcAvg();
	  
	  //Calc Var
	  ag_stats.calcVar();
	  
	  //Calc Std_err
	  ag_stats.calcStdErr();

	  //Ouput
	  ag_stats.printStats(ag_out,pop[j].getName(),g);

	  if(full_rich)
	    {
	      ag_stats.printData(ag_full_out,pop[j].getName(),g);
	    }

	  //Delete

	  delete [] ag;	  
	}
      ag_out << endl;
      ag_full_out << endl;
    }

  if(param.pp.val) bar.done();

  if(ag_full_out.is_open())
    {
      ag_full_out.close();
    }

  ag_out.close();
  

  return;
}




double calcPg(Population pop[],int j,int locus,int g,int numDivs)
{
  int Nj;
  
  for(int p = 0; p < numDivs; p++)
    {
      Nj = pop[p].getNj(locus);

      if(Nj < g)
	{
	  return -9;
	}      
    }

  int maxI = pop[j].getNjiColLength(locus);
  double pg = 0;
  double Q = 1;
  double P;

  /*
   * Calculate the private allelic richness
   *            m              J
   *            _             ___ 
   * __(j)     \   /        / | |       \ \
   * ||g  =    /_  \ Pijg * \ | | Qij'g / /
   *           i=1           j'=1
   *                         j'!=j
   */

  for(int i = 0; i < maxI; i++)
    {
      for(int p = 0; p < numDivs; p++)
	{
	  if(p != j)
	    {
	      Q *= pop[p].calcQjig(i,g,locus);
	    }
	}

      P = 1 - pop[j].calcQjig(i,g,locus);
      pg += P*Q;
      Q = 1; //reset Q
    }

  return pg;
}


/*
 * Calculates all Nji's for every locus
 * Stores in Population objects
 *
 */
void calcNji(Population pop[],int numDivs,string missing)
{
  int foundAt, count;
  string data, current;
  int numLoci = pop[0].getNumLoci();
  bool good, seen = 0;
  int NjiColLength;
  int** Nji = new int*[numDivs];
  int tally = 0;

  IntStrMap alleleName;
  IntStrMap::iterator pos;

  //At each locus
  for(int locus = 0; locus < numLoci; locus++)
    {
      for(int j = 0; j < numDivs; j++)
	{
	  for(int row = 0; row < pop[j].getNumRows(); row++)
	    {
	      data = pop[j].getDataElement(row,locus);
	      if(data.compare(missing) == 0)
		{
		  //do nothing
		}
	      else if(alleleName.size() == 0)
		{
		  alleleName[1] = data;
		}
	      else
		{
		  for(pos = alleleName.begin(); pos != alleleName.end();pos++)
		    {
		      if(data.compare(pos->second) == 0)
			{
			  seen = 1;
			}
		      count = pos->first;
		    }
		  if(!seen)
		    {
		      alleleName[count+1] = data;
		    }
		  seen = 0;
		}
	    }
	}

      pos = alleleName.end();
      pos--;
      
      NjiColLength = pos->first;

      for (int j = 0; j < numDivs; j++)
	{
	  pop[j].setNjiColLength(NjiColLength,locus);
	}


      for(int j = 0; j < numDivs; j++)
	{
	  Nji[j] = new int[NjiColLength];
	  for(int i = 0; i < NjiColLength;i++)
	    {
	      Nji[j][i] = 0;
	    }
	}
      
      //Look in each division and bin up all the alleles
      for(int j = 0; j < numDivs; j++)
	{
	  //look at each allele
	  for(int row = 0; row < pop[j].getNumRows(); row++)
	    {
	      current = pop[j].getDataElement(row,locus);
	      
	      if(current.compare(missing) == 0)
		{
		  //Missing data, do nothing
		}
	      //Otherwise, we have to check the other bins
	      else
		{  
		  //Look to see if current has been seen before and where
		  for(pos=alleleName.begin(); pos!=alleleName.end();pos++)
		    {
		      count = pos->first;
		      if(current.compare(pos->second) == 0)
			{
			  Nji[j][count-1]++;
			}
		    }	  
		}
	    }
	}

      //store Nji's in population objects
      for (int j = 0; j < numDivs;j++)
	{
	  //good = pop[j].setNjiColLength(Nji[j].size(),locus);
	  for(int i = 0; i < NjiColLength;i++)
	    {
	      good = pop[j].putNji(Nji[j][i],i,locus);
	    }
	}

      //Clear out Nji for use at next locus
      for (int j = 0; j < numDivs; j++)
	{
	  delete [] Nji[j];
	}

      //clear names for next locus
      alleleName.clear();
    }
  
  delete [] Nji;
  return;  
}

/* CALCULATE Nj
 *      m
 *      __
 * Nj = \  Nji
 *      /_
 *     i = 0
 *
 * where m = Nji[j].size()
 */
void calcNj(Population pop[],int numDivs)
{
  int loci = pop[0].getNumLoci();
  int Nj, tally = 0;
  bool good;

  for(int locus = 0; locus < loci; locus++)
    {
      for (int j = 0; j < numDivs; j++)
	{
	  Nj = 0;
	  for (int i = 0; i < pop[j].getNjiColLength(locus); i++)
	    {
	      Nj += pop[j].getNji(i,locus);
	    }
	  good = pop[j].putNj(Nj,locus);
	}
    }
  return;
}

int seenBefore(vector<int>& allele, int current)
{
  for (int i = 0; i < allele.size(); i++)
    {
      if (allele[i] == current)
	{
	  return i;
	}
    }
  
  return -9;
}

void readData(/*ifstream& data,*/ Population pop[],
	      vector<string>& sortLabel,int numDivs, const ParamSet &param)
{
  ifstream data;
  data.open(param.dfile.val.c_str());
  string junk;

  for(int i = 0; i < param.nd_rows.val;i++)
    {
      getline(data,junk);
    }

  string *sortOption;
  sortOption = new string[param.nd_cols.val];
  //cout << "param.nd_cols.val " << param.nd_cols.val << endl;

  //track which line of the division data block to write to
  int *divLine = new int[numDivs]; 
  
  for (int i = 0; i < numDivs; i++)
    {
      divLine[i] = 0;
    }
  
  string element; // int element;
  bool goodStore;
  
  //Initial read of data line 1
  for (int j = 0; j < param.nd_cols.val; j++)
    {
      data >> sortOption[j];
    }
  
  int p = 0; //population index to write in
  
  for (int l = 0; l < param.loci.val; l++)
    {
      data >> element;
      goodStore = pop[p].putDataElement(element,divLine[p],l);
    }
  
  divLine[p]++;
  
  //Go line by line through the data file
  //First read the preceding columns before the data
  //and store as a string array
  //split the data into chunks based on requested parameter
  //store in the class population
  for (int i = 1; i < param.dlines.val; i++)
    {
      
      for (int j = 0; j < param.nd_cols.val; j++)
        {
	  data >> sortOption[j];
	  //cout << "sortOption[j] " << sortOption[j] << endl;
	}
      //Check to see if the current data line falls in the current category
      if (!sameStr(sortLabel[p],sortOption[(param.sort_by.val - 1)]))
	{
	  //It does not, so check to see which division the data line
	  //fits into
	  /*cout << "sortLabel,sortOption[(param.sort_by.val - 1) = " 
	       << sortLabel[p] << " " << sortOption[(param.sort_by.val - 1)] 
	       << endl;*/
	  p = seenLabel(sortLabel,sortOption[(param.sort_by.val - 1)]);
	  //cout << p << endl;
	}
      
      //Read in the data to the appropriate object indexed by p
      for (int l = 0; l < param.loci.val; l++)
	{
	  data >> element;
	  //cout << i << " " << l << " " << element << endl;
	  goodStore = pop[p].putDataElement(element,divLine[p],l);
	}
      
      divLine[p]++;
    }
  
  delete [] divLine;
  delete [] sortOption;
  
  return;
}

void checkDatafile(const ParamSet &p)
{

  ifstream data;
  data.open(p.dfile.val.c_str());

  string line;
  int count;
 
  getline(data,line);
  count = countCols(line);
  
  if(count != p.loci.val)
    {
      cout << "ERROR: Expected " << p.loci.val; 
      if(p.loci.val == 1) cout << " locus ";
      else cout << " loci ";
      cout << "in " << p.dfile.val << " but found " << count << ".\n";
      data.close();
      BAD_PARAM x;
      throw x;
    }
 
  for(int i = 1; i < p.nd_rows.val; i++)
    {
      getline(data,line);
    }

  count = 0;

  int count2 = 0;
  int blank = 0;
  //getline(data,line);
  while(!data.eof())
    {
      getline(data,line);
      count2++;
      count = countCols(line);
      if(count == 0) blank++;
      if(count2-blank > p.dlines.val)
	{
	  cout << "ERROR: Expected " << p.dlines.val << " data lines "
	       << "in " << p.dfile.val << "\nbut found at least "
	       << count2-blank << ".\n";
	  data.close();
	  BAD_PARAM x;
	  throw x;
	}
      if(count != (p.nd_cols.val + p.loci.val) && count != 0)
	{
	  cout << "ERROR: Expected " << p.nd_cols.val + p.loci.val
	       << " columns at data line " << count2 << " in " 
	       << p.dfile.val << "\nbut found " << count
	       << ". Possible bad NON_DATA_COLS, NON_DATA_ROWS, "
	       << "or LOCI value.\n";
	  data.close();
	  BAD_PARAM x;
	  throw x;
	}
    }

  if(count2-blank != p.dlines.val)
    {
      cout << "ERROR: Expected " << p.dlines.val << " data lines "
	   << "in " << p.dfile.val << "\nbut found " << count2-blank << ".\n";
      data.close();
      BAD_PARAM x;
      throw x;
    }

  data.close();
  return;
}

int countCols(string s)
{
  char c;
  int count = 0;
  bool countThis = 1;
  string::iterator i;
  for(i = s.begin(); i != s.end(); i++)
    {
      c = *i;
      if(isgraph(c) && countThis)
	{
	  count++;
	  countThis = 0;
	}
      else if(!isgraph(c)) countThis = 1;
    }
  return count;
}

void getLociNames(ifstream& data, string names[], const ParamSet &p)
{
  string junk;
  //Read in all the names of the loci
  for (int l = 0; l < p.loci.val; l++)
    {
      data >> names[l];
    }
  
  getline(data,junk);

  return;
}

void getDivLines(ifstream& data,vector<string>& divisionNames,
		 vector<int>& lines,const ParamSet &p)
{
  string junk;
  
  string *sortOption;
  sortOption = new string[p.nd_cols.val];
  
  //Initial read of before data strings
  for (int j = 0; j < p.nd_cols.val; j++)
    {
      data >> sortOption[j];
    }
  getline(data,junk); //throw away the rest of the line
  
  int index, foundAt;
  
  divisionNames.push_back(sortOption[(p.sort_by.val - 1)]);
  lines.push_back(1);
  index = 0;
  
  //Read through data and count lines in each division, and keep their names
  for (int i = 1; i < p.dlines.val; i++)
    {
      //Read in before data strings
      for (int j = 0; j < p.nd_cols.val; j++)
	{
	  data >> sortOption[j];
	}
      getline(data,junk); //Throw away rest of line
      
      //Look to see if the current line division has been seen before
      foundAt = seenLabel(divisionNames,sortOption[(p.sort_by.val - 1)]);
      
      if (foundAt < 0) //Its a new division
	{
	  divisionNames.push_back(sortOption[(p.sort_by.val - 1)]);
	  lines.push_back(1);
	}
      else // foundAt  >= 0, its a previously seen divison
	{
	  lines[foundAt] = lines[foundAt]+1;
	}
    }
  
  delete [] sortOption;
  
  return;
}

/*
 * seenLabel
 * INPUT:
 *       a string vector in which to search
 *       a string to search for
 * OUTPUT:
 *       returns the index at which the string was found in the vector
 *       otherwise returns -9
 * FUNCTION:
 *       uses sameStr to compare the string of interest to every
 *       element in the vector
 */
int seenLabel(vector<string>& seen, string current)
{
  for (int i = 0; i < seen.size(); i++)
    {
      //cout << "seen " << seen[i] << " current " << current << endl;
      if (sameStr(seen[i],current))
	{
	  return i;
	}
    }
  
  return -9;
}


/*
 * sameStr
 * INPUT:
 *	two strings
 * OUTPUT:
 *	TRUE of the strings are the same
 *	FALSE if the strings are different
 * FUNCTION:
 *	uses the compare() class function to compare strings
 */
bool sameStr(string str1, string str2)
{
  int outcome = str1.compare(str2);
  
  if (outcome < 0 || outcome > 0)
    {
      return 0;
    }
  else
    {
      return 1;
    }
}

string nameCreate(string name,string toPut)
{
  name += toPut;
  return name;

  /*
  int len = name.length();
  string temp_front, temp_back;
  int put_here = name.length();

  for(int l = 0; l < len; l++)
    {
      if(name[l] == '.')
	{
	  put_here = l;
	}
    }

  name.insert(put_here,toPut);

  return name;
  */
}
