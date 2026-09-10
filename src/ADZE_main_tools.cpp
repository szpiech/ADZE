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


/*
 * Fill q[(p*numAlleles + i)*gStride + g] with Qpig -- the probability that
 * allele i is absent from a sample of g gene copies drawn from grouping p --
 * for one locus, every grouping, every allele, and every g in 1..gMax.
 *
 *           / Nj - Nji \        g-1
 *           \    g     /       ----   Nj - Nji - u
 *  Qjig =  --------------  =   |  |  --------------
 *             /  Nj \          |  |      Nj - u
 *             \  g  /          u = 0
 *
 * ADZE 1.0 evaluated that product from scratch inside the g sweep, so the
 * sweep cost sum_g g = G^2/2 divisions per (grouping, allele, locus) instead
 * of G.  The recurrence
 *
 *      Q(g) = Q(g-1) * (Nj - Nji - (g-1)) / (Nj - (g-1))
 *
 * performs exactly the same multiplications in the same order as the loop it
 * replaces, so every cached value is bit-identical to what 1.0 computed.
 *
 * Entries for g > Nj are left at zero and never read: the callers guard those
 * cases with the -9 sentinel, exactly as 1.0 did.  Not writing them also keeps
 * the recurrence away from a zero denominator.
 */
void buildQTable(Population pop[], int numDivs, int locus, int numAlleles,
		 int gMax, int gStride, vector<double>& q)
{
  q.assign(size_t(numDivs) * numAlleles * gStride, 0.0);

  for(int p = 0; p < numDivs; p++)
    {
      const int Nj = pop[p].getNj(locus);
      const int gTop = (gMax < Nj) ? gMax : Nj;

      for(int i = 0; i < numAlleles; i++)
	{
	  const int Nji = pop[p].getNji(i,locus);
	  double* qpi = &q[(size_t(p) * numAlleles + i) * gStride];
	  double Q = 1;

	  for(int g = 1; g <= gTop; g++)
	    {
	      Q *= double(Nj - Nji - (g-1))/double(Nj - (g-1));
	      qpi[g] = Q;
	    }
	}
    }

  return;
}

void calcAllPgComb(Population pop[], int numDivs, int k, const ParamSet &param,
		   bool full_comb, string comb_out)
{
  const int numLoci = param.loci.val;
  const int tot_m = int(nCk(numDivs,k));

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

  //Largest g the sweep reaches; see calcAllPgs for why this replaces breakG.
  int minNjAll = 0;
  for(int locus = 0; locus < numLoci; locus++)
    {
      for(int p = 0; p < numDivs; p++)
	{
	  int Nj = pop[p].getNj(locus);
	  if((locus == 0 && p == 0) || Nj < minNjAll) minNjAll = Nj;
	}
    }

  int gLast = param.g.val;
  if(minNjAll < gLast) gLast = minNjAll;
  if(gLast < 2) gLast = 2;
  const int gStride = gLast + 1;

  gsl_combination* c = gsl_combination_calloc(numDivs,k);
  gsl_combination_init_first(c);

  ProgressBar bar(&cout,double(tot_m)*(gLast-1)*numLoci,BARLEN[min(k-1,3)]);
  if(param.pp.val)
    {
      bar.init();
    }

  vector<string> names(k);
  vector<char> inComb(numDivs,0);
  vector<double> pgcomb(size_t(gStride) * numLoci, 0.0); //[g][locus]
  vector<double> q;

  do
    {
      /*
       * Both of these are properties of the combination, not of a locus, an
       * allele or a sample size.  1.0 rebuilt the name array inside the
       * innermost allele loop -- one std::string copy per allele per g per
       * locus -- and asked isIn() to rescan the combination for every
       * non-member grouping at every allele.  Hoisting them out of three loop
       * levels leaves a membership mask that answers in constant time.
       */
      for(int j = 0; j < numDivs; j++) inComb[j] = 0;
      for(size_t j = 0; j < gsl_combination_k(c); j++)
	{
	  int member = int(gsl_combination_get(c,j));
	  names[j] = pop[member].getName();
	  inComb[member] = 1;
	}
      const string all_names = combineNames(&names[0],k);

      for(int locus = 0; locus < numLoci; locus++)
	{
	  const int numAlleles = pop[0].getNjiColLength(locus);
	  buildQTable(pop,numDivs,locus,numAlleles,gLast,gStride,q);

	  for(int g = 2; g <= gLast; g++)
	    {
	      double pg = 0;

	      //Sum over all alleles
	      for(int i = 0; i < numAlleles; i++)
		{
		  double P = 1, Q = 1;

		  //Calc P's over the groupings in the combination
		  for(size_t j = 0; j < gsl_combination_k(c); j++)
		    {
		      int member = int(gsl_combination_get(c,j));
		      P *= (1-q[(size_t(member)*numAlleles + i)*gStride + g]);
		    }

		  //Calc Q's over the groupings outside it
		  for(int j_p = 0; j_p < numDivs; j_p++)
		    {
		      if(!inComb[j_p])
			{
			  Q *= q[(size_t(j_p)*numAlleles + i)*gStride + g];
			}
		    }

		  //Multiply together and add to total
		  pg += (P*Q);
		}

	      pgcomb[size_t(g)*numLoci + locus] = pg;
	    }

	  if(param.pp.val) bar.adv(gLast-1);
	}

      for(int g = 2; g <= gLast; g++)
	{
	  Stats comb_stats;
	  comb_stats.putData(&pgcomb[size_t(g)*numLoci],numLoci);
	  comb_stats.calcAvg();
	  comb_stats.calcVar();
	  comb_stats.calcStdErr();

	  comb_stats.printStats(reg_out,all_names,g);

	  if(full_comb)
	    {
	      comb_stats.printData(full_out,all_names,g);
	    }
	}

      reg_out << endl;
      full_out << endl;
    }while(gsl_combination_next(c) == GSL_SUCCESS);
  
  if(param.pp.val) bar.done();

  gsl_combination_free(c);

  return;
}

void calcAllPgs(Population pop[],int numDivs,const ParamSet &param,
		bool full_priv,string private_out)
{
  const int numLoci = param.loci.val;
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

  /*
   * Smallest Nj per locus (for the -9 sentinel) and over the whole dataset.
   *
   * 1.0 walked g upward and set a breakG flag as soon as any grouping at any
   * locus had Nj < g+1, leaving the sweep after finishing that g -- so the
   * largest g it ever evaluated was the smallest Nj in the dataset.  Computing
   * that bound up front gives the same set of g values without the flag, and
   * lets the per-locus result buffer be sized to the g values that are
   * actually reached.
   */
  vector<int> minNjLocus(numLoci,0);
  int minNjAll = 0;
  for(int locus = 0; locus < numLoci; locus++)
    {
      int m = pop[0].getNj(locus);
      for(int p = 1; p < numDivs; p++)
	{
	  int Nj = pop[p].getNj(locus);
	  if(Nj < m) m = Nj;
	}
      minNjLocus[locus] = m;
      if(locus == 0 || m < minNjAll) minNjAll = m;
    }

  int gLast = param.g.val;
  if(minNjAll < gLast) gLast = minNjAll;
  if(gLast < 2) gLast = 2; //1.0 always evaluated g = 2 at least once
  const int gStride = gLast + 1;

  vector<double> pg(size_t(gStride) * numLoci, 0.0); //[g][locus]
  vector<double> q;

  ProgressBar bar(&cout,double(numDivs)*(gLast-1)*numLoci,BARLEN[0]);
  if(param.pp.val)
    {
      bar.init();
    }

  /*
   * Calculate the private allelic richness
   *            m              J
   *            _             ___ 
   * __(j)     \   /        / | |       \ \
   * ||g  =    /_  \ Pijg * \ | | Qij'g / /
   *           i=1           j'=1
   *                         j'!=j
   */
  for(int j = 0; j < numDivs; j++)
    {
      for(int locus = 0; locus < numLoci; locus++)
	{
	  const int numAlleles = pop[j].getNjiColLength(locus);
	  buildQTable(pop,numDivs,locus,numAlleles,gLast,gStride,q);

	  for(int g = 2; g <= gLast; g++)
	    {
	      if(minNjLocus[locus] < g)
		{
		  pg[size_t(g)*numLoci + locus] = -9;
		  continue;
		}

	      double total = 0;
	      for(int i = 0; i < numAlleles; i++)
		{
		  double Q = 1;
		  for(int p = 0; p < numDivs; p++)
		    {
		      if(p != j) Q *= q[(size_t(p)*numAlleles + i)*gStride + g];
		    }

		  double P = 1 - q[(size_t(j)*numAlleles + i)*gStride + g];
		  total += P*Q;
		}

	      pg[size_t(g)*numLoci + locus] = total;
	    }

	  if(param.pp.val) bar.adv(gLast-1);
	}

      for(int g = 2; g <= gLast; g++)
	{
	  Stats pg_stats;
	  pg_stats.putData(&pg[size_t(g)*numLoci],numLoci);
	  pg_stats.calcAvg();
	  pg_stats.calcVar();
	  pg_stats.calcStdErr();

	  pg_stats.printStats(pg_out,pop[j].getName(),g);

	  if(full_priv)
	    {
	      pg_stats.printData(pg_full_out,pop[j].getName(),g);
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
  const int numLoci = param.loci.val;
  const int gTop = param.g.val;
  const int gStride = gTop + 1;
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

  ProgressBar bar(&cout,double(numDivs)*(gTop-1)*numLoci,BARLEN[0]);
  if(param.pp.val)
    {
      bar.init();
    }

  vector<double> ag(size_t(gStride) * numLoci, 0.0); //[g][locus]
  vector<double> q;

  /*
   * Calculate the allelic richness
   *            m
   *            _
   *   (j)     \
   *  A g  =   /_  Pijg
   *           i=1
   *
   * One grouping and one locus at a time, so the Qjig table for that locus is
   * built once (by the recurrence in g) and read by every g, instead of being
   * recomputed from scratch for each g as in 1.0.  Loci where g exceeds the
   * grouping's sample size keep the -9 sentinel that calcAg returned.
   */
  for(int j = 0; j < numDivs; j++)
    {
      for(int locus = 0; locus < numLoci; locus++)
	{
	  const int numAlleles = pop[j].getNjiColLength(locus);
	  const int Nj = pop[j].getNj(locus);

	  buildQTable(&pop[j],1,locus,numAlleles,gTop,gStride,q);

	  for(int g = 2; g <= gTop; g++)
	    {
	      if(g > Nj)
		{
		  ag[size_t(g)*numLoci + locus] = -9;
		  continue;
		}

	      double total = 0;
	      for(int i = 0; i < numAlleles; i++)
		{
		  total += 1 - q[size_t(i)*gStride + g];
		}

	      ag[size_t(g)*numLoci + locus] = total;
	    }

	  if(param.pp.val) bar.adv(gTop-1);
	}

      for(int g = 2; g <= gTop; g++)
	{
	  Stats ag_stats;
	  ag_stats.putData(&ag[size_t(g)*numLoci],numLoci);
	  ag_stats.calcAvg();
	  ag_stats.calcVar();
	  ag_stats.calcStdErr();

	  ag_stats.printStats(ag_out,pop[j].getName(),g);

	  if(full_rich)
	    {
	      ag_stats.printData(ag_full_out,pop[j].getName(),g);
	    }
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


/*
 * Calculates all Nji's for every locus
 * Stores in Population objects
 *
 */
void calcNji(Population pop[],int numDivs,const string& missing)
{
  int numLoci = pop[0].getNumLoci();
  int emptyLoci = 0;

  /*
   * ADZE 1.0 kept a map<int,string> of index -> allele label and, to ask
   * whether a label had been seen, walked the whole map (it could not even
   * break early, because the same loop supplied the next free index).  It then
   * walked the map again for every genotype to find the bin to increment, so
   * binning cost O(loci * rows * alleles) string comparisons where O(loci *
   * rows) suffices.  A hash index from label to dense bin does the same work
   * in one pass; allele indices are still assigned in order of first
   * appearance, scanning groupings then rows, so Nji column order -- and
   * therefore the summation order of every downstream statistic -- is
   * unchanged.
   */
  unordered_map<string,int> alleleIndex;
  vector< vector<int> > Nji(numDivs);

  for(int locus = 0; locus < numLoci; locus++)
    {
      alleleIndex.clear();
      for(int j = 0; j < numDivs; j++) Nji[j].clear();

      for(int j = 0; j < numDivs; j++)
	{
	  int numRows = pop[j].getNumRows();
	  for(int row = 0; row < numRows; row++)
	    {
	      const string& allele = pop[j].getDataElement(row,locus);
	      if(allele.compare(missing) == 0) continue; //missing data

	      pair<unordered_map<string,int>::iterator,bool> found =
		alleleIndex.insert(make_pair(allele,int(alleleIndex.size())));

	      if(found.second) //first sighting of this allele at this locus
		{
		  for(int jj = 0; jj < numDivs; jj++) Nji[jj].push_back(0);
		}

	      Nji[j][found.first->second]++;
	    }
	}

      int NjiColLength = int(alleleIndex.size());

      /*
       * A locus with no observed allele in any grouping is legal input (it
       * happens in merged panels) and 1.0 crashed on it: it decremented the
       * end() iterator of the empty allele map.  Record zero alleles instead.
       * Nj then comes out 0, so the locus reports -9 like any locus whose
       * sample size is smaller than g.
       */
      if(NjiColLength == 0) emptyLoci++;

      for(int j = 0; j < numDivs; j++)
	{
	  pop[j].setNjiColLength(NjiColLength,locus);
	  for(int i = 0; i < NjiColLength; i++)
	    {
	      pop[j].putNji(Nji[j][i],i,locus);
	    }
	}
    }

  if(emptyLoci > 0)
    {
      cout << "WARNING: " << emptyLoci
	   << ((emptyLoci == 1) ? " locus has" : " loci have")
	   << " no observed alleles in any grouping.\n"
	   << "         Such loci make every statistic undefined at every g; "
	   << "set TOLERANCE < 1 to drop them.\n";
    }

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
