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




/*
 * ------------------------------------------------------------------------
 * Reading the data file
 * ------------------------------------------------------------------------
 *
 * ADZE 1.0 read the file three times: checkDatafile validated the column
 * counts, getLociNames plus getDivLines discovered the groupings and counted
 * their rows, and readData stored every genotype as a std::string in a
 * numLoci x rows matrix per grouping.  Each field was extracted with
 * `stream >> string`, at roughly 230 ns per genotype, and the resident matrix
 * cost 34 bytes per gene copy -- 8.5x the size of the input file.
 *
 * One pass does all of it.  Fields are tokenized in place out of a reused line
 * buffer; allele labels are interned per locus to a dense slot index; and the
 * only things retained are the per-locus allele counts, the missing-data
 * tallies, and the row counts -- which is all the estimators ever read.  The
 * genotype matrix is never materialized.
 *
 * Allele slots are renumbered at the end of the pass into the order 1.0
 * assigned them (first appearance scanning groupings in discovery order, then
 * rows within a grouping), so Nji column order -- and therefore the summation
 * order of every statistic -- is unchanged.
 */


//Orders allele slots by where each was first seen, so that Nji columns come out
//in the same order ADZE 1.0 produced.
struct FirstSeenLess
{
  const vector<long long>& key;
  FirstSeenLess(const vector<long long>& k) : key(k) {}
  bool operator()(int a, int b) const { return key[a] < key[b]; }
};

typedef pair<const char*,size_t> Field;

//Split a line on whitespace, into spans pointing back into the line itself.
static void tokenize(const string& line, vector<Field>& out)
{
  out.clear();
  const char* p = line.c_str();
  const char* end = p + line.size();

  while(p < end)
    {
      while(p < end && isspace((unsigned char)*p)) p++;
      if(p >= end) break;
      const char* start = p;
      while(p < end && !isspace((unsigned char)*p)) p++;
      out.push_back(Field(start,size_t(p-start)));
    }

  return;
}

namespace {

  //Per-locus allele bookkeeping while the file streams past.
  struct LocusTally
  {
    unordered_map<string,int> slotOf;  //allele label -> dense slot
    vector<long long> firstSeen;       //per slot: (grouping, row) of first sighting
    vector<int> count;                 //slot-major, stride groupCap
    vector<int> missing;               //per grouping
  };

  struct Accumulator
  {
    vector<string> groupName;
    vector<int> groupRows;
    unordered_map<string,int> groupOf;
    vector<string> locusName;
    vector<LocusTally> locus;
    int groupCap;

    Accumulator() : groupCap(4) {}

    int group(const string& label)
    {
      unordered_map<string,int>::iterator it = groupOf.find(label);
      if(it != groupOf.end()) return it->second;

      int index = int(groupName.size());
      groupOf.insert(make_pair(label,index));
      groupName.push_back(label);
      groupRows.push_back(0);

      if(index >= groupCap)
	{
	  //Widen every locus's count block. Groupings are discovered early in
	  //practice, and doubling bounds the total reshuffling to O(loci *
	  //alleles * groupings).
	  int newCap = groupCap * 2;
	  for(size_t l = 0; l < locus.size(); l++)
	    {
	      LocusTally& t = locus[l];
	      int slots = int(t.firstSeen.size());
	      vector<int> wider(size_t(slots) * newCap, 0);
	      for(int sl = 0; sl < slots; sl++)
		{
		  for(int g = 0; g < groupCap; g++)
		    {
		      wider[size_t(sl)*newCap + g] = t.count[size_t(sl)*groupCap + g];
		    }
		}
	      t.count.swap(wider);
	    }
	  groupCap = newCap;
	}

      for(size_t l = 0; l < locus.size(); l++) locus[l].missing.push_back(0);

      return index;
    }
  };

} //anonymous namespace

static void badData(const string& msg)
{
  cout << "ERROR: " << msg << "\n";
  BAD_PARAM x;
  throw x;
}

/*
 * Read the data file and return a freshly allocated array of numDivs
 * Population objects, one per grouping, with allele counts, sample sizes and
 * missing-data tallies already filled in.
 *
 * Throws BAD_FILE if the file cannot be opened and BAD_PARAM if its shape
 * disagrees with the declared LOCI, DATA_LINES or NON_DATA_COLS.
 */
Population* readDataset(ParamSet& p, vector<string>& groupNames, int& numDivs)
{
  ifstream in(p.dfile.val.c_str());
  if(in.fail())
    {
      cout << "ERROR: Could not open " << p.dfile.val << "\n";
      BAD_FILE x;
      throw x;
    }

  const int declaredLoci = p.loci.val;
  const int ndCols = p.nd_cols.val;
  const int groupCol = p.sort_by.val - 1;
  const string& missingLabel = p.miss.val;

  string line;
  vector<Field> fields;

  //Row 1 of the non-data rows carries the locus names.
  if(!getline(in,line)) badData("no locus-name row in " + p.dfile.val + ".");
  tokenize(line,fields);

  if(int(fields.size()) != declaredLoci)
    {
      ostringstream m;
      m << "Expected " << declaredLoci << (declaredLoci == 1 ? " locus " : " loci ")
	<< "in " << p.dfile.val << " but found " << fields.size() << ".";
      badData(m.str());
    }

  Accumulator acc;
  acc.locusName.reserve(fields.size());
  for(size_t l = 0; l < fields.size(); l++)
    {
      acc.locusName.push_back(string(fields[l].first,fields[l].second));
    }
  acc.locus.resize(acc.locusName.size());

  //Remaining non-data rows are ignored, exactly as in 1.0.
  for(int skip = 1; skip < p.nd_rows.val; skip++) getline(in,line);

  const int expected = ndCols + declaredLoci;
  string token;
  long long dataRows = 0, physicalRows = 0;

  while(getline(in,line))
    {
      physicalRows++;
      tokenize(line,fields);
      if(fields.empty()) continue; //blank separator line

      if(int(fields.size()) != expected)
	{
	  ostringstream m;
	  m << "Expected " << expected << " columns at data line " << physicalRows
	    << " in " << p.dfile.val << " but found " << fields.size()
	    << ". Possible bad NON_DATA_COLS, NON_DATA_ROWS, or LOCI value.";
	  badData(m.str());
	}

      if(dataRows >= p.dlines.val)
	{
	  ostringstream m;
	  m << "Expected " << p.dlines.val << " data lines in " << p.dfile.val
	    << " but found at least " << dataRows+1 << ".";
	  badData(m.str());
	}

      token.assign(fields[groupCol].first,fields[groupCol].second);
      const int g = acc.group(token);
      const int rowInGroup = acc.groupRows[g]++;
      const long long firstKey = (long long)(g) * 4294967296LL + rowInGroup;
      dataRows++;

      for(int l = 0; l < declaredLoci; l++)
	{
	  LocusTally& t = acc.locus[l];
	  const Field& f = fields[ndCols + l];

	  token.assign(f.first,f.second);
	  if(token.compare(missingLabel) == 0)
	    {
	      t.missing[g]++;
	      continue;
	    }

	  pair<unordered_map<string,int>::iterator,bool> found =
	    t.slotOf.insert(make_pair(token,int(t.firstSeen.size())));

	  if(found.second)
	    {
	      t.firstSeen.push_back(firstKey);
	      t.count.resize(t.count.size() + acc.groupCap, 0);
	    }

	  t.count[size_t(found.first->second)*acc.groupCap + g]++;
	}
    }

  in.close();

  if(dataRows != p.dlines.val)
    {
      ostringstream m;
      m << "Expected " << p.dlines.val << " data lines in " << p.dfile.val
	<< " but found " << dataRows << ".";
      badData(m.str());
    }

  numDivs = int(acc.groupName.size());
  groupNames = acc.groupName;

  Population* pop = new Population[numDivs];
  for(int j = 0; j < numDivs; j++)
    {
      pop[j].setLoci(declaredLoci);
      pop[j].setName(acc.groupName[j]);
      pop[j].setRows(acc.groupRows[j]);
      for(int l = 0; l < declaredLoci; l++)
	{
	  pop[j].setLocusName(acc.locusName[l],l);
	}
    }

  //Renumber slots into 1.0's first-appearance order and publish the counts.
  vector<int> order;
  int emptyLoci = 0;

  for(int l = 0; l < declaredLoci; l++)
    {
      LocusTally& t = acc.locus[l];
      const int slots = int(t.firstSeen.size());
      if(slots == 0) emptyLoci++;

      order.resize(slots);
      for(int sl = 0; sl < slots; sl++) order[sl] = sl;
      sort(order.begin(),order.end(),FirstSeenLess(t.firstSeen));

      for(int j = 0; j < numDivs; j++)
	{
	  pop[j].setNjiColLength(slots,l);
	  for(int i = 0; i < slots; i++)
	    {
	      pop[j].putNji(t.count[size_t(order[i])*acc.groupCap + j],i,l);
	    }
	  pop[j].putMissing(t.missing[j],l);
	}

      //Release this locus's bookkeeping as soon as it is published.
      unordered_map<string,int>().swap(t.slotOf);
      vector<int>().swap(t.count);
      vector<long long>().swap(t.firstSeen);
    }

  for(int j = 0; j < numDivs; j++) pop[j].sumNj();

  if(emptyLoci > 0)
    {
      cout << "WARNING: " << emptyLoci
	   << ((emptyLoci == 1) ? " locus has" : " loci have")
	   << " no observed alleles in any grouping.\n"
	   << "         Such loci make every statistic undefined at every g; "
	   << "set TOLERANCE < 1 to drop them.\n";
    }

  return pop;
}

void filterLoci(Population pop[],int numDivs, double tol, string file,
		bool pp)
{
  vector<char> toDelete(pop[0].getNumLoci(),0);

  for(int n = 0; n < numDivs; n++)
    {
      pop[n].recLociDelete(tol,toDelete);
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
