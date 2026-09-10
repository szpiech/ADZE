#include "ADZE_main_tools.h"
#include <unistd.h>

/*
 * Progress and informational messages go to stderr, so that redirecting
 * stdout captures results only.  1.0 wrote errors, the progress bar and the
 * results to the same stream, which made "adze ... > log" put error text into
 * the log and the backspace-driven progress bar into it as binary noise.
 */
bool ADZE_QUIET = false;

ostream& adzelog()
{
  static ofstream sink; //never opened: writes are discarded
  return ADZE_QUIET ? static_cast<ostream&>(sink) : cerr;
}

bool stderrIsTerminal()
{
  return isatty(fileno(stderr)) ? true : false;
}


using namespace std;

/*
 * Parse a K_RANGE value: a comma- or space-separated list of tuple sizes and
 * inclusive ranges, e.g. "2", "1-3", "1,3,5-7".  Returns the sizes sorted and
 * deduplicated.
 *
 * 1.0 hand-rolled this with a fixed 10-byte buffer, a range convention that
 * pushed negative sentinels to be expanded afterwards against an
 * uninitialized `prev`, and single-character whitespace skipping that could
 * not handle "1, 2".
 */
list<int> parseKVals(string str)
{
  string err = "ERROR: \"";
  err += str;
  err += "\" not a valid K_RANGE definition.\n";

  list<int> nums;

  size_t i = 0;
  const size_t n = str.size();

  while(i < n)
    {
      if(str[i] == '#') break; //trailing comment
      if(isspace((unsigned char)str[i]) || str[i] == ',') { i++; continue; }
      if(!isdigit((unsigned char)str[i])) throw err;

      long lo = 0;
      while(i < n && isdigit((unsigned char)str[i]))
	{
	  lo = lo*10 + (str[i]-'0');
	  if(lo > 1000000) throw err;
	  i++;
	}

      long hi = lo;
      if(i < n && str[i] == '-')
	{
	  i++;
	  if(i >= n || !isdigit((unsigned char)str[i])) throw err;
	  hi = 0;
	  while(i < n && isdigit((unsigned char)str[i]))
	    {
	      hi = hi*10 + (str[i]-'0');
	      if(hi > 1000000) throw err;
	      i++;
	    }
	  if(hi < lo) throw err;
	}

      for(long v = lo; v <= hi; v++) nums.push_back(int(v));

      if(i < n && !(isspace((unsigned char)str[i]) || str[i] == ',' || str[i] == '#'))
	{
	  throw err;
	}
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
  cerr << "ERROR: " << msg << "\n";
  BAD_PARAM x;
  throw x;
}

//Split a comma/space-separated list into names.
static void splitList(const string& text, vector<string>& out)
{
  out.clear();
  string cur;
  for(size_t i = 0; i < text.size(); i++)
    {
      const char c = text[i];
      if(c == ',' || isspace((unsigned char)c))
	{
	  if(!cur.empty()) { out.push_back(cur); cur.clear(); }
	}
      else cur += c;
    }
  if(!cur.empty()) out.push_back(cur);
  return;
}

/*
 * Read the data file and return a freshly allocated array of numDivs
 * Population objects, one per grouping, with allele counts, sample sizes and
 * missing-data tallies already filled in.
 *
 * LOCI, NON_DATA_COLS and DATA_LINES are measured from the file rather than
 * demanded from the user: the locus-name row gives the locus count, the width
 * of the first data row gives the number of label columns, and the data rows
 * count themselves.  A value the user did declare is checked against what is
 * there, and a disagreement is reported as a warning against the measurement,
 * not as a fatal error about the declaration.
 *
 * Throws BAD_FILE if the file cannot be opened and BAD_PARAM if the file's
 * shape is internally inconsistent.
 */
Population* readDataset(ParamSet& p, vector<string>& groupNames, int& numDivs)
{
  ifstream in(p.dfile.val.c_str());
  if(in.fail())
    {
      cerr << "ERROR: could not open " << p.dfile.val << "\n";
      BAD_FILE x;
      throw x;
    }

  vector<string> keepList, dropList;
  splitList(p.pops.val,keepList);
  splitList(p.expops.val,dropList);

  string line;
  vector<Field> fields;

  //Row 1 of the non-data rows carries the locus names.
  if(!getline(in,line)) badData("no locus-name row in " + p.dfile.val + ".");
  tokenize(line,fields);

  const int foundLoci = int(fields.size());
  if(foundLoci < 1) badData("no locus names in the first row of " + p.dfile.val + ".");

  if(p.loci.set && p.loci.val != foundLoci)
    {
      adzelog() << "WARNING: LOCI says " << p.loci.val << " but "
		<< p.dfile.val << " has " << foundLoci
		<< " locus names; using " << foundLoci << ".\n";
    }
  p.loci.val = foundLoci;
  const int declaredLoci = foundLoci;

  Accumulator acc;
  acc.locusName.reserve(fields.size());
  for(size_t l = 0; l < fields.size(); l++)
    {
      acc.locusName.push_back(string(fields[l].first,fields[l].second));
    }
  acc.locus.resize(acc.locusName.size());

  //Remaining non-data rows are ignored, exactly as in 1.0.
  for(int skip = 1; skip < p.nd_rows.val; skip++) getline(in,line);

  int ndCols = p.nd_cols.set ? p.nd_cols.val : 0;
  int groupCol = -1;
  int expected = 0;

  string token;
  long long dataRows = 0, physicalRows = 0, keptRows = 0;

  while(getline(in,line))
    {
      physicalRows++;
      tokenize(line,fields);
      if(fields.empty()) continue; //blank separator line

      if(expected == 0)
	{
	  //First data row fixes the layout.
	  const int found = int(fields.size()) - declaredLoci;
	  if(found < 1)
	    {
	      ostringstream m;
	      m << "the first data row of " << p.dfile.val << " has "
		<< fields.size() << " columns, which leaves no room for "
		<< declaredLoci << " loci plus at least one label column.";
	      badData(m.str());
	    }

	  if(p.nd_cols.set && ndCols != found)
	    {
	      adzelog() << "WARNING: NON_DATA_COLS says " << ndCols
			<< " but the data rows leave room for " << found
			<< "; using " << found << ".\n";
	    }
	  ndCols = found;
	  p.nd_cols.val = ndCols;

	  if(!p.sort_by.set) p.sort_by.val = ndCols;
	  groupCol = p.sort_by.val - 1;

	  if(groupCol < 0 || groupCol >= ndCols)
	    {
	      ostringstream m;
	      m << "GROUP_BY_COL " << p.sort_by.val << " is not one of the "
		<< ndCols << " label columns in " << p.dfile.val << ".";
	      badData(m.str());
	    }

	  expected = ndCols + declaredLoci;
	}

      if(int(fields.size()) != expected)
	{
	  ostringstream m;
	  m << "expected " << expected << " columns at data line " << physicalRows
	    << " in " << p.dfile.val << " but found " << fields.size()
	    << ". Check NON_DATA_ROWS, or whether the row is truncated.";
	  badData(m.str());
	}

      dataRows++;

      token.assign(fields[groupCol].first,fields[groupCol].second);

      if(!keepList.empty() &&
	 find(keepList.begin(),keepList.end(),token) == keepList.end()) continue;
      if(!dropList.empty() &&
	 find(dropList.begin(),dropList.end(),token) != dropList.end()) continue;

      keptRows++;
      const int g = acc.group(token);
      const int rowInGroup = acc.groupRows[g]++;
      const long long firstKey = (long long)(g) * 4294967296LL + rowInGroup;

      for(int l = 0; l < declaredLoci; l++)
	{
	  LocusTally& t = acc.locus[l];
	  const Field& f = fields[ndCols + l];

	  token.assign(f.first,f.second);
	  if(token.compare(p.miss.val) == 0)
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

  if(p.dlines.set && p.dlines.val != int(dataRows))
    {
      adzelog() << "WARNING: DATA_LINES says " << p.dlines.val << " but "
		<< p.dfile.val << " has " << dataRows << " data rows; using "
		<< dataRows << ".\n";
    }
  p.dlines.val = int(dataRows);

  if(dataRows == 0) badData("no data rows in " + p.dfile.val + ".");

  numDivs = int(acc.groupName.size());
  if(numDivs == 0)
    {
      badData("no grouping in " + p.dfile.val + " survived --pops/--exclude-pops.");
    }

  if(keptRows != dataRows)
    {
      adzelog() << "Using " << keptRows << " of " << dataRows
		<< " gene copies in " << numDivs
		<< (numDivs == 1 ? " grouping.\n" : " groupings.\n");
    }

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
      adzelog() << "WARNING: " << emptyLoci
		<< ((emptyLoci == 1) ? " locus has" : " loci have")
		<< " no observed alleles in any grouping.\n"
		<< "         Such loci make every statistic undefined at "
		<< "every g; set TOLERANCE < 1 to drop them.\n";
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

  ProgressBar bar(&adzelog(),size*numDivs,BARLEN[0]);
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

  const string report = pop[0].deletedSummary();
  adzelog() << report;
  lout << report;
  pop[0].printDeleted(lout);

  return;
}




/*
 * Elapsed wall-clock time since the program started, as d:h:m:s.
 *
 * 1.0 reported clock()/CLOCKS_PER_SEC, an integer division that truncated
 * everything to whole seconds -- which is why every phase in the shipped
 * example reports 0:0:0:0 -- and measured CPU rather than elapsed time, so it
 * would also have counted every thread separately once the loops are
 * parallelized.
 */
double displayTime(ostream& out)
{
  static const chrono::steady_clock::time_point start =
    chrono::steady_clock::now();

  const double sec =
    chrono::duration<double>(chrono::steady_clock::now() - start).count();

  double rest = sec;
  const ios::fmtflags flags = out.flags();
  const streamsize prec = out.precision();

  out.unsetf(ios::floatfield);
  out.unsetf(ios::showpoint);

  out << long(rest/86400.0) << ":";
  rest -= floor(rest/86400.0)*86400.0;
  out << long(rest/3600.0) << ":";
  rest -= floor(rest/3600.0)*3600.0;
  out << long(rest/60.0) << ":";
  rest -= floor(rest/60.0)*60.0;
  out << fixed << setprecision(2) << rest;

  out.flags(flags);
  out.precision(prec);

  return sec;
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


/*
 * Label for a tuple of groupings.  1.0 joined the names with a space, which
 * makes the label indistinguishable from the surrounding space-separated
 * fields; in tab-separated output the names are joined with a comma so the
 * label stays one field.
 */
string combineNames(string names[],int k,char sep)
{
  string tmp = "";

  for(int i = 0; i < k; i++)
    {
      tmp += names[i];
      if(i != k-1) tmp += sep;
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

/*
 * Every k-subset of the groupings, in lexicographic order -- the order 1.0's
 * gsl_combination walk produced.
 *
 * This is the only thing ADZE 1.0 used GSL for, so generating the subsets here
 * removes the program's sole external dependency.
 */
void buildKTuples(int numDivs, int k, vector< vector<int> >& out)
{
  out.clear();
  if(k < 1 || k > numDivs) return;

  vector<int> c(k);
  for(int i = 0; i < k; i++) c[i] = i;

  while(1)
    {
      out.push_back(c);

      //Advance the rightmost index that has room, then repack those after it.
      int i = k-1;
      while(i >= 0 && c[i] == numDivs-k+i) i--;
      if(i < 0) break;

      c[i]++;
      for(int j = i+1; j < k; j++) c[j] = c[j-1]+1;
    }

  return;
}

/*
 * Read named tuples: one tuple per line, grouping names separated by commas or
 * whitespace.  Lets an analysis ask for the handful of groupings it cares
 * about instead of enumerating all 2^J subsets to read three of them.
 */
bool readTupleFile(const string& file, Population pop[], int numDivs,
		   vector< vector<int> >& out)
{
  ifstream in(file.c_str());
  if(in.fail())
    {
      cerr << "ERROR: could not open tuple file " << file << "\n";
      return 0;
    }

  out.clear();
  string line;
  long lineNo = 0;
  bool ok = 1;

  while(getline(in,line))
    {
      lineNo++;
      size_t hash = line.find('#');
      if(hash != string::npos) line = line.substr(0,hash);

      vector<string> names;
      splitList(line,names);
      if(names.empty()) continue;

      vector<int> tuple;
      for(size_t i = 0; i < names.size(); i++)
	{
	  int found = -1;
	  for(int j = 0; j < numDivs; j++)
	    {
	      if(pop[j].getName().compare(names[i]) == 0) found = j;
	    }

	  if(found < 0)
	    {
	      cerr << "ERROR: " << file << ":" << lineNo << ": no grouping named \""
		   << names[i] << "\" in the data.\n";
	      ok = 0;
	    }
	  else tuple.push_back(found);
	}

      sort(tuple.begin(),tuple.end());
      tuple.erase(unique(tuple.begin(),tuple.end()),tuple.end());
      if(!tuple.empty()) out.push_back(tuple);
    }

  in.close();

  if(out.empty())
    {
      cerr << "ERROR: no tuples found in " << file << "\n";
      ok = 0;
    }

  return ok;
}

/*
 * Private allelic richness of grouping tuples.
 *
 *  __(T)      m /                          \
 *  ||        __ |  | |                | |  |
 *    g   =   \  |  | | Pijg  *  | |   Qij'g |
 *            /_ |  |j in T      |  |j' not  |
 *           i=1 \                          /
 *
 * 1.0 recomputed every Qjig from scratch for each of the C(J,k) combinations,
 * rebuilt the grouping-name array inside the innermost allele loop, and
 * rescanned the combination to test membership for every non-member grouping
 * at every allele.  The Qjig table is now built once per locus per tuple, and
 * the names and a constant-time membership mask once per tuple.
 */
void calcPgTuples(Population pop[], int numDivs,
		  const vector< vector<int> >& tuples, const ParamSet &param,
		  bool full_comb, string comb_out, bool namedTuples)
{
  const int numLoci = param.loci.val;
  const int tot_m = int(tuples.size());
  if(tot_m == 0) return;

  size_t widest = 0;
  for(size_t t = 0; t < tuples.size(); t++)
    {
      if(tuples[t].size() > widest) widest = tuples[t].size();
    }

  ofstream full_out, reg_out;

  if(full_comb)
    {
      string name;
      name = nameCreate(comb_out,"_fulldata");
      full_out.open(name.c_str());

      if(namedTuples) full_out << "TUPLE ";
      else
	{
	  for(size_t i = 1; i <= widest; i++) full_out << "POP_GROUPING" << i << " ";
	}

      full_out << "G NUM_LOCI ";
      for(int l = 0; l < numLoci; l++)
	{
	  full_out << pop[0].getLocusName(l) << " ";
	}
      
      full_out << "MEAN VAR STD_ERR\n";
    }

  reg_out.open(comb_out.c_str());
  if(param.tsv.val)
    {
      reg_out << "TUPLE\tG\tNUM_LOCI\tMEAN\tVAR\tSTD_ERR\n";
    }

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

  ProgressBar bar(&adzelog(),double(tot_m)*(gLast-1)*numLoci,
		  BARLEN[min(int(widest)-1,3)]);
  if(param.pp.val)
    {
      bar.init();
    }

  vector<string> names(widest);
  vector<char> inTuple(numDivs,0);
  vector<double> pgcomb(size_t(gStride) * numLoci, 0.0); //[g][locus]

  for(int m = 0; m < tot_m; m++)
    {
      const vector<int>& tuple = tuples[m];
      const int k = int(tuple.size());

      /*
       * Properties of the tuple, not of a locus, an allele or a sample size:
       * hoisted out of three loop levels.
       */
      for(int j = 0; j < numDivs; j++) inTuple[j] = 0;
      for(int j = 0; j < k; j++)
	{
	  names[j] = pop[tuple[j]].getName();
	  inTuple[tuple[j]] = 1;
	}
      const string all_names =
	combineNames(&names[0],k,param.tsv.val ? ',' : ' ');

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
      vector<double> q; //one per thread

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
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

		  //Calc P's over the groupings in the tuple
		  for(int j = 0; j < k; j++)
		    {
		      P *= (1-q[(size_t(tuple[j])*numAlleles + i)*gStride + g]);
		    }

		  //Calc Q's over the groupings outside it
		  for(int j_p = 0; j_p < numDivs; j_p++)
		    {
		      if(!inTuple[j_p])
			{
			  Q *= q[(size_t(j_p)*numAlleles + i)*gStride + g];
			}
		    }

		  //Multiply together and add to total
		  pg += (P*Q);
		}

	      pgcomb[size_t(g)*numLoci + locus] = pg;
	    }

	  if(param.pp.val)
	    {
#ifdef _OPENMP
#pragma omp critical(progress)
#endif
	      bar.adv(gLast-1);
	    }
	}
      } //end parallel region

      for(int g = 2; g <= gLast; g++)
	{
	  Stats comb_stats;
	  comb_stats.putData(&pgcomb[size_t(g)*numLoci],numLoci);
	  comb_stats.calcAvg();
	  comb_stats.calcVar();
	  comb_stats.calcStdErr();

	  comb_stats.printStats(reg_out,all_names,g,param.tsv.val);

	  if(full_comb)
	    {
	      comb_stats.printData(full_out,all_names,g);
	    }
	}

      if(!param.tsv.val) reg_out << endl;
      full_out << endl;
    }

  if(param.pp.val) bar.done();

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
  if(param.tsv.val) pg_out << "POP_GROUPING\tG\tNUM_LOCI\tMEAN\tVAR\tSTD_ERR\n";

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

  ProgressBar bar(&adzelog(),double(numDivs)*(gLast-1)*numLoci,BARLEN[0]);
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
      //Independent per locus; see calcAllAgs on why threading cannot move a
      //reported value.
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
      vector<double> q; //one per thread

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
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

	  if(param.pp.val)
	    {
#ifdef _OPENMP
#pragma omp critical(progress)
#endif
	      bar.adv(gLast-1);
	    }
	}
      } //end parallel region

      for(int g = 2; g <= gLast; g++)
	{
	  Stats pg_stats;
	  pg_stats.putData(&pg[size_t(g)*numLoci],numLoci);
	  pg_stats.calcAvg();
	  pg_stats.calcVar();
	  pg_stats.calcStdErr();

	  pg_stats.printStats(pg_out,pop[j].getName(),g,param.tsv.val);

	  if(full_priv)
	    {
	      pg_stats.printData(pg_full_out,pop[j].getName(),g);
	    }
	}

      if(!param.tsv.val) pg_out << endl;
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
  if(param.tsv.val) ag_out << "POP_GROUPING\tG\tNUM_LOCI\tMEAN\tVAR\tSTD_ERR\n";

  ProgressBar bar(&adzelog(),double(numDivs)*(gTop-1)*numLoci,BARLEN[0]);
  if(param.pp.val)
    {
      bar.init();
    }

  vector<double> ag(size_t(gStride) * numLoci, 0.0); //[g][locus]

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
      /*
       * Loci are independent, and each writes only its own column of ag, so
       * the results do not depend on the thread count: the reductions over
       * loci happen afterwards, in ascending locus order, inside Stats.
       */
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
      vector<double> q; //one per thread

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
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

	  if(param.pp.val)
	    {
#ifdef _OPENMP
#pragma omp critical(progress)
#endif
	      bar.adv(gTop-1);
	    }
	}
      } //end parallel region

      for(int g = 2; g <= gTop; g++)
	{
	  Stats ag_stats;
	  ag_stats.putData(&ag[size_t(g)*numLoci],numLoci);
	  ag_stats.calcAvg();
	  ag_stats.calcVar();
	  ag_stats.calcStdErr();

	  ag_stats.printStats(ag_out,pop[j].getName(),g,param.tsv.val);

	  if(full_rich)
	    {
	      ag_stats.printData(ag_full_out,pop[j].getName(),g);
	    }
	}

      if(!param.tsv.val) ag_out << endl;
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
