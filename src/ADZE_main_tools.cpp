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

/*
 * ------------------------------------------------------------------------
 * Locus coordinates
 * ------------------------------------------------------------------------
 */

int LocusMap::chromIndex(const string& name)
{
  for(size_t i = 0; i < chromName.size(); i++)
    {
      if(chromName[i] == name) return int(i);
    }
  chromName.push_back(name);
  return int(chromName.size()) - 1;
}

//Drop the flagged loci, keeping the table aligned with the surviving loci.
void LocusMap::compact(const vector<char>& del)
{
  if(chrom.empty()) return;

  size_t keep = 0;
  for(size_t l = 0; l < chrom.size() && l < del.size(); l++)
    {
      if(del[l]) continue;
      chrom[keep] = chrom[l];
      pos[keep] = pos[l];
      keep++;
    }

  chrom.resize(keep);
  pos.resize(keep);
  return;
}

/*
 * Windows are runs of consecutive loci, so each chromosome's loci must form
 * one block and positions must increase inside it.  Rather than sorting --
 * which would mean holding the whole dataset to reorder it -- say what is
 * wrong and let the caller sort the input, which is one bcftools call.
 */
string LocusMap::checkOrder(const vector<string>& locusName) const
{
  vector<char> seen(chromName.size(),0);
  ostringstream m;

  for(size_t l = 0; l < chrom.size(); l++)
    {
      const string name = (l < locusName.size()) ? locusName[l] : string("?");

      if(l == 0 || chrom[l] != chrom[l-1])
	{
	  if(seen[chrom[l]])
	    {
	      m << "chromosome " << chromName[chrom[l]] << " appears in more "
		<< "than one block, starting again at locus " << name
		<< ". Sort the input by chromosome and position.";
	      return m.str();
	    }
	  seen[chrom[l]] = 1;
	  continue;
	}

      if(pos[l] < pos[l-1])
	{
	  m << "locus " << name << " is at " << chromName[chrom[l]] << ":"
	    << pos[l] << ", behind the previous locus at " << pos[l-1]
	    << ". Sort the input by chromosome and position.";
	  return m.str();
	}
    }

  return "";
}

/*
 * Lay the windows out over the surviving loci, one chromosome at a time so no
 * window spans a boundary.
 *
 * Basepair windows are anchored at position 1 and advance by the step
 * regardless of where the loci fall, so the same region gives the same
 * intervals in every run and two datasets can be compared window by window.
 * Locus windows count surviving loci instead, so their width in basepairs
 * varies with locus density; their reported interval is the span of the loci
 * they hold.
 */
long long buildWindows(const LocusMap& lmap, const ParamSet& p,
		       vector<Window>& out)
{
  out.clear();

  const int n = int(lmap.size());
  long long sparse = 0;
  if(n == 0) return 0;

  int b = 0;
  while(b < n)
    {
      int e = b;
      while(e < n && lmap.chrom[e] == lmap.chrom[b]) e++;   //[b,e) = one chromosome

      if(p.win_bp.set)
	{
	  const long long width = p.win_bp.val;
	  const long long step = p.step_bp.val;
	  const long long firstPos = lmap.pos[b];
	  const long long lastPos = lmap.pos[e-1];

	  //First window that can reach the first locus, counting from 1.
	  long long k = (firstPos - 1) / step;
	  if(k < 0) k = 0;

	  int lo = b;
	  for(long long start = 1 + k*step; start <= lastPos; start += step)
	    {
	      const long long stop = start + width - 1;

	      while(lo < e && lmap.pos[lo] < start) lo++;
	      int hi = lo;
	      while(hi < e && lmap.pos[hi] <= stop) hi++;

	      if(hi == lo) continue;                //empty window: not reported
	      if(hi - lo < p.min_win_loci.val) { sparse++; continue; }

	      Window w;
	      w.chrom = lmap.chrom[b];
	      w.start = start;
	      w.end = stop;
	      w.first = lo;
	      w.last = hi;
	      out.push_back(w);
	    }
	}
      else
	{
	  const int width = int(p.win_loci.val);
	  const int step = int(p.step_loci.val);

	  for(int i = b; i < e; i += step)
	    {
	      int hi = i + width;
	      if(hi > e) hi = e;

	      if(hi - i < p.min_win_loci.val) { sparse++; continue; }

	      Window w;
	      w.chrom = lmap.chrom[b];
	      w.start = lmap.pos[i];
	      w.end = lmap.pos[hi-1];
	      w.first = i;
	      w.last = hi;
	      out.push_back(w);
	    }
	}

      b = e;
    }

  return sparse;
}

namespace {

  /*
   * A line source that reads plain text or gzip-compressed text.
   *
   * When the build has zlib, everything goes through gzFile: zlib reads an
   * uncompressed file transparently, so one path serves .stru, .vcf and
   * .vcf.gz alike, and bgzip output is just gzip with extra block structure.
   * Without zlib the plain path is used and a compressed file is refused with
   * a message rather than misparsed.
   */
  class LineSource
  {
  public:
    LineSource() : ok(false)
#ifdef ADZE_HAVE_ZLIB
      , gz(0)
#endif
    {}

    ~LineSource() { close(); }

    bool open(const string& path)
    {
      name = path;
#ifdef ADZE_HAVE_ZLIB
      gz = gzopen(path.c_str(),"rb");
      ok = (gz != 0);
      if(ok) gzbuffer(gz,1<<20);
#else
      plain.open(path.c_str());
      ok = !plain.fail();
#endif
      return ok;
    }

    bool next(string& line)
    {
#ifdef ADZE_HAVE_ZLIB
      line.clear();
      char buf[65536];
      while(gzgets(gz,buf,sizeof(buf)) != 0)
	{
	  line += buf;
	  if(!line.empty() && line[line.size()-1] == '\n')
	    {
	      line.erase(line.size()-1);
	      if(!line.empty() && line[line.size()-1] == '\r') line.erase(line.size()-1);
	      return true;
	    }
	}
      return !line.empty();
#else
      if(!getline(plain,line)) return false;
      if(!line.empty() && line[line.size()-1] == '\r') line.erase(line.size()-1);
      return true;
#endif
    }

    void close()
    {
#ifdef ADZE_HAVE_ZLIB
      if(gz) { gzclose(gz); gz = 0; }
#else
      if(plain.is_open()) plain.close();
#endif
      ok = false;
    }

  private:
    string name;
    bool ok;
#ifdef ADZE_HAVE_ZLIB
    gzFile gz;
#else
    ifstream plain;
#endif
  };

} //anonymous namespace

//True if the name ends in .gz or .bgz, case-insensitively.
static bool looksCompressed(const string& path)
{
  string low = path;
  for(size_t i = 0; i < low.size(); i++) low[i] = char(tolower((unsigned char)low[i]));
  return (low.size() > 3 && low.compare(low.size()-3,3,".gz") == 0) ||
         (low.size() > 4 && low.compare(low.size()-4,4,".bgz") == 0);
}

/*
 * Resolve FORMAT.  "auto" looks at the file name: anything ending .vcf, after
 * an optional compression suffix, is VCF; everything else is the
 * STRUCTURE-like layout ADZE has always read.
 */
bool wantsVCF(const string& format, const string& path)
{
  if(format == "vcf") return true;
  if(format == "structure") return false;

  string low = path;
  for(size_t i = 0; i < low.size(); i++) low[i] = char(tolower((unsigned char)low[i]));
  if(looksCompressed(low)) low.erase(low.rfind('.'));
  return low.size() > 4 && low.compare(low.size()-4,4,".vcf") == 0;
}

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
 * Read the STRUCTURE-like layout: a row of locus names, then one row per gene
 * copy, each beginning with label columns.
 *
 * LOCI, NON_DATA_COLS and DATA_LINES are measured from the file rather than
 * demanded from the user: the locus-name row gives the locus count, the width
 * of the first data row gives the number of label columns, and the data rows
 * count themselves.  A value the user did declare is checked against what is
 * there, and a disagreement is reported as a warning against the measurement,
 * not as a fatal error about the declaration.
 */
static void readStructureInto(ParamSet& p, Accumulator& acc, LineSource& in,
			      const vector<string>& keepList,
			      const vector<string>& dropList,
			      long long& dataRows, long long& keptRows)
{
  string line;
  vector<Field> fields;

  //Row 1 of the non-data rows carries the locus names.
  if(!in.next(line)) badData("no locus-name row in " + p.dfile.val + ".");
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

  acc.locusName.reserve(fields.size());
  for(size_t l = 0; l < fields.size(); l++)
    {
      acc.locusName.push_back(string(fields[l].first,fields[l].second));
    }
  acc.locus.resize(acc.locusName.size());

  //Remaining non-data rows are ignored, exactly as in 1.0.
  for(int skip = 1; skip < p.nd_rows.val; skip++) in.next(line);

  int ndCols = p.nd_cols.set ? p.nd_cols.val : 0;
  int groupCol = -1;
  int expected = 0;

  string token;
  long long physicalRows = 0;
  dataRows = 0;
  keptRows = 0;

  while(in.next(line))
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
  return;
}

/*
 * Read the sample-to-grouping map that VCF input needs, since a VCF carries
 * sample names but no population labels.  Two whitespace-separated columns,
 * sample then grouping; '#' starts a comment; further columns are ignored so
 * a fuller sample sheet can be used as it stands.
 *
 * Groupings are numbered in the order they first appear in this file, which
 * is therefore the order they appear in the output.
 */
static void readSampleMap(const string& file,
			  unordered_map<string,string>& groupOf,
			  vector<string>& order)
{
  LineSource in;
  if(!in.open(file))
    {
      cerr << "ERROR: could not open sample file " << file << "\n";
      BAD_FILE x;
      throw x;
    }

  string line;
  vector<Field> fields;
  long long lineNo = 0;

  while(in.next(line))
    {
      lineNo++;
      const size_t hash = line.find('#');
      if(hash != string::npos) line.erase(hash);
      tokenize(line,fields);
      if(fields.empty()) continue;

      if(fields.size() < 2)
	{
	  ostringstream m;
	  m << "line " << lineNo << " of " << file << " has one column; "
	    << "each line needs a sample name and a grouping name.";
	  badData(m.str());
	}

      const string sample(fields[0].first,fields[0].second);
      const string group(fields[1].first,fields[1].second);

      if(groupOf.find(sample) != groupOf.end())
	{
	  badData("sample " + sample + " appears twice in " + file + ".");
	}

      groupOf.insert(make_pair(sample,group));
      if(find(order.begin(),order.end(),group) == order.end()) order.push_back(group);
    }

  if(groupOf.empty()) badData("no sample assignments in " + file + ".");

  return;
}

/*
 * Read the locus-to-coordinate map that the STRUCTURE layout needs, since it
 * carries locus names but no positions.  Three whitespace-separated columns,
 * locus then chromosome then position; '#' starts a comment; further columns
 * are ignored.
 *
 * Every locus in the data must appear: a windowed statistic over loci whose
 * positions are unknown would be a window in name only, so a missing
 * coordinate is an error rather than a warning.
 */
static void readLocusMap(const string& file, const vector<string>& locusName,
			 LocusMap& lmap)
{
  LineSource in;
  if(!in.open(file))
    {
      cerr << "ERROR: could not open locus map " << file << "\n";
      BAD_FILE x;
      throw x;
    }

  unordered_map<string,size_t> row;      //locus name -> line in the map
  vector<string> chrom;
  vector<long long> pos;

  string line;
  vector<Field> fields;
  long long lineNo = 0;

  while(in.next(line))
    {
      lineNo++;
      const size_t hash = line.find('#');
      if(hash != string::npos) line.erase(hash);
      tokenize(line,fields);
      if(fields.empty()) continue;

      if(fields.size() < 3)
	{
	  ostringstream m;
	  m << "line " << lineNo << " of " << file << " has "
	    << fields.size() << (fields.size() == 1 ? " column" : " columns")
	    << "; each line needs a locus name, a chromosome and a position.";
	  badData(m.str());
	}

      const string name(fields[0].first,fields[0].second);
      if(row.find(name) != row.end())
	{
	  badData("locus " + name + " appears twice in " + file + ".");
	}

      const string posText(fields[2].first,fields[2].second);
      for(size_t i = 0; i < posText.size(); i++)
	{
	  if(!isdigit((unsigned char)posText[i]))
	    {
	      ostringstream m;
	      m << "position \"" << posText << "\" for locus " << name
		<< " in " << file << " is not a whole number.";
	      badData(m.str());
	    }
	}

      row.insert(make_pair(name,chrom.size()));
      chrom.push_back(string(fields[1].first,fields[1].second));
      pos.push_back(atoll(posText.c_str()));
    }

  //Every locus must be placed; name the first few that are not.
  vector<string> unplaced;
  for(size_t l = 0; l < locusName.size(); l++)
    {
      if(row.find(locusName[l]) == row.end()) unplaced.push_back(locusName[l]);
    }

  if(!unplaced.empty())
    {
      ostringstream m;
      m << unplaced.size() << " of " << locusName.size()
	<< " loci have no coordinate in " << file << ":";
      for(size_t i = 0; i < unplaced.size() && i < 5; i++) m << " " << unplaced[i];
      if(unplaced.size() > 5) m << " ...";
      m << ".";
      badData(m.str());
    }

  if(row.size() > locusName.size())
    {
      adzelog() << "WARNING: " << (row.size() - locusName.size())
		<< " loci named in " << file << " are not in the data.\n";
    }

  lmap.chrom.resize(locusName.size());
  lmap.pos.resize(locusName.size());
  for(size_t l = 0; l < locusName.size(); l++)
    {
      const size_t r = row[locusName[l]];
      lmap.chrom[l] = lmap.chromIndex(chrom[r]);
      lmap.pos[l] = pos[r];
    }

  return;
}

/*
 * Read VCF (optionally gzip- or bgzip-compressed).
 *
 * Each record is one locus and each allele index is one allele type, so REF is
 * allele 0 and the ALT alleles follow.  A sample contributes as many gene
 * copies as its GT field has alleles, which lets haploid and diploid records
 * mix; '.' is an uncalled copy.
 *
 * A sample's ploidy is the largest number of alleles any of its GT fields
 * holds, because that is how many gene copies it could have contributed.  A
 * record where it carries fewer -- a haploid call in a file that is otherwise
 * diploid -- therefore counts the difference as missing, which is what
 * TOLERANCE is asked about.
 *
 * LocusTally::missing holds observed calls during the pass and is converted to
 * missing counts at the end, once every sample's ploidy -- and therefore each
 * grouping's total gene copies -- is known.
 *
 * Allele slots are keyed the same way as in the STRUCTURE reader: by grouping,
 * then by position within the grouping, then by gene copy.  The two readers
 * therefore lay out Nji identically for the same genotypes, so every statistic
 * sums its terms in the same order and the results are bit-identical rather
 * than merely equal to within rounding (test/formats.py checks this).
 */
static void readVCFInto(ParamSet& p, Accumulator& acc, LineSource& in,
			const vector<string>& keepList,
			const vector<string>& dropList,
			long long& geneCopies, LocusMap& lmap)
{
  unordered_map<string,string> groupOfSample;
  vector<string> groupOrder;
  readSampleMap(p.samples.val,groupOfSample,groupOrder);

  //Groupings are registered up front, in sample-file order.
  for(size_t i = 0; i < groupOrder.size(); i++)
    {
      const string& g = groupOrder[i];
      if(!keepList.empty() &&
	 find(keepList.begin(),keepList.end(),g) == keepList.end()) continue;
      if(!dropList.empty() &&
	 find(dropList.begin(),dropList.end(),g) != dropList.end()) continue;
      acc.group(g);
    }

  if(acc.groupName.empty())
    {
      badData("no grouping in " + p.samples.val + " survived --pops/--exclude-pops.");
    }

  string line;
  vector<Field> fields;
  bool haveHeader = false;

  vector<int> sampleGroup;    //per VCF sample column: grouping index, or -1
  vector<int> samplePloidy;   //per VCF sample column: gene copies, 0 until seen
  vector<int> sampleRank;     //per VCF sample column: position within grouping
  int usedSamples = 0;

  long long records = 0;
  string token;

  while(in.next(line))
    {
      if(line.empty()) continue;

      if(line.compare(0,2,"##") == 0) continue;

      if(!haveHeader)
	{
	  if(line[0] != '#')
	    {
	      badData("no #CHROM header line before the records in " + p.dfile.val + ".");
	    }

	  tokenize(line,fields);
	  if(fields.size() < 10)
	    {
	      badData("the #CHROM line of " + p.dfile.val +
		      " names no samples; there is nothing to count.");
	    }

	  sampleGroup.assign(fields.size()-9,-1);
	  samplePloidy.assign(fields.size()-9,0);
	  sampleRank.assign(fields.size()-9,0);
	  vector<int> seenInGroup(acc.groupName.size(),0);

	  int unmapped = 0;
	  for(size_t c = 9; c < fields.size(); c++)
	    {
	      token.assign(fields[c].first,fields[c].second);
	      unordered_map<string,string>::const_iterator it = groupOfSample.find(token);
	      if(it == groupOfSample.end()) { unmapped++; continue; }

	      unordered_map<string,int>::const_iterator gi = acc.groupOf.find(it->second);
	      if(gi == acc.groupOf.end()) continue;   //grouping filtered out

	      sampleGroup[c-9] = gi->second;
	      sampleRank[c-9] = seenInGroup[gi->second]++;
	      usedSamples++;
	    }

	  if(usedSamples == 0)
	    {
	      badData("none of the samples in " + p.dfile.val + " appears in " +
		      p.samples.val + " under a grouping being analysed.");
	    }

	  if(unmapped > 0)
	    {
	      adzelog() << "WARNING: " << unmapped << " of " << (fields.size()-9)
			<< " samples in " << p.dfile.val << " are absent from "
			<< p.samples.val << " and were skipped.\n";
	    }

	  const int mapOnly = int(groupOfSample.size()) - (int(fields.size()-9) - unmapped);
	  if(mapOnly > 0)
	    {
	      adzelog() << "WARNING: " << mapOnly << " samples named in "
			<< p.samples.val << " are not in " << p.dfile.val << ".\n";
	    }

	  haveHeader = true;
	  continue;
	}

      tokenize(line,fields);
      if(fields.empty()) continue;

      if(fields.size() != sampleGroup.size() + 9)
	{
	  ostringstream m;
	  m << "record " << (records+1) << " of " << p.dfile.val << " has "
	    << fields.size() << " columns but the header declares "
	    << (sampleGroup.size() + 9) << ".";
	  badData(m.str());
	}

      //Locus name: the ID column when it carries one, else CHROM:POS.
      string locusName(fields[2].first,fields[2].second);
      if(locusName == ".")
	{
	  locusName.assign(fields[0].first,fields[0].second);
	  locusName += ":";
	  locusName.append(fields[1].first,fields[1].second);
	}

      const int l = int(acc.locusName.size());
      acc.locusName.push_back(locusName);

      //CHROM and POS give this format its coordinates for free.
      {
	const string chrom(fields[0].first,fields[0].second);
	const string posText(fields[1].first,fields[1].second);
	lmap.chrom.push_back(lmap.chromIndex(chrom));
	lmap.pos.push_back(atoll(posText.c_str()));
      }
      acc.locus.resize(acc.locusName.size());
      LocusTally& t = acc.locus[l];
      t.missing.assign(acc.groupName.size(),0);   //observed calls for now

      //Where GT sits in the colon-separated FORMAT field.
      int gtField = -1;
      {
	const char* f = fields[8].first;
	const size_t n = fields[8].second;
	int index = 0;
	size_t i = 0;
	while(i <= n)
	  {
	    size_t j = i;
	    while(j < n && f[j] != ':') j++;
	    if(j - i == 2 && f[i] == 'G' && f[i+1] == 'T') { gtField = index; break; }
	    if(j >= n) break;
	    i = j + 1;
	    index++;
	  }
      }
      if(gtField < 0)
	{
	  ostringstream m;
	  m << "record " << (records+1) << " of " << p.dfile.val
	    << " has no GT in its FORMAT column.";
	  badData(m.str());
	}

      records++;

      for(size_t c = 0; c < sampleGroup.size(); c++)
	{
	  const int g = sampleGroup[c];
	  if(g < 0) continue;

	  const char* f = fields[c+9].first;
	  const size_t n = fields[c+9].second;

	  //Walk to the GT subfield.
	  size_t i = 0;
	  for(int skip = 0; skip < gtField && i < n; skip++)
	    {
	      while(i < n && f[i] != ':') i++;
	      if(i < n) i++;
	    }
	  size_t stop = i;
	  while(stop < n && f[stop] != ':') stop++;

	  if(i >= stop) continue;             //empty GT: no gene copies here
	  if(stop - i == 1 && f[i] == '.') { if(samplePloidy[c] == 0) samplePloidy[c] = 1; continue; }

	  int copies = 0;
	  size_t a = i;
	  const long long sampleKey =
	    (long long)(g) * 4294967296LL + (long long)(sampleRank[c]) * 64LL;
	  while(a < stop)
	    {
	      size_t b = a;
	      while(b < stop && f[b] != '/' && f[b] != '|') b++;

	      copies++;
	      if(!(b - a == 1 && f[a] == '.'))
		{
		  token.assign(f + a, b - a);

		  pair<unordered_map<string,int>::iterator,bool> found =
		    t.slotOf.insert(make_pair(token,int(t.firstSeen.size())));

		  if(found.second)
		    {
		      t.firstSeen.push_back(sampleKey + (copies - 1));
		      t.count.resize(t.count.size() + acc.groupCap, 0);
		    }

		  t.count[size_t(found.first->second)*acc.groupCap + g]++;
		  t.missing[g]++;               //observed for now
		}

	      a = (b < stop) ? b + 1 : stop;
	    }

	  if(samplePloidy[c] < copies) samplePloidy[c] = copies;
	}
    }

  if(!haveHeader) badData("no #CHROM header line in " + p.dfile.val + ".");
  if(records == 0) badData("no variant records in " + p.dfile.val + ".");

  //Gene copies per grouping, and hence the missing-data denominator.
  geneCopies = 0;
  for(size_t c = 0; c < sampleGroup.size(); c++)
    {
      const int g = sampleGroup[c];
      if(g < 0) continue;
      const int ploidy = (samplePloidy[c] > 0) ? samplePloidy[c] : 2;
      acc.groupRows[g] += ploidy;
      geneCopies += ploidy;
    }

  //Convert observed calls into missing counts now that the totals are known.
  for(size_t l = 0; l < acc.locus.size(); l++)
    {
      LocusTally& t = acc.locus[l];
      for(size_t g = 0; g < acc.groupName.size(); g++)
	{
	  const int seen = t.missing[g];
	  t.missing[g] = (acc.groupRows[g] > seen) ? acc.groupRows[g] - seen : 0;
	}
    }

  if(p.loci.set && p.loci.val != int(records))
    {
      adzelog() << "WARNING: LOCI says " << p.loci.val << " but " << p.dfile.val
		<< " has " << records << " records; using " << records << ".\n";
    }
  p.loci.val = int(records);

  adzelog() << "Read " << records << (records == 1 ? " record" : " records")
	    << " for " << usedSamples
	    << (usedSamples == 1 ? " sample" : " samples") << ".\n";

  return;
}

/*
 * Read the data file and return a freshly allocated array of numDivs
 * Population objects, one per grouping, with allele counts, sample sizes and
 * missing-data tallies already filled in.
 *
 * Throws BAD_FILE if the file cannot be opened and BAD_PARAM if its shape is
 * internally inconsistent.
 */
Population* readDataset(ParamSet& p, vector<string>& groupNames, int& numDivs,
			LocusMap& lmap)
{
  const bool vcf = wantsVCF(p.format.val,p.dfile.val);

#ifndef ADZE_HAVE_ZLIB
  if(looksCompressed(p.dfile.val))
    {
      cerr << "ERROR: " << p.dfile.val << " looks compressed, but this build "
	   << "has no zlib support.\n       Rebuild with zlib, or decompress "
	   << "the file first.\n";
      BAD_FILE x;
      throw x;
    }
#endif

  LineSource in;
  if(!in.open(p.dfile.val))
    {
      cerr << "ERROR: could not open " << p.dfile.val << "\n";
      BAD_FILE x;
      throw x;
    }

  vector<string> keepList, dropList;
  splitList(p.pops.val,keepList);
  splitList(p.expops.val,dropList);

  Accumulator acc;
  long long dataRows = 0, keptRows = 0;

  if(vcf) readVCFInto(p,acc,in,keepList,dropList,dataRows,lmap);
  else readStructureInto(p,acc,in,keepList,dropList,dataRows,keptRows);

  in.close();

  if(!vcf && p.loci_map.set)
    {
      readLocusMap(p.loci_map.val,acc.locusName,lmap);
    }
  else if(vcf && p.loci_map.set)
    {
      adzelog() << "WARNING: --loci-map applies to STRUCTURE input only; "
		<< "ignoring it. VCF coordinates come from CHROM and POS.\n";
    }

  const int declaredLoci = int(acc.locusName.size());

  if(!vcf)
    {
      if(p.dlines.set && p.dlines.val != int(dataRows))
	{
	  adzelog() << "WARNING: DATA_LINES says " << p.dlines.val << " but "
		    << p.dfile.val << " has " << dataRows << " data rows; using "
		    << dataRows << ".\n";
	}
      if(dataRows == 0) badData("no data rows in " + p.dfile.val + ".");
    }
  p.dlines.val = int(dataRows);

  numDivs = int(acc.groupName.size());
  if(numDivs == 0)
    {
      badData("no grouping in " + p.dfile.val + " survived --pops/--exclude-pops.");
    }

  if(!vcf && keptRows != dataRows)
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

  /*
   * Publish the counts, ordering each locus's allele slots by where the allele
   * was first seen -- grouping, then row (or sample) within it, then gene copy
   * -- which reproduces 1.0's Nji column order for STRUCTURE input and gives
   * VCF input the same layout for the same genotypes.
   */
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
		bool pp, LocusMap& lmap)
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

  //Coordinates are indexed by locus, so they follow the same compaction.
  lmap.compact(toDelete);
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
		  bool full_comb, string comb_out, bool namedTuples,
		  const vector<Window>& windows, const LocusMap& lmap)
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
  if(gLast < 1) gLast = 1;
  /*
   * --at-g: the recurrence in g is sequential, so it still has to climb to the
   * g being reported -- but no higher, and the buffer need hold no more than
   * that. Asking for g = 2 where MAX_G is 33 therefore builds a table of 3
   * rows instead of 34 and stops the climb 31 steps early.
   */
  const int gCeil = param.at_g_val ? param.at_g_val : gLast;
  const int gStride = gCeil + 1;

  ProgressBar bar(&adzelog(),double(tot_m)*(param.at_g_val ? 1 : gLast)*numLoci,
		  BARLEN[min(int(widest)-1,3)]);
  if(param.pp.val)
    {
      bar.init();
    }

  vector<string> names(widest);
  vector<char> inTuple(numDivs,0);
  ofstream comb_win_out;
  if(!windows.empty())
    {
      string name = nameCreate(comb_out,"_windows");
      comb_win_out.open(name.c_str());
      writeWindowHeader(comb_win_out,"TUPLE");
    }

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
	  buildQTable(pop,numDivs,locus,numAlleles,gCeil,gStride,q);

	  for(int g = 1; g <= gCeil; g++)
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
	      bar.adv(param.at_g_val ? 1 : gLast);
	    }
	}
      } //end parallel region

      for(int g = 1; g <= gCeil; g++)
	{
	  //--at-g: one row, not the ladder.  The sweep above still had to
	  //climb to it, the recurrence in g being sequential.
	  if(param.at_g_val && g != param.at_g_val) continue;

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

      if(!windows.empty())
	{
	  /*
	   * The tuple label is comma-joined here whatever the genome-wide
	   * format, since window output is always tab-separated and the label
	   * has to stay one field.
	   */
	  const string win_names = combineNames(&names[0],k,',');
	  writeWindowStats(comb_win_out,pgcomb,numLoci,windows,lmap,
			   win_names,
			   param.at_g_val ? param.at_g_val : 1,
			   param.at_g_val ? param.at_g_val : gCeil);
	}

      if(!param.tsv.val) reg_out << endl;
      full_out << endl;
    }

  if(param.pp.val) bar.done();

  if(comb_win_out.is_open()) comb_win_out.close();

  return;
}

/*
 * ------------------------------------------------------------------------
 * Windowed output
 * ------------------------------------------------------------------------
 *
 * A window is a contiguous run of surviving loci, and each statistics pass
 * already holds its per-locus values as a [g][locus] buffer -- the same array
 * the _fulldata files print.  A windowed row is therefore the existing
 * summary over a subrange of that array: same mean, variance and standard
 * error code, applied to the loci inside the window instead of all of them.
 * Nothing about the estimators changes.
 *
 * Window files always carry a header and are always tab-separated. They are
 * new output with no 1.0 consumer to keep compatible, and a genome scan is
 * read by machine.
 */
void writeWindowHeader(ostream& out, const char* groupColumn)
{
  out << "CHROM\tSTART\tEND\t" << groupColumn
      << "\tG\tNUM_LOCI\tMEAN\tVAR\tSTD_ERR\n";
  return;
}

void writeWindowStats(ostream& out, vector<double>& perLocus, int numLoci,
		      const vector<Window>& windows, const LocusMap& lmap,
		      const string& label, int gFirst, int gLast)
{
  for(int g = gFirst; g <= gLast; g++)
    {
      for(size_t w = 0; w < windows.size(); w++)
	{
	  const Window& win = windows[w];

	  Stats st;
	  st.putData(&perLocus[size_t(g)*numLoci + win.first],win.numLoci());
	  st.calcAvg();
	  st.calcVar();
	  st.calcStdErr();

	  ostringstream head;
	  head << lmap.chromName[win.chrom] << '\t' << win.start << '\t'
	       << win.end << '\t' << label;

	  //tsv = true: a window undefined at this g is reported as NA rather
	  //than dropped, so a scan keeps one row per window per g.
	  st.printStats(out,head.str(),g,1);
	}
    }

  return;
}

void calcAllPgs(Population pop[],int numDivs,const ParamSet &param,
		bool full_priv,string private_out,
		const vector<Window>& windows,const LocusMap& lmap)
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

  ofstream pg_win_out;
  if(!windows.empty())
    {
      string name = nameCreate(private_out,"_windows");
      pg_win_out.open(name.c_str());
      writeWindowHeader(pg_win_out,"POP_GROUPING");
    }

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
  if(gLast < 1) gLast = 1;
  /*
   * --at-g: the recurrence in g is sequential, so it still has to climb to the
   * g being reported -- but no higher, and the buffer need hold no more than
   * that. Asking for g = 2 where MAX_G is 33 therefore builds a table of 3
   * rows instead of 34 and stops the climb 31 steps early.
   */
  const int gCeil = param.at_g_val ? param.at_g_val : gLast;
  const int gStride = gCeil + 1;

  vector<double> pg(size_t(gStride) * numLoci, 0.0); //[g][locus]

  ProgressBar bar(&adzelog(),double(numDivs)*(param.at_g_val ? 1 : gLast)*numLoci,BARLEN[0]);
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
	  buildQTable(pop,numDivs,locus,numAlleles,gCeil,gStride,q);

	  for(int g = 1; g <= gCeil; g++)
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
	      bar.adv(param.at_g_val ? 1 : gLast);
	    }
	}
      } //end parallel region

      for(int g = 1; g <= gCeil; g++)
	{
	  //--at-g: one row, not the ladder.  The sweep above still had to
	  //climb to it, the recurrence in g being sequential.
	  if(param.at_g_val && g != param.at_g_val) continue;

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

      if(!windows.empty())
	{
	  writeWindowStats(pg_win_out,pg,numLoci,windows,lmap,
			   pop[j].getName(),
			   param.at_g_val ? param.at_g_val : 1,
			   param.at_g_val ? param.at_g_val : gCeil);
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
  if(pg_win_out.is_open()) pg_win_out.close();

  return;
}


void calcAllAgs(Population pop[],int numDivs,const ParamSet &param,
		bool full_rich,string richness_out,
		const vector<Window>& windows,const LocusMap& lmap)
{
  const int numLoci = param.loci.val;
  const int gTop = param.g.val;
  /*
   * --at-g: the recurrence in g is sequential, so it still has to climb to the
   * g being reported -- but no higher, and the buffer need hold no more than
   * that. Asking for g = 2 where MAX_G is 33 therefore builds a table of 3
   * rows instead of 34 and stops the climb 31 steps early.
   */
  const int gCeil = param.at_g_val ? param.at_g_val : gTop;
  const int gStride = gCeil + 1;
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

  ofstream ag_win_out;
  if(!windows.empty())
    {
      string name = nameCreate(richness_out,"_windows");
      ag_win_out.open(name.c_str());
      writeWindowHeader(ag_win_out,"POP_GROUPING");
    }

  ProgressBar bar(&adzelog(),double(numDivs)*(param.at_g_val ? 1 : gTop)*numLoci,BARLEN[0]);
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

	  buildQTable(&pop[j],1,locus,numAlleles,gCeil,gStride,q);

	  for(int g = 1; g <= gCeil; g++)
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
	      bar.adv(param.at_g_val ? 1 : gTop);
	    }
	}
      } //end parallel region

      for(int g = 1; g <= gCeil; g++)
	{
	  //--at-g: one row, not the ladder.  The sweep above still had to
	  //climb to it, the recurrence in g being sequential.
	  if(param.at_g_val && g != param.at_g_val) continue;

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

      if(!windows.empty())
	{
	  writeWindowStats(ag_win_out,ag,numLoci,windows,lmap,
			   pop[j].getName(),
			   param.at_g_val ? param.at_g_val : 1,
			   param.at_g_val ? param.at_g_val : gCeil);
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
  if(ag_win_out.is_open()) ag_win_out.close();

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
