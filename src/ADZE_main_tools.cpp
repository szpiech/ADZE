#include "ADZE_main_tools.h"
#include <limits>
#include <unistd.h>
#include <sys/stat.h>

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
    vector<string> label;              //allele label per slot, first-seen order
    vector<long long> firstSeen;       //per slot: (grouping, row) of first sighting
    vector<int> count;                 //slot-major, stride groupCap
    vector<int> missing;               //per grouping

    /*
     * The slot for an allele label, appended if this is the first sighting.
     *
     * A scan, not a hash. A locus carries two alleles for a SNP and tens for
     * a microsatellite, and an unordered_map per locus cost far more in
     * bucket array and node headers than the counts it indexed -- about 200
     * bytes per locus, which at genome scale was the largest single thing
     * the reader held. The scan also orders itself: slots are created in
     * first-seen order, so a locus's common alleles sit at the front, where
     * the scan reaches them first.
     */
    int slot(const string& allele, long long firstKey, int groupCap)
    {
      for(size_t i = 0; i < label.size(); i++)
	{
	  if(label[i] == allele) return int(i);
	}

      /*
       * Exact growth, not geometric: a locus holds two alleles for a SNP and
       * tens for a microsatellite, so the copying is trivial, while the
       * doubling these vectors do by default left each locus holding up to
       * twice the counts it needed -- across a genome, the largest single
       * thing the reader carried.
       */
      label.reserve(label.size() + 1);
      label.push_back(allele);

      firstSeen.reserve(firstSeen.size() + 1);
      firstSeen.push_back(firstKey);

      count.reserve(count.size() + groupCap);
      count.resize(count.size() + groupCap, 0);

      return int(label.size()) - 1;
    }
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

    /*
     * Fix the count stride before reading begins, when the number of
     * groupings is already known -- from the sample map for VCF input, or
     * from --pops. Nothing has been allocated yet at that point, so this
     * costs nothing and leaves no padding at all. Grouping indices are still
     * assigned in first-seen order; only the layout is affected.
     */
    void setGroupCap(int n)
    {
      if(locus.empty() && n > 0) groupCap = n;
      return;
    }

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
	  /*
	   * Widen every locus's count block, doubling rather than growing to
	   * exactly the groupings seen.
	   *
	   * Doubling leaves the stride at the next power of two -- 16 columns
	   * for 10 groupings -- which is waste. Growing by one removes it but
	   * re-strides every locus J-1 times instead of log2(J), and freeing
	   * 80 000 blocks of one size while allocating 80 000 of the next leaves
	   * the pages mapped: measured on 80 000 biallelic loci it cost 10 MB
	   * more than the padding it saved, while saving 6 MB on 20 alleles per
	   * locus. Neither is needed when the groupings are known before the
	   * first locus exists, which is the usual case -- see setGroupCap.
	   */
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

	  const int sl = t.slot(token,firstKey,acc.groupCap);
	  t.count[size_t(sl)*acc.groupCap + g]++;
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

  //Every grouping is named in the sample map, and no locus exists yet.
  acc.setGroupCap(int(acc.groupName.size()));

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

		  const int sl = t.slot(token,sampleKey + (copies - 1),acc.groupCap);
		  t.count[size_t(sl)*acc.groupCap + g]++;
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
 * A VCF read one record at a time.
 *
 * The parsing is readVCFInto's, with one difference that is the point of the
 * rewrite: nothing accumulates. A record's alleles are interned into a tally
 * that is reused for the next record, the slots are put into the order the
 * output has always used, and the locus is handed to the caller. The sample
 * map, the grouping order, the ploidies and the per-grouping totals all come
 * from the scan, which has already read the file -- so this pass needs no
 * end-of-file fixups and can hand out the first locus immediately.
 */
namespace {

  class VCFSource : public LocusSource
  {
  public:
    VCFSource(ParamSet& p,const ScanResult& scan)
      : param(p), rows(scan.groupRows), name(scan.groupName), records(0)
    {
      if(!in.open(p.dfile.val))
	{
	  cerr << "ERROR: could not open " << p.dfile.val << "\n";
	  throw BAD_FILE();
	}
      readHeader();
    }

    const vector<string>& groupNames() const { return name; }
    long long geneCopies(int grouping) const { return rows[grouping]; }

    bool next(LocusCounts& out)
    {
      const int J = int(name.size());
      string line;

      while(in.next(line))
	{
	  if(line.empty()) continue;
	  if(line[0] == '#') continue;

	  tokenize(line,fields);
	  if(fields.empty()) continue;

	  out.name.assign(fields[2].first,fields[2].second);
	  if(out.name == ".")
	    {
	      out.name.assign(fields[0].first,fields[0].second);
	      out.name += ":";
	      out.name.append(fields[1].first,fields[1].second);
	    }

	  {
	    const string chrom(fields[0].first,fields[0].second);
	    const string posText(fields[1].first,fields[1].second);
	    out.chrom = chromIndex(chrom);
	    out.pos = atoll(posText.c_str());
	  }

	  //One tally, reused: a locus's alleles are interned and forgotten.
	  tally.label.clear();
	  tally.firstSeen.clear();
	  tally.count.clear();

	  const int gtField = gtPosition(fields[8]);

	  for(size_t c = 0; c < sampleGroup.size(); c++)
	    {
	      const int g = sampleGroup[c];
	      if(g < 0) continue;

	      const char* f = fields[c+9].first;
	      const size_t n = fields[c+9].second;

	      size_t i = 0;
	      for(int skip = 0; skip < gtField && i < n; skip++)
		{
		  while(i < n && f[i] != ':') i++;
		  if(i < n) i++;
		}
	      size_t stop = i;
	      while(stop < n && f[stop] != ':') stop++;

	      if(i >= stop) continue;
	      if(stop - i == 1 && f[i] == '.') continue;

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
		      const int sl = tally.slot(token,sampleKey + (copies - 1),J);
		      tally.count[size_t(sl)*J + g]++;
		    }

		  a = (b < stop) ? b + 1 : stop;
		}
	    }

	  records++;
	  finish(out,J);
	  return true;
	}

      return false;
    }

  private:
    ParamSet& param;
    const vector<long long>& rows;
    const vector<string>& name;
    LineSource in;
    LocusTally tally;
    vector<Field> fields;
    vector<int> sampleGroup, sampleRank;
    vector<string> chromName;
    unordered_map<string,int> chromOf;
    string token;
    long long records;

    int chromIndex(const string& c)
    {
      unordered_map<string,int>::iterator it = chromOf.find(c);
      if(it != chromOf.end()) return it->second;
      const int i = int(chromName.size());
      chromName.push_back(c);
      chromOf.insert(make_pair(c,i));
      return i;
    }

    int gtPosition(const Field& format)
    {
      const char* f = format.first;
      const size_t n = format.second;
      int index = 0;
      size_t i = 0;
      while(i <= n)
	{
	  size_t j = i;
	  while(j < n && f[j] != ':') j++;
	  if(j - i == 2 && f[i] == 'G' && f[i+1] == 'T') return index;
	  if(j >= n) break;
	  i = j + 1;
	  index++;
	}

      ostringstream m;
      m << "record " << (records+1) << " of " << param.dfile.val
	<< " has no GT in its FORMAT column.";
      badData(m.str());
      return -1;
    }

    //Slots into output order, and Nj from the counts.
    void finish(LocusCounts& out,int J)
    {
      const int slots = int(tally.label.size());
      out.slots = slots;
      out.count.assign(size_t(slots)*J,0);
      out.nj.assign(J,0);

      vector<int> order(slots);
      for(int sl = 0; sl < slots; sl++) order[sl] = sl;
      sort(order.begin(),order.end(),FirstSeenLess(tally.firstSeen));

      for(int i = 0; i < slots; i++)
	{
	  for(int g = 0; g < J; g++)
	    {
	      const int c = tally.count[size_t(order[i])*J + g];
	      out.count[size_t(i)*J + g] = c;
	      out.nj[g] += c;
	    }
	}
      return;
    }

    void readHeader()
    {
      unordered_map<string,string> groupOfSample;
      vector<string> groupOrder;
      readSampleMap(param.samples.val,groupOfSample,groupOrder);

      unordered_map<string,int> indexOf;
      for(size_t g = 0; g < name.size(); g++) indexOf.insert(make_pair(name[g],int(g)));

      string line;
      while(in.next(line))
	{
	  if(line.empty()) continue;
	  if(line.compare(0,2,"##") == 0) continue;
	  if(line[0] != '#')
	    {
	      badData("no #CHROM header line before the records in " + param.dfile.val + ".");
	    }

	  tokenize(line,fields);
	  sampleGroup.assign(fields.size()-9,-1);
	  sampleRank.assign(fields.size()-9,0);
	  vector<int> seenInGroup(name.size(),0);

	  for(size_t c = 9; c < fields.size(); c++)
	    {
	      token.assign(fields[c].first,fields[c].second);
	      unordered_map<string,string>::const_iterator it = groupOfSample.find(token);
	      if(it == groupOfSample.end()) continue;

	      unordered_map<string,int>::const_iterator gi = indexOf.find(it->second);
	      if(gi == indexOf.end()) continue;

	      sampleGroup[c-9] = gi->second;
	      sampleRank[c-9] = seenInGroup[gi->second]++;
	    }
	  return;
	}

      badData("no #CHROM header line in " + param.dfile.val + ".");
      return;
    }
  };

} //anonymous namespace

LocusSource* openVCFSource(ParamSet& p,const ScanResult& scan)
{
  return new VCFSource(p,scan);
}

/*
 * STRUCTURE input, converted into the counts the engine reads.
 *
 * A STRUCTURE file is individual-major: one row per gene copy, one column
 * per locus. There is no way to finish a locus without visiting every row,
 * so a locus cannot simply be streamed off it. Rather than teach the engine
 * a second input shape, the file is converted: a chunk of locus columns is
 * counted, written out, and the next chunk begun. The chunk is sized from a
 * memory budget, so a microsatellite dataset of a few thousand loci is one
 * pass and only a very large file needs more.
 *
 * (A strictly single-pass conversion is possible -- write each row's tokens
 * to per-chunk spill files, then count each spill in turn -- but it trades
 * the re-reads for writing and re-reading the genotypes themselves, which is
 * more I/O than the re-reads unless the chunk count is large. Worth doing
 * only if the pass count ever becomes the measured problem.)
 *
 * The file format is transient and host-native: a header naming the
 * groupings and their gene copies, then one record per locus -- name, allele
 * count, then the counts, slot-major with stride J, in the same order
 * VCFSource yields and the reader published.
 */
/*
 * Size and last-modified time of a file, or -1 when it cannot be read. Both
 * are recorded in a converted file so it can be refused when the data has
 * moved on. Neither is a hash: a same-size edit made within the filesystem's
 * timestamp resolution would go unnoticed, which is why doc/streaming.md
 * says to name a converted file only when the input is settled.
 */
static long long fileSize(const string& path)
{
  ifstream f(path.c_str(),ios::binary|ios::ate);
  if(!f.is_open()) return -1;
  return (long long)(f.tellg());
}

static long long fileStamp(const string& path)
{
  struct stat st;
  if(stat(path.c_str(),&st) != 0) return -1;
  return (long long)(st.st_mtime);
}

namespace {

  const char COUNT_MAGIC[8] = {'A','D','Z','E','C','N','T','1'};

  void writeInt(ostream& out,int v) { out.write((const char*)&v,sizeof(int)); }
  void writeLL(ostream& out,long long v) { out.write((const char*)&v,sizeof(long long)); }
  void writeStr(ostream& out,const string& s)
  {
    writeInt(out,int(s.size()));
    out.write(s.data(),s.size());
  }

  bool readInt(istream& in,int& v)
  {
    in.read((char*)&v,sizeof(int));
    return in.good();
  }
  bool readLL(istream& in,long long& v)
  {
    in.read((char*)&v,sizeof(long long));
    return in.good();
  }
  bool readStr(istream& in,string& s)
  {
    int n = 0;
    if(!readInt(in,n) || n < 0) return false;
    s.resize(size_t(n));
    if(n) in.read(&s[0],n);
    return in.good();
  }

  class CountFileSource : public LocusSource
  {
  public:
    CountFileSource(const string& path,const ScanResult& scan)
      : name(scan.groupName), rows(scan.groupRows), in(path.c_str(),ios::binary),
	lmap(&scan.lmap), index(0)
    {
      if(!in.is_open())
	{
	  cerr << "ERROR: could not read the converted counts at " << path << "\n";
	  throw BAD_FILE();
	}

      char magic[8];
      in.read(magic,8);
      if(!in.good() || memcmp(magic,COUNT_MAGIC,8) != 0)
	{
	  cerr << "ERROR: " << path << " is not an adze count file\n";
	  throw BAD_FILE();
	}

      int J = 0;
      long long loci = 0;
      readInt(in,J);
      readLL(in,loci);

      string s1;
      long long l1 = 0;
      int i1 = 0;
      readStr(in,s1); readLL(in,l1); readLL(in,l1);  //source, size, timestamp
      readStr(in,s1); readStr(in,s1);          //--pops, --exclude-pops
      readStr(in,s1);                          //missing code
      readInt(in,i1); readInt(in,i1);          //group column, header rows

      for(int g = 0; g < J; g++)
	{
	  string nm;
	  long long copies = 0;
	  readStr(in,nm);
	  readLL(in,copies);
	}
    }

    const vector<string>& groupNames() const { return name; }
    long long geneCopies(int grouping) const { return rows[grouping]; }

    bool next(LocusCounts& out)
    {
      const int J = int(name.size());
      int slots = 0;

      if(!readStr(in,out.name)) return false;
      if(!readInt(in,slots)) return false;

      out.slots = slots;
      out.count.assign(size_t(slots)*J,0);
      out.nj.assign(J,0);

      if(slots)
	{
	  in.read((char*)&out.count[0],streamsize(sizeof(int))*slots*J);
	  if(!in.good()) return false;
	}

      for(int i = 0; i < slots; i++)
	{
	  for(int g = 0; g < J; g++) out.nj[g] += out.count[size_t(i)*J + g];
	}

      //Coordinates stay with the scan: the conversion need not repeat them.
      if(lmap && lmap->size() > size_t(index))
	{
	  out.chrom = lmap->chrom[size_t(index)];
	  out.pos = lmap->pos[size_t(index)];
	}
      else
	{
	  out.chrom = -1;
	  out.pos = -1;
	}

      index++;
      return true;
    }

  private:
    const vector<string>& name;
    const vector<long long>& rows;
    ifstream in;
    const LocusMap* lmap;
    long long index;
  };

} //anonymous namespace

long long transposeStructure(ParamSet& p,const ScanResult& scan,
			     const string& path,long long lociPerPass)
{
  const int J = int(scan.groupName.size());
  const long long numLoci = scan.numLoci;
  if(lociPerPass < 1) lociPerPass = 1;

  unordered_map<string,int> indexOf;
  for(int g = 0; g < J; g++) indexOf.insert(make_pair(scan.groupName[g],g));

  vector<string> keepList, dropList;
  splitList(p.pops.val,keepList);
  splitList(p.expops.val,dropList);

  ofstream out(path.c_str(),ios::binary);
  if(!out.is_open())
    {
      cerr << "ERROR: could not write the converted counts to " << path << "\n";
      throw BAD_FILE();
    }

  out.write(COUNT_MAGIC,8);
  writeInt(out,J);
  writeLL(out,numLoci);

  /*
   * Provenance, so a converted file cannot be used against data it did not
   * come from. Everything here changes what the counts would be: the source
   * file's identity and size, the grouping filters, the missing code, and
   * which column names the grouping. A mismatch is reported and the file
   * converted again -- never used with a warning.
   */
  writeStr(out,p.dfile.val);
  writeLL(out,fileSize(p.dfile.val));
  writeLL(out,fileStamp(p.dfile.val));
  writeStr(out,p.pops.val);
  writeStr(out,p.expops.val);
  writeStr(out,p.miss.val);
  writeInt(out,p.sort_by.val);
  writeInt(out,p.nd_rows.val);

  for(int g = 0; g < J; g++)
    {
      writeStr(out,scan.groupName[g]);
      writeLL(out,scan.groupRows[g]);
    }

  vector<LocusTally> chunk;
  vector<Field> fields;
  vector<int> order;
  string line, token;
  long long passes = 0;

  for(long long base = 0; base < numLoci; base += lociPerPass)
    {
      const long long stop = (base + lociPerPass < numLoci) ? base + lociPerPass : numLoci;
      const int width = int(stop - base);

      LineSource in;
      if(!in.open(p.dfile.val))
	{
	  cerr << "ERROR: could not open " << p.dfile.val << "\n";
	  throw BAD_FILE();
	}
      passes++;

      //The locus-name row, then any further non-data rows.
      if(!in.next(line)) badData("no locus-name row in " + p.dfile.val + ".");
      tokenize(line,fields);
      vector<string> names;
      names.reserve(size_t(width));
      for(int i = 0; i < width; i++)
	{
	  const Field& f = fields[size_t(base) + i];
	  names.push_back(string(f.first,f.second));
	}
      for(int skip = 1; skip < p.nd_rows.val; skip++) in.next(line);

      chunk.clear();
      chunk.resize(size_t(width));

      const int ndCols = p.nd_cols.val;
      const int groupCol = p.sort_by.val - 1;
      vector<int> rowInGroup(J,0);

      while(in.next(line))
	{
	  tokenize(line,fields);
	  if(fields.empty()) continue;

	  token.assign(fields[groupCol].first,fields[groupCol].second);
	  if(!keepList.empty() &&
	     find(keepList.begin(),keepList.end(),token) == keepList.end()) continue;
	  if(!dropList.empty() &&
	     find(dropList.begin(),dropList.end(),token) != dropList.end()) continue;

	  unordered_map<string,int>::const_iterator gi = indexOf.find(token);
	  if(gi == indexOf.end()) continue;

	  const int g = gi->second;
	  const long long firstKey = (long long)(g) * 4294967296LL + rowInGroup[g];
	  rowInGroup[g]++;

	  for(int i = 0; i < width; i++)
	    {
	      const Field& f = fields[ndCols + size_t(base) + i];
	      token.assign(f.first,f.second);
	      if(token.compare(p.miss.val) == 0) continue;

	      LocusTally& t = chunk[size_t(i)];
	      const int sl = t.slot(token,firstKey,J);
	      t.count[size_t(sl)*J + g]++;
	    }
	}

      for(int i = 0; i < width; i++)
	{
	  LocusTally& t = chunk[size_t(i)];
	  const int slots = int(t.label.size());

	  order.resize(size_t(slots));
	  for(int sl = 0; sl < slots; sl++) order[size_t(sl)] = sl;
	  sort(order.begin(),order.end(),FirstSeenLess(t.firstSeen));

	  writeStr(out,names[size_t(i)]);
	  writeInt(out,slots);
	  for(int sl = 0; sl < slots; sl++)
	    {
	      for(int g = 0; g < J; g++)
		{
		  writeInt(out,t.count[size_t(order[size_t(sl)])*J + g]);
		}
	    }
	}
    }

  out.close();
  return passes;
}


/*
 * Where a conversion goes when the user has not named one: beside the
 * output, not in a system temporary directory, because it can be the size of
 * the counts and the output directory is the one the user chose for large
 * files. It is removed when the run ends, successfully or not; only a file
 * named with --counts is kept.
 */
string countFilePath(const ParamSet& p)
{
  return p.out_prefix.val + ".counts.tmp";
}

//How much disk a conversion of this dataset needs, near enough to warn with.
long long countFileBytes(const ScanResult& scan)
{
  const long long J = (long long)(scan.groupName.size());
  long long slots = 0;
  for(long long l = 0; l < scan.numLoci; l++) slots += 1;   //at least one each
  return scan.numLoci*(64 + 8) + slots*J*(long long)(sizeof(int));
}

/*
 * Loci per conversion pass. A pass holds one tally per locus in the chunk --
 * labels, first-sighting keys and a count per grouping per allele -- so the
 * budget buys loci in inverse proportion to the groupings. 64 MB is enough
 * for a microsatellite dataset many times over, and for a million SNPs at
 * ten groupings it means a handful of passes rather than one enormous one.
 */
long long convertChunk(const ParamSet& p,const ScanResult& scan)
{
  const long long budget = 64LL*1024*1024;
  const long long J = (long long)(scan.groupName.size());
  const long long perLocus = 4*(long long)(sizeof(int))*(J > 0 ? J : 1) + 64;
  long long chunk = budget/perLocus;
  if(chunk < 1) chunk = 1;
  if(chunk > scan.numLoci) chunk = scan.numLoci;
  return chunk;
}

/*
 * Whether a converted file on disk describes this run's data. Every field
 * compared changes what the counts would be, so a mismatch means convert
 * again -- the alternative, using it with a warning, is how a user ends up
 * with yesterday's genotypes in today's table.
 */
bool countFileUsable(const string& path,const ParamSet& p,
		     const ScanResult& scan,string& why)
{
  ifstream in(path.c_str(),ios::binary);
  if(!in.is_open()) { why = "it does not exist"; return false; }

  char magic[8];
  in.read(magic,8);
  if(!in.good() || memcmp(magic,COUNT_MAGIC,8) != 0)
    {
      why = "it is not an adze count file";
      return false;
    }

  int J = 0;
  long long loci = 0;
  if(!readInt(in,J) || !readLL(in,loci)) { why = "its header is truncated"; return false; }

  string src, pops, expops, miss;
  long long bytes = 0, stamp = 0;
  int groupCol = 0, ndRows = 0;
  if(!readStr(in,src) || !readLL(in,bytes) || !readLL(in,stamp)
     || !readStr(in,pops) || !readStr(in,expops)
     || !readStr(in,miss) || !readInt(in,groupCol) || !readInt(in,ndRows))
    {
      why = "its header is truncated";
      return false;
    }

  if(src != p.dfile.val) { why = "it was converted from " + src; return false; }
  if(bytes != fileSize(p.dfile.val)) { why = "the data file has changed size since"; return false; }
  if(stamp != fileStamp(p.dfile.val)) { why = "the data file has been modified since"; return false; }
  if(pops != p.pops.val || expops != p.expops.val) { why = "the grouping filters differ"; return false; }
  if(miss != p.miss.val) { why = "the missing-data code differs"; return false; }
  if(groupCol != p.sort_by.val) { why = "a different column names the grouping"; return false; }
  if(ndRows != p.nd_rows.val) { why = "a different number of header rows"; return false; }
  if(loci != scan.numLoci || J != int(scan.groupName.size()))
    {
      why = "it holds a different number of loci or groupings";
      return false;
    }

  for(int g = 0; g < J; g++)
    {
      string nm;
      long long copies = 0;
      if(!readStr(in,nm) || !readLL(in,copies)) { why = "its header is truncated"; return false; }
      if(nm != scan.groupName[g] || copies != scan.groupRows[g])
	{
	  why = "its groupings differ from this run's";
	  return false;
	}
    }

  why = "";
  return true;
}

LocusSource* openCountFileSource(const string& path,const ScanResult& scan)
{
  return new CountFileSource(path,scan);
}

/*
 * The dry run: what the program would do with this dataset, from the counts
 * alone. Every figure here -- the dimensions, each grouping's sample size and
 * the range of Nj it scored, the largest feasible MAX_G, the window layout
 * and the row counts -- follows from the scan, so a dry run never reads the
 * dataset into memory.
 */
void printDryRun(const ParamSet& p,const ScanResult& scan,
		 const vector< vector<int> >& tuples,list<int>& k,
		 bool do_rich,bool do_priv,bool do_tuple)
{
  const int numDivs = int(scan.groupName.size());

  const bool vcf = wantsVCF(p.format.val,p.dfile.val);

  cout << "Data file:        " << p.dfile.val
         << (vcf ? "  (VCF)\n" : "  (STRUCTURE)\n")
         << "Loci:             " << p.loci.val << "\n"
         << "Gene copies:      " << p.dlines.val << "\n";

  if(vcf)
      {
        cout << "Sample map:       " << p.samples.val << "\n";
      }
  else
      {
        cout << "Label columns:    " << p.nd_cols.val
             << " (grouping from column " << p.sort_by.val << ")\n"
             << "Header rows:      " << p.nd_rows.val << "\n"
             << "Missing code:     " << p.miss.val << "\n";
      }

  cout << "Groupings:        " << numDivs << "\n";

  for(int j = 0; j < numDivs; j++)
      {
        int minNj = 0, maxNj = 0;
        for(long long l = 0; l < scan.numLoci; l++)
          {
            const int Nj = scan.nj(j,l);
            if(l == 0 || Nj < minNj) minNj = Nj;
            if(l == 0 || Nj > maxNj) maxNj = Nj;
          }
        cout << "  " << scan.groupName[j] << ": " << scan.groupRows[j]
             << " gene copies, Nj per locus " << minNj << "-" << maxNj << "\n";
      }

  cout << "Largest feasible MAX_G: " << scan.feasibleGAll << "\n";
  if(p.g.set) cout << "MAX_G to be used:       " << p.g.val << "\n";
  else cout << "MAX_G to be used:       " << "the largest feasible\n";
  if(p.tol.val != 1)
      {
        cout << "  (TOLERANCE " << p.tol.val << " will drop loci first, "
             << "which can raise both numbers)\n";
      }

  /*
   * The g that will actually be reported, against the ceilings known
   * before filtering: the same rule the run itself applies, so the dry
   * run says what the run will do rather than what was asked for.
   */
  if(p.at_g.set)
      {
        const int sweepTop = p.g.set ? p.g.val : ((scan.feasibleGAll < 1) ? 1 : scan.feasibleGAll);
        int ceiling = (scan.feasibleGAll < sweepTop) ? scan.feasibleGAll : sweepTop;
        if(ceiling < 1) ceiling = 1;

        const int asked = (p.at_g.val == "max") ? sweepTop : p.at_g_val;
        cout << "Reporting at g:         " << ((asked > ceiling) ? ceiling : asked);
        if(asked > ceiling) cout << "  (--at-g " << p.at_g.val
      			   << " clamped to the ceiling)";
        cout << "\n";
      }

  if(p.windowed())
      {
        /*
         * Windows are laid over the loci that survive filtering, so work
         * out which those are -- without writing the report file, this
         * being a dry run -- and lay the windows over a copy.
         */
        LocusMap trial = scan.lmap;
        trial.compact(scan.dropped);

        vector<Window> windows;
        const long long sparse = buildWindows(trial,p,windows);

        int lo = 0, hi = 0;
        for(size_t w = 0; w < windows.size(); w++)
          {
            const int n = windows[w].numLoci();
            if(w == 0 || n < lo) lo = n;
            if(w == 0 || n > hi) hi = n;
          }

        cout << "Windows:          " << windows.size();
        if(p.win_bp.set)
          {
            cout << " of " << p.win_bp.val << " bp, step " << p.step_bp.val;
          }
        else
          {
            cout << " of " << p.win_loci.val << " loci, step "
      	   << p.step_loci.val;
          }
        cout << "\n";

        if(!windows.empty())
          {
            cout << "  loci per window: " << lo << "-" << hi << "\n";
          }
        if(sparse > 0)
          {
            cout << "  " << sparse << " window" << (sparse == 1 ? "" : "s")
      	   << " below --min-window-loci " << p.min_win_loci.val
      	   << ", not reported\n";
          }

        /*
         * Windowed output is one row per window per g per grouping, and the
         * tuple statistics multiply that by the number of tuples, so say how
         * large the files will be before anyone waits for them.
         */
        //Rows per window: the whole ladder from g = 1, or one row for --at-g.
        const int gRange = p.at_g.set
          ? 1 : (p.g.set ? p.g.val : ((scan.feasibleGAll < 1) ? 1 : scan.feasibleGAll));
        if(gRange > 0)
          {
            long long rows = (long long)(windows.size()) * gRange;
            if(do_rich) cout << "  richness rows:   " << rows*numDivs << "\n";
            if(do_priv) cout << "  private rows:    " << rows*numDivs << "\n";
            if(do_tuple)
      	{
      	  double nTuples = double(tuples.size());
      	  if(!p.tuple_file.set)
      	    {
      	      nTuples = 0;
      	      for(list<int>::iterator i = k.begin(); i != k.end(); i++)
      		{
      		  nTuples += nCk(numDivs,*i);
      		}
      	    }
      	  cout << "  tuple rows:      " << long(double(rows)*nTuples)
      	       << " (over " << long(nTuples) << " tuple"
      	       << (nTuples == 1 ? "" : "s") << ")\n";
      	}
          }
      }

  if(do_tuple)
      {
        if(p.tuple_file.set) cout << "Named tuples:           " << tuples.size() << "\n";
        else
          {
            double total = 0;
            for(list<int>::iterator i = k.begin(); i != k.end(); i++)
      	{
      	  total += nCk(numDivs,*i);
      	}
            cout << "Tuples to evaluate:     " << long(total) << "\n";
          }
      }

  cout << "Statistics:            ";
  if(do_rich) cout << " richness";
  if(do_priv) cout << " private";
  if(do_tuple) cout << " tuples";
  cout << "\n";

    return;
}

/*
 * Everything the sweep needs to know about the dataset as a whole, derived
 * from the counts alone.
 *
 * The filter is the one in Population::recLociDelete: a locus goes if any
 * grouping failed to score more than a fraction tol of its gene copies
 * there, and a grouping that scored nothing anywhere has no say. The
 * ceilings are smallestNj's, taken over the loci that survive, and the count
 * of loci sitting at a grouping's ceiling is kept because it is what tells a
 * user whether one locus is costing them the sweep or a third of the genome.
 */
void ScanResult::resolve(double tol)
{
  const int J = int(groupName.size());

  dropped.assign(size_t(numLoci),0);
  for(int g = 0; g < J; g++)
    {
      if(groupRows[g] == 0) continue;   //no gene copies: no opinion

      const double copies = double(groupRows[g]);
      for(long long l = 0; l < numLoci; l++)
	{
	  if(dropped[size_t(l)]) continue;
	  if(double(missing(g,l))/copies > tol) dropped[size_t(l)] = 1;
	}
    }

  survivors = 0;
  for(long long l = 0; l < numLoci; l++)
    {
      if(!dropped[size_t(l)]) survivors++;
    }

  ceiling.assign(J,0);
  binding.assign(J,0);
  emptyAt.assign(J,0);
  feasibleG = 0;
  feasibleGAll = 0;

  for(int g = 0; g < J; g++)
    {
      int least = 0, atLeast = 0, none = 0;
      bool first = true;

      for(long long l = 0; l < numLoci; l++)
	{
	  if(dropped[size_t(l)]) continue;

	  const int Nj = nj(g,l);
	  if(first || Nj < least) { least = Nj; atLeast = 0; first = false; }
	  if(Nj == least) atLeast++;
	  if(Nj == 0) none++;
	}

      ceiling[g] = least;
      binding[g] = atLeast;
      emptyAt[g] = none;
      if(g == 0 || least < feasibleG) feasibleG = least;

      int all = 0;
      for(long long l = 0; l < numLoci; l++)
	{
	  const int Nj = nj(g,l);
	  if(l == 0 || Nj < all) all = Nj;
	}
      if(g == 0 || all < feasibleGAll) feasibleGAll = all;
    }

  return;
}

/*
 * One pass over the input that keeps no alleles.
 *
 * See ScanResult. The work is deliberately the reader's minus everything
 * expensive: no allele labels are interned, no counts per allele are kept,
 * and nothing per locus survives except one integer per grouping. What comes
 * out is enough to decide which loci survive, what MAX_G resolves to and
 * where the windows fall -- the three things the sweep cannot start without.
 *
 * The parsing mirrors readStructureInto and readVCFInto exactly, including
 * their dimension detection, their filtering by --pops/--exclude-pops and
 * their notion of a gene copy, because a scan that disagreed with the reader
 * about any of that would filter one set of loci and compute over another.
 * ADZE_DUMP_SCAN is how that agreement is checked against ADZE_DUMP_COUNTS.
 */
static int scanGroup(ScanResult& out,const string& label,
		     unordered_map<string,int>& groupOf,long long loci)
{
  unordered_map<string,int>::iterator it = groupOf.find(label);
  if(it != groupOf.end()) return it->second;

  const int index = int(out.groupName.size());
  groupOf.insert(make_pair(label,index));
  out.groupName.push_back(label);
  out.groupRows.push_back(0);
  //A grouping first seen at locus l scored nothing at the loci before it.
  out.observed.push_back(vector<int>(size_t(loci),0));
  return index;
}

static void scanStructure(ParamSet& p,ScanResult& out,LineSource& in,
			  const vector<string>& keepList,
			  const vector<string>& dropList,
			  bool needNames,long long& dataRows)
{
  string line;
  vector<Field> fields;

  if(!in.next(line)) badData("no locus-name row in " + p.dfile.val + ".");
  tokenize(line,fields);

  const int declaredLoci = int(fields.size());
  if(declaredLoci < 1) badData("no locus names in the first row of " + p.dfile.val + ".");

  if(p.loci.set && p.loci.val != declaredLoci)
    {
      adzelog() << "WARNING: LOCI says " << p.loci.val << " but "
		<< p.dfile.val << " has " << declaredLoci
		<< " locus names; using " << declaredLoci << ".\n";
    }
  p.loci.val = declaredLoci;
  out.numLoci = declaredLoci;

  if(needNames)
    {
      out.locusName.reserve(fields.size());
      for(size_t l = 0; l < fields.size(); l++)
	{
	  out.locusName.push_back(string(fields[l].first,fields[l].second));
	}
    }

  for(int skip = 1; skip < p.nd_rows.val; skip++) in.next(line);

  unordered_map<string,int> groupOf;
  int ndCols = p.nd_cols.set ? p.nd_cols.val : 0;
  int groupCol = -1;
  int expected = 0;
  string token;
  long long physicalRows = 0;
  dataRows = 0;

  while(in.next(line))
    {
      physicalRows++;
      tokenize(line,fields);
      if(fields.empty()) continue;

      if(expected == 0)
	{
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

      const int g = scanGroup(out,token,groupOf,out.numLoci);
      out.groupRows[g]++;
      out.geneCopies++;

      vector<int>& seen = out.observed[g];
      for(int l = 0; l < declaredLoci; l++)
	{
	  const Field& f = fields[ndCols + l];
	  token.assign(f.first,f.second);
	  if(token.compare(p.miss.val) != 0) seen[size_t(l)]++;
	}
    }

  return;
}

static void scanVCF(ParamSet& p,ScanResult& out,LineSource& in,
		    const vector<string>& keepList,
		    const vector<string>& dropList,
		    bool announce)
{
  unordered_map<string,string> groupOfSample;
  vector<string> groupOrder;
  readSampleMap(p.samples.val,groupOfSample,groupOrder);

  unordered_map<string,int> groupOf;
  for(size_t i = 0; i < groupOrder.size(); i++)
    {
      const string& g = groupOrder[i];
      if(!keepList.empty() &&
	 find(keepList.begin(),keepList.end(),g) == keepList.end()) continue;
      if(!dropList.empty() &&
	 find(dropList.begin(),dropList.end(),g) != dropList.end()) continue;
      scanGroup(out,g,groupOf,0);
    }

  if(out.groupName.empty())
    {
      badData("no grouping in " + p.samples.val + " survived --pops/--exclude-pops.");
    }

  string line, token;
  vector<Field> fields;
  bool haveHeader = false;
  vector<int> sampleGroup;
  vector<int> samplePloidy;
  int usedSamples = 0;
  long long records = 0;

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

	  for(size_t c = 9; c < fields.size(); c++)
	    {
	      token.assign(fields[c].first,fields[c].second);
	      unordered_map<string,string>::const_iterator it = groupOfSample.find(token);
	      if(it == groupOfSample.end()) continue;

	      unordered_map<string,int>::const_iterator gi = groupOf.find(it->second);
	      if(gi == groupOf.end()) continue;

	      sampleGroup[c-9] = gi->second;
	      usedSamples++;
	    }

	  if(usedSamples == 0)
	    {
	      badData("none of the samples in " + p.dfile.val + " appears in " +
		      p.samples.val + " under a grouping being analysed.");
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

      {
	const string chrom(fields[0].first,fields[0].second);
	const string posText(fields[1].first,fields[1].second);
	out.lmap.chrom.push_back(out.lmap.chromIndex(chrom));
	out.lmap.pos.push_back(atoll(posText.c_str()));
      }

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

      const size_t l = size_t(records);
      records++;
      for(size_t g = 0; g < out.observed.size(); g++) out.observed[g].push_back(0);

      for(size_t c = 0; c < sampleGroup.size(); c++)
	{
	  const int g = sampleGroup[c];
	  if(g < 0) continue;

	  const char* f = fields[c+9].first;
	  const size_t n = fields[c+9].second;

	  size_t i = 0;
	  for(int skip = 0; skip < gtField && i < n; skip++)
	    {
	      while(i < n && f[i] != ':') i++;
	      if(i < n) i++;
	    }
	  size_t stop = i;
	  while(stop < n && f[stop] != ':') stop++;

	  if(i >= stop) continue;
	  if(stop - i == 1 && f[i] == '.') { if(samplePloidy[c] == 0) samplePloidy[c] = 1; continue; }

	  int copies = 0;
	  size_t a = i;
	  while(a < stop)
	    {
	      size_t b = a;
	      while(b < stop && f[b] != '/' && f[b] != '|') b++;

	      copies++;
	      if(!(b - a == 1 && f[a] == '.')) out.observed[g][l]++;

	      a = (b < stop) ? b + 1 : stop;
	    }

	  if(samplePloidy[c] < copies) samplePloidy[c] = copies;
	}
    }

  if(!haveHeader) badData("no #CHROM header line in " + p.dfile.val + ".");
  if(records == 0) badData("no variant records in " + p.dfile.val + ".");

  //A sample's ploidy is the largest GT it carries anywhere, so the totals are
  //only knowable here. A sample with no call at any record counts as diploid,
  //as the reader assumes.
  for(size_t c = 0; c < sampleGroup.size(); c++)
    {
      const int g = sampleGroup[c];
      if(g < 0) continue;
      const int ploidy = (samplePloidy[c] > 0) ? samplePloidy[c] : 2;
      out.groupRows[g] += ploidy;
      out.geneCopies += ploidy;
    }

  if(p.loci.set && p.loci.val != int(records))
    {
      adzelog() << "WARNING: LOCI says " << p.loci.val << " but " << p.dfile.val
		<< " has " << records << " records; using " << records << ".\n";
    }
  p.loci.val = int(records);
  out.numLoci = records;

  //The reader's words, since the scan can replace it for a dry run and the
  //user should not be told something different about the same file.
  if(announce)
    {
      adzelog() << "Read " << records << (records == 1 ? " record" : " records")
		<< " for " << usedSamples
		<< (usedSamples == 1 ? " sample" : " samples") << ".\n";
    }

  return;
}

void scanDataset(ParamSet& p,ScanResult& out,bool announce)
{
  const bool vcf = wantsVCF(p.format.val,p.dfile.val);

  if(vcf && !p.samples.set)
    {
      badData("VCF input needs --samples FILE, a sample-to-grouping map.");
    }
  if(!vcf && p.samples.set)
    {
      adzelog() << "WARNING: --samples is for VCF input and is ignored here.\n";
    }

  LineSource in;
  if(!in.open(p.dfile.val))
    {
      cerr << "ERROR: could not open " << p.dfile.val << "\n";
      throw BAD_FILE();
    }
  if(!vcf && looksCompressed(p.dfile.val))
    {
      badData("compressed input is supported for VCF only.");
    }

  vector<string> keepList, dropList;
  splitList(p.pops.val,keepList);
  splitList(p.expops.val,dropList);

  long long dataRows = 0;
  if(vcf) scanVCF(p,out,in,keepList,dropList,announce);
  else scanStructure(p,out,in,keepList,dropList,p.loci_map.set,dataRows);

  if(!vcf)
    {
      if(p.dlines.set && p.dlines.val != int(dataRows))
	{
	  adzelog() << "WARNING: DATA_LINES says " << p.dlines.val << " but "
		    << p.dfile.val << " has " << dataRows << " data rows; using "
		    << dataRows << ".\n";
	}
      p.dlines.val = int(dataRows);
      if(dataRows == 0) badData("no data rows in " + p.dfile.val + ".");
    }
  else
    {
      p.dlines.val = int(out.geneCopies);
    }

  if(out.groupName.empty()) badData("no grouping survived the filters.");

  if(announce && !vcf && (!keepList.empty() || !dropList.empty()))
    {
      long long kept = 0;
      for(size_t g = 0; g < out.groupRows.size(); g++) kept += out.groupRows[g];
      adzelog() << "Using " << kept << " of " << dataRows
		<< " gene copies in " << out.groupName.size()
		<< ((out.groupName.size() == 1) ? " grouping.\n" : " groupings.\n");
    }

  if(!vcf && p.loci_map.set) readLocusMap(p.loci_map.val,out.locusName,out.lmap);

  out.resolve(p.tol.val);

  //Development facility, the counterpart of ADZE_DUMP_COUNTS: the same two
  //numbers per locus and grouping from a pass that interned nothing, plus
  //what the scan concluded from them, so the filter, the ceilings and the
  //window layout can be checked against the program's own output.
  if(const char* path = getenv("ADZE_DUMP_SCAN"))
    {
      ofstream d(path);
      d << "LOCUS_INDEX\tGROUPING\tNJ\tMISSING\n";
      for(long long l = 0; l < out.numLoci; l++)
	{
	  for(size_t g = 0; g < out.groupName.size(); g++)
	    {
	      d << l << "\t" << out.groupName[g] << "\t" << out.nj(int(g),l)
		<< "\t" << out.missing(int(g),l) << "\n";
	    }
	}

      d << "#SUMMARY\tloci=" << out.numLoci << "\tsurvivors=" << out.survivors
	<< "\tfeasibleG=" << out.feasibleG
	<< "\tfeasibleGAll=" << out.feasibleGAll << "\n";

      for(size_t g = 0; g < out.groupName.size(); g++)
	{
	  d << "#GROUP\t" << out.groupName[g] << "\trows=" << out.groupRows[g]
	    << "\tceiling=" << out.ceiling[g] << "\tbinding=" << out.binding[g]
	    << "\tempty=" << out.emptyAt[g] << "\n";
	}

      //Dropped loci, highest index first: the order 1.0 wrote _deletedloci.
      for(long long l = out.numLoci - 1; l >= 0; l--)
	{
	  if(!out.dropped[size_t(l)]) continue;
	  d << "#DROPPED\t" << l << "\t"
	    << (out.locusName.empty() ? string(".") : out.locusName[size_t(l)]) << "\n";
	}

      //The window layout follows the survivors, as it does after filtering.
      if(out.lmap.size() == size_t(out.numLoci) &&
	 (p.win_bp.set || p.win_loci.set))
	{
	  LocusMap kept = out.lmap;
	  kept.compact(out.dropped);

	  vector<Window> windows;
	  buildWindows(kept,p,windows);
	  for(size_t w = 0; w < windows.size(); w++)
	    {
	      d << "#WINDOW\t" << kept.chromName[windows[w].chrom] << "\t"
		<< windows[w].start << "\t" << windows[w].end << "\t"
		<< windows[w].numLoci() << "\n";
	    }
	}

      d.close();
    }

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
			LocusTable& loci, LocusMap& lmap)
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

  /*
   * With --pops the groupings are listed, so the count stride can be set
   * before anything is read. Without it, STRUCTURE input discovers them row
   * by row and the doubling cap stands; VCF input sets it from the sample
   * map, inside the reader.
   */
  if(!keepList.empty()) acc.setGroupCap(int(keepList.size()));

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

  /*
   * The allele labels are not needed once the file is read: they exist to
   * give an allele its slot, and the slot is what gets published. Released
   * before any grouping's count block is allocated, so the two do not have
   * to coexist.
   */
  for(size_t l = 0; l < acc.locus.size(); l++)
    {
      vector<string>().swap(acc.locus[l].label);
    }

  /*
   * The shared per-locus table: names moved rather than copied, and the
   * allele-slot offsets, which every grouping indexes its counts through.
   */
  loci.name.swap(acc.locusName);
  loci.deleted.clear();
  loci.offset.assign(size_t(declaredLoci)+1,0);
  for(int l = 0; l < declaredLoci; l++)
    {
      loci.offset[l+1] = loci.offset[l] + int(acc.locus[l].firstSeen.size());
    }

  Population* pop = new Population[numDivs];
  for(int j = 0; j < numDivs; j++)
    {
      pop[j].setLoci(declaredLoci);
      pop[j].setName(acc.groupName[j]);
      pop[j].setRows(acc.groupRows[j]);
      pop[j].setTable(&loci);
      pop[j].allocNji();
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
	  for(int i = 0; i < slots; i++)
	    {
	      pop[j].putNji(t.count[size_t(order[i])*acc.groupCap + j],i,l);
	    }
	  pop[j].putMissing(t.missing[j],l);
	}

      //Release this locus's bookkeeping as soon as it is published.
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
		<< "every g; lower --tolerance to drop them.\n";
    }

  return pop;
}

void filterLoci(Population pop[],int numDivs, double tol, string file,
		bool pp, LocusMap& lmap, LocusTable& loci)
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
  
  //Each grouping moves its counts down within its own block, reading from
  //where the loci used to start; the shared offsets are rewritten afterwards.
  const vector<int> oldOffset = loci.offset;

  for(int n = 0; n < numDivs; n++)
    {
      pop[n].deleteLoci(toDelete,oldOffset);
      if(pp) bar.adv(size);
    }
  if(pp) bar.done();

  //Names, offsets and coordinates are indexed by locus, so they follow the
  //same compaction -- once for the run, the filter having condemned each
  //locus in every grouping at once.
  loci.compact(toDelete);
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
  vector<string> groupNames;
  for(int j = 0; j < numDivs; j++) groupNames.push_back(pop[j].getName());
  return readTupleFile(file,groupNames,out);
}

/*
 * The same, from the grouping names alone: the dry run knows them from the
 * scan and has no Population to ask.
 */
bool readTupleFile(const string& file, const vector<string>& groupName,
		   vector< vector<int> >& out)
{
  const int numDivs = int(groupName.size());
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
	      if(groupName[j].compare(names[i]) == 0) found = j;
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
 * A running mean and sum of squared deviations, one pair per g (and per
 * window per g, for the window files).
 *
 * The passes used to fill a (gCeil+1) x numLoci array of doubles and read it
 * twice -- once for the mean, once for the deviations -- which is 8 bytes per
 * locus per g of memory for arithmetic that needs none: 168 MB for a million
 * loci swept to g = 20. Welford's recurrence forms both in one pass, so
 * nothing per-locus is stored. The cost is that the last bits of a variance
 * no longer match version 1.0's two-pass sum; see test/regress.py and the
 * manual's Verification section for what is compared instead.
 *
 * add() is called in ascending locus order, from the ordered region of the
 * pass, so the result does not depend on the thread count.
 */
struct Running
{
  vector<long long> n;
  vector<double> mean, m2;
  vector<char> bad;   //some value was undefined at this g

  void init(size_t slots)
  {
    n.assign(slots,0);
    mean.assign(slots,0.0);
    m2.assign(slots,0.0);
    bad.assign(slots,0);
  }

  void add(size_t i, double x)
  {
    if(x == -9) { bad[i] = 1; return; }

    n[i]++;
    const double d = x - mean[i];
    mean[i] += d/double(n[i]);
    m2[i] += d*(x - mean[i]);
  }

  /*
   * The divisor is the number of loci the statistic covers, not the number
   * of values that turned out to be defined: one undefined locus makes the
   * whole grouping undefined at that g, which is what bad records.
   */
  void into(Stats& st, size_t i, int numLoci) const
  {
    if(bad[i] || numLoci < 1)
      {
	st.putSummary(-9,-9,-9,numLoci);
	return;
      }

    const double var = (numLoci < 2)
      ? numeric_limits<double>::quiet_NaN()
      : m2[i]/double(numLoci-1);

    st.putSummary(mean[i],var,sqrt(var/double(numLoci)),numLoci);
  }
};

/*
 * Window rows, in the order the buffered writer produced them: g outer,
 * window inner. Slot (w,g) of the accumulator is w*gStride + g.
 */
void writeWindowRunning(ostream& out,const Running& acc,int gStride,
			const vector<Window>& windows,const LocusMap& lmap,
			const string& label,int gFrom,int gTo)
{
  for(int g = gFrom; g <= gTo; g++)
    {
      for(size_t w = 0; w < windows.size(); w++)
	{
	  const Window& win = windows[w];

	  Stats st;
	  acc.into(st,w*size_t(gStride) + size_t(g),win.numLoci());

	  ostringstream head;
	  head << lmap.chromName[win.chrom] << '\t' << win.start << '\t'
	       << win.end << '\t' << label;

	  //Tabbed: a window undefined at this g is reported as NA rather
	  //than dropped, so a scan keeps one row per window per g.
	  st.printStats(out,head.str(),g,1);
	}
    }

  return;
}

/*
 * The windows covering each locus, as the sweep reaches it. Windows are laid
 * out in ascending locus order and may overlap, so the set is maintained by
 * opening windows whose first locus has arrived and closing those whose last
 * has passed.
 */
struct WindowCursor
{
  const vector<Window>* windows;
  size_t next;
  vector<size_t> open;

  void init(const vector<Window>& w) { windows = &w; next = 0; open.clear(); }

  const vector<size_t>& at(int locus)
  {
    while(next < windows->size() && (*windows)[next].first <= locus)
      {
	open.push_back(next);
	next++;
      }

    size_t keep = 0;
    for(size_t i = 0; i < open.size(); i++)
      {
	if((*windows)[open[i]].last > locus) open[keep++] = open[i];
      }
    open.resize(keep);

    return open;
  }
};

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

      //Header once the g ceiling is known, below.
    }

  reg_out.open(comb_out.c_str());
  if(param.tabbed())
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

  /*
   * One row per g, or a single row under --at-g: the recurrence in g still
   * has to climb to the requested g, but the values below it are read by
   * nothing. gBase is the g that row 0 holds.
   */
  //See calcAllAgs: nothing per-locus is stored.
  const int gFrom = param.at_g_val ? param.at_g_val : 1;
  const int gTo = param.at_g_val ? param.at_g_val : gCeil;

  if(full_comb) writeFullDataHeader(full_out,"TUPLE",gFrom,gTo);

  Running acc, wacc;
  const int BLOCK = 4096;
  const int nBlocks = (numLoci + BLOCK - 1)/BLOCK;

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
	combineNames(&names[0],k,param.tabbed() ? ',' : ' ');

      acc.init(size_t(gStride));
      if(!windows.empty()) wacc.init(windows.size()*size_t(gStride));

      WindowCursor cursor;
      cursor.init(windows);

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
      vector<double> q;
      vector<double> vals(size_t(BLOCK)*gStride);

#ifdef _OPENMP
#pragma omp for ordered schedule(static,1)
#endif
      for(int b = 0; b < nBlocks; b++)
	{
	  const int lo = b*BLOCK;
	  const int hi = (lo + BLOCK < numLoci) ? lo + BLOCK : numLoci;

	  for(int locus = lo; locus < hi; locus++)
	    {
	      const int numAlleles = pop[0].getNjiColLength(locus);
	      buildQTable(pop,numDivs,locus,numAlleles,gCeil,gStride,q);

	      for(int g = gFrom; g <= gTo; g++)
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

		  vals[size_t(locus-lo)*gStride + g] = pg;
		}
	    }

#ifdef _OPENMP
#pragma omp ordered
#endif
	  {
	  for(int locus = lo; locus < hi; locus++)
	    {
	      for(int g = gFrom; g <= gTo; g++)
		{
		  acc.add(size_t(g),vals[size_t(locus-lo)*gStride + g]);
		}

	      if(!windows.empty())
		{
		  const vector<size_t>& here = cursor.at(locus);
		  for(size_t i = 0; i < here.size(); i++)
		    {
		      for(int g = gFrom; g <= gTo; g++)
			{
			  wacc.add(here[i]*size_t(gStride) + size_t(g),
				   vals[size_t(locus-lo)*gStride + g]);
			}
		    }
		}
	    }

	  if(full_comb)
	    {
	      //The tuple label is comma-joined here so it stays one field.
	      writeFullDataRows(full_out,combineNames(&names[0],k,','),vals,lo,hi,
				gFrom,gTo,gStride,pop[0]);
	    }

	  if(param.pp.val) bar.adv((hi-lo)*(param.at_g_val ? 1 : gLast));
	  }
	}
      } //end parallel region

      for(int g = gFrom; g <= gTo; g++)
	{
	  Stats comb_stats;
	  acc.into(comb_stats,size_t(g),numLoci);
	  comb_stats.printStats(reg_out,all_names,g,param.tabbed());
	}

      if(!windows.empty())
	{
	  /*
	   * The tuple label is comma-joined here whatever the genome-wide
	   * format, since window output is always tab-separated and the label
	   * has to stay one field.
	   */
	  writeWindowRunning(comb_win_out,wacc,gStride,windows,lmap,
			     combineNames(&names[0],k,','),gFrom,gTo);
	}

      if(!param.tabbed()) reg_out << endl;
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

/*
 * Every count the reader produced, as read: one row per locus and grouping,
 * with that grouping's total, its missing copies, and the per-allele counts
 * in slot order.
 *
 * A development facility, switched on by ADZE_DUMP_COUNTS=<file> rather than
 * by a flag, because it exists to compare two builds' reading and not to be
 * part of the interface. Two readers that produce the same dump produce the
 * same statistics by construction; without it, a counting difference only
 * surfaces as a difference in a mean, with no way to say which locus.
 */
void dumpCounts(const char* path, Population pop[], int numDivs, int numLoci,
		const LocusMap& lmap)
{
  ofstream out(path);
  if(!out.is_open())
    {
      adzelog() << "WARNING: could not write the count dump to " << path << "\n";
      return;
    }

  out << "LOCUS\tCHROM\tPOS\tGROUPING\tNJ\tMISSING\tSLOTS\tNJI\n";

  for(int l = 0; l < numLoci; l++)
    {
      const bool placed = (lmap.size() == size_t(numLoci));

      for(int j = 0; j < numDivs; j++)
	{
	  const int slots = pop[j].getNjiColLength(l);

	  out << pop[0].getLocusName(l) << "\t"
	      << (placed ? lmap.chromName[lmap.chrom[l]] : ".") << "\t"
	      << (placed ? lmap.pos[l] : -1) << "\t"
	      << pop[j].getName() << "\t"
	      << pop[j].getNj(l) << "\t"
	      << pop[j].getMissing(l) << "\t"
	      << slots << "\t";

	  for(int i = 0; i < slots; i++)
	    {
	      if(i) out << ",";
	      out << pop[j].getNji(i,l);
	    }

	  out << "\n";
	}
    }

  out.close();
  return;
}

/*
 * The _fulldata layout: one row per grouping and locus, one column per g.
 *
 * Version 1.0 wrote the transpose of this -- a row per (grouping, g) holding
 * every locus across it -- which put the locus names in the header and made
 * the file as wide as the dataset: a million-locus run produced a handful of
 * rows of a million fields, which no ordinary tool will read, and which
 * cannot be written until every locus has been computed. One row per locus
 * streams, stays within column limits, and joins to a locus table the way a
 * reader expects. See the manual, "Changes from version 1.0".
 *
 * A value undefined at that g is NA, in both output modes; the row is kept
 * either way, so a locus that is undefined everywhere is still visible.
 */
/*
 * The _fulldata header. locusMajor puts the locus first, matching the row
 * order the streaming sweep writes: a locus's values for every grouping
 * arrive together and nothing comes back to them (doc/streaming.md,
 * decision 1).
 */
void writeFullDataHeader(ostream& out,const string& labelCols,int gFrom,int gTo,
			 bool locusMajor)
{
  if(locusMajor) out << "LOCUS\t" << labelCols;
  else out << labelCols << "\tLOCUS";
  for(int g = gFrom; g <= gTo; g++) out << "\tG" << g;
  out << "\n";
  return;
}

void writeFullDataRows(ostream& out,const string& label,const vector<double>& vals,
		       int lo,int hi,int gFrom,int gTo,int gStride,
		       const Population& pop)
{
  for(int l = lo; l < hi; l++)
    {
      out << label << "\t" << pop.getLocusName(l);

      for(int g = gFrom; g <= gTo; g++)
	{
	  const double v = vals[size_t(l-lo)*gStride + g];
	  out << "\t";
	  if(v == -9) out << "NA";
	  else out << v;
	}

      out << "\n";
    }

  return;
}

void writeWindowStats(ostream& out, vector<double>& perLocus, int numLoci,
		      const vector<Window>& windows, const LocusMap& lmap,
		      const string& label, int gFirst, int gLast, int gBase)
{
  //gBase is the g held by row 0: 0 for a full ladder, and the requested g
  //when --at-g means only that row was kept.
  for(int g = gFirst; g <= gLast; g++)
    {
      for(size_t w = 0; w < windows.size(); w++)
	{
	  const Window& win = windows[w];

	  Stats st;
	  st.putData(&perLocus[size_t(g-gBase)*numLoci + win.first],win.numLoci());
	  st.calcAvg();
	  st.calcVar();
	  st.calcStdErr();

	  ostringstream head;
	  head << lmap.chromName[win.chrom] << '\t' << win.start << '\t'
	       << win.end << '\t' << label;

	  //Tabbed: a window undefined at this g is reported as NA rather
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
      
      //gLast is not known yet; the header is written once the ceiling is.
    }

  pg_out.open(private_out.c_str());
  if(param.tabbed()) pg_out << "POP_GROUPING\tG\tNUM_LOCI\tMEAN\tVAR\tSTD_ERR\n";

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

  //See calcAllAgs: nothing per-locus is stored.
  const int gFrom = param.at_g_val ? param.at_g_val : 1;
  const int gTo = param.at_g_val ? param.at_g_val : gCeil;

  if(full_priv) writeFullDataHeader(pg_full_out,"POP_GROUPING",gFrom,gTo);

  Running acc, wacc;
  const int BLOCK = 4096;
  const int nBlocks = (numLoci + BLOCK - 1)/BLOCK;

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
      acc.init(size_t(gStride));
      if(!windows.empty()) wacc.init(windows.size()*size_t(gStride));

      WindowCursor cursor;
      cursor.init(windows);

#ifdef _OPENMP
#pragma omp parallel
#endif
      {
      vector<double> q;
      vector<double> vals(size_t(BLOCK)*gStride);

#ifdef _OPENMP
#pragma omp for ordered schedule(static,1)
#endif
      for(int b = 0; b < nBlocks; b++)
	{
	  const int lo = b*BLOCK;
	  const int hi = (lo + BLOCK < numLoci) ? lo + BLOCK : numLoci;

	  for(int locus = lo; locus < hi; locus++)
	    {
	      const int numAlleles = pop[j].getNjiColLength(locus);
	      buildQTable(pop,numDivs,locus,numAlleles,gCeil,gStride,q);

	      for(int g = gFrom; g <= gTo; g++)
		{
		  if(minNjLocus[locus] < g)
		    {
		      vals[size_t(locus-lo)*gStride + g] = -9;
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

		  vals[size_t(locus-lo)*gStride + g] = total;
		}
	    }

#ifdef _OPENMP
#pragma omp ordered
#endif
	  {
	  for(int locus = lo; locus < hi; locus++)
	    {
	      for(int g = gFrom; g <= gTo; g++)
		{
		  acc.add(size_t(g),vals[size_t(locus-lo)*gStride + g]);
		}

	      if(!windows.empty())
		{
		  const vector<size_t>& here = cursor.at(locus);
		  for(size_t i = 0; i < here.size(); i++)
		    {
		      for(int g = gFrom; g <= gTo; g++)
			{
			  wacc.add(here[i]*size_t(gStride) + size_t(g),
				   vals[size_t(locus-lo)*gStride + g]);
			}
		    }
		}
	    }

	  if(full_priv)
	    {
	      writeFullDataRows(pg_full_out,pop[j].getName(),vals,lo,hi,
				gFrom,gTo,gStride,pop[0]);
	    }

	  if(param.pp.val) bar.adv((hi-lo)*(param.at_g_val ? 1 : gLast));
	  }
	}
      } //end parallel region

      for(int g = gFrom; g <= gTo; g++)
	{
	  Stats pg_stats;
	  acc.into(pg_stats,size_t(g),numLoci);
	  pg_stats.printStats(pg_out,pop[j].getName(),g,param.tabbed());
	}

      if(!windows.empty())
	{
	  writeWindowRunning(pg_win_out,wacc,gStride,windows,lmap,
			     pop[j].getName(),gFrom,gTo);
	}

      if(!param.tabbed()) pg_out << endl;
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

      writeFullDataHeader(ag_full_out,"POP_GROUPING",
			  param.at_g_val ? param.at_g_val : 1,
			  param.at_g_val ? param.at_g_val : gTop);
    }

  ag_out.open(richness_out.c_str());
  if(param.tabbed()) ag_out << "POP_GROUPING\tG\tNUM_LOCI\tMEAN\tVAR\tSTD_ERR\n";

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

  /*
   * Nothing per-locus is stored: each g carries a running mean and sum of
   * squared deviations, updated as the sweep passes each locus, and the
   * windows carry one such pair each. --at-g narrows the sweep to the one g
   * asked for; the recurrence in g still has to climb to it.
   */
  const int gFrom = param.at_g_val ? param.at_g_val : 1;
  const int gTo = param.at_g_val ? param.at_g_val : gCeil;

  Running acc, wacc;
  const int BLOCK = 4096;
  const int nBlocks = (numLoci + BLOCK - 1)/BLOCK;

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
   * grouping's sample size are undefined, and one such locus makes the
   * grouping undefined at that g.
   */
  for(int j = 0; j < numDivs; j++)
    {
      acc.init(size_t(gStride));
      if(!windows.empty()) wacc.init(windows.size()*size_t(gStride));

      WindowCursor cursor;
      cursor.init(windows);

      /*
       * Loci are computed in parallel in blocks, and every accumulation
       * happens in the ordered region below, in ascending locus order -- so
       * a reported value does not depend on the thread count.
       */
#ifdef _OPENMP
#pragma omp parallel
#endif
      {
      vector<double> q;                       //one per thread
      vector<double> vals(size_t(BLOCK)*gStride);

#ifdef _OPENMP
#pragma omp for ordered schedule(static,1)
#endif
      for(int b = 0; b < nBlocks; b++)
	{
	  const int lo = b*BLOCK;
	  const int hi = (lo + BLOCK < numLoci) ? lo + BLOCK : numLoci;

	  for(int locus = lo; locus < hi; locus++)
	    {
	      const int numAlleles = pop[j].getNjiColLength(locus);
	      const int Nj = pop[j].getNj(locus);

	      buildQTable(&pop[j],1,locus,numAlleles,gCeil,gStride,q);

	      for(int g = gFrom; g <= gTo; g++)
		{
		  if(g > Nj)
		    {
		      vals[size_t(locus-lo)*gStride + g] = -9;
		      continue;
		    }

		  double total = 0;
		  for(int i = 0; i < numAlleles; i++)
		    {
		      total += 1 - q[size_t(i)*gStride + g];
		    }

		  vals[size_t(locus-lo)*gStride + g] = total;
		}
	    }

#ifdef _OPENMP
#pragma omp ordered
#endif
	  {
	  for(int locus = lo; locus < hi; locus++)
	    {
	      for(int g = gFrom; g <= gTo; g++)
		{
		  acc.add(size_t(g),vals[size_t(locus-lo)*gStride + g]);
		}

	      if(!windows.empty())
		{
		  const vector<size_t>& here = cursor.at(locus);
		  for(size_t i = 0; i < here.size(); i++)
		    {
		      for(int g = gFrom; g <= gTo; g++)
			{
			  wacc.add(here[i]*size_t(gStride) + size_t(g),
				   vals[size_t(locus-lo)*gStride + g]);
			}
		    }
		}
	    }

	  if(full_rich)
	    {
	      writeFullDataRows(ag_full_out,pop[j].getName(),vals,lo,hi,
				gFrom,gTo,gStride,pop[0]);
	    }

	  if(param.pp.val) bar.adv((hi-lo)*(param.at_g_val ? 1 : gTop));
	  }
	}
      } //end parallel region

      for(int g = gFrom; g <= gTo; g++)
	{
	  Stats ag_stats;
	  acc.into(ag_stats,size_t(g),numLoci);
	  ag_stats.printStats(ag_out,pop[j].getName(),g,param.tabbed());
	}

      if(!windows.empty())
	{
	  writeWindowRunning(ag_win_out,wacc,gStride,windows,lmap,
			     pop[j].getName(),gFrom,gTo);
	}

      if(!param.tabbed()) ag_out << endl;
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

string deletedHeader(long long gone,double tol)
{
  ostringstream out;
  out << gone << ((gone == 1) ? " locus has " : " loci have ")
      << "at least one grouping with";
  if(tol > 0) out << " more than " << 100*tol << "%";
  out << " missing data.\n";
  return out.str();
}

/*
 * The deleted-loci report when there is nothing left to sweep.
 *
 * Normally the sweep writes this file as the dropped loci go past, but a run
 * where every locus is dropped stops before the sweep -- and that is exactly
 * the run whose user most needs the list. The names are read again rather
 * than kept: a STRUCTURE file names its loci in the header row, and a VCF is
 * re-read, which is affordable on a path that is about to exit.
 */
void writeDeletedLoci(ParamSet& p,const ScanResult& scan,const string& out)
{
  const string name = nameCreate(out,"_deletedloci");
  ofstream d(name.c_str());
  if(!d.is_open()) return;

  d << deletedHeader(scan.numLoci - scan.survivors,p.tol.val);

  if(wantsVCF(p.format.val,p.dfile.val))
    {
      LocusSource* src = openVCFSource(p,scan);
      LocusCounts locus;
      long long l = 0;
      while(src->next(locus))
	{
	  if(size_t(l) < scan.dropped.size() && scan.dropped[size_t(l)])
	    {
	      d << locus.name << "\n";
	    }
	  l++;
	}
      delete src;
    }
  else
    {
      LineSource in;
      if(!in.open(p.dfile.val)) return;

      string line;
      if(!in.next(line)) return;

      vector<Field> fields;
      tokenize(line,fields);
      for(long long l = 0; l < scan.numLoci; l++)
	{
	  if(size_t(l) >= scan.dropped.size() || !scan.dropped[size_t(l)]) continue;
	  if(size_t(l) < fields.size())
	    {
	      d << string(fields[size_t(l)].first,fields[size_t(l)].second) << "\n";
	    }
	}
    }

  d.close();
  return;
}

/*
 * One locus's rows in the _fulldata file: a row per label, the locus first.
 * The streaming sweep has a locus's values for every grouping at once and
 * nothing afterwards, so this is where the file's row order comes from --
 * doc/streaming.md, decision 1.
 */
void writeFullDataRow(ostream& out,const string& locus,
		      const vector<string>& label,const vector<double>& vals,
		      int gFrom,int gTo,int gStride)
{
  for(size_t x = 0; x < label.size(); x++)
    {
      out << locus << '\t' << label[x];
      for(int g = gFrom; g <= gTo; g++)
	{
	  const double v = vals[x*size_t(gStride) + size_t(g)];
	  out << '\t';
	  if(v == -9) out << "NA";
	  else out << v;
	}
      out << '\n';
    }
  return;
}


/*
 * How many tuples one sweep can carry.
 *
 * Every tuple needs a running accumulator per g, and another per window per
 * g, all live at once -- so a k range over many groupings can ask for more
 * memory in accumulators than the whole rest of the run uses: all k at 20
 * groupings is a million tuples, and at g <= 20 that is 550 MB, against the
 * 97 MB the old one-k-at-a-time passes peaked at. The tuple list is
 * therefore swept in batches that fit a budget. One batch is the common
 * case; each extra batch costs one more pass over the loci, which for
 * STRUCTURE input is a re-read of the converted counts rather than of the
 * data.
 */
long long tupleBatch(const ParamSet& param,const ScanResult& scan,
		     size_t numWindows,long long budget)
{
  int gTop = param.g.val;
  if(scan.feasibleG < gTop) gTop = scan.feasibleG;
  if(gTop < 1) gTop = 1;
  const long long gStride = (param.at_g_val ? param.at_g_val : gTop) + 1;

  //Running holds n, mean and m2 per slot, plus a byte of "some value was
  //undefined here"; one set per tuple and one per tuple per window.
  const long long perSlot = 8 + 8 + 8 + 1;
  const long long perTuple = perSlot*gStride*(1 + (long long)(numWindows));

  long long n = budget/(perTuple > 0 ? perTuple : 1);
  if(n < 1) n = 1;
  return n;
}

/*
 * The sweep, over a stream of loci.
 *
 * Where calcAllAgs, calcAllPgs and calcPgTuples each walked the whole
 * dataset once -- and each rebuilt the per-locus Q table as it went, the
 * tuple pass rebuilding it once per tuple -- this walks it once in total and
 * builds that table once per locus. Everything reads from the same table:
 * allelic richness per grouping, private allelic richness per grouping, and
 * each requested tuple.
 *
 * Nothing per-locus is kept. Each statistic carries a running mean and sum
 * of squared deviations per g, and per window per g, exactly as the passes
 * do now; the per-locus files stream as the loci go by, which is why their
 * rows are locus-major (doc/streaming.md, decision 1).
 *
 * The arithmetic is the passes' arithmetic, term for term and in the same
 * order, because the point of the rewrite is that no number changes:
 *
 *   richness   undefined where g exceeds that grouping's Nj at the locus
 *   private    undefined where g exceeds the smallest Nj at the locus
 *   tuples     never marked undefined; their ladder is capped instead, at
 *              the smallest Nj anywhere
 */
void sweepLoci(LocusSource& src,const ScanResult& scan,const ParamSet& param,
	       const vector<Window>& windows,const LocusMap& lmap,
	       const vector< vector<int> >& tuples,
	       const vector<int>& tupleFile,const vector<string>& tupleOut,
	       bool do_rich,bool do_priv,bool do_tuple,
	       const string& richness_out,const string& private_out,
	       bool full_rich,bool full_priv,bool full_comb,
	       const string& deleted_out,bool appendTuples)
{
  const int J = int(scan.groupName.size());
  const int numLoci = int(scan.survivors);
  const int T = int(tuples.size());

  //Each statistic's ladder, as each pass computes it today.
  const int gRich = param.g.val;
  int gPair = param.g.val;
  if(scan.feasibleG < gPair) gPair = scan.feasibleG;
  if(gPair < 1) gPair = 1;

  const int gCeilRich = param.at_g_val ? param.at_g_val : gRich;
  const int gCeilPair = param.at_g_val ? param.at_g_val : gPair;
  const int gCeil = (gCeilRich > gCeilPair) ? gCeilRich : gCeilPair;
  const int gStride = gCeil + 1;

  const int gFrom = param.at_g_val ? param.at_g_val : 1;
  const int gToRich = param.at_g_val ? param.at_g_val : gCeilRich;
  const int gToPair = param.at_g_val ? param.at_g_val : gCeilPair;

  //Output files, opened in the order the passes open them.
  ofstream ag_out,ag_full_out,ag_win_out;
  ofstream pg_out,pg_full_out,pg_win_out;

  if(do_rich)
    {
      ag_out.open(richness_out.c_str());
      if(param.tabbed()) ag_out << "POP_GROUPING\tG\tNUM_LOCI\tMEAN\tVAR\tSTD_ERR\n";
      if(full_rich)
	{
	  const string name = nameCreate(richness_out,"_fulldata");
	  ag_full_out.open(name.c_str());
	  writeFullDataHeader(ag_full_out,"POP_GROUPING",gFrom,gToRich,true);
	}
      if(!windows.empty())
	{
	  const string name = nameCreate(richness_out,"_windows");
	  ag_win_out.open(name.c_str());
	  writeWindowHeader(ag_win_out,"POP_GROUPING");
	}
    }

  if(do_priv)
    {
      pg_out.open(private_out.c_str());
      if(param.tabbed()) pg_out << "POP_GROUPING\tG\tNUM_LOCI\tMEAN\tVAR\tSTD_ERR\n";
      if(full_priv)
	{
	  const string name = nameCreate(private_out,"_fulldata");
	  pg_full_out.open(name.c_str());
	  writeFullDataHeader(pg_full_out,"POP_GROUPING",gFrom,gToPair,true);
	}
      if(!windows.empty())
	{
	  const string name = nameCreate(private_out,"_windows");
	  pg_win_out.open(name.c_str());
	  writeWindowHeader(pg_win_out,"POP_GROUPING");
	}
    }

  /*
   * A k range writes one file per k, so the tuples of a run belong to
   * several output files. They are all accumulated in the same pass and
   * sorted into their files when the sweep ends.
   */
  const int F = int(tupleOut.size());
  vector<ofstream*> cbOut(F,(ofstream*)0), cbFull(F,(ofstream*)0), cbWin(F,(ofstream*)0);

  if(do_tuple && T > 0)
    {
      /*
       * A later batch of tuples appends to the files the first batch opened
       * and writes no second header: the batches are consecutive slices of
       * the same tuple list, so the rows land in the order they would have
       * had in one pass.
       */
      const ios::openmode mode = appendTuples ? (ios::out|ios::app) : ios::out;

      for(int f = 0; f < F; f++)
	{
	  cbOut[f] = new ofstream(tupleOut[f].c_str(),mode);
	  if(param.tabbed() && !appendTuples)
	    {
	      *cbOut[f] << "TUPLE\tG\tNUM_LOCI\tMEAN\tVAR\tSTD_ERR\n";
	    }
	  if(full_comb)
	    {
	      const string name = nameCreate(tupleOut[f],"_fulldata");
	      cbFull[f] = new ofstream(name.c_str(),mode);
	      if(!appendTuples) writeFullDataHeader(*cbFull[f],"TUPLE",gFrom,gToPair,true);
	    }
	  if(!windows.empty())
	    {
	      const string name = nameCreate(tupleOut[f],"_windows");
	      cbWin[f] = new ofstream(name.c_str(),mode);
	      if(!appendTuples) writeWindowHeader(*cbWin[f],"TUPLE");
	    }
	}
    }

  //One accumulator per grouping, per grouping again, and per tuple.
  vector<Running> rich(J), priv(J), comb(T);
  vector<Running> richWin(J), privWin(J), combWin(T);
  for(int j = 0; j < J; j++)
    {
      rich[j].init(size_t(gStride));
      priv[j].init(size_t(gStride));
      if(!windows.empty())
	{
	  richWin[j].init(windows.size()*size_t(gStride));
	  privWin[j].init(windows.size()*size_t(gStride));
	}
    }
  for(int m = 0; m < T; m++)
    {
      comb[m].init(size_t(gStride));
      if(!windows.empty()) combWin[m].init(windows.size()*size_t(gStride));
    }

  //Tuple membership, hoisted out of the locus loop as the tuple pass does.
  vector< vector<char> > inTuple(T,vector<char>(J,0));
  vector<string> tupleLabel(T), tupleLabelComma(T);
  {
    vector<string> names;
    for(int m = 0; m < T; m++)
      {
	names.assign(tuples[m].size(),string());
	for(size_t x = 0; x < tuples[m].size(); x++)
	  {
	    names[x] = scan.groupName[tuples[m][x]];
	    inTuple[m][tuples[m][x]] = 1;
	  }

	/*
	 * Two spellings, as the passes wrote them: the statistics file
	 * follows the layout's separator, so a legacy space-separated file
	 * gets a space-joined label; the _fulldata and window files are
	 * tab-separated whatever the layout, so their label is always
	 * comma-joined and stays one field.
	 */
	tupleLabel[m] = combineNames(&names[0],int(names.size()),
				     param.tabbed() ? ',' : ' ');
	tupleLabelComma[m] = combineNames(&names[0],int(names.size()),',');
      }
  }

  /*
   * The dropped loci are named as they go past, not collected and written at
   * the end: the scan knows how many there are but not what they are called,
   * and holding a name per dropped locus would put back a per-locus cost
   * this rewrite exists to remove. They therefore appear in file order --
   * doc/streaming.md, decision 4.
   */
  ofstream del_out;
  if(!deleted_out.empty())
    {
      const string name = nameCreate(deleted_out,"_deletedloci");
      del_out.open(name.c_str());
      del_out << deletedHeader(scan.numLoci - scan.survivors,param.tol.val);
    }

  WindowCursor richCursor, privCursor, combCursor;
  richCursor.init(windows);
  privCursor.init(windows);
  combCursor.init(windows);

  ProgressBar bar(&adzelog(),double(numLoci),BARLEN[0]);
  if(param.pp.val) bar.init();

  vector<double> q;
  vector<double> vRich(size_t(J)*gStride), vPriv(size_t(J)*gStride),
    vComb(size_t(T > 0 ? T : 1)*gStride);
  vector<double> rowRich(size_t(J)*gStride), rowPriv(size_t(J)*gStride);

  LocusCounts locus;
  long long index = 0;        //position among the surviving loci
  long long seen = 0;         //position in the file, dropped loci included

  while(src.next(locus))
    {
      const bool keep = (size_t(seen) >= scan.dropped.size()) || !scan.dropped[size_t(seen)];
      seen++;
      if(!keep)
	{
	  if(del_out.is_open()) del_out << locus.name << "\n";
	  continue;
	}

      const int slots = locus.slots;
      const int here = int(index);
      index++;

      //Q for every grouping at this locus, once.
      q.assign(size_t(J)*slots*gStride,0.0);
      int minNj = 0;
      for(int j = 0; j < J; j++)
	{
	  const int Nj = locus.nj[j];
	  if(j == 0 || Nj < minNj) minNj = Nj;

	  const int gTop = (gCeil < Nj) ? gCeil : Nj;
	  for(int i = 0; i < slots; i++)
	    {
	      const int Nji = locus.count[size_t(i)*J + j];
	      double* qpi = &q[(size_t(j)*slots + i)*gStride];
	      double Q = 1;
	      for(int g = 1; g <= gTop; g++)
		{
		  Q *= double(Nj - Nji - (g-1))/double(Nj - (g-1));
		  qpi[g] = Q;
		}
	    }
	}

      if(do_rich)
	{
	  for(int j = 0; j < J; j++)
	    {
	      const int Nj = locus.nj[j];
	      for(int g = gFrom; g <= gToRich; g++)
		{
		  double v = -9;
		  if(g <= Nj)
		    {
		      double total = 0;
		      for(int i = 0; i < slots; i++)
			{
			  total += 1 - q[(size_t(j)*slots + i)*gStride + g];
			}
		      v = total;
		    }
		  rowRich[size_t(j)*gStride + g] = v;
		  rich[j].add(size_t(g),v);
		}
	    }
	}

      if(do_priv)
	{
	  for(int j = 0; j < J; j++)
	    {
	      for(int g = gFrom; g <= gToPair; g++)
		{
		  double v = -9;
		  if(minNj >= g)
		    {
		      double total = 0;
		      for(int i = 0; i < slots; i++)
			{
			  double Q = 1;
			  for(int pj = 0; pj < J; pj++)
			    {
			      if(pj != j) Q *= q[(size_t(pj)*slots + i)*gStride + g];
			    }
			  const double P = 1 - q[(size_t(j)*slots + i)*gStride + g];
			  total += P*Q;
			}
		      v = total;
		    }
		  rowPriv[size_t(j)*gStride + g] = v;
		  priv[j].add(size_t(g),v);
		}
	    }
	}

      if(do_tuple && T > 0)
	{
	  for(int m = 0; m < T; m++)
	    {
	      const vector<int>& tuple = tuples[m];
	      const int k = int(tuple.size());

	      for(int g = gFrom; g <= gToPair; g++)
		{
		  double pg = 0;
		  for(int i = 0; i < slots; i++)
		    {
		      double P = 1, Q = 1;
		      for(int x = 0; x < k; x++)
			{
			  P *= (1 - q[(size_t(tuple[x])*slots + i)*gStride + g]);
			}
		      for(int pj = 0; pj < J; pj++)
			{
			  if(!inTuple[m][pj]) Q *= q[(size_t(pj)*slots + i)*gStride + g];
			}
		      pg += (P*Q);
		    }
		  vComb[size_t(m)*gStride + g] = pg;
		  comb[m].add(size_t(g),pg);
		}
	    }
	}

      //Windows: the same value, into whichever windows cover this locus.
      if(!windows.empty())
	{
	  if(do_rich)
	    {
	      const vector<size_t>& open = richCursor.at(here);
	      for(size_t w = 0; w < open.size(); w++)
		{
		  for(int j = 0; j < J; j++)
		    {
		      for(int g = gFrom; g <= gToRich; g++)
			{
			  richWin[j].add(open[w]*size_t(gStride) + size_t(g),
					 rowRich[size_t(j)*gStride + g]);
			}
		    }
		}
	    }
	  if(do_priv)
	    {
	      const vector<size_t>& open = privCursor.at(here);
	      for(size_t w = 0; w < open.size(); w++)
		{
		  for(int j = 0; j < J; j++)
		    {
		      for(int g = gFrom; g <= gToPair; g++)
			{
			  privWin[j].add(open[w]*size_t(gStride) + size_t(g),
					 rowPriv[size_t(j)*gStride + g]);
			}
		    }
		}
	    }
	  if(do_tuple && T > 0)
	    {
	      const vector<size_t>& open = combCursor.at(here);
	      for(size_t w = 0; w < open.size(); w++)
		{
		  for(int m = 0; m < T; m++)
		    {
		      for(int g = gFrom; g <= gToPair; g++)
			{
			  combWin[m].add(open[w]*size_t(gStride) + size_t(g),
					 vComb[size_t(m)*gStride + g]);
			}
		    }
		}
	    }
	}

      //Per-locus rows, as the loci go by.
      if(do_rich && full_rich)
	{
	  writeFullDataRow(ag_full_out,locus.name,scan.groupName,rowRich,
			   gFrom,gToRich,gStride);
	}
      if(do_priv && full_priv)
	{
	  writeFullDataRow(pg_full_out,locus.name,scan.groupName,rowPriv,
			   gFrom,gToPair,gStride);
	}
      if(do_tuple && full_comb && T > 0)
	{
	  for(int f = 0; f < F; f++)
	    {
	      for(int m = 0; m < T; m++)
		{
		  if(tupleFile[m] != f) continue;
		  *cbFull[f] << locus.name << '\t' << tupleLabelComma[m];
		  for(int g = gFrom; g <= gToPair; g++)
		    {
		      const double v = vComb[size_t(m)*gStride + g];
		      *cbFull[f] << '\t';
		      if(v == -9) *cbFull[f] << "NA";
		      else *cbFull[f] << v;
		    }
		  *cbFull[f] << '\n';
		}
	    }
	}

      if(param.pp.val) bar.adv(1);
    }

  if(param.pp.val) bar.done();

  //Summaries, in the order the passes wrote them: grouping outer, g inner.
  for(int j = 0; do_rich && j < J; j++)
    {
      for(int g = gFrom; g <= gToRich; g++)
	{
	  Stats st;
	  rich[j].into(st,size_t(g),numLoci);
	  st.printStats(ag_out,scan.groupName[j],g,param.tabbed());
	}
      if(!windows.empty())
	{
	  writeWindowRunning(ag_win_out,richWin[j],gStride,windows,lmap,
			     scan.groupName[j],gFrom,gToRich);
	}
      if(!param.tabbed()) ag_out << endl;
    }

  for(int j = 0; do_priv && j < J; j++)
    {
      for(int g = gFrom; g <= gToPair; g++)
	{
	  Stats st;
	  priv[j].into(st,size_t(g),numLoci);
	  st.printStats(pg_out,scan.groupName[j],g,param.tabbed());
	}
      if(!windows.empty())
	{
	  writeWindowRunning(pg_win_out,privWin[j],gStride,windows,lmap,
			     scan.groupName[j],gFrom,gToPair);
	}
      if(!param.tabbed()) pg_out << endl;
    }

  for(int f = 0; do_tuple && f < F; f++)
    {
      for(int m = 0; m < T; m++)
	{
	  if(tupleFile[m] != f) continue;

	  for(int g = gFrom; g <= gToPair; g++)
	    {
	      Stats st;
	      comb[m].into(st,size_t(g),numLoci);
	      st.printStats(*cbOut[f],tupleLabel[m],g,param.tabbed());
	    }
	  if(!windows.empty())
	    {
	      writeWindowRunning(*cbWin[f],combWin[m],gStride,windows,lmap,
				 tupleLabelComma[m],gFrom,gToPair);
	    }
	  if(!param.tabbed()) *cbOut[f] << endl;
	}
    }

  for(int f = 0; f < F; f++)
    {
      delete cbOut[f];
      delete cbFull[f];
      delete cbWin[f];
    }

  return;
}
