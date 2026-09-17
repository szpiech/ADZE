#include "ADZE_pfile.h"
#include <sstream>
#include <cstring>

//Defined in ADZE_main_tools.cpp; declared here to avoid an include cycle.
bool wantsVCF(const string& format,const string& path);

using namespace std;

const char* ADZE_VERSION = "2.0-dev";

/*
 * The single description of every parameter.  The command-line parser, the
 * paramfile parser, --help and the generated template all read this table, so
 * a flag cannot exist in one and be missing (or take a different number of
 * arguments) in another -- which is what happened in 1.0, whose template
 * documented -tnocalc and -skipchk as valueless while the argv scanner walked
 * strict pairs and rejected them without a value.
 *
 * legacy holds the 1.0 flag spelling, still accepted so existing scripts and
 * paramfiles keep working.
 */
const OptSpec OPTIONS[] = {
  {DFILE, "DATA_FILE", "--data", "-f", OPT_STRING, "FILE", 0,
   "Input", "genotype file: STRUCTURE layout or VCF, optionally gzipped (required)"},
  {FORMAT, "FORMAT", "--format", 0, OPT_STRING, "FMT", "auto",
   "Input", "input format: auto, structure or vcf; auto reads the file extension"},
  {SAMPLES, "SAMPLE_FILE", "--samples", 0, OPT_STRING, "FILE", 0,
   "Input", "sample-to-grouping map, required for VCF input (two columns: sample grouping)"},
  {LOCI_MAP, "LOCI_MAP", "--loci-map", 0, OPT_STRING, "FILE", 0,
   "Input", "locus coordinates for STRUCTURE input (three columns: locus chromosome position)"},
  {G, "MAX_G", "--max-g", "-g", OPT_INT, "N", "the largest the data supports",
   "Input", "largest standardized sample size"},
  {ND_ROWS, "NON_DATA_ROWS", "--non-data-rows", "-nr", OPT_INT, "N", "1",
   "Input", "header rows before the genotypes"},
  {ND_COLS, "NON_DATA_COLS", "--non-data-cols", "-nc", OPT_INT, "N", "detected",
   "Input", "label columns before the genotypes"},
  {LOCI, "LOCI", "--loci", "-l", OPT_INT, "N", "detected",
   "Input", "number of loci"},
  {DLINES, "DATA_LINES", "--data-lines", "-d", OPT_INT, "N", "detected",
   "Input", "number of data rows"},
  {SORT_BY, "GROUP_BY_COL", "--group-col", "-s", OPT_INT, "N", "the last label column",
   "Input", "which label column names the grouping"},
  {MISS, "MISSING", "--missing", "-m", OPT_STRING, "STR", "-9",
   "Input", "code for a missing allele"},
  {POPS, 0, "--pops", 0, OPT_STRING, "LIST", "every grouping",
   "Input", "analyse only these groupings (comma-separated)"},
  {EXPOPS, 0, "--exclude-pops", 0, OPT_STRING, "LIST", "none excluded",
   "Input", "analyse everything except these groupings"},

  {WIN_BP, "WINDOW_BP", "--window-bp", 0, OPT_LONG, "N", "off",
   "Analysis", "sliding windows N basepairs wide (needs locus coordinates)"},
  {WIN_LOCI, "WINDOW_LOCI", "--window-loci", 0, OPT_LONG, "N", "off",
   "Analysis", "sliding windows N loci wide (surviving loci; needs coordinates for chromosomes)"},
  {STEP_BP, "STEP_BP", "--step-bp", 0, OPT_LONG, "N", "the window width",
   "Analysis", "advance basepair windows by N; a step equal to the width means no overlap"},
  {STEP_LOCI, "STEP_LOCI", "--step-loci", 0, OPT_LONG, "N", "the window width",
   "Analysis", "advance locus windows by N loci"},
  {MIN_WIN_LOCI, "MIN_WINDOW_LOCI", "--min-window-loci", 0, OPT_INT, "N", "1",
   "Analysis", "do not report a window holding fewer than N loci"},
  {AT_G, "AT_G", "--at-g", 0, OPT_STRING, "N|max", "every g",
   "Analysis", "report only this g, or max for whatever MAX_G resolves to"},
  {STAT, 0, "--stat", 0, OPT_STRING, "LIST", "richness,private",
   "Analysis", "which statistics to compute: richness,private,tuples"},
  {TOL, "TOLERANCE", "--tolerance", "-t", OPT_DOUBLE, "X", "0.1",
   "Analysis", "drop a locus if any grouping exceeds this missing fraction; 1 keeps every locus"},
  {COMB, "COMB", "--combinations", "-c", OPT_BOOL, 0, "off",
   "Analysis", "also compute private alleles of grouping tuples"},
  {K, "K_RANGE", "--tuples-k", "-k", OPT_STRING, "LIST", 0,
   "Analysis", "tuple sizes for --combinations, e.g. 2 or 1-3 or 1,3,5-7"},
  {TUPLE_FILE, 0, "--tuples", 0, OPT_STRING, "FILE", "every k-subset",
   "Analysis", "compute only the named tuples in FILE (one per line) instead of every k-subset"},
  {THREADS, 0, "--threads", 0, OPT_INT, "N", "1",
   "Analysis", "worker threads; requires an OpenMP build"},

  {OUT_PREFIX, 0, "--out-prefix", 0, OPT_STRING, "PREFIX", "adze",
   "Output", "compose output names from PREFIX"},
  {R_OUT, "R_OUT", "--out-richness", "-r", OPT_STRING, "FILE", "from --out-prefix",
   "Output", "allelic richness output file"},
  {P_OUT, "P_OUT", "--out-private", "-p", OPT_STRING, "FILE", "from --out-prefix",
   "Output", "private allelic richness output file"},
  {C_OUT, "C_OUT", "--out-tuples", "-o", OPT_STRING, "FILE", "from --out-prefix",
   "Output", "tuple output file (one per tuple size, suffixed _k)"},
  {FULL_R, "FULL_R", "--full-richness", "-fr", OPT_BOOL, 0, "off",
   "Output", "also write per-locus allelic richness"},
  {FULL_P, "FULL_P", "--full-private", "-fp", OPT_BOOL, 0, "off",
   "Output", "also write per-locus private allelic richness"},
  {FULL_C, "FULL_C", "--full-tuples", "-fc", OPT_BOOL, 0, "off",
   "Output", "also write per-locus tuple values"},
  {LEGACY, "LEGACY_FORMAT", "--legacy", 0, OPT_BOOL, 0, "off",
   "Output", "write version 1.0's layout: space-separated, no header, undefined rows omitted"},
  {TSV, 0, "--tsv", 0, OPT_BOOL, 0, "on",
   "Output", "tab-separated output with a header row and NA for undefined values"},

  {PARAMS, 0, "--params", 0, OPT_STRING, "FILE", 0,
   "Other", "read parameters from FILE (also accepted as the first argument)"},
  {DRY_RUN, 0, "--dry-run", 0, OPT_BOOL, 0, "off",
   "Other", "report the detected layout, groupings and feasible MAX_G, then stop"},
  {PP, "PRINT_PROGRESS", "--progress", "-pp", OPT_BOOL, 0, "on when stderr is a terminal",
   "Other", "progress bars"},
  {QUIET, 0, "--quiet", 0, OPT_BOOL, 0, "off",
   "Other", "suppress progress and informational messages"},
  {TNC, "TNC", "--tolerance-only", "-tnocalc", OPT_BOOL, 0, "off",
   "Other", "apply the missing-data filter, report it, and stop"},
  {SKIP_CHK, "SKIP_CHK", "--no-check", "-skipchk", OPT_BOOL, 0, "off",
   "Other", "accepted for compatibility; validation is now free and always on"}
};

const int NUM_OPTIONS = int(sizeof(OPTIONS)/sizeof(OPTIONS[0]));

/*----------------------------------------------------------------------------*/

ParamSet::ParamSet()
{
  for(int i = 0; i < LABEL_SIZE; i++) { LABEL_SEEN[i] = 0; LABEL_CL[i] = 0; }

  g.val = 0;              //0 => resolve from the data
  loci.val = 0;           //0 => detect
  nd_rows.val = 1;
  nd_cols.val = 0;        //0 => detect
  dlines.val = 0;         //0 => detect
  sort_by.val = 0;        //0 => last label column
  tol.val = 0.1;
  k.val = "none";
  dfile.val = "none";
  r_out.val = "none";
  p_out.val = "none";
  c_out.val = "none";
  miss.val = "-9";
  comb.val = 0;
  full_r.val = 0;
  full_p.val = 0;
  full_c.val = 0;
  pp.val = 0;
  tnc.val = 0;
  skip_chk.val = 0;
  out_prefix.val = "adze";
  stat.val = "";
  pops.val = "";
  expops.val = "";
  tuple_file.val = "";
  format.val = "auto";
  samples.val = "";
  loci_map.val = "";
  win_bp.val = 0;
  win_loci.val = 0;
  step_bp.val = 0;
  step_loci.val = 0;
  min_win_loci.val = 1;
  at_g.val = "";
  at_g_val = 0;
  tsv.val = 1;
  legacy.val = 0;
  dry_run.val = 0;
  quiet.val = 0;
  threads.val = 1;
  params.val = "";
}

const OptSpec* ParamSet::byKey(const string& key) const
{
  for(int i = 0; i < NUM_OPTIONS; i++)
    {
      if(OPTIONS[i].key && key.compare(OPTIONS[i].key) == 0) return &OPTIONS[i];
    }
  return 0;
}

const OptSpec* ParamSet::byFlag(const string& flag) const
{
  for(int i = 0; i < NUM_OPTIONS; i++)
    {
      if(OPTIONS[i].lng && flag.compare(OPTIONS[i].lng) == 0) return &OPTIONS[i];
      if(OPTIONS[i].legacy && flag.compare(OPTIONS[i].legacy) == 0) return &OPTIONS[i];
    }
  return 0;
}

string ParamSet::trim(const string& s)
{
  size_t a = 0, b = s.size();
  while(a < b && isspace((unsigned char)s[a])) a++;
  while(b > a && isspace((unsigned char)s[b-1])) b--;
  return s.substr(a,b-a);
}

bool ParamSet::isint(string s)
{
  s = trim(s);
  if(s.empty()) return 0;
  size_t i = (s[0] == '+' || s[0] == '-') ? 1 : 0;
  if(i >= s.size()) return 0;
  for(; i < s.size(); i++) if(!isdigit((unsigned char)s[i])) return 0;
  return 1;
}

bool ParamSet::isdouble(string s)
{
  s = trim(s);
  if(s.empty()) return 0;
  char* end = 0;
  strtod(s.c_str(),&end);
  return (end && *end == '\0');
}

bool ParamSet::isbool(string s)
{
  s = trim(s);
  return (s == "0" || s == "1" || s == "true" || s == "false" ||
	  s == "TRUE" || s == "FALSE" || s == "yes" || s == "no");
}

static bool boolValue(const string& s)
{
  return !(s == "0" || s == "false" || s == "FALSE" || s == "no");
}

/*
 * Store one value.  cmd distinguishes the command line, which wins over the
 * paramfile no matter which is read first.
 */
void ParamSet::storeVal(int id,const string& raw,bool cmd)
{
  const OptSpec* spec = 0;
  for(int i = 0; i < NUM_OPTIONS; i++) if(OPTIONS[i].id == id) spec = &OPTIONS[i];
  if(!spec) return;

  const string val = trim(raw);
  const string where = cmd ? string(spec->lng) : string(spec->key ? spec->key : spec->lng);

  //The command line wins over the paramfile, whichever is parsed first.
  if(!cmd && LABEL_CL[id]) return;

  bool bad = 0;
  switch(spec->type)
    {
    case OPT_INT:    bad = !isint(val);    break;
    case OPT_LONG:   bad = !isint(val);    break;
    case OPT_DOUBLE: bad = !isdouble(val); break;
    case OPT_BOOL:   bad = !isbool(val);   break;
    case OPT_STRING: bad = val.empty();    break;
    }

  if(bad)
    {
      cerr << "ERROR: " << where << " needs ";
      switch(spec->type)
	{
	case OPT_INT:    cerr << "an integer";             break;
	case OPT_LONG:   cerr << "a whole number";         break;
	case OPT_DOUBLE: cerr << "a number";               break;
	case OPT_BOOL:   cerr << "0 or 1";                 break;
	case OPT_STRING: cerr << "a value";                break;
	}
      cerr << ", got \"" << val << "\".\n";
      BAD_PARAM x;
      throw x;
    }

#define SETP(field)							\
  do {									\
    (field).cl = cmd;							\
    (field).set = 1;							\
  } while(0)

  switch(id)
    {
    case G:          SETP(g);          g.val = atoi(val.c_str());        break;
    case LOCI:       SETP(loci);       loci.val = atoi(val.c_str());     break;
    case ND_ROWS:    SETP(nd_rows);    nd_rows.val = atoi(val.c_str());  break;
    case ND_COLS:    SETP(nd_cols);    nd_cols.val = atoi(val.c_str());  break;
    case DLINES:     SETP(dlines);     dlines.val = atoi(val.c_str());   break;
    case SORT_BY:    SETP(sort_by);    sort_by.val = atoi(val.c_str());  break;
    case THREADS:    SETP(threads);    threads.val = atoi(val.c_str());  break;
    case TOL:        SETP(tol);        tol.val = atof(val.c_str());      break;
    case K:          SETP(k);          k.val = val;                      break;
    case DFILE:      SETP(dfile);      dfile.val = val;                  break;
    case R_OUT:      SETP(r_out);      r_out.val = val;                  break;
    case P_OUT:      SETP(p_out);      p_out.val = val;                  break;
    case C_OUT:      SETP(c_out);      c_out.val = val;                  break;
    case MISS:       SETP(miss);       miss.val = val;                   break;
    case OUT_PREFIX: SETP(out_prefix); out_prefix.val = val;             break;
    case STAT:       SETP(stat);       stat.val = val;                   break;
    case POPS:       SETP(pops);       pops.val = val;                   break;
    case EXPOPS:     SETP(expops);     expops.val = val;                 break;
    case TUPLE_FILE: SETP(tuple_file); tuple_file.val = val;             break;
    case FORMAT:     SETP(format);     format.val = val;                 break;
    case SAMPLES:    SETP(samples);    samples.val = val;                break;
    case LOCI_MAP:   SETP(loci_map);   loci_map.val = val;               break;
    case WIN_BP:     SETP(win_bp);     win_bp.val = atoll(val.c_str());   break;
    case WIN_LOCI:   SETP(win_loci);   win_loci.val = atoll(val.c_str()); break;
    case STEP_BP:    SETP(step_bp);    step_bp.val = atoll(val.c_str());  break;
    case STEP_LOCI:  SETP(step_loci);  step_loci.val = atoll(val.c_str());break;
    case MIN_WIN_LOCI: SETP(min_win_loci); min_win_loci.val = atoi(val.c_str()); break;
    case AT_G:       SETP(at_g);       at_g.val = val;                   break;
    case PARAMS:     SETP(params);     params.val = val;                 break;
    case COMB:       SETP(comb);       comb.val = boolValue(val);        break;
    case FULL_R:     SETP(full_r);     full_r.val = boolValue(val);      break;
    case FULL_P:     SETP(full_p);     full_p.val = boolValue(val);      break;
    case FULL_C:     SETP(full_c);     full_c.val = boolValue(val);      break;
    case PP:         SETP(pp);         pp.val = boolValue(val);          break;
    case TNC:        SETP(tnc);        tnc.val = boolValue(val);         break;
    case SKIP_CHK:   SETP(skip_chk);   skip_chk.val = boolValue(val);    break;
    case TSV:        SETP(tsv);        tsv.val = boolValue(val);         break;
    case LEGACY:     SETP(legacy);     legacy.val = boolValue(val);      break;
    case DRY_RUN:    SETP(dry_run);    dry_run.val = boolValue(val);     break;
    case QUIET:      SETP(quiet);      quiet.val = boolValue(val);       break;
    default: break;
    }

#undef SETP

  LABEL_SEEN[id] = 1;
  if(cmd) LABEL_CL[id] = 1;
  return;
}

/*
 * Parse the command line.  Long options, the 1.0 short flags, and switches
 * that take no value: --combinations and -c both work, and so does the 1.0
 * spelling '-c 1'.  A value is consumed for a switch only when it looks like
 * one, so '--progress --quiet' is not misread.
 */
void ParamSet::CMDread(int argc, char* argv[])
{
  int i = 1;

  //A bare first argument is a paramfile, as in 1.0.
  if(argc > 1 && argv[1][0] != '-')
    {
      params.val = argv[1];
      params.set = 1;
      params.cl = 1;
      i = 2;
    }

  for(; i < argc; i++)
    {
      string flag = argv[i];

      //--flag=value
      string inlineVal;
      bool hasInline = 0;
      size_t eq = flag.find('=');
      if(flag.size() > 2 && flag[0] == '-' && flag[1] == '-' && eq != string::npos)
	{
	  inlineVal = flag.substr(eq+1);
	  flag = flag.substr(0,eq);
	  hasInline = 1;
	}

      const OptSpec* spec = byFlag(flag);
      if(!spec)
	{
	  cerr << "ERROR: unrecognized option \"" << flag << "\".\n"
	       << "Try 'adze --help'.\n";
	  BAD_PARAM x;
	  throw x;
	}

      string val;
      if(hasInline) val = inlineVal;
      else if(spec->arg)                      //takes a value
	{
	  if(i+1 >= argc)
	    {
	      cerr << "ERROR: " << flag << " needs a value.\n";
	      BAD_PARAM x;
	      throw x;
	    }
	  val = argv[++i];
	}
      else                                    //switch
	{
	  if(i+1 < argc && isbool(argv[i+1])) val = argv[++i];
	  else val = "1";
	}

      storeVal(spec->id,val,1);
    }

  return;
}

/*
 * Read a paramfile.  Each line is a keyword followed by the rest of the line.
 *
 * 1.0 instead searched every line for any keyword as a substring, in enum
 * order, and then assumed the keyword it found was the line's leading label:
 * 'R_OUT out_LOCI_richness' matched LOCI first, sliced the first four
 * characters as the label, and rejected a perfectly valid line.  It also
 * deleted all whitespace from the line, so no value could contain a space.
 */
void ParamSet::read(const string& file)
{
  ifstream in(file.c_str());
  if(in.fail())
    {
      cerr << "ERROR: Could not find paramfile \"" << file << "\"\n";
      BAD_FILE x;
      throw x;
    }

  string line;
  long lineNo = 0;

  while(getline(in,line))
    {
      lineNo++;

      size_t hash = line.find('#');
      if(hash != string::npos) line = line.substr(0,hash);

      line = trim(line);
      if(line.empty()) continue;

      size_t sep = 0;
      while(sep < line.size() && !isspace((unsigned char)line[sep])) sep++;

      const string key = line.substr(0,sep);
      const string val = trim(line.substr(sep));

      const OptSpec* spec = byKey(key);
      if(!spec)
	{
	  cerr << "ERROR: " << file << ":" << lineNo
	       << ": unrecognized keyword \"" << key << "\".\n";
	  BAD_PARAM x;
	  throw x;
	}

      if(!spec->arg && val.empty()) storeVal(spec->id,"1",0);
      else storeVal(spec->id,val,0);
    }

  in.close();
  return;
}

/*
 * Apply defaults that do not depend on the data file and check the parameters
 * that can be checked before reading it.  Dimensions left unset are detected
 * while the file is read; MAX_G is resolved from the sample sizes afterwards.
 */
bool ParamSet::finish()
{
  bool ok = 1;

  if(!r_out.set) r_out.val = out_prefix.val + ".richness";
  if(!p_out.set) p_out.val = out_prefix.val + ".private";
  if(!c_out.set) c_out.val = out_prefix.val + ".tuples";

  if(tuple_file.set) comb.val = 1;

  if(stat.set)
    {
      //--stat tuples implies the tuple pass
      if(stat.val.find("tuple") != string::npos) comb.val = 1;
    }

  if(dfile.val.compare("none") == 0)
    {
      cerr << "ERROR: no data file given (--data FILE).\n";
      ok = 0;
    }

  if(g.set && g.val < 2)
    {
      cerr << "ERROR: --max-g must be an integer greater than 1.\n";
      ok = 0;
    }

  if(nd_rows.val < 1)
    {
      cerr << "ERROR: --non-data-rows must be a positive integer.\n";
      ok = 0;
    }

  if(loci.set && loci.val < 1)
    {
      cerr << "ERROR: --loci must be a positive integer.\n";
      ok = 0;
    }

  if(nd_cols.set && nd_cols.val < 1)
    {
      cerr << "ERROR: --non-data-cols must be a positive integer.\n";
      ok = 0;
    }

  if(dlines.set && dlines.val < 1)
    {
      cerr << "ERROR: --data-lines must be a positive integer.\n";
      ok = 0;
    }

  if(sort_by.set && sort_by.val < 1)
    {
      cerr << "ERROR: --group-col must be a positive integer.\n";
      ok = 0;
    }

  if(tol.val < 0 || tol.val > 1)
    {
      cerr << "ERROR: --tolerance must be between 0 and 1, inclusive.\n";
      ok = 0;
    }

  if(threads.val < 1)
    {
      cerr << "ERROR: --threads must be a positive integer.\n";
      ok = 0;
    }

  if(pops.set && expops.set)
    {
      cerr << "ERROR: --pops and --exclude-pops are mutually exclusive.\n";
      ok = 0;
    }

  /*
   * Windows are defined one way or the other, never both, and each step
   * belongs to its own unit.  A step smaller than the window is the sliding
   * case; the default step is the window itself, i.e. a plain tiling.
   */
  if(win_bp.set && win_loci.set)
    {
      cerr << "ERROR: --window-bp and --window-loci are mutually exclusive; "
	   << "a window is measured one way or the other.\n";
      ok = 0;
    }

  if((win_bp.set && win_bp.val < 1) || (win_loci.set && win_loci.val < 1))
    {
      cerr << "ERROR: the window width must be a positive number.\n";
      ok = 0;
    }

  if(step_bp.set && !win_bp.set)
    {
      cerr << "ERROR: --step-bp needs --window-bp.\n";
      ok = 0;
    }

  if(step_loci.set && !win_loci.set)
    {
      cerr << "ERROR: --step-loci needs --window-loci.\n";
      ok = 0;
    }

  if((step_bp.set && step_bp.val < 1) || (step_loci.set && step_loci.val < 1))
    {
      cerr << "ERROR: the window step must be a positive number.\n";
      ok = 0;
    }

  if(min_win_loci.val < 1)
    {
      cerr << "ERROR: --min-window-loci must be at least 1.\n";
      ok = 0;
    }

  if(!step_bp.set) step_bp.val = win_bp.val;
  if(!step_loci.set) step_loci.val = win_loci.val;

  if(min_win_loci.set && !windowed())
    {
      cerr << "WARNING: --min-window-loci applies to a windowed run; "
	   << "ignoring it.\n";
    }

  /*
   * Output layout: tab-separated with a header by default, 1.0's layout on
   * request.  Asking for both is a contradiction rather than a precedence
   * puzzle, so it is refused.
   */
  if(legacy.set && legacy.val)
    {
      if(tsv.set && tsv.val)
	{
	  cerr << "ERROR: --tsv and --legacy ask for different layouts; "
	       << "choose one.\n";
	  ok = 0;
	}
      else tsv.val = 0;
    }

  /*
   * --at-g takes a number or the word max.  A value below 1 is refused rather
   * than rounded up: g is a sample size, and 0 gene copies is not a request
   * the program can meet halfway.
   */
  if(at_g.set && at_g.val != "max")
    {
      if(!isint(at_g.val) || atoi(at_g.val.c_str()) < 1)
	{
	  cerr << "ERROR: --at-g takes a whole number of 1 or more, or the "
	       << "word max; got \"" << at_g.val << "\".\n";
	  ok = 0;
	}
      else at_g_val = atoi(at_g.val.c_str());
    }

  if(format.set && format.val != "auto" && format.val != "structure" &&
     format.val != "vcf")
    {
      cerr << "ERROR: --format must be auto, structure or vcf, not \""
	   << format.val << "\".\n";
      ok = 0;
    }
  else
    {
      /*
       * A VCF names its samples but says nothing about which population each
       * belongs to, so the map is required rather than optional.
       */
      const bool vcf = wantsVCF(format.val,dfile.val);

      if(vcf && !samples.set)
	{
	  cerr << "ERROR: VCF input needs --samples FILE, a two-column map "
	       << "from sample name to grouping.\n";
	  ok = 0;
	}

      /*
       * Windows are ranges of the genome, so every locus has to be placed on
       * it.  VCF input always is; the STRUCTURE layout needs the map.
       */
      if(windowed() && !vcf && !loci_map.set)
	{
	  cerr << "ERROR: a windowed run needs locus coordinates. STRUCTURE "
	       << "input carries none,\n       so supply --loci-map FILE "
	       << "(locus, chromosome, position).\n";
	  ok = 0;
	}

      if(!vcf && samples.set)
	{
	  cerr << "WARNING: --samples applies to VCF input only; ignoring it. "
	       << "Groupings come from GROUP_BY_COL.\n";
	}

      if(vcf && (nd_rows.set || nd_cols.set || sort_by.set))
	{
	  cerr << "WARNING: NON_DATA_ROWS, NON_DATA_COLS and GROUP_BY_COL "
	       << "describe the STRUCTURE layout and are ignored for VCF "
	       << "input.\n";
	}
    }

  if(comb.val && !tuple_file.set)
    {
      if(k.val.compare("none") == 0)
	{
	  cerr << "ERROR: --combinations needs --tuples-k (or --tuples FILE).\n";
	  ok = 0;
	}
      else if(!isvalidk(k.val))
	{
	  cerr << "ERROR: \"" << k.val << "\" not a valid K_RANGE definition.\n";
	  ok = 0;
	}
    }

  return ok;
}

/*
 * 1.0 appended suffixes to the end of the given name, so "-r results.txt"
 * produced results.txt, results.txt_summary and results.txt_fulldata. Explicit
 * names keep that behaviour for compatibility; names composed from
 * --out-prefix get the suffix before the extension.
 */
string ParamSet::summaryName() const
{
  if(r_out.set) return r_out.val + "_summary";
  return out_prefix.val + ".summary.txt";
}

void ParamSet::echo(ostream& out)
{
  out << "###----Main Parameters----###\n"
      << "MAX_G ";
  if(g.set) out << g.val << endl;
  else out << "auto" << endl;
  out
      << "DATA_LINES " << dlines.val << endl
      << "LOCI " << loci.val << endl
      << "NON_DATA_ROWS " << nd_rows.val << endl
      << "NON_DATA_COLS " << nd_cols.val << endl
      << "GROUP_BY_COL " << sort_by.val << endl
      << "DATA_FILE " << dfile.val << endl
      << "R_OUT " << r_out.val << endl
      << "P_OUT " << p_out.val << endl
      << "\n###----Combination Parameters----###\n"
      << "COMB " << comb.val << endl
      << "K_RANGE " << k.val << endl
      << "C_OUT " << c_out.val << endl;
  if(tuple_file.set) out << "TUPLE_FILE " << tuple_file.val << endl;
  if(format.set) out << "FORMAT " << format.val << endl;
  if(samples.set) out << "SAMPLE_FILE " << samples.val << endl;
  if(loci_map.set) out << "LOCI_MAP " << loci_map.val << endl;
  if(win_bp.set)
    {
      out << "WINDOW_BP " << win_bp.val << endl
	  << "STEP_BP " << step_bp.val << endl;
    }
  if(win_loci.set)
    {
      out << "WINDOW_LOCI " << win_loci.val << endl
	  << "STEP_LOCI " << step_loci.val << endl;
    }
  if(windowed()) out << "MIN_WINDOW_LOCI " << min_win_loci.val << endl;
  if(at_g.set) out << "AT_G " << at_g.val << endl;
  out << "\n###-----------Advanced Options-----------###\n"
      << "MISSING " << miss.val << endl
      << "TOLERANCE " << tol.val << endl
      << "FULL_R " << full_r.val << endl
      << "FULL_P " << full_p.val << endl
      << "FULL_C " << full_c.val << endl
      << "PRINT_PROGRESS " << pp.val;
  if(pops.set)   out << endl << "POPS " << pops.val;
  if(expops.set) out << endl << "EXCLUDE_POPS " << expops.val;
  if(threads.val != 1) out << endl << "THREADS " << threads.val;
}

void ParamSet::version(ostream& out)
{
  out << "adze " << ADZE_VERSION << "\n"
      << "Allelic Diversity Analyzer -- rarefaction estimators of allelic\n"
      << "richness and private allelic richness.\n"
      << "Method: Szpiech ZA, Jakobsson M, Rosenberg NA (2008) Bioinformatics\n"
      << "24:2498-2504. doi:10.1093/bioinformatics/btn478\n";
}

void ParamSet::usage(ostream& out)
{
  out << "adze " << ADZE_VERSION
      << " -- allelic richness and private allelic richness by rarefaction\n\n"
      << "Usage:\n"
      << "  adze --data FILE [options]\n"
      << "  adze PARAMFILE [options]        (1.0-style; options override the file)\n"
      << "  adze --write-template [FILE]    write a commented paramfile template\n\n";

  const char* sections[] = {"Input","Analysis","Output","Other"};

  for(int s = 0; s < 4; s++)
    {
      out << sections[s] << " options:\n";
      for(int i = 0; i < NUM_OPTIONS; i++)
	{
	  const OptSpec& o = OPTIONS[i];
	  if(strcmp(o.section,sections[s]) != 0) continue;

	  ostringstream left;
	  left << "  " << o.lng;
	  if(o.arg) left << " " << o.arg;
	  if(o.legacy) left << ", " << o.legacy;

	  string col = left.str();
	  out << col;
	  if(col.size() < 34) out << string(34-col.size(),' ');
	  else out << "\n" << string(34,' ');
	  out << o.help;

	  /*
	   * The default comes from the table rather than from the help text,
	   * so --help, the manual and the README cannot disagree about it.
	   */
	  if(o.dflt) out << "  [default: " << o.dflt << "]";
	  out << "\n";
	}
      out << "\n";
    }

  out << "Other:\n"
      << "  --help, -h                        this message\n"
      << "  --version, -V                     version and citation\n\n"
      << "Exit status: 0 success, " << EXIT_USAGE << " bad usage, "
      << EXIT_IO << " I/O error, " << EXIT_DATA << " data validation error.\n";
}

/*
 * Write a paramfile template generated from the option table, so it can no
 * longer disagree with what the parser accepts.
 */
void ParamSet::makeParamFile(const string& file)
{
  ofstream out(file.c_str());
  if(out.fail())
    {
      cerr << "ERROR: could not write " << file << "\n";
      BAD_FILE x;
      throw x;
    }

  out << "# adze " << ADZE_VERSION << " parameter file\n"
      << "#\n"
      << "# One KEYWORD per line, followed by its value; blank lines and\n"
      << "# text after # are ignored. Every keyword below has an equivalent\n"
      << "# command-line option, which overrides the value set here.\n"
      << "#\n";

  const char* sections[] = {"Input","Analysis","Output","Other"};

  for(int s = 0; s < 4; s++)
    {
      out << "\n###----" << sections[s] << "----###\n";
      for(int i = 0; i < NUM_OPTIONS; i++)
	{
	  const OptSpec& o = OPTIONS[i];
	  if(!o.key) continue;
	  if(strcmp(o.section,sections[s]) != 0) continue;

	  out << "# " << o.help << "\n"
	      << "#   command line: " << o.lng;
	  if(o.arg) out << " " << o.arg;
	  out << "\n";

	  /*
	   * Keywords with a usable default are written with it, so a template
	   * whose DATA_FILE has been filled in runs as it stands.  Keywords
	   * that are detected from the data, or that have no default, are
	   * written commented out with a placeholder.
	   */
	  switch(o.id)
	    {
	    case DFILE:    out << "DATA_FILE your_data.stru\n";     break;
	    case ND_ROWS:  out << "NON_DATA_ROWS 1\n";              break;
	    case MISS:     out << "MISSING -9\n";                   break;
	    case TOL:      out << "TOLERANCE 0.1\n";                break;
	    case COMB:     out << "COMB 0\n";                       break;
	    case R_OUT:    out << "#R_OUT richness.txt\n";          break;
	    case P_OUT:    out << "#P_OUT private.txt\n";           break;
	    case C_OUT:    out << "#C_OUT tuples.txt\n";            break;
	    case K:        out << "#K_RANGE 2\n";                   break;
	    case FULL_R:   out << "FULL_R 0\n";                     break;
	    case FULL_P:   out << "FULL_P 0\n";                     break;
	    case FULL_C:   out << "FULL_C 0\n";                     break;
	    case TNC:      out << "TNC 0\n";                        break;
	    case SKIP_CHK: out << "#SKIP_CHK 0\n";                  break;
	    default:
	      out << "#" << o.key << " " << (o.arg ? o.arg : "0") << "\n";
	      break;
	    }
	}
    }

  out.close();
  return;
}

/*
 * Syntax check for a K_RANGE value: comma- or space-separated tuple sizes and
 * inclusive ranges, e.g. "2", "1-3", "1,3,5-7".
 *
 * 1.0's version compared characters against the first and last character of
 * the string to decide whether it was safe to look ahead, then read *(i+1) --
 * which is one past the last character when the guard misfired.
 */
bool ParamSet::isvalidk(string s)
{
  size_t i = 0;
  const size_t n = s.size();
  bool any = false;

  while(i < n)
    {
      if(s[i] == '#') break;
      if(isspace((unsigned char)s[i]) || s[i] == ',') { i++; continue; }

      if(!isdigit((unsigned char)s[i])) return 0;
      while(i < n && isdigit((unsigned char)s[i])) i++;
      any = true;

      if(i < n && s[i] == '-')
	{
	  i++;
	  if(i >= n || !isdigit((unsigned char)s[i])) return 0;
	  while(i < n && isdigit((unsigned char)s[i])) i++;
	}

      if(i < n && !(isspace((unsigned char)s[i]) || s[i] == ',' || s[i] == '#'))
	{
	  return 0;
	}
    }

  return any;
}

