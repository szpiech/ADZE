#include <cctype>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

enum{FALSE,TRUE};

/*
 * Every parameter is described exactly once, in the table in ADZE_pfile.cpp:
 * its paramfile keyword, its long option, the 1.0 short flag it replaces, its
 * type, and its help text.  The parser, the --help text and the paramfile
 * template are all generated from that one table, so they cannot drift apart
 * the way 1.0's hand-written template drifted from its argv scanner.
 */
enum{G,LOCI,ND_ROWS,ND_COLS,DLINES,SORT_BY,TOL,K,DFILE,R_OUT,P_OUT,C_OUT,
     MISS,COMB,FULL_R,FULL_P,FULL_C,PP,TNC,SKIP_CHK,
     OUT_PREFIX,STAT,POPS,EXPOPS,TUPLE_FILE,TSV,DRY_RUN,QUIET,THREADS,PARAMS,
     FORMAT,SAMPLES,
     LABEL_SIZE};

enum OptType{OPT_INT,OPT_DOUBLE,OPT_BOOL,OPT_STRING};

struct OptSpec
{
  int id;
  const char* key;     //paramfile keyword, 0 if command line only
  const char* lng;     //long option
  const char* legacy;  //1.0 short flag, 0 if none
  OptType type;
  const char* arg;     //argument placeholder, 0 for a switch
  const char* section;
  const char* help;
};

extern const char* ADZE_VERSION;
extern const OptSpec OPTIONS[];
extern const int NUM_OPTIONS;

//Exit statuses, so a wrapper script can tell the failures apart.
const int EXIT_OK = 0;
const int EXIT_USAGE = 2;   //bad command line or paramfile
const int EXIT_IO = 3;      //could not read or write a file
const int EXIT_DATA = 4;    //the data file disagrees with itself or the params

class BAD_FILE{}; //exception class
class BAD_PARAM{}; 

class ParamSet
{
  template<class T> class Param
    {
      public:
      bool cl;   //set on the command line (overrides the paramfile)
      bool set;  //given at all (as opposed to defaulted)
      T val;
      Param() : cl(0), set(0) {}
    };

  bool LABEL_SEEN[LABEL_SIZE];
  bool LABEL_CL[LABEL_SIZE];   //given on the command line

  bool isint(string);
  bool isdouble(string);
  bool isbool(string);
  bool valid();
  void storeVal(int id,const string& val,bool cmd);
  bool isvalidk(string);
  static string trim(const string& s);
  const OptSpec* byKey(const string& key) const;
  const OptSpec* byFlag(const string& flag) const;

 public:
  
  Param<int> g;
  Param<int> loci;
  Param<int> nd_rows;
  Param<int> nd_cols;
  Param<int> dlines;
  Param<int> sort_by;
  Param<double> tol;
  Param<string> k;
  Param<string> dfile;
  Param<string> r_out;
  Param<string> p_out;
  Param<string> c_out;
  Param<string> miss;
  Param<bool> comb;
  Param<bool> full_r;
  Param<bool> full_p;
  Param<bool> full_c;
  Param<bool> pp;
  Param<bool> tnc;
  Param<bool> skip_chk;
  Param<string> out_prefix;
  Param<string> stat;
  Param<string> pops;
  Param<string> expops;
  Param<string> tuple_file;
  Param<string> format;   //auto | structure | vcf
  Param<string> samples;  //sample -> grouping map, for VCF input
  Param<bool> tsv;
  Param<bool> dry_run;
  Param<bool> quiet;
  Param<int> threads;
  Param<string> params;

  ParamSet();
  void read(const string& file);      //read a paramfile
  void echo(ostream& out);
  void CMDread(int,char**);           //parse the command line
  bool finish();                      //apply defaults, validate
  void makeParamFile(const string& file);
  string summaryName() const;
  static void usage(ostream& out);
  static void version(ostream& out);
};
