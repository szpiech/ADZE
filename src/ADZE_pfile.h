#include <cctype>
#include <cstdlib>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

enum{FALSE,TRUE};

enum{G,LOCI,ND_ROWS,ND_COLS,DLINES,SORT_BY,TOL,K,DFILE,R_OUT,P_OUT,C_OUT,
       MISS,COMB,FULL_R,FULL_P,FULL_C,PP,TNC,SKIP_CHK};

const int LABEL_SIZE = 20;
const string CMD_LABEL[] = {"-g","-l","-nr","-nc","-d","-s","-t","-k","-f",
			    "-r","-p","-o","-m","-c","-fr","-fp","-fc","-pp",
			    "-tnocalc","-skipchk"};
const string LABEL[] = {"MAX_G","LOCI","NON_DATA_ROWS","NON_DATA_COLS",
			"DATA_LINES","GROUP_BY_COL","TOLERANCE","K_RANGE",
			"DATA_FILE","R_OUT","P_OUT","C_OUT","MISSING","COMB",
			"FULL_R","FULL_P","FULL_C","PRINT_PROGRESS","TNC",
			"SKIP_CHK"};

class BAD_FILE{}; //exception class
class BAD_PARAM{}; 

class ParamSet
{
  template<class T> class Param
    {
      public:
      bool cl;
      T val;
    };

  ifstream pin;

  bool SKIP_LABEL[LABEL_SIZE];
  bool LABEL_SEEN[LABEL_SIZE];

  bool allWhiteSpace(string);
  bool isint(string);
  bool isdouble(string);
  bool isbool(string);
  string replaceWhite(string);
  bool valid();
  bool isgoodstr(string);
  void storeVal(string,string,bool cmd = 0);
  bool validCMD(string);
  int strtoen(string);
  bool a_label(string);
  bool isvalidk(string);
  string eraseWhite(string);
  //string eraseBackWhite(string);

 public:
  
  /*0*/  Param<int> g;
  /*1*/ Param<int> loci;
  /*2*/  Param<int> nd_rows;
  /*3*/ Param<int> nd_cols;
  /*4*/ Param<int> dlines;
  /*5*/ Param<int> sort_by;
  /*6*/ Param<double> tol;
  /*7*/ Param<string> k;
  /*8*/ Param<string> dfile;
  /*9*/ Param<string> r_out;
  /*10*/ Param<string> p_out;
  /*11*/ Param<string> c_out;
  /*12*/ Param<string> miss;
  /*13*/ Param<bool> comb;
  /*14*/ Param<bool> full_r;
  /*15*/ Param<bool> full_p;
  /*16*/ Param<bool> full_c;
  /*17*/ Param<bool> pp;
  /*18*/ Param<bool> tnc;
  /*19*/ Param<bool> skip_chk;
  
  ParamSet();
  void read();
  void echo(ostream& out);
  void CMDread(int,char**);
  void open(char*);
  void close();
  void makeParamFile();
};
