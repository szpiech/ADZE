#include "ADZE_pfile.h"

using namespace std;

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

bool ParamSet::a_label(string s)
{
  for(int i = 0; i < LABEL_SIZE; i ++)
    {
      if(s.compare(LABEL[i]) == 0) return 1;
    }

  return 0;
}


void ParamSet::close()
{
  try
    {
      pin.close();
    }
  catch(...)
    {
      BAD_FILE x;
      throw x;
    }
  return;
}

void ParamSet::makeParamFile()
{
  ofstream pout;
  pout.open("paramfile.txt");

  pout << "# This is a parameter file for ADZE.\n"
       << "# The entire line after a \'#\' will be ignored.\n\n"
       << "###----Main Parameters----###\n"
       << "MAX_G                #Max standardized sample size\n\n"
       << "DATA_LINES           #Number of lines of data\n\n"
       << "LOCI                 #Number of loci\n\n"
       << "NON_DATA_ROWS        #Number of rows preceeding data\n"
       << "                     #including at least the locus names\n\n"
       << "NON_DATA_COLS        #Number of classifier columns\n"
       << "                     #at the beginning of each data line\n\n"
       << "GROUP_BY_COL         #The column number by which to\n"
       << "                     #group the data\n\n"
       << "DATA_FILE none       #Name of the datafile\n\n"
       << "R_OUT none           #Name of the allelic richness output file\n\n"
       << "P_OUT none           #Name of the private allelic richness\n"
       << "                     #output file\n\n"
       << "###----Combination Parameters----###\n"
       << "COMB 0               #Calculate private allelic richness for\n"
       << "                     #combinations of groupings?\n\n"
       << "K_RANGE              #A listing of combinations to calculate\n\n"
       << "C_OUT none           #Name of the private allelic richness of\n"
       << "                     #combinations output file\n\n"
       << "###----Advanced Options----###\n"
       << "MISSING -9           #Missing data representation\n\n"
       << "TOLERANCE 1          #Filter loci with a grouping having more\n"
       << "                     #than this fraction of missing data\n\n"
       << "FULL_R 0             #Output allelic richness results for\n"
       << "                     #all loci?\n\n"
       << "FULL_P 0             #Output private allelic richness results\n"
       << "                     #for all loci?\n\n"
       << "FULL_C 0             #Output private allelic richness for\n"
       << "                     #combinations results for all loci?\n\n"
       << "PRINT_PROGRESS 1     #Track calculation progress on screen?\n\n"
       << "###----Command line arguments----###\n"
       << "# -g MAX_G" << endl
       << "# -d DATA_LINES" << endl
       << "# -l LOCI" << endl
       << "# -nr NON_DATA_ROWS" << endl
       << "# -nc NON_DATA_COLS" << endl
       << "# -s GROUP_BY_COL" << endl
       << "# -f DATA_FILE" << endl
       << "# -r R_OUT" << endl
       << "# -p P_OUT" << endl
       << "# -c COMB" << endl
       << "# -k K_RANGE" << endl
       << "# -o COUT" << endl
       << "# -m MISSING" << endl
       << "# -t TOLERANCE" << endl
       << "# -tnocalc" << endl
       << "# -fr FULL_R" << endl
       << "# -fp FULL_P" << endl
       << "# -fc FULL_C" << endl
       << "# -pp PRINT_PROGRESS" << endl
       << "###----End of file----###\n";

  pout.close();
  return;
}

void ParamSet::open(char* s)
{
  pin.open(s);
  if(pin.fail())
    {
      cout << "ERROR: Could not find paramfile \"" << s << "\"\n";
      BAD_FILE x;
      throw x;
    }
  return;
}

void ParamSet::CMDread(int argc, char* argv[])
{
  //cout << argc << endl;
  bool CMD_READ[LABEL_SIZE];
  for(int i = 0; i < LABEL_SIZE;i++)
    {
      CMD_READ[i] = FALSE;
    }
  
  int index;
  
  for(int i = 2; i <= argc-1; i+=2)
    {
      if(validCMD(argv[i]) && i+1 <= argc-1)
	{
	  try
	    {
	      index = strtoen(argv[i]);
	    }
	  catch (BAD_PARAM x)
	    {
	      cout << "Undefined ERROR 1 in ParamSet::CMDread().\n";
	      throw x;
	    }
	  
	  if(validCMD(argv[i+1]))
	    {		  
	      cout << "ERROR: " << argv[i] << " needs a value.\n";
	      BAD_PARAM x;
	      throw x;
	    }
	  else if(CMD_READ[index])
	    {
	      string line = argv[i];
	      line += ' ';
	      line += argv[i+1];
	      cout << "ERROR: Duplicate parameter definition:\n\t\"" 
		   << line << "\"\n ";
	      BAD_PARAM x;
	      throw x;
	    }
	  else
	    {
	      try
		{
		  storeVal(argv[i+1],argv[i],TRUE);
		  CMD_READ[index] = TRUE;
		}
	      catch(BAD_PARAM x)
		{
		  //cout << "Program terminated.\n";
		  throw x;
		}
	    }
	}
      else if(validCMD(argv[i]) && i+1 > argc-1)
	{
	   try
	    {
	      index = strtoen(argv[i]);
	    }
	  catch (BAD_PARAM x)
	    {
	      cout << "Undefined ERROR 2 in ParamSet::CMDread().\n";
	      throw x;
	    }
	  if(CMD_READ[index])
	    {
	      cout << "ERROR: Duplicate parameter definition:\n\t\"" 
		   << argv[i] << "\"\n ";
	    }
	  else
	    {
	      cout << "ERROR: " << argv[i] << " needs a value.\n";
	    }
	  BAD_PARAM x;
	  throw x;
	}
      else
	{
	  cout << "ERROR: " << argv[i] << " not a valid command line "
	       << "argument.\n";
	  BAD_PARAM x;
	  throw x;
	}
    }
  return;
}

bool ParamSet::validCMD(string cmd)
{
  for(int i = 0; i < LABEL_SIZE;i++)
    {
      if(cmd.compare(CMD_LABEL[i]) == 0) return 1;
    }

  return 0;
}

int ParamSet::strtoen(string s)
{
  for(int i = 0; i < LABEL_SIZE; i++)
    {
      if(s.compare(CMD_LABEL[i]) == 0) return i;
    }
  
  BAD_PARAM x;
  throw x;
}

bool ParamSet::isgoodstr(string s)
{
  char c;
  string::iterator i;
  for(i = s.begin(); i != s.end(); i++)
    {
      c = *i;
      if(isgraph(c)) return 1;
    }

  return 0;
 }

bool ParamSet::valid()
{
  bool passAll = 1;
  
  if(g.val < 2)
    {
      cout << "ERROR: Maximum g value must be an integer > 1.\n";
      passAll = 0;
    }
  if(loci.val < 1)
    {
      cout << "ERROR: Number of loci must be an integer > 0.\n";
      passAll = 0;
    }
  if(nd_rows.val < 1)
    {
      cout << "ERROR: Number of non-data rows must be an integer > 0.\n";
      passAll = 0;
    }
  if(nd_cols.val < 1)
    {
      cout << "ERROR: Number of non-data cols must be an integer > 0.\n";
      passAll = 0;
    }
  if(dlines.val < 1)
    {
      cout << "ERROR: Number of data lines must be an integer > 0.\n";
      passAll = 0;
    }
  if(sort_by.val < 1 || sort_by.val > nd_cols.val)
    {
      cout << "ERROR: Must sort by a positive integer less than or\n"
	   << "equal to the number of non-data columns.\n";
      passAll = 0;
    }
  if(tol.val < 0 || tol.val > 1)
    {
      cout << "ERROR: Tolerance must be between 0 and 1, inclusive.\n";
      passAll = 0;
    }
  if(r_out.val.compare("none") == 0)
    {
      cout << "ERROR: Must specify an allelic richness output file.\n";
      passAll = 0;
    }
  if(p_out.val.compare("none") == 0)
    {
      cout << "ERROR: Must specify a private richness output file.\n";
      passAll = 0;
    }
  if(dfile.val.compare("none") == 0)
    {
      cout << "ERROR: Must specify a data file.\n";
      passAll = 0;
    }
  if(comb.val)
    {
      if(c_out.val.compare("none") == 0)
	{
	  cout << "ERROR: Must specify a combination output file.\n";
	  passAll = 0;
	}
      if(k.val.compare("none") == 0)
	{
	  cout << "ERROR: Must specify a range of k values to " 
	       << "calculate combinations.\n";
	  passAll = 0;
	}
      else if(!isvalidk(k.val))
	{
	  cout << "Error: \"" << k.val 
	       << "\" not a valid K_RANGE definition.\n";
	  passAll = 0;
	}
      
    }
  return passAll;
}

void ParamSet::echo(ostream& out)
{
  out << "###----Main Parameters----###\n"
      << LABEL[G] << " " << g.val << endl
      << LABEL[DLINES] << " " << dlines.val << endl
      << LABEL[LOCI] << " " << loci.val << endl
      << LABEL[ND_ROWS] << " " << nd_rows.val << endl
      << LABEL[ND_COLS] << " " << nd_cols.val << endl
      << LABEL[SORT_BY] << " " << sort_by.val << endl
      << LABEL[DFILE] << " " << dfile.val << endl
      << LABEL[R_OUT] << " " << r_out.val << endl
      << LABEL[P_OUT] << " " << p_out.val << endl
      << "\n###----Combination Parameters----###\n"
      << LABEL[COMB] << " " << comb.val << endl
      << LABEL[K] << " " << k.val << endl
      << LABEL[C_OUT] << " " << c_out.val << endl
      << "\n###-----------Advanced Options-----------###\n"
      << LABEL[MISS] << " " << miss.val << endl
      << LABEL[TOL] << " " << tol.val << endl
      << LABEL[FULL_R] << " " << full_r.val << endl
      << LABEL[FULL_P] << " " << full_p.val << endl
      << LABEL[FULL_C] << " " << full_c.val << endl
      << LABEL[PP] << " " << pp.val;
}

bool ParamSet::isbool(string s)
{
  //  if(s.size() == 0) return 0;
  string::iterator i;
  char c;
  for(i = s.begin();i != s.end(); i++)
    {
      c = *i;
      if(!(c == '0' || c == '1')) return 0;
    }
  
  int num = atoi(s.c_str());
  
  if(num < 0 || num > 1) return 0;
  
  return 1;
}

bool ParamSet::isdouble(string s)
{
  //  if(s.size() == 0) return 0;
  string::iterator i;
  char c;
  for(i = s.begin();i != s.end();i++)
    {
      c = *i;
      if(!(isdigit(c) || c == '.')) return 0;
    }
  return 1;
}

string ParamSet::replaceWhite(string s)
{
  string::iterator i;
  char c;
  for(i = s.begin();i != s.end();i++)
    {
      c = *i;
      if(isspace(c)) *i = '#';
    }
  return s;
}

bool ParamSet::isint(string s)
{
  char c;
  string::iterator i;
  for(i = s.begin();i != s.end();i++)
    {
      c = *i;
      
      if(i == s.begin() && !(isdigit(c) || c == '-')) return 0;
      if(!isdigit(c)) return 0;
    }

  return 1;
}

bool ParamSet::allWhiteSpace(string s)
{
  string::iterator i;
  char c;
  for(i = s.begin();i != s.end();i++)
    {
      c = *i;
      if(!isspace(c)) return 0;
    }

  return 1;
}

void ParamSet::storeVal(string val, string label, bool cmd)
{
  string L[LABEL_SIZE];
  bool err = 0;

  if(cmd)
    {
      for(int i = 0; i < LABEL_SIZE; i++)
	{
	  L[i] = CMD_LABEL[i];
	}
    }
  else
    {
      for(int i = 0; i < LABEL_SIZE; i++)
	{
	  L[i] = LABEL[i];
	}
    }

  
  if(label.compare(L[G]) == 0 ||
     label.compare(L[LOCI]) == 0 ||
     label.compare(L[ND_ROWS]) == 0 ||
     label.compare(L[ND_COLS]) == 0 ||
     label.compare(L[DLINES]) == 0 ||
     label.compare(L[SORT_BY]) == 0)
    {
      if(label.compare(L[G]) == 0 && !g.cl) 
	{
	  if(isint(val))
	    {
	      int i = atoi(val.c_str());
	      g.cl = cmd;
	      g.val = i;
	    }
	  else err = 1;
	}
      if(label.compare(L[LOCI]) == 0 && !loci.cl)
	{
	  if(isint(val))
	    {
	      int i = atoi(val.c_str());
	      loci.cl = cmd;
	      loci.val = i;
	    }
	  else err = 1;
	}
      if(label.compare(L[ND_ROWS]) == 0 && !nd_rows.cl)
	{
	  if(isint(val))
	    {
	      int i = atoi(val.c_str());
	      nd_rows.cl = cmd;
	      nd_rows.val = i;
	    }
	   else err = 1;
	}
      if(label.compare(L[ND_COLS]) == 0 && !nd_cols.cl)
	{
	  if(isint(val))
	    {
	      int i = atoi(val.c_str());
	      nd_cols.cl = cmd;
	      nd_cols.val = i;
	    }
	   else err = 1;
	}
      if(label.compare(L[DLINES]) == 0 && !dlines.cl) 
	{
	  if(isint(val))
	    {
	      int i = atoi(val.c_str());
	      dlines.cl = cmd;
	      dlines.val = i;
	    }
	   else err = 1;
	}
      if(label.compare(L[SORT_BY]) == 0 && !sort_by.cl) 
	{
	  if(isint(val))
	    {
	      int i = atoi(val.c_str());
	      sort_by.cl = cmd;
	      sort_by.val = i;
	    }
	   else err = 1;
	}
    
      if(err) 
	{
	  cout << "ERROR: " << label << " must specify a positive integer.\n";
	  BAD_PARAM x;
	  throw x;
	}
    }
  else if(label.compare(L[TOL]) == 0)
    {
      if(!tol.cl && val.size() != 0)
	{
	  if(isdouble(val))
	    {
	      tol.cl = cmd;
	      tol.val = atof(val.c_str());
	    }
	  else err = 1;
	}
	
      if(err)
	{
	  cout << "ERROR: " << label << " must specify a double.\n";
	  BAD_PARAM x;
	  throw x;
	}
    }
  else if(label.compare(L[K]) == 0 ||
	  label.compare(L[DFILE]) == 0 ||
	  label.compare(L[R_OUT]) == 0 ||
	  label.compare(L[P_OUT]) == 0 ||
	  label.compare(L[C_OUT]) == 0 ||
	  label.compare(L[MISS]) == 0)
    {     
      if(isgoodstr(val))
	{
	  if(label.compare(L[K]) == 0 && !k.cl /*&& (isvalidk(val) || 
						    val.compare("none") == 0)*/)
	    {
	      k.cl = cmd;
	      k.val = val;
	    }
	  /*
	  else if(label.compare(L[K]) == 0 && !k.cl && !isvalidk(val))
	    {
	      cout << "Error: \"" << val 
		   << "\" not a valid K_RANGE definition.\n";
	      BAD_PARAM x;
	      throw x;
	    }
	  */
	  if(label.compare(L[DFILE]) == 0 && !dfile.cl)
	    {
	      dfile.cl = cmd;
	      dfile.val = val;
	    }
	  if(label.compare(L[R_OUT]) == 0 && !r_out.cl)
	    {
	      r_out.cl = cmd;
	      r_out.val = val;
	    }
	  if(label.compare(L[P_OUT]) == 0 && !p_out.cl)
	    {
	      p_out.cl = cmd;
	      p_out.val = val;
	    }
	  if(label.compare(L[C_OUT]) == 0 && !c_out.cl)
	    {
	      c_out.cl = cmd;
	      c_out.val = val;
	    }
	  if(label.compare(L[MISS]) == 0 && !miss.cl)
	    {
	      miss.cl = cmd;
	      miss.val = val;
	    }
	}
    }
  else if(label.compare(L[COMB]) == 0 ||
	  label.compare(L[FULL_R]) == 0 ||
	  label.compare(L[FULL_P]) == 0 ||
	  label.compare(L[FULL_C]) == 0 ||
	  label.compare(L[TNC]) == 0 ||
	  label.compare(L[PP]) == 0 ||
	  label.compare(L[SKIP_CHK]) == 0)
    {
      if(label.compare(L[COMB]) == 0 && !comb.cl)
	{
	  if(isbool(val))
	    {
	      bool b = atoi(val.c_str());
	      comb.cl = cmd;
	      comb.val = b;
	    }
	  else err = 1;
	}
      if(label.compare(L[FULL_R]) == 0 && !full_r.cl)
	{
	  if(isbool(val))
	    {
	      bool b = atoi(val.c_str());
	      full_r.cl = cmd;
	      full_r.val = b;
	    }
	  else err = 1;
	}
      if(label.compare(L[FULL_P]) == 0 && !full_p.cl)
	{
	  if(isbool(val))
	    {
	      bool b = atoi(val.c_str());
	      full_p.cl = cmd;
	      full_p.val = b;
	    }
	  else err = 1;
	}
      if(label.compare(L[FULL_C]) == 0 && !full_c.cl)
	{
	  if(isbool(val))
	    {
	      bool b = atoi(val.c_str());
	      full_c.cl = cmd;
	      full_c.val = b;
	    }
	  else err = 1;
	}
      if(label.compare(L[TNC]) == 0 && !tnc.cl)
	{
	  if(isbool(val))
	    {
	      bool b = atoi(val.c_str());
	      tnc.cl = cmd;
	      tnc.val = b;
	    }
	  else err = 1;
	}
      if(label.compare(L[PP]) == 0 && !pp.cl)
	{
	  if(isbool(val))
	    {
	      bool b = atoi(val.c_str());
	      pp.cl = cmd;
	      pp.val = b;
	    }
	  else err = 1;
	}
      if(label.compare(L[SKIP_CHK]) == 0 && !skip_chk.cl)
	{
	  if(isbool(val))
	    {
	      bool b = atoi(val.c_str());
	      skip_chk.cl = cmd;
	      skip_chk.val = b;
	    }
	  else err = 1;
	}

      if(err)
	{
	  cout << "ERROR: " << label << " must equal 1 or 0.\n";
	  BAD_PARAM x;
	  throw x;
	}
    }
  else
    {
      cout << "Undefined ERROR in ParamSet::storeVal().\n";
      BAD_PARAM x;
      throw x;
    }
 
  return;
}


string ParamSet::eraseWhite(string s)
{
  char c;
  string::iterator i;
  i = s.begin();
  c = *i;
  while(i != s.end())
    {
      if(!isgraph(c)) i = s.erase(i);
      else i++;
      c = *i;
    }
  return s;
}

void ParamSet::read()
{
  string line1, line, candidate;
  int found;
  
  if(!pin)
    {
      cout << "ERROR: ifstream not open.\n";
      BAD_FILE x;
      throw x;
    }

  while(!pin.eof())
    {
      getline(pin,line1);
   
      found = line1.find_first_of('#');
      if(found != string::npos) line = line1.erase(found);
      
      line = eraseWhite(line1);

      for(int i=0;i<LABEL_SIZE;i++)
	{
	  string find_this = LABEL[i];
	  //find_this += ' ';
	  
	  found = line.find(find_this);
	  if(found != string::npos)
	    {
	      //line = eraseLeadWhite(line);
	      candidate = line.substr(0,LABEL[i].size());

	      if(!a_label(candidate))
		{
		  cout << "Line not recognized:\n\t\"" << line1 << "\"\n";
		  BAD_PARAM x;
		  throw x;
		}

	      line = line.erase(0,LABEL[i].size());
	      
	      if(SKIP_LABEL[i])
		{
		  cout << "Duplicate parameter definition:\n\t\"" 
		       << line1 << "\"\n";
		  BAD_PARAM x;
		  throw x;
		}
	      else
		{
		  SKIP_LABEL[i] = TRUE;
		  try
		    {
		      storeVal(line,/*LABEL[i]*/candidate);
		    }
		  catch(BAD_PARAM x)
		    {
		      throw x;
		    }
		  break;
		}
	    }
	  else if(i == LABEL_SIZE-1 && !allWhiteSpace(line) && !a_label(line))
	    {
	      cout << "Line not recognized:\n\t\"" << line1 << "\"\n";
	      BAD_PARAM x;
	      throw x;
	    }
	}      
    }

  if(!valid())
    {
      BAD_PARAM x;
      throw x;
    }

  return;
}




ParamSet::ParamSet()
{
  for(int i = 0; i < LABEL_SIZE; i++)
    {
      SKIP_LABEL[i] = FALSE;
      LABEL_SEEN[i] = FALSE;
    }
  
  g.cl = FALSE;
  loci.cl = FALSE;
  nd_rows.cl = FALSE;
  nd_cols.cl = FALSE;
  dlines.cl = FALSE;
  sort_by.cl = FALSE;
  tol.cl = FALSE;
  k.cl = FALSE;
  dfile.cl = FALSE;
  r_out.cl = FALSE;
  p_out.cl = FALSE;
  c_out.cl = FALSE;
  miss.cl = FALSE;
  comb.cl = FALSE;
  full_r.cl = FALSE;
  full_p.cl = FALSE;
  full_c.cl = FALSE;
  tnc.cl = FALSE;
  pp.cl = FALSE;
  skip_chk.cl = FALSE;

  g.val = -9;
  loci.val = -9;
  nd_rows.val = -9;
  nd_cols.val = -9;
  dlines.val = -9;
  sort_by.val = -9;
  tol.val = 1;
  k.val = "none";
  dfile.val = "none";
  r_out.val = "none";
  p_out.val = "none";
  c_out.val = "none";
  miss.val = "-9";
  comb.val = FALSE;
  full_r.val = FALSE;
  full_p.val = FALSE;
  full_c.val = FALSE;
  tnc.val = FALSE;
  pp.val = TRUE;
  skip_chk.val = FALSE;
}
