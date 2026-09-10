#include "ADZE_stats.h"

using namespace std;

void Stats::printData(ostream& out,string name,int g)
{
  for(int l = 0; l < numLoci; l++)
    {
      if(data[l] == -9) return;
    }

  out << name << " "
      << g << " "
      << numLoci << " ";
  
  for(int l = 0; l < numLoci;l++)
    {
      out << data[l] << " ";
      
    }
  
  out << avg << " "
      << var << " "
      << std_err << endl;
	
  return;
}


void Stats::calcStdErr()
{
  if(var_flag && !std_err_flag)
    {
      std_err = sqrt(var/double(numLoci));
    }

  std_err_flag = 1;
  return;
}

void Stats::calcVar()
{
  if(avg_flag && !var_flag)
    {
      for(int l = 0; l < numLoci; l++)
	{
	  if(data[l] == -9)
	    {
	      var = -9;
	      l = numLoci; //break out of locus
	    }
	  else
	    {
	      var += (data[l] - avg)*(data[l] - avg);
	    }
	}
      
      if(var != -9)
	{
	  /*
	   * With one locus there is no across-locus variance to estimate.
	   * 1.0 divided by numLoci-1 regardless, so the reported value was
	   * whatever 0.0/0.0 happened to produce; make it an explicit NaN.
	   */
	  if(numLoci < 2) var = numeric_limits<double>::quiet_NaN();
	  else var /= double(numLoci-1);
	}
    }
  
  var_flag = 1;
  return;
}

void Stats::calcAvg()
{
  if(!avg_flag && data != NULL)
    {
      for(int l = 0; l < numLoci; l++)
	{
	  if(data[l] == -9)
	    {
	      avg = -9;
	      l = numLoci; //break out of locus
	    }
	  else
	    {
	      avg += data[l];
	    }
	}
      
      if(avg != -9)
	{
	  avg /= double(numLoci);
	}
    }
  
  avg_flag = 1;
  return;
}

Stats::Stats()
{
  avg_flag = 0;
  var_flag = 0;
  std_err_flag = 0;
  avg = 0;
  var = 0;
  numLoci = -9;
  data = NULL;
}


void Stats::putData(double* d,int l)
{
  numLoci = l;
  data = d;
  return;
}

/*
 * One result row.  In legacy format a row whose value is undefined at this g
 * is omitted entirely, so a consumer cannot tell "not computable here" from
 * "not run"; with tsv, the row is written with NA in the value columns and the
 * fields are tab-separated under a header.
 */
void Stats::printStats(ostream& out,string name,int g,bool tsv)
{
  const bool undefined = (avg == -9 || var == -9 || std_err == -9);

  if(undefined && !tsv) return;

  const char sep = tsv ? '\t' : ' ';

  out << name << sep << g << sep << numLoci << sep;

  if(undefined) out << "NA" << sep << "NA" << sep << "NA" << endl;
  else out << avg << sep << var << sep << std_err << endl;

  return;
}

Stats::~Stats()
{
  data = NULL;
}
