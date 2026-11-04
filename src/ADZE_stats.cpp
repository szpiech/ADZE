/*
 * ADZE: Allelic Diversity Analyzer
 * Copyright (C) 2008-2026 Zachary A Szpiech
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "ADZE_stats.h"

using namespace std;



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
 * "not run"; when tabbed, the row is written with NA in the value columns and the
 * fields are tab-separated under a header.
 */
void Stats::printStats(ostream& out,const string& name,int g,bool tabbed)
{
  const bool undefined = (avg == -9 || var == -9 || std_err == -9);

  if(undefined && !tabbed) return;

  const char sep = tabbed ? '\t' : ' ';

  out << name << sep << g << sep << numLoci << sep;

  //'\n' rather than endl: endl flushes the stream, and a window file is
  //millions of these rows.
  if(undefined) out << "NA" << sep << "NA" << sep << "NA" << '\n';
  else out << avg << sep << var << sep << std_err << '\n';

  return;
}

Stats::~Stats()
{
  data = NULL;
}
