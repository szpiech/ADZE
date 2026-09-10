#include "ADZE_main_tools.h"

using namespace std;

int main(int argc, char* argv[])
{
  cout << "Allelic Diversity Analyzer v1.0\n";

  ParamSet p;
  
  if(argc == 1)
    {
      try
	{
	  p.open("paramfile.txt");
	}
      catch(BAD_FILE x)
	{
	  p.makeParamFile(); 
	  cout << "\tA template has been created in paramfile.txt\n";
	  cout << "Program terminated.\n";
	  return -1;
	}
    }
  else if(argc == 2)
    {
      if(argv[1][0] == '-')
	{
	  cout << "ERROR: A paramfile must be the first argument and "
	       << "must not begin with '-'.\n"
	       << "Program terminated.\n";
	  return -1;
	}

      try
	{
	  p.open(argv[1]);
	}
      catch(BAD_FILE x)
	{
	  cout << "Program terminated.\n";
	  return -1;
	}
    }
  else
    {
      if(argv[1][0] == '-')
	{
	  cout << "ERROR: A paramfile must be the first argument and "
	       << "must not begin with '-'.\n"
	       << "Program terminated.\n";
	  return -1;
	}

      try
	{
	  p.open(argv[1]);
	}
      catch(BAD_FILE x)
	{
	  cout << "Program terminated.\n";
	  return -1;
	}
      
      try
	{
	  p.CMDread(argc,argv);
	}
      catch (BAD_PARAM x)
	{
	  cout << "Program terminated.\n";
	  return -1;
	}
    }

  try
    {
      p.read();
    }
  catch(BAD_FILE x)
    {
      cout << "Program terminated.\n";
      return -1;
    }
  catch (BAD_PARAM x)
    {
      cout << "Program terminated.\n";
      return -1;
    }
  
  p.close();

  cout << "Parameters read at (d:h:m:s) ";
  displayTime(cout);
  cout << endl;

  /*-------------------------------------------------------------------------*/
  ofstream summary;
  string sum_out = nameCreate(p.r_out.val,"_summary");
  summary.open(sum_out.c_str());
  p.echo(summary);
  summary << endl;


  list<int> k;
  try
    {
      if(p.comb.val) k = parseKVals(p.k.val);
    }
  catch (string err)
    {
      cout << err << "Program terminated.\n";
      return -1;
    }

  cout << "Reading " << p.dfile.val << "...\n";

  /*
   * One pass over the data file yields the groupings, the allele counts, the
   * sample sizes and the missing-data tallies.  1.0 read the file three times
   * (validate, discover groupings, load) and kept every genotype in memory.
   */
  vector<string> divisionNames;
  int numDivs = 0;
  Population* pop = NULL;

  try
    {
      pop = readDataset(p,divisionNames,numDivs);
    }
  catch(BAD_FILE x)
    {
      cout << "Program terminated.\n";
      return -1;
    }
  catch(BAD_PARAM x)
    {
      cout << "Program terminated.\n";
      return -1;
    }

  cout << "Done\n";
  cout << "Completed at (d:h:m:s) ";
  displayTime(cout);
  cout << endl;

  if(p.comb.val)
    if(!validK(numDivs,k))
      {
	delete [] pop;
	return -1;
      }

  if(p.tol.val != 1)
    {
      cout << "Throwing out bad loci...\n";
      filterLoci(pop,numDivs,p.tol.val,p.p_out.val,p.pp.val);
      //cout << "Done\n";
      
      cout << "Completed at (d:h:m:s) ";
      displayTime(cout);
      cout << endl;

      if(p.tnc.val)
	{
	  if(pop) delete [] pop;
	  return 0;
	} 
    }
  p.loci.val = pop[0].getNumLoci();

  /*
   * With no surviving locus every statistic is undefined.  1.0 carried on and
   * printed rows of "nan -0 nan"; say so and stop instead.
   */
  if(p.loci.val == 0)
    {
      cout << "ERROR: no locus survived filtering at TOLERANCE " << p.tol.val
	   << ".\n       Every statistic would be undefined; "
	   << "raise TOLERANCE or check MISSING.\n"
	   << "Program terminated.\n";
      summary << "ERROR: no locus survived filtering at TOLERANCE "
	      << p.tol.val << ".\n";
      summary.close();
      if(pop) delete [] pop;
      return -1;
    }

  int numLoci = p.loci.val;
  

  cout << "Calculating total alleles...\n";
  calcAllAgs(pop,numDivs,p,p.full_r.val,p.r_out.val);
  if(p.pp.val) cout << endl;
  cout << "Completed at (d:h:m:s) ";
  displayTime(cout);
  cout << endl;

  summary << "Total alleles completed at (d:h:m:s) ";
  displayTime(summary);
  summary << endl;
  
  cout << "Calculating private alleles...\n";
  calcAllPgs(pop,numDivs,p,p.full_p.val,p.p_out.val);
  if(p.pp.val) cout << endl;
  cout << "Completed at (d:h:m:s) ";
  displayTime(cout);
  cout << endl;

  summary << "Private alleles completed at (d:h:m:s) ";
  displayTime(summary);
  summary << endl;


  /*----------------------------------------------------------------
  
  
  ----------------------------------------------------------------*/


  if (p.comb.val)
    {
      for(list<int>::iterator i = k.begin();
	  i != k.end(); i++)
	{
	  cout << "Calculating private alleles for all possible " 
	       << *i << "-tuples...\n";
	  char suffix[50];
	  int junk;
	  junk = sprintf(suffix,"_%i",*i);
	  string new_comb_out;
	  new_comb_out = nameCreate(p.c_out.val,suffix);
	  calcAllPgComb(pop,numDivs,*i,p,p.full_c.val,new_comb_out);
	  if(p.pp.val) cout << endl;
	  cout << "Completed at (d:h:m:s) ";
	  displayTime(cout);
	  cout << endl;
	  /*
	  estimate(eout,numDivs,*i,pop,p,time-off);
	  off = time;
	  */
	  
	  summary << "All " << *i << "-tuples completed at (d:h:m:s) ";
	  displayTime(summary);
	  summary << endl;
	}
    }

 
  summary.close();
  /*
  cout << "Cleaning up...\n";
  //cout.flush();
  if(pop) delete [] pop;
  cout << "Completed at (d:h:m:s) ";
  displayTime(cout);
  */

  if(pop) delete [] pop;

  cout << "\nADZE finished in (d:h:m:s) ";
  summary << "\nADZE finished in (d:h:m:s) ";
  displayTime(cout);
  displayTime(summary);
  cout << endl;
   
  return 0;
}

