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


  ifstream data;

  data.open(p.dfile.val.c_str());

  if(data.fail())
    {
      cout << "ERROR: Could not open " << p.dfile.val 
	   << "\nProgram terminated.\n";
      return -1;
    }
  else
    {
      cout << p.dfile.val << " sucessfully opened.\n";
    }
     
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


  int begin;
    
  cout << "Reading and sorting...\n";
  if(!p.skip_chk.val)
    {
      try
	{
	  checkDatafile(p);
	}
      catch(BAD_PARAM x)
	{
	  cout << "Program terminated.\n";
	  return -1;
	}
    }

  //allocate an array to hold lociNames
  string *lociNames = new string[p.loci.val];
  
  getLociNames(data,lociNames,p);
  
  //Get rid of non_data_rows
  string junk;
  for(int skip = 1; skip < p.nd_rows.val; skip++)
    {
      getline(data,junk);
      //cout << junk << endl;
    }


 
  //note the beginning of the data
  //begin = data.tellg();
  
  //Vectors to store the name of each division
  //and the number of data lines per division
  vector<string> divisionNames;
  vector<int> lines;
  
  //fill up the above vectors
  try
    {
      getDivLines(data,divisionNames,lines,p);
    }
  catch(BAD_PARAM x)
    {
      delete [] lociNames;
      cout << "Program terminated.\n";
      return -1;
    }
  
  int numDivs = divisionNames.size();
  
  if(p.comb.val)
    if(!validK(numDivs,k))
      {
	delete [] lociNames;
	return -1;
      }
  
  //Allocate and initialize population objects
  Population *pop = NULL;
  pop = new Population[numDivs];

  for (int i = 0; i < numDivs; i++)
    {
      pop[i].setRowsLoci(lines[i],p.loci.val);
      pop[i].setName(divisionNames[i]);
      
      for (int j = 0; j < p.loci.val; j++)
	{
	  pop[i].setLocusName(lociNames[j],j);
	}
    }
  
  delete [] lociNames;	
  
  //data.seekg(begin); //go back to beginning of data
  data.close();
  readData(/*data,*/pop,divisionNames,numDivs,p);
  //data.close();
 
  cout << "Done\n"; //done reading file

  if(p.tol.val != 1)
    {
      cout << "Throwing out bad loci...\n";
      filterLoci(pop,numDivs,p.tol.val,p.p_out.val,p.miss.val,p.pp.val);
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
  

  //CALCULATE Nji's
  cout << "Calculating Nji's...\n";
  //cout.flush();
  calcNji(pop,numDivs,p.miss.val);
  cout << "Completed at (d:h:m:s) ";
  displayTime(cout);
  cout << endl;

  //CALCULATE Nj's
  cout << "Calculating Nj's...\n";
  //cout.flush();
  calcNj(pop,numDivs);
  cout << "Completed at (d:h:m:s) ";
  displayTime(cout);
  cout << endl;
  
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

  cout << "\nADZE finished in (d:h:m:s) ";
  summary << "\nADZE finished in (d:h:m:s) ";
  displayTime(cout);
  displayTime(summary);
  cout << endl;
   
  return 0;
}

