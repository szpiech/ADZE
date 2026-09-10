#include "ADZE_main_tools.h"

using namespace std;

//Smallest number of gene copies scored anywhere: the largest feasible MAX_G.
static int smallestNj(Population pop[], int numDivs, int numLoci)
{
  int smallest = 0;
  for(int j = 0; j < numDivs; j++)
    {
      for(int l = 0; l < numLoci; l++)
	{
	  const int Nj = pop[j].getNj(l);
	  if((j == 0 && l == 0) || Nj < smallest) smallest = Nj;
	}
    }
  return smallest;
}

/*
 * Interface, in brief:
 *
 *   adze --data FILE [options]        run from flags alone
 *   adze PARAMFILE [options]          1.0-style; flags override the file
 *   adze --help                       every option, generated from one table
 *
 * 1.0 required a paramfile as argv[1], rejected anything starting with '-'
 * there, offered no --help (it printed "A paramfile must be the first
 * argument" and exited 255), and answered a bare invocation by writing
 * paramfile.txt into the working directory. Flags could only override a
 * paramfile, never stand on their own.
 */
int main(int argc, char* argv[])
{
  /*
   * Answered before anything else is parsed, so they work with no data file,
   * no paramfile and no valid parameters.
   */
  for(int i = 1; i < argc; i++)
    {
      const string a = argv[i];

      if(a == "--help" || a == "-h")
	{
	  ParamSet::usage(cout);
	  return EXIT_OK;
	}

      if(a == "--version" || a == "-V")
	{
	  ParamSet::version(cout);
	  return EXIT_OK;
	}

      if(a == "--write-template")
	{
	  const string file = (i+1 < argc && argv[i+1][0] != '-')
	    ? argv[i+1] : "paramfile.txt";
	  ParamSet p;
	  try
	    {
	      p.makeParamFile(file);
	    }
	  catch(BAD_FILE x)
	    {
	      return EXIT_IO;
	    }
	  cout << "Wrote a parameter-file template to " << file << "\n";
	  return EXIT_OK;
	}
    }

  if(argc == 1)
    {
      ParamSet::usage(cerr);
      return EXIT_USAGE;
    }

  ParamSet p;

  try
    {
      p.CMDread(argc,argv);
      if(p.params.set) p.read(p.params.val);
    }
  catch(BAD_FILE x)
    {
      return EXIT_IO;
    }
  catch(BAD_PARAM x)
    {
      return EXIT_USAGE;
    }

  ADZE_QUIET = p.quiet.val;

#ifdef _OPENMP
  omp_set_num_threads(p.threads.val);
#else
  if(p.threads.set && p.threads.val > 1)
    {
      cerr << "WARNING: this build has no OpenMP support; --threads "
	   << p.threads.val << " ignored.\n";
      p.threads.val = 1;
    }
#endif

  //Progress bars are for a terminal; a redirected log gets none by default.
  if(!p.pp.set) p.pp.val = stderrIsTerminal() && !p.quiet.val;
  if(p.quiet.val) p.pp.val = 0;

  if(!p.finish()) return EXIT_USAGE;

  adzelog() << "Allelic Diversity Analyzer v" << ADZE_VERSION << "\n";

  //Which statistics to run: everything named, or 1.0's defaults.
  bool do_rich = 1, do_priv = 1, do_tuple = p.comb.val;

  if(p.stat.set)
    {
      do_rich = (p.stat.val.find("rich") != string::npos);
      do_priv = (p.stat.val.find("priv") != string::npos);
      do_tuple = (p.stat.val.find("tuple") != string::npos);

      if(!do_rich && !do_priv && !do_tuple)
	{
	  cerr << "ERROR: --stat \"" << p.stat.val
	       << "\" names no statistic; use richness, private and/or tuples.\n";
	  return EXIT_USAGE;
	}
    }

  list<int> k;
  try
    {
      if(do_tuple && !p.tuple_file.set) k = parseKVals(p.k.val);
    }
  catch (string err)
    {
      cerr << err;
      return EXIT_USAGE;
    }

  adzelog() << "Reading " << p.dfile.val << "...\n";

  /*
   * One pass over the data file yields the groupings, the allele counts, the
   * sample sizes and the missing-data tallies. 1.0 read the file three times
   * (validate, discover groupings, load) and kept every genotype in memory.
   */
  vector<string> divisionNames;
  int numDivs = 0;
  Population* pop = NULL;
  LocusMap lmap;

  try
    {
      pop = readDataset(p,divisionNames,numDivs,lmap);
    }
  catch(BAD_FILE x)
    {
      return EXIT_IO;
    }
  catch(BAD_PARAM x)
    {
      return EXIT_DATA;
    }

  adzelog() << "Read " << p.dlines.val << " gene copies at " << p.loci.val
	    << (p.loci.val == 1 ? " locus in " : " loci in ") << numDivs
	    << (numDivs == 1 ? " grouping" : " groupings") << " (d:h:m:s ";
  displayTime(adzelog());
  adzelog() << ")\n";

  int feasibleG = smallestNj(pop,numDivs,p.loci.val);

  vector< vector<int> > tuples;

  if(do_tuple)
    {
      if(p.tuple_file.set)
	{
	  if(!readTupleFile(p.tuple_file.val,pop,numDivs,tuples))
	    {
	      delete [] pop;
	      return EXIT_USAGE;
	    }
	}
      else if(!validK(numDivs,k))
	{
	  delete [] pop;
	  return EXIT_USAGE;
	}
    }

  if(p.dry_run.val)
    {
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
	  for(int l = 0; l < p.loci.val; l++)
	    {
	      const int Nj = pop[j].getNj(l);
	      if(l == 0 || Nj < minNj) minNj = Nj;
	      if(l == 0 || Nj > maxNj) maxNj = Nj;
	    }
	  cout << "  " << pop[j].getName() << ": " << pop[j].getNumRows()
	       << " gene copies, Nj per locus " << minNj << "-" << maxNj << "\n";
	}

      cout << "Largest feasible MAX_G: " << feasibleG << "\n";
      if(p.g.set) cout << "MAX_G to be used:       " << p.g.val << "\n";
      else cout << "MAX_G to be used:       " << "the largest feasible\n";
      if(p.tol.val != 1)
	{
	  cout << "  (TOLERANCE " << p.tol.val << " will drop loci first, "
	       << "which can raise both numbers)\n";
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

      delete [] pop;
      return EXIT_OK;
    }

  ofstream summary;
  const string sum_out = p.summaryName();
  summary.open(sum_out.c_str());
  if(summary.fail())
    {
      cerr << "ERROR: could not write " << sum_out << "\n";
      delete [] pop;
      return EXIT_IO;
    }
  p.echo(summary);
  summary << endl;

  if(p.tol.val != 1)
    {
      adzelog() << "Applying the missing-data filter...\n";
      filterLoci(pop,numDivs,p.tol.val,p.p_out.val,p.pp.val,lmap);

      adzelog() << "Completed at (d:h:m:s) ";
      displayTime(adzelog());
      adzelog() << endl;

      if(p.tnc.val)
	{
	  summary.close();
	  delete [] pop;
	  return EXIT_OK;
	} 
    }

  p.loci.val = pop[0].getNumLoci();

  /*
   * With no surviving locus every statistic is undefined. 1.0 carried on and
   * printed rows of "nan -0 nan"; say so and stop instead.
   */
  if(p.loci.val == 0)
    {
      cerr << "ERROR: no locus survived filtering at TOLERANCE " << p.tol.val
	   << ".\n       Every statistic would be undefined; "
	   << "raise TOLERANCE or check MISSING.\n";
      summary << "ERROR: no locus survived filtering at TOLERANCE "
	      << p.tol.val << ".\n";
      summary.close();
      delete [] pop;
      return EXIT_DATA;
    }

  /*
   * MAX_G is resolved here, not before filtering: dropping loci with heavy
   * missing data raises the smallest sample size, so the largest usable g is a
   * property of the surviving loci.
   */
  feasibleG = smallestNj(pop,numDivs,p.loci.val);

  if(!p.g.set)
    {
      p.g.val = (feasibleG < 2) ? 2 : feasibleG;
      adzelog() << "MAX_G not given; using " << p.g.val << ", the largest the "
		<< p.loci.val << " surviving loci support.\n";
      summary << "MAX_G resolved to " << p.g.val << " over " << p.loci.val
	      << (p.loci.val == 1 ? " locus\n" : " loci\n");
    }
  else if(p.g.val > feasibleG)
    {
      adzelog() << "WARNING: MAX_G " << p.g.val << " exceeds the smallest "
		<< "sample size in the data (" << feasibleG << ").\n"
		<< "         Rows above g = " << feasibleG
		<< " will be undefined.\n";
    }

  if(do_rich)
    {
      adzelog() << "Calculating allelic richness...\n";
      calcAllAgs(pop,numDivs,p,p.full_r.val,p.r_out.val);
      if(p.pp.val) adzelog() << endl;
      adzelog() << "Completed at (d:h:m:s) ";
      displayTime(adzelog());
      adzelog() << endl;

      summary << "Total alleles completed at (d:h:m:s) ";
      displayTime(summary);
      summary << endl;
    }

  if(do_priv)
    {
      adzelog() << "Calculating private allelic richness...\n";
      calcAllPgs(pop,numDivs,p,p.full_p.val,p.p_out.val);
      if(p.pp.val) adzelog() << endl;
      adzelog() << "Completed at (d:h:m:s) ";
      displayTime(adzelog());
      adzelog() << endl;

      summary << "Private alleles completed at (d:h:m:s) ";
      displayTime(summary);
      summary << endl;
    }

  if(do_tuple)
    {
      if(p.tuple_file.set)
	{
	  adzelog() << "Calculating private alleles for " << tuples.size()
		    << " named tuples...\n";
	  calcPgTuples(pop,numDivs,tuples,p,p.full_c.val,p.c_out.val,1);
	  if(p.pp.val) adzelog() << endl;
	  adzelog() << "Completed at (d:h:m:s) ";
	  displayTime(adzelog());
	  adzelog() << endl;

	  summary << "Named tuples completed at (d:h:m:s) ";
	  displayTime(summary);
	  summary << endl;
	}
      else
	{
	  for(list<int>::iterator i = k.begin(); i != k.end(); i++)
	    {
	      adzelog() << "Calculating private alleles for all possible "
			<< *i << "-tuples...\n";

	      ostringstream suffix;
	      suffix << "_" << *i;
	      const string new_comb_out = nameCreate(p.c_out.val,suffix.str());

	      buildKTuples(numDivs,*i,tuples);
	      calcPgTuples(pop,numDivs,tuples,p,p.full_c.val,new_comb_out,0);

	      if(p.pp.val) adzelog() << endl;
	      adzelog() << "Completed at (d:h:m:s) ";
	      displayTime(adzelog());
	      adzelog() << endl;

	      summary << *i << "-tuples completed at (d:h:m:s) ";
	      displayTime(summary);
	      summary << endl;
	    }
	}
    }

  delete [] pop;

  adzelog() << "\nadze finished in (d:h:m:s) ";
  const double total = displayTime(adzelog());
  adzelog() << endl;

  summary << "\nadze finished in (d:h:m:s) ";
  displayTime(summary);
  summary << endl;
  summary.close();

  (void)total;

  return EXIT_OK;
}
