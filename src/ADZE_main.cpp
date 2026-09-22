#include "ADZE_main_tools.h"
#include <cstdlib>
#include <cstdio>

using namespace std;

/*
 * Per-grouping ceiling on g: the fewest gene copies a grouping scored at any
 * surviving locus.  A locus where it scored fewer than g copies makes the
 * statistic undefined at that g, and Stats propagates that sentinel through
 * the average, so a single such locus takes the whole grouping out at that g.
 * The count of loci sitting at the ceiling therefore matters as much as the
 * ceiling itself: one locus is worth dropping, a third of the genome is not.
 *
 * ceiling[j] is that minimum, binding[j] how many loci sit on it, and empty[j]
 * how many loci the grouping did not score at all. The smallest ceiling is the
 * largest feasible MAX_G.
 */
static int smallestNj(Population pop[], int numDivs, int numLoci,
		      vector<int>* ceiling = 0, vector<int>* binding = 0,
		      vector<int>* empty = 0)
{
  int smallest = 0;

  if(ceiling) ceiling->assign(numDivs,0);
  if(binding) binding->assign(numDivs,0);
  if(empty) empty->assign(numDivs,0);

  for(int j = 0; j < numDivs; j++)
    {
      int least = 0;
      int atLeast = 0;
      int none = 0;

      for(int l = 0; l < numLoci; l++)
	{
	  const int Nj = pop[j].getNj(l);
	  if(l == 0 || Nj < least) { least = Nj; atLeast = 0; }
	  if(Nj == least) atLeast++;
	  if(Nj == 0) none++;
	}

      if(ceiling) (*ceiling)[j] = least;
      if(binding) (*binding)[j] = atLeast;
      if(empty) (*empty)[j] = none;

      if(j == 0 || least < smallest) smallest = least;
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
   * Development facility: walk the locus source and write the counts it
   * yields in the same layout the reader's dump uses (see dumpCounts), so
   * the two can be diffed directly. The sweep will read exactly these.
   */
  if(const char* sourcePath = getenv("ADZE_DUMP_SOURCE"))
    {
      try
	{
	  ParamSet sourceParams = p;
	  ScanResult scan;
	  scanDataset(sourceParams,scan,false);

	  /*
	   * STRUCTURE is converted first; VCF is read as it stands. Either way
	   * the engine sees the same interface, which is the point.
	   */
	  LocusSource* src = 0;
	  string tmp;
	  long long passes = 1;
	  if(wantsVCF(sourceParams.format.val,sourceParams.dfile.val))
	    {
	      src = openVCFSource(sourceParams,scan);
	    }
	  else
	    {
	      /*
	       * The lifecycle the sweep will use: a converted file that matches
	       * this run's data is read as it stands, anything else is converted
	       * again, and a conversion nobody asked to keep is removed.
	       */
	      const string kept = getenv("ADZE_COUNT_FILE")
		? string(getenv("ADZE_COUNT_FILE")) : string();
	      const string path = kept.empty() ? countFilePath(sourceParams) : kept;

	      string why;
	      if(countFileUsable(path,sourceParams,scan,why))
		{
		  adzelog() << "Reusing the converted counts in " << path << "\n";
		  passes = 0;
		}
	      else
		{
		  if(!kept.empty() && why != "it does not exist")
		    {
		      adzelog() << "Converting again: " << path << " cannot be used, "
				<< why << ".\n";
		    }
		  long long chunk = convertChunk(sourceParams,scan);
		  if(const char* forced = getenv("ADZE_CONVERT_CHUNK")) chunk = atoll(forced);
		  passes = transposeStructure(sourceParams,scan,path,chunk);
		}

	      if(kept.empty()) tmp = path;     //ours to remove
	      src = openCountFileSource(path,scan);
	    }
	  adzelog() << "converted in " << passes
		    << (passes == 1 ? " pass\n" : " passes\n");
	  ofstream d(sourcePath);
	  d << "LOCUS\tCHROM\tPOS\tGROUPING\tNJ\tMISSING\tSLOTS\tNJI\n";

	  LocusCounts locus;
	  const int J = int(src->groupNames().size());
	  while(src->next(locus))
	    {
	      for(int j = 0; j < J; j++)
		{
		  d << locus.name << "\t"
		    << (locus.chrom < 0 ? string(".") : scan.lmap.chromName[locus.chrom]) << "\t"
		    << locus.pos << "\t"
		    << src->groupNames()[j] << "\t"
		    << locus.nj[j] << "\t"
		    << (src->geneCopies(j) - locus.nj[j]) << "\t"
		    << locus.slots << "\t";

		  for(int i = 0; i < locus.slots; i++)
		    {
		      if(i) d << ",";
		      d << locus.count[size_t(i)*J + j];
		    }
		  d << "\n";
		}
	    }
	  d.close();
	  delete src;
	  if(!tmp.empty()) remove(tmp.c_str());
	}
      catch(BAD_FILE x) { return EXIT_IO; }
      catch(BAD_PARAM x) { return EXIT_DATA; }
    }

  /*
   * A dry run answers questions about the dataset -- its dimensions, each
   * grouping's sample size, the largest feasible MAX_G, how many windows and
   * rows the run would produce -- and every one of them follows from the
   * counts. So it is served by the scan, and the dataset is never read into
   * memory: asking what a 40 GB file would do no longer costs what doing it
   * would.
   */
  if(p.dry_run.val)
    {
      ScanResult scan;

      try
	{
	  scanDataset(p,scan);
	}
      catch(BAD_FILE x)
	{
	  return EXIT_IO;
	}
      catch(BAD_PARAM x)
	{
	  return EXIT_DATA;
	}

      const int numDivs = int(scan.groupName.size());

      vector< vector<int> > tuples;
      if(do_tuple)
	{
	  if(p.tuple_file.set)
	    {
	      if(!readTupleFile(p.tuple_file.val,scan.groupName,tuples)) return EXIT_USAGE;
	    }
	  else if(!validK(numDivs,k)) return EXIT_USAGE;
	}

      printDryRun(p,scan,tuples,k,do_rich,do_priv,do_tuple);
      return EXIT_OK;
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
  LocusTable locusTable;
  LocusMap lmap;

  /*
   * Development facility: run the scan alongside the reader so the two can be
   * compared where they should agree (see ADZE_DUMP_SCAN and
   * ADZE_DUMP_COUNTS). It costs an extra pass and is off unless the variable
   * is set. The scan resolves the same dimensions the reader does, so it runs
   * on a copy of the parameters rather than on the ones the reader will use.
   */
  if(getenv("ADZE_DUMP_SCAN"))
    {
      try
	{
	  ParamSet scanParams = p;
	  ScanResult scan;
	  scanDataset(scanParams,scan,false);
	}
      catch(BAD_FILE x) { return EXIT_IO; }
      catch(BAD_PARAM x) { return EXIT_DATA; }
    }

  try
    {
      pop = readDataset(p,divisionNames,numDivs,locusTable,lmap);
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

  //Development facility: see dumpCounts. Written before the filter runs, so
  //it is what the reader produced rather than what survived.
  if(const char* dumpPath = getenv("ADZE_DUMP_COUNTS"))
    {
      dumpCounts(dumpPath,pop,numDivs,p.loci.val,lmap);
    }

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
      filterLoci(pop,numDivs,p.tol.val,p.p_out.val,p.pp.val,lmap,locusTable);
    }

  /*
   * Windows are laid out over the loci that survived, so this has to follow
   * filtering.  A window is a run of consecutive loci, which only means
   * anything if the loci arrive in genome order.
   */
  vector<Window> windows;

  if(p.windowed())
    {
      vector<string> names(p.loci.val);
      for(int l = 0; l < p.loci.val; l++) names[l] = pop[0].getLocusName(l);

      const string bad = lmap.checkOrder(names);
      if(!bad.empty())
	{
	  cerr << "ERROR: " << bad << "\n";
	  delete [] pop;
	  return EXIT_DATA;
	}

      buildWindows(lmap,p,windows);

      if(windows.empty())
	{
	  cerr << "ERROR: no window holds at least --min-window-loci "
	       << p.min_win_loci.val << " loci.\n";
	  delete [] pop;
	  return EXIT_DATA;
	}

      adzelog() << "Laid out " << windows.size() << " window"
		<< (windows.size() == 1 ? "" : "s") << " over "
		<< p.loci.val << " loci.\n";

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
  vector<int> ceiling, binding, emptyIn;
  feasibleG = smallestNj(pop,numDivs,p.loci.val,&ceiling,&binding,&emptyIn);

  if(!p.g.set)
    {
      p.g.val = (feasibleG < 1) ? 1 : feasibleG;
      adzelog() << "MAX_G not given; using " << p.g.val << ", the largest the "
		<< p.loci.val << " surviving loci support.\n";
      summary << "MAX_G resolved to " << p.g.val << " over " << p.loci.val
	      << (p.loci.val == 1 ? " locus\n" : " loci\n");
    }

  /*
   * Which groupings cannot reach MAX_G, and what holds them back.  Naming them
   * is the whole point: a run whose ceiling is 1 produces nothing usable, and
   * in the legacy format it produces nothing visible either, since undefined
   * rows are omitted rather than written as NA. Reporting only the global
   * minimum -- which is what 1.0 did, and what this did until now -- leaves a
   * user staring at blank result files with no idea which grouping, or how
   * many loci, caused it.
   */
  for(int j = 0; j < numDivs; j++)
    {
      if(ceiling[j] >= p.g.val) continue;

      adzelog() << "WARNING: " << pop[j].getName() << " scored only "
		<< ceiling[j] << " gene cop" << ((ceiling[j] == 1) ? "y" : "ies")
		<< " at " << binding[j]
		<< ((binding[j] == 1) ? " locus" : " loci")
		<< " of " << p.loci.val;
      if(emptyIn[j] > 0)
	{
	  adzelog() << " (" << emptyIn[j]
		    << ((emptyIn[j] == 1) ? " locus is" : " loci are")
		    << " not scored at all)";
	}
      if(ceiling[j] == 0)
	adzelog() << ",\n         so its statistics are undefined at every g.\n";
      else
	adzelog() << ",\n         so its statistics are undefined above g = "
		  << ceiling[j] << ".\n";

      if(ceiling[j] < 2)
	{
	  adzelog() << "         That leaves nothing usable for "
		    << pop[j].getName() << ": lower --tolerance to drop those "
		    << "loci\n         (--dry-run reports the ceiling without "
		    << "computing anything).\n";
	  if(!p.tabbed())
	    {
	      adzelog() << "         --legacy omits undefined rows entirely, "
			<< "so they will be missing from the\n"
			<< "         result file rather than marked NA.\n";
	    }
	}
    }

  /*
   * --at-g reports one g instead of the ladder.  Two ceilings bind: MAX_G,
   * which is what the run was asked to sweep, and the smallest number of gene
   * copies scored anywhere, which is the largest g with defined values.  A
   * request above either is met at the ceiling with a warning -- unlike an
   * over-large MAX_G, which keeps its shape and reports the unreachable rows
   * as undefined, because a single-g report of nothing but NA would be no
   * answer at all.
   */
  if(p.at_g.set)
    {
      const bool wantsMax = (p.at_g.val == "max");
      int target = wantsMax ? p.g.val : p.at_g_val;

      int ceiling = p.g.val;
      string bound = "MAX_G";
      if(feasibleG < ceiling)
	{
	  ceiling = feasibleG;
	  bound = "the smallest sample size in the data";
	}
      if(ceiling < 1) ceiling = 1;

      if(target > ceiling)
	{
	  adzelog() << "WARNING: --at-g " << p.at_g.val << " exceeds " << bound
		    << " (" << ceiling << "); reporting at g = " << ceiling
		    << ".\n";
	  target = ceiling;
	}

      p.at_g_val = target;

      adzelog() << "Reporting at g = " << p.at_g_val << " only.\n";
      summary << "AT_G resolved to " << p.at_g_val << endl;
    }

  if(do_rich)
    {
      adzelog() << "Calculating allelic richness...\n";
      calcAllAgs(pop,numDivs,p,p.full_r.val,p.r_out.val,windows,lmap);
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
      calcAllPgs(pop,numDivs,p,p.full_p.val,p.p_out.val,windows,lmap);
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
	  calcPgTuples(pop,numDivs,tuples,p,p.full_c.val,p.c_out.val,1,windows,lmap);
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
	      calcPgTuples(pop,numDivs,tuples,p,p.full_c.val,new_comb_out,0,windows,lmap);

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
