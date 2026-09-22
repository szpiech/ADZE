#include "ADZE_main_tools.h"
#include <cstdlib>
#include <cstdio>

using namespace std;


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
  /*
   * The run: a scan, then one sweep over the loci. Nothing between them
   * holds the dataset -- the scan keeps a count per locus per grouping, the
   * sweep keeps one locus at a time -- so peak memory no longer follows the
   * genotypes. See doc/streaming.md.
   */
  ScanResult scan;
  string countPath;

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

  numDivs = int(scan.groupName.size());

  adzelog() << "Read " << p.dlines.val << " gene copies at " << p.loci.val
	    << " loci in " << numDivs
	    << (numDivs == 1 ? " grouping" : " groupings") << " (d:h:m:s ";
  displayTime(adzelog());
  adzelog() << ")" << endl;

  vector< vector<int> > tuples;
  vector<int> tupleFile;
  vector<string> tupleOut;

  if(do_tuple)
    {
      if(p.tuple_file.set)
	{
	  if(!readTupleFile(p.tuple_file.val,scan.groupName,tuples)) return EXIT_USAGE;
	  tupleOut.push_back(p.c_out.val);
	  tupleFile.assign(tuples.size(),0);
	}
      else if(!validK(numDivs,k)) return EXIT_USAGE;
    }

  ofstream summary;
  const string sum_out = p.summaryName();
  summary.open(sum_out.c_str());
  if(summary.fail())
    {
      cerr << "ERROR: could not write " << sum_out << "\n";
      return EXIT_IO;
    }
  p.echo(summary);
  summary << endl;

  /*
   * The filter is decided here, from the counts, and applied by the sweep as
   * the loci go past. Which loci were dropped is reported by the sweep too:
   * it knows their names, and the scan deliberately does not keep them.
   */
  if(p.tol.val != 1) adzelog() << "Applying the missing-data filter...\n";
  scan.resolve(p.tol.val);

  if(p.tol.val != 1)
    {
      adzelog() << "\n"
		<< deletedHeader(scan.numLoci - scan.survivors,p.tol.val);
    }

  vector<Window> windows;

  if(p.windowed())
    {
      /*
       * Names are only there to point at the offending locus; the scan keeps
       * them for a locus map and not for a VCF, where the message falls back
       * to the position.
       */
      const string bad = scan.lmap.checkOrder(scan.locusName);
      if(!bad.empty())
	{
	  cerr << "ERROR: " << bad << "\n";
	  return EXIT_DATA;
	}

      LocusMap kept = scan.lmap;
      kept.compact(scan.dropped);
      buildWindows(kept,p,windows);

      if(windows.empty())
	{
	  cerr << "ERROR: no window holds at least --min-window-loci "
	       << p.min_win_loci.val << " loci.\n";
	  return EXIT_DATA;
	}

      adzelog() << "Laid out " << windows.size() << " window"
		<< (windows.size() == 1 ? "" : "s") << " over "
		<< scan.survivors << " loci.\n";

      adzelog() << "Completed at (d:h:m:s) ";
      displayTime(adzelog());
      adzelog() << endl;

      if(p.tnc.val)
	{
	  summary.close();
	  return EXIT_OK;
	}
    }

  p.loci.val = int(scan.survivors);

  if(scan.survivors == 0)
    {
      //The list is the diagnosis here, so it is written before giving up.
      writeDeletedLoci(p,scan,p.p_out.val);

      cerr << "ERROR: no locus survived filtering at TOLERANCE " << p.tol.val
	   << ".\n       Every statistic would be undefined; "
	   << "raise TOLERANCE or check MISSING.\n";
      summary << "ERROR: no locus survived filtering at TOLERANCE "
	      << p.tol.val << ".\n";
      summary.close();
      return EXIT_DATA;
    }

  /*
   * MAX_G is resolved over the surviving loci, not over all of them:
   * dropping loci with heavy missing data raises the smallest sample size.
   * The scan has both, having applied the filter to its own counts.
   */
  const int feasibleG = scan.feasibleG;

  if(!p.g.set)
    {
      p.g.val = (feasibleG < 1) ? 1 : feasibleG;
      adzelog() << "MAX_G not given; using " << p.g.val << ", the largest the "
		<< p.loci.val << " surviving loci support.\n";
      summary << "MAX_G resolved to " << p.g.val << " over " << p.loci.val
	      << (p.loci.val == 1 ? " locus\n" : " loci\n");
    }

  /*
   * Which groupings cannot reach MAX_G, and what holds them back. A run
   * whose ceiling is 1 produces nothing usable, and reporting only the
   * global minimum leaves a user staring at blank files with no idea which
   * grouping, or how many loci, caused it.
   */
  for(int j = 0; j < numDivs; j++)
    {
      if(scan.ceiling[j] >= p.g.val) continue;

      adzelog() << "WARNING: " << scan.groupName[j] << " scored only "
		<< scan.ceiling[j] << " gene cop"
		<< ((scan.ceiling[j] == 1) ? "y" : "ies")
		<< " at " << scan.binding[j]
		<< ((scan.binding[j] == 1) ? " locus" : " loci")
		<< " of " << p.loci.val;

      if(scan.emptyAt[j] > 0)
	{
	  adzelog() << " (" << scan.emptyAt[j]
		    << (scan.emptyAt[j] == 1 ? " locus is" : " loci are")
		    << " not scored at all)";
	}

      if(scan.ceiling[j] == 0)
	adzelog() << ",\n         so its statistics are undefined at every g.\n";
      else
	adzelog() << ",\n         so its statistics are undefined above g = "
		  << scan.ceiling[j] << ".\n";

      if(scan.ceiling[j] == 0)
	{
	  adzelog() << "         That leaves nothing usable for "
		    << scan.groupName[j] << ": lower --tolerance to drop "
		    << "those loci,\n         or --exclude-pops to leave the "
		    << "grouping out.\n";
	}
    }

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

  //Every requested tuple, whichever file it belongs to.
  if(do_tuple && !p.tuple_file.set)
    {
      for(list<int>::iterator i = k.begin(); i != k.end(); i++)
	{
	  ostringstream suffix;
	  suffix << "_" << *i;
	  tupleOut.push_back(nameCreate(p.c_out.val,suffix.str()));

	  vector< vector<int> > ofK;
	  buildKTuples(numDivs,*i,ofK);
	  for(size_t x = 0; x < ofK.size(); x++)
	    {
	      tuples.push_back(ofK[x]);
	      tupleFile.push_back(int(tupleOut.size())-1);
	    }
	}
    }

  /*
   * One announcement per statistic, then one pass computing all of them: a
   * locus's Q table serves every statistic, so they finish together and
   * there is one completion line rather than three.
   */
  if(do_rich) adzelog() << "Calculating allelic richness...\n";
  if(do_priv) adzelog() << "Calculating private allelic richness...\n";
  if(do_tuple)
    {
      if(p.tuple_file.set)
	{
	  adzelog() << "Calculating private alleles for " << tuples.size()
		    << " named tuples...\n";
	}
      else
	{
	  for(list<int>::iterator i = k.begin(); i != k.end(); i++)
	    {
	      adzelog() << "Calculating private alleles for all possible "
			<< *i << "-tuples...\n";
	    }
	}
    }

  try
    {
      /*
       * One sweep, unless the tuple accumulators would not fit: every tuple
       * carries a running accumulator per g and per window per g, all live
       * together, so a large k range is swept in batches. Each batch after
       * the first costs another pass over the loci and appends to the files
       * the first opened.
       */
      long long budget = 256LL*1024*1024;
      //A lever for the tests: batching is otherwise unreachable without a
      //dataset large enough to need hundreds of megabytes of accumulators.
      if(const char* forced = getenv("ADZE_TUPLE_BUDGET")) budget = atoll(forced);
      const long long batch = do_tuple
	? tupleBatch(p,scan,windows.size(),budget) : (long long)(tuples.size());
      const long long batches = tuples.empty()
	? 1 : (long long)((tuples.size() + batch - 1)/batch);

      if(batches > 1)
	{
	  adzelog() << "Sweeping " << tuples.size() << " tuples in " << batches
		    << " batches of " << batch << ", to keep their "
		    << "accumulators inside the memory budget.\n";
	  if(p.full_c.val)
	    {
	      adzelog() << "         With --full-tuples the per-locus rows "
			<< "arrive one batch at a time, so a locus's rows are\n"
			<< "         contiguous within a batch rather than "
			<< "across the whole file.\n";
	    }
	}

      if(!wantsVCF(p.format.val,p.dfile.val))
	{
	  countPath = countFilePath(p);
	  string why;
	  if(!countFileUsable(countPath,p,scan,why))
	    {
	      transposeStructure(p,scan,countPath,convertChunk(p,scan));
	    }
	}

      for(long long b = 0; b < batches; b++)
	{
	  const size_t lo = size_t(b*batch);
	  const size_t hi = (lo + size_t(batch) < tuples.size())
	    ? lo + size_t(batch) : tuples.size();

	  vector< vector<int> > slice(tuples.begin() + lo, tuples.begin() + hi);
	  vector<int> sliceFile(tupleFile.begin() + lo, tupleFile.begin() + hi);

	  LocusSource* src = wantsVCF(p.format.val,p.dfile.val)
	    ? openVCFSource(p,scan) : openCountFileSource(countPath,scan);

	  //Everything but the tuples is computed once, in the first batch.
	  sweepLoci(*src,scan,p,windows,scan.lmap,slice,sliceFile,tupleOut,
		    do_rich && b == 0,do_priv && b == 0,do_tuple,
		    p.r_out.val,p.p_out.val,
		    p.full_r.val,p.full_p.val,p.full_c.val,
		    (b == 0 && p.tol.val != 1) ? p.p_out.val : string(),
		    b > 0);

	  delete src;
	}
    }
  catch(BAD_FILE x)
    {
      if(!countPath.empty()) remove(countPath.c_str());
      return EXIT_IO;
    }
  catch(BAD_PARAM x)
    {
      if(!countPath.empty()) remove(countPath.c_str());
      return EXIT_DATA;
    }

  if(!countPath.empty()) remove(countPath.c_str());

  if(p.pp.val) adzelog() << endl;
  adzelog() << "Completed at (d:h:m:s) ";
  displayTime(adzelog());
  adzelog() << endl;

  summary << "Statistics completed at (d:h:m:s) ";
  displayTime(summary);
  summary << endl;

  adzelog() << "\nadze finished in (d:h:m:s) ";
  displayTime(adzelog());
  adzelog() << endl;

  summary << "\nadze finished in (d:h:m:s) ";
  displayTime(summary);
  summary << endl;
  summary.close();

  return EXIT_OK;
}
