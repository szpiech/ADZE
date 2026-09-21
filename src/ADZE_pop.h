#include <climits>
#include <string>
//#include "binom.h"
#include <list>
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace std;


/*
 * What every grouping knows identically about the loci, held once.
 *
 * Names: each Population used to own a string[numLoci], so J groupings held J
 * copies of identical names -- 32 bytes each even for a name short enough to
 * live inside the string object. At a million loci and ten groupings that was
 * some 290 MB of duplicates.
 *
 * Offsets: locus l occupies allele slots [offset[l],offset[l+1]) of each
 * grouping's count block. The slot count is a property of the locus, not of
 * the grouping -- the reader assigns one slot per allele observed anywhere --
 * so it too belongs here rather than in a per-grouping array.
 */
struct LocusTable
{
  vector<string> name;     //one per locus, compacted with the data
  vector<string> deleted;  //what the missing-data filter dropped, in report order
  vector<int> offset;      //size numLoci+1, prefix sums of the slot counts

  void compact(const vector<char>& del);
};

/*
  A Class to handle diversity data in structure format
*/
class Population
{
 private:
  
  double tol;
  int numLoci; //The number of loci in the data
  int rows; //The number of gene copies (data rows) in this grouping
  string name; //The name of the population
  const LocusTable* loci; //Shared across groupings; not owned here
  void fillNj(int);

  /*
   * The allele counts, one contiguous block rather than a row per locus.
   * A row per locus meant L allocations per grouping, each paying an
   * allocator header of 16-24 bytes and an 8-byte pointer to reach it --
   * more bookkeeping than counts for a biallelic locus. The block is indexed
   * through the shared offsets.
   */
  vector<int> nji;
  int* Nj; //A vector holding the total number of alleles in the jth population
  int* missing; //Number of missing gene copies at each locus

 public:

  int minG;
  Population();
  ~Population();

  void recLociDelete(double tolerance, vector<char>& del);

  void setLoci(int numLoci);
  void setName(string str) {name = str;};
  void setRows(int r) {rows = r;};
  void setTable(const LocusTable* table) {loci = table;};
  void allocNji();
  bool putNji(int,int,int);
  bool putNj(int,int);
  void putMissing(int count, int locus);
  void sumNj();

  int getNj(int);
  int getNji(int,int);
  int getNjiColLength(int);
  string getName() {return name;};
  int getNumLoci() {return numLoci;};
  int getNumRows() {return rows;};
  const string& getLocusName(int) const;

  string deletedSummary() const;
  void printDeleted(ostream&);
  void deleteLoci(const vector<char>& del, const vector<int>& oldOffset);

};
