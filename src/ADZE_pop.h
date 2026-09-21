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
 * The locus names, held once for the whole run.
 *
 * Every Population used to own a string[numLoci] of its own, so J groupings
 * held J copies of identical names -- 32 bytes each even for a name short
 * enough to live inside the string object, and more when it is not. At a
 * million loci and ten groupings that was some 290 MB of duplicates. One
 * table, pointed at by every grouping, costs the same as one grouping used
 * to.
 */
struct LocusNames
{
  vector<string> name;     //one per locus, compacted with the data
  vector<string> deleted;  //what the missing-data filter dropped, in report order

  void compact(const vector<char>& del);
};

/*
  A Class to handle diversity data in structure format
*/
class Population
{
 private:
  
  double tol;
  int numLociOrig;
  int numLoci; //The number of loci in the data
  int rows; //The number of gene copies (data rows) in this grouping
  string name; //The name of the population
  LocusNames* names; //Shared across groupings; not owned here
  void fillNj(int);

  int** Nji; //A matrix whose entries correspond to the number of i alleles in the jth population
  int* NjiColLength; //The length of the columns (loci) in the matrix Nji
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
  void setNames(LocusNames* table) {names = table;};
  bool setNjiColLength(int,int);
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
  void deleteLoci(const vector<char>& del);

};
