#include <string>
//#include "binom.h"
#include <list>
#include <vector>
#include <iomanip>
#include <iostream>

using namespace std;


/*
  A Class to handle diversity data in structure format
*/
class Population
{
 private:
  
  double tol;
  int numLociOrig;
  int numLoci; //The number of loci in the data
  int rows; //The number of individuals in the data
  string name; //The name of the population
  string** data; //The data block to be allocated with dimentions numLoci, rows
  string* locusName; //A vector to hold all the loci names
  vector<string> deletedLocus;//to hold the names of the deleted loci
  void fillData(string);
  void fillRichness(int);
  void fillNj(int);

  int** Nji; //A matrix whose entries correspond to the number of i alleles in the jth population
  int* NjiColLength; //The length of the columns (loci) in the matrix Nji
  int* Nj; //A vector holding the total number of alleles in the jth population
  
 public:

  int minG;
  Population();
  Population(string n, int numLoci, int rows);
  ~Population();

  void recLociDelete(double tolerance, const string& missing, vector<char>& del);

  void setRowsLoci(int,int);
  void setName(string str) {name = str;};
  bool setLocusName(string str, int pos);
  bool setAllelicRichness(double val, int locus);
  bool putDataElement(string dataElem, int row, int col);
  bool setNjiColLength(int,int);
  bool putNji(int,int,int);
  bool putNj(int,int);

  int getNj(int);
  int getNji(int,int);
  int getNjiColLength(int);
  string getName() {return name;};
  int getNumLoci() {return numLoci;};
  int getNumRows() {return rows;};
  const string& getDataElement(int row, int col) const;
  string getLocusName(int);

  void printDeleted(ostream&);
  void deleteLoci(const vector<char>& del);

  double calcAg(int,int);
  double calcQjig(int,int,int);
};
