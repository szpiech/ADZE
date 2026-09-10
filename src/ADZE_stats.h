#include <string>
#include <cmath>
#include <limits>
#include <iostream>

using namespace std;

class Stats
{
private:
  double* data;
  double avg;
  double var;
  double std_err;
  int numLoci;
  bool avg_flag;
  bool var_flag;
  bool std_err_flag;

public:
  //  void allocate(int,int);
  void printStats(ostream&,string,int);
  void printData(ostream&,string,int);
  void calcAvg();
  void calcVar();
  void calcStdErr();
  void putData(double*,int);
  ~Stats();
  Stats();
};
