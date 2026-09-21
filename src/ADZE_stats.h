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
  void printStats(ostream&,string,int,bool tabbed = 0);
  void calcAvg();
  void calcVar();
  void calcStdErr();
  void putData(double*,int);
  /*
   * Set the summary directly. The passes accumulate mean and variance as
   * they sweep rather than storing the per-locus values, so there is nothing
   * for calcAvg/calcVar to read; -9 in all three marks a statistic that is
   * undefined at this g, exactly as the two-pass path produced it.
   */
  void putSummary(double a,double v,double s,int n)
  {avg = a; var = v; std_err = s; numLoci = n; data = NULL;
   avg_flag = 1; var_flag = 1; std_err_flag = 1;}
  ~Stats();
  Stats();
};
