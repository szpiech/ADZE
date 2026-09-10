#include <iostream>
#include <time.h>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

const char BKSPC = '\b';
const char SPC = ' ';
const int SCRLEN = 80;
const string est = " Finish in ~";

class ProgressBar
{
 private:
  int width;
  double diff;
  bool timer;
  double progress;
  double frac;
  double BARLEN;
  int filled;
  int blank;
  bool full;
  double hits;
  ostream* out;
  void outchar(int n, char c); //Outputs n of char c to out
  void clr(int n = 80);
  void update();
  void estimate();

 public:
  ProgressBar(ostream* o,double n, double blen = 76, bool t = 0);
  void init();
  void operator++();
  void done();
  void adv(int n);
};

