/*
 * ADZE: Allelic Diversity Analyzer
 * Copyright (C) 2008-2026 Zachary A Szpiech
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
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
  void printStats(ostream&,const string&,int,bool tabbed = 0);
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
