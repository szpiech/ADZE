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

