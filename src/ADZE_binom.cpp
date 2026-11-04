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
#include "ADZE_binom.h"

#include <climits>

using namespace std;

static unsigned long long gcd(unsigned long long a,unsigned long long b)
{
  while(b != 0)
    {
      const unsigned long long t = a % b;
      a = b;
      b = t;
    }
  return a;
}

/*
 * The number of ways to choose k of n items.
 *
 * Built up one factor at a time: after step i the running value is
 * C(n-k+i, i), always an integer, so the arithmetic is exact for as long as
 * the coefficient fits in 64 bits. Each step divides out the common factor of
 * the running value and i before multiplying, which keeps the intermediate
 * product no larger than the next coefficient. Should a coefficient pass
 * 2^64, the remaining steps continue in long double, so a count too large to
 * be exact is still a close approximation rather than a wrapped integer. The
 * result is returned as a double, which holds it exactly up to 2^53.
 *
 * Used for the number of combinations of groupings that a run would
 * evaluate, which --dry-run reports.
 */
double choose(int n,int k)
{
  if(k < 0 || n < 0 || k > n) return 0;
  if(k > n - k) k = n - k;

  unsigned long long r = 1;
  int i = 1;
  for(; i <= k; i++)
    {
      const unsigned long long g = gcd(r,(unsigned long long)(i));
      const unsigned long long d = (unsigned long long)(i) / g;
      const unsigned long long m = (unsigned long long)(n - k + i) / d;
      const unsigned long long q = r / g;
      if(m != 0 && q > ULLONG_MAX / m) break;   //r is still C(n-k+i-1, i-1)
      r = q * m;
    }
  if(i > k) return double(r);

  long double x = (long double)(r);
  for(; i <= k; i++) x = x * (long double)(n - k + i) / (long double)(i);
  return double(x);
}
