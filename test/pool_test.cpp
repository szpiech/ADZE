// The thread pool in src/ADZE_threads.h, on its own.
//
// Checks, at several pool sizes, that every item runs exactly once, that the
// worker index is in range, that results written per item do not depend on
// the pool size, that the pool can be reused for many batches, and that an
// exception thrown by one item reaches the caller. Exit status 0 on success.
#include <cstdio>
#include <stdexcept>
#include <vector>
#include "ADZE_threads.h"

static int failures = 0;

static void check(bool ok, const char* what, int threads)
{
  if(!ok)
    {
      std::printf("FAIL  %s (threads %d)\n", what, threads);
      failures++;
    }
}

int main()
{
  const int sizes[] = {1, 2, 3, 8, 64};
  std::vector<double> reference;

  for(int t : sizes)
    {
      WorkerPool pool(t);
      check(pool.size() == t, "size", t);

      for(size_t n : {size_t(0), size_t(1), size_t(7), size_t(1000)})
	{
	  std::vector<int> hits(n, 0);
	  std::vector<int> badWorker(1, 0);
	  pool.run(n, [&](size_t i, int w) {
	    hits[i]++;
	    if(w < 0 || w >= t) badWorker[0] = 1;
	  });
	  bool once = true;
	  for(size_t i = 0; i < n; i++) if(hits[i] != 1) once = false;
	  check(once, "every item exactly once", t);
	  check(!badWorker[0], "worker index in range", t);
	}

      //Per-item results, then a serial fold in item order: the pattern the
      //sweep uses. The fold must give the same bits at every pool size.
      std::vector<double> v(5000);
      double sum = 0;
      for(int batch = 0; batch < 200; batch++)
	{
	  pool.run(v.size(), [&](size_t i, int) {
	    double x = 1.0 / double(i + 1 + batch);
	    for(int k = 0; k < 20; k++) x = x * 0.999 + 1e-7 * double(k);
	    v[i] = x;
	  });
	  for(size_t i = 0; i < v.size(); i++) sum += v[i];
	}
      if(reference.empty()) reference.push_back(sum);
      check(sum == reference[0], "ordered fold identical", t);

      bool caught = false;
      try
	{
	  pool.run(100, [&](size_t i, int) {
	    if(i == 37) throw std::runtime_error("item 37");
	  });
	}
      catch(const std::runtime_error&) { caught = true; }
      check(caught, "exception reaches the caller", t);

      //and the pool still works afterwards
      std::vector<int> after(50, 0);
      pool.run(after.size(), [&](size_t i, int) { after[i] = 1; });
      bool all = true;
      for(int x : after) if(!x) all = false;
      check(all, "usable after an exception", t);
    }

  std::printf("%s\n", failures ? "pool: FAILED" : "pool: all checks passed");
  return failures ? 1 : 0;
}
