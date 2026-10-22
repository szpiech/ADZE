#ifndef __ADZE_THREADS_H__
#define __ADZE_THREADS_H__

/*
 * A fixed set of worker threads that runs a batch of independent items and
 * waits for all of them.
 *
 * This is the only parallel pattern adze needs: a batch of loci, each worked
 * on independently, followed by a serial step that consumes the batch in
 * order. It is written on std::thread rather than OpenMP so that every build,
 * including the macOS release binary, can use it without a runtime library
 * the user would have to install.
 *
 * run(n, f) calls f(item, worker) once for each item in [0, n). worker is in
 * [0, size()) and identifies the thread, so f can use per-thread scratch
 * without locking. Items are handed out one at a time from a shared counter,
 * so the assignment of items to threads varies from run to run; f must write
 * only to storage belonging to its item or its worker. With one thread no
 * worker is started and run() is a plain loop on the calling thread.
 *
 * An exception thrown by f is caught, the remaining items are skipped, and
 * the first exception is rethrown by run() on the calling thread.
 */

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

class WorkerPool
{
 public:
  explicit WorkerPool(int threads)
    : nThreads(threads < 1 ? 1 : threads), job(0), items(0), next(0),
      generation(0), busy(0), stopping(false)
  {
    //The calling thread is worker 0, so n threads means n-1 started here.
    for(int w = 1; w < nThreads; w++)
      {
	workers.push_back(std::thread(&WorkerPool::loop, this, w));
      }
  }

  ~WorkerPool()
  {
    {
      std::lock_guard<std::mutex> lock(m);
      stopping = true;
    }
    wake.notify_all();
    for(size_t i = 0; i < workers.size(); i++) workers[i].join();
  }

  int size() const { return nThreads; }

  void run(size_t n, const std::function<void(size_t,int)>& f)
  {
    if(n == 0) return;

    if(nThreads == 1 || n == 1)
      {
	for(size_t i = 0; i < n; i++) f(i, 0);
	return;
      }

    {
      std::lock_guard<std::mutex> lock(m);
      job = &f;
      items = n;
      next.store(0);
      failure = std::exception_ptr();
      busy = int(workers.size());
      generation++;
    }
    wake.notify_all();

    work(0);                              //the caller takes items too

    std::unique_lock<std::mutex> lock(m);
    done.wait(lock, [this]{ return busy == 0; });
    job = 0;

    if(failure) std::rethrow_exception(failure);
  }

 private:
  WorkerPool(const WorkerPool&);
  WorkerPool& operator=(const WorkerPool&);

  void work(int worker)
  {
    for(;;)
      {
	const size_t i = next.fetch_add(1);
	if(i >= items) return;
	try
	  {
	    (*job)(i, worker);
	  }
	catch(...)
	  {
	    std::lock_guard<std::mutex> lock(m);
	    if(!failure) failure = std::current_exception();
	    next.store(items);                //skip what is left
	  }
      }
  }

  void loop(int worker)
  {
    unsigned long seen = 0;
    for(;;)
      {
	{
	  std::unique_lock<std::mutex> lock(m);
	  wake.wait(lock, [&]{ return stopping || generation != seen; });
	  if(stopping) return;
	  seen = generation;
	}

	work(worker);

	{
	  std::lock_guard<std::mutex> lock(m);
	  busy--;
	}
	done.notify_one();
      }
  }

  const int nThreads;
  std::vector<std::thread> workers;
  std::mutex m;
  std::condition_variable wake, done;
  const std::function<void(size_t,int)>* job;
  size_t items;
  std::atomic<size_t> next;
  unsigned long generation;
  int busy;
  bool stopping;
  std::exception_ptr failure;
};

#endif
