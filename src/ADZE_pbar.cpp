#include "ADZE_pbar.h"

using namespace std;

void ProgressBar::adv(int n)
{
  //progress += n;
  for(int i = 0; i < n; i++) this->operator++();
  return;
}

void ProgressBar::operator++()
{
  if(full) return;
  //static double i = 0;
  //static time_t start,end;
  progress++;
  /*
  if(progress == 1) time(&start);
  if(progress == 100)
    {
      time(&end);
      double diff = difftime(end,start);
      double t = (diff/100)*hits;
      *out << t;
    }
  */

  if(progress/hits <= double(filled)*frac) return;

  while(progress/hits > double(filled)*frac)
    {
      update();
      if(progress/hits >= 1)
	{
	  done();
	  return;
	}    
    }
  
  /*
  if(i/hits > filled*frac) update();
  i++;
  */
  if(progress/hits >= 1)
    {
      done();
    }

  return;
}

void ProgressBar::estimate()
{
  static time_t start, end;
  double avg = 0;
  
  if(filled == 1) time(&start);   
  time(&end);
  /*
  cout.unsetf(ios::fixed);
  cout.unsetf(ios::showpoint);
  */
  *out << est;

  diff += difftime(end,start);
  avg = diff/filled;
  double t = avg*(BARLEN-filled);
  out->setf(ios::left,ios::adjustfield);
  *out << setw(2);

  if(floor(t/86400.0) > 0)
    {
      *out << t/86400.0 << " days.";
    }
  else if(floor(t/3600.0))
    {
      *out << t/3600.0 << " hours.";
    }
  else if(floor(t/60.0))
    {
      *out << t/60.0 << " minutes.";
    }
  else if(filled == 1)
    {
      *out << "??.??";
    }
  else
    {
      *out << t << " seconds.";
    }
  
  outchar(10,' ');

  /*
  *out << setw(5) << floor(t/86400.0) << ":";
  t = t - floor(t/86400.0)*86400.0;
  *out << setw(2) << floor(t/3600.0) << ":";
  t = t - floor(t/3600.0)*3600.0;
  *out << setw(2) <<floor(t/60.0) << ":";
  t = t - floor(t/60.0)*60.0;
  */  
  //*out << left << setw(10) << t; 
  out->flush();

  time(&start);
  
  return;
}


void ProgressBar::update()
{
  if(full) return;
  ++filled;

  clr();
  out->setf(ios::fixed,ios::floatfield);
  out->setf(ios::showpoint);
  *out << setprecision(2) << (double(filled)/BARLEN)*100 << "%";

  estimate();

  /*
  *out << "=|";
  outchar(++filled,'#');
  outchar(--blank,' ');
  *out << "|=";
  //if(timer) estimate();
  */  


  out->flush();
  
  if(filled == BARLEN) full = 1;
  
  return;
}

void ProgressBar::clr(int n)
{
    for(int i = 0; i < n; i ++)
    {
      *out << '\b'; //<< '\b';
      out->flush();
      /*
       *out << '\b' << '\b';
      out->flush();
      */
    }
}

void ProgressBar::outchar(int n, char c)
{
  for(int i = 0; i < n; i++) *out << c;
  return;
}

ProgressBar::ProgressBar(ostream* o, double n, double blen, bool t)
{
  diff = 0;
  timer = t;
  progress = 0;
  full = 0;
  out = o;
  hits = n;
  filled = 0;
  BARLEN = blen;
  frac = 1/BARLEN;
  blank = int(BARLEN);
  width = SCRLEN - (int(BARLEN)+est.size()+4);
}

void ProgressBar::init()
{
  if(full) return;
  out->setf(ios::fixed,ios::floatfield);
  out->setf(ios::showpoint);

  *out << setprecision(2)<< 0.0 << "%";

  /*
  *out << "=|";
  outchar(int(BARLEN),' ');
  *out << "|=";
  if(timer)
    {
      *out << est;
    }
  */
  out->flush();

  return;
}

void ProgressBar::done()
{
  clr();
  out->setf(ios::fixed,ios::floatfield);
  out->setf(ios::showpoint);
  *out << setprecision(2) << 100.0 << "% Finished.";
  outchar(20,' ');

  /*
  *out << "=|";
  outchar(int(BARLEN),'#');
  *out << "|=";
  */
  
  out->flush();
  full = 1;  
  
  cout.unsetf(ios::fixed);
  cout.unsetf(ios::showpoint);
  //cout.unsetf(cout.precision);
  
  //*out << fixed;
  return;
}
