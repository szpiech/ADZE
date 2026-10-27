/*
 * Moving a finished file into place.
 *
 * This is in a file of its own because the Windows version needs
 * <windows.h>, which defines macros (min, max, near, far, ERROR and others)
 * that could change the meaning of ordinary names in the rest of the code.
 */
#include "ADZE_fileops.h"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <sstream>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#endif

using namespace std;

/*
 * Puts a finished conversion at the name the user gave with --count-file.
 * It is written under a temporary name first and moved when complete, so a
 * run that is reading the file, or one that checks it for reuse, never sees
 * it half written. On POSIX systems rename() replaces the old file in one
 * step, and a run that already has it open keeps reading the old one.
 *
 * Windows needs MoveFileEx to replace a file in one step. Removing the old
 * file and then renaming leaves a gap in which another run can put its own
 * file there, and the rename then fails because the name exists. Windows
 * also refuses to replace a file that another run has open. Both cases mean
 * another run has just published a conversion, so this returns false with
 * the reason in err, and the caller decides whether that file will do.
 * Either way the temporary file is removed.
 */
bool publishFile(const string& tmp,const string& path,string& err)
{
#ifdef _WIN32
  const bool moved = MoveFileExA(tmp.c_str(),path.c_str(),
				 MOVEFILE_REPLACE_EXISTING) != 0;
  if(!moved)
    {
      ostringstream e;
      e << "Windows error " << GetLastError();
      err = e.str();
    }
#else
  const bool moved = rename(tmp.c_str(),path.c_str()) == 0;
  if(!moved) err = strerror(errno);
#endif
  if(!moved) remove(tmp.c_str());
  return moved;
}
