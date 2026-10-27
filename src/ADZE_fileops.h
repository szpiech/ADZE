#ifndef __ADZE_FILEOPS_H__
#define __ADZE_FILEOPS_H__

#include <string>

//Moves tmp to path, replacing any file there, in one step. On failure,
//removes tmp, puts the reason in err and returns false. See ADZE_fileops.cpp.
bool publishFile(const std::string& tmp,const std::string& path,
		 std::string& err);

#endif
