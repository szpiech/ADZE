/*
 * ADZE: Allelic Diversity Analyzer
 * Copyright (C) 2026 Zachary A Szpiech
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
#ifndef __ADZE_FILEOPS_H__
#define __ADZE_FILEOPS_H__

#include <string>

//Moves tmp to path, replacing any file there, in one step. On failure,
//removes tmp, puts the reason in err and returns false. See ADZE_fileops.cpp.
bool publishFile(const std::string& tmp,const std::string& path,
		 std::string& err);

#endif
