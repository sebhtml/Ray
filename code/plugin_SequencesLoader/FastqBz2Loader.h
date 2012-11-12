/*
 	Ray
    Copyright (C) 2010, 2011, 2012 Sébastien Boisvert

	http://DeNovoAssembler.SourceForge.Net/

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You have received a copy of the GNU General Public License
    along with this program (gpl-3.0.txt).  
	see <http://www.gnu.org/licenses/>

*/

#ifndef _FastqBz2Loader
#define _FastqBz2Loader

#ifdef CONFIG_HAVE_LIBBZ2

#include "ArrayOfReads.h"
#include "Read.h"
#include "BzReader.h"

#include <RayPlatform/memory/MyAllocator.h>

#include <string>
#include <vector>
using namespace std;

/**
 * \author Sébastien Boisvert
 */
class FastqBz2Loader{
	int m_loaded;
	int m_size;
	BzReader m_reader;
public:
	int open(string file,int period);
	int getSize();
	void load(int maxToLoad,ArrayOfReads*reads,MyAllocator*seqMyAllocator,int period);
};

#endif
#endif

