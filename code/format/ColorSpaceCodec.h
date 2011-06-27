/*
 	Ray
    Copyright (C)  2010  Sébastien Boisvert

	http://DeNovoAssembler.SourceForge.Net/

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You have received a copy of the GNU General Public License
    along with this program (COPYING).  
	see <http://www.gnu.org/licenses/>

*/

#ifndef _ColorSpaceDecoder
#define _ColorSpaceDecoder

#include<string>
using namespace std;

class ColorSpaceCodec{
  static const char csColours[5];
  static const char bsBases[5];
 public:
	ColorSpaceCodec();
	int csChrToInt(char tChr);
	int bsChrToInt(char tChr);
	string decodeCStoBS(string csInput);
	string decodeCStoBS(string csInput, bool reverseComplement);
	string encodeBStoCS(string bsInput);
  bool check();
};

#endif
