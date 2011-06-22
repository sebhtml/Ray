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

class ColorSpaceDecoder{
  static const char csBases[5];
  static const char csValues[5];
  static const char csMap[25];
 public:
  // could be condensed into a char[5], but a 5x5 array is easier to understand
	ColorSpaceDecoder();
  int csChrToInt(char tChr);
	string decode(string x);
	string decodeRC(string csInputRC, bool reverseComplement);
	string encode(string x);
	string encodeRC(string bsInput);
  bool check();
};

#endif
