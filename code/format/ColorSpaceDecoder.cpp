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

#include<string>
#include<iostream>
#include<format/ColorSpaceDecoder.h>

/*
 * see http://www.ploscompbiol.org/article/slideshow.action?uri=info:doi/10.1371/journal.pcbi.1000386&imageURI=info:doi/10.1371/journal.pcbi.1000386.g002
 */

/* This decoder is approximately equivalent to the python code from Galaxy library lib/galaxy_utils/sequence/transform.py */

/* Declare constant array values */
const char ColorSpaceDecoder::csBases[5] = {'A','C','G','T','N'};
const char ColorSpaceDecoder::csValues[5] = {'0','1','2','3','4'};

ColorSpaceDecoder::ColorSpaceDecoder(){
}

int ColorSpaceDecoder::csChrToInt(char tChr){
  switch(toupper(tChr)){
  case 'A':
  case '0':
    return(0);
    break;
  case 'C':
  case '1':
    return(1);
    break;
  case 'G':
  case '2':
    return(2);
    break;
  case 'T':
  case '3':
    return(3);
    break;
  }
  return(4);
}

/*
 * decode color-space read.
 */
string ColorSpaceDecoder::decode(string x){
  if(x.length() == 0){
    return "";
  }
  string output("");
  // convert first base to A/C/G/T/N
  output += csBases[csChrToInt(x.at(0))];
  for(int pos = 1; pos < x.length(); pos++){
    int mapX = csChrToInt(output.at(pos-1));
    int mapY = csChrToInt(x.at(pos));
    //cout << "pos=" << pos << ",mapX=" << mapX << ",mapY=" << mapY << endl;
    if((mapX >= 4) || (mapY >= 4)){
      output += 'N';
    } else {
      // add 4 to avoid modulus becoming negative
      output += csBases[((mapX - mapY) + 4) % 4];
    }
  }
  return(output);
}

/*
 * Convert base-pair sequence to colour space.
 */
string ColorSpaceDecoder::encode(string x){
  if(x.length() == 0){
    return "";
  }
  // convert to A/C/G/T/N
  string output("");
  output += csBases[csChrToInt(x.at(0))];
  for(int pos = 0; pos < (x.length()-1); pos++){
    int mapX = csChrToInt(x.at(pos));
    int mapY = csChrToInt(x.at(pos+1));
    if((mapX >= 4) || (mapY >= 4)){
      output += 'N';
    } else {
      // add 4 to avoid modulus becoming negative
      output += csValues[((mapX - mapY) + 4) % 4];
    }
  }
  return(output);
}

bool ColorSpaceDecoder::check(){
  ColorSpaceDecoder cd;
  string csSeq1("T3.020000223213003122002213101232303020002301033000");
  string cConv1("TANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
  string csSeq2("T32002333220000303130320033020032123301032223033002");
  string cConv2("TAGGGACGTCTTTTTAACACCGAAACGGAAACTGACGGCCGAGACCGTTTC");
  return ((cd.decode(csSeq1).compare(cConv1) == 0) &&
          (cd.decode(csSeq2).compare(cConv2) == 0) &&
          (cd.encode(cd.decode(csSeq2)).compare(csSeq2) == 0));
}

// int main(){
//   ColorSpaceDecoder cd;
//   cout << "The result of the check is " << (cd.check()?"true":"false") << endl;
// }
