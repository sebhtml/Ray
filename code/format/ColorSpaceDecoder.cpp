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
#include<core/common_functions.h>
#include<format/ColorSpaceDecoder.h>

/*
 * see http://www.ploscompbiol.org/article/slideshow.action?uri=info:doi/10.1371/journal.pcbi.1000386&imageURI=info:doi/10.1371/journal.pcbi.1000386.g002
 */

/* Note: by default, this will not trim off the first base + first colour space read */

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
 * decode color-space read, assuming it has a starting base.
 */
string ColorSpaceDecoder::decode(string x){
	if(x.length() == 0){
		return "";
	}
	string output("");
	// convert first base to A/C/G/T/N
	output += csBases[csChrToInt(x.at(0))];
	for(unsigned int pos = 1; pos < x.length(); pos++){
		// get next bases, first from output, second from input
		int mapX = csChrToInt(output.at(pos-1));
		int mapY = csChrToInt(x.at(pos));
		if((mapX >= 4) || (mapY >= 4)){
			output += 'N';
		} else {
			// convert dimer into colour space
			// add 4 to avoid modulus becoming negative
			output += csBases[((mapX - mapY) + 4) % 4];
		}
	}
	return(output);
}

/*
 * decode color-space read, assuming an sequence has a starting base,
 * and a reverse-complement starting base. The assumed input format is as follows:
 * [ACGT][0-3]+[ACGT]
 * i.e. <starting base><colour-space encoding><reverse complement starting base>
 *
 */
string ColorSpaceDecoder::decodeRC(string csInputRC, bool reverseComplement){
	string basicCS("");
	basicCS.reserve(csInputRC.length()); // reserve space for input
	if(reverseComplement){
		// reverse input string
		for(string::reverse_iterator rit = csInputRC.rbegin();
				rit < csInputRC.rend(); rit++){
			basicCS+= *rit;
		}
		// trim off forward start base (which is at the end of the reversed string)
		basicCS.erase(basicCS.length() - 1, 1);
	} else {
		basicCS = csInputRC.substr(0,csInputRC.length()-1); // trim off RC start base
	}
	return(decode(basicCS));
}

/*
 * Convert base-pair sequence to colour space, with a starting base.
 */
string ColorSpaceDecoder::encode(string x){
	if(x.length() == 0){
		return "";
	}
	// convert to A/C/G/T/N
	string output("");
	output.reserve(x.length());
	output += csBases[csChrToInt(x.at(0))];
	for(int pos = 0; pos < (x.length()-1); pos++){
		// get next colours, both from input
		int mapX = csChrToInt(x.at(pos));
		int mapY = csChrToInt(x.at(pos+1));
		if((mapX >= 4) || (mapY >= 4)){
			output += 'N';
		} else {
			// convert dimer from colour space
			// add 4 to avoid modulus becoming negative
			output += csValues[((mapX - mapY) + 4) % 4];
		}
	}
	return(output);
}

/*
 * Convert base-pair sequence to colour space, with a starting base and ending reverse-complement base.
 */
string ColorSpaceDecoder::encodeRC(string bsInput){
	string output = encode(bsInput);
	output.reserve(bsInput.length()+1);
	output += complementNucleotide(bsInput.at(bsInput.length()-1));
	return(output);
}


bool ColorSpaceDecoder::check(){
	ColorSpaceDecoder cd;
	string csSeq1("T3.020000223213003122002213101232303020002301033000");
	string cConv1("TANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
	string csSeq2("T32002333220000303130320033020032123301032223033002");
	string cConv2("TAGGGACGTCTTTTTAACACCGAAACGGAAACTGACGGCCGAGACCGTTTC");
	if(cd.decode(csSeq1).compare(cConv1) != 0){
		cout << "Decode failed:"
				<< csSeq1 << endl << "->\n"
				<< cd.decode(csSeq1) << endl
				<< cConv1 << " expected";
	}
	if(cd.decode(csSeq2).compare(cConv2) != 0){
		cout << "Decode failed:"
				<< csSeq2 << endl << "->\n"
				<< cd.decode(csSeq2) << endl
				<< cConv2 << " expected";
	}
	if(cd.encode(cd.decode(csSeq2)).compare(csSeq2) != 0){
		cout << "Encode+Decode failed... encoding is not the inverse of decoding "
				<< "for fully informative input:"
				<< csSeq2 << endl << "->\n"
				<< cd.decode(csSeq2) << "->\n"
				<< cd.encode(cd.decode(csSeq2)) << endl
				<< csSeq2 << " expected";
	}
	return ((cd.decode(csSeq1).compare(cConv1) == 0) &&
			(cd.decode(csSeq2).compare(cConv2) == 0) &&
			(cd.encode(cd.decode(csSeq2)).compare(csSeq2) == 0));
}
