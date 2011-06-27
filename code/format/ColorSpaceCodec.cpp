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
#include<algorithm>
#include<core/common_functions.h>
#include<format/ColorSpaceCodec.h>

/*
 * see http://www.ploscompbiol.org/article/slideshow.action?uri=info:doi/10.1371/journal.pcbi.1000386&imageURI=info:doi/10.1371/journal.pcbi.1000386.g002
 */

/* Note: by default, this will not trim off the first base + first colour space read */

/* Declare constant array values */
const char ColorSpaceCodec::bsBases[5] = {'A','C','G','T','N'};
const char ColorSpaceCodec::csColours[5] = {'0','1','2','3','4'};

ColorSpaceCodec::ColorSpaceCodec(){
}

/*
 * Convert colour-space character to a number. Anything other than [ACGT]
 * (i.e. double encoded) or [0123] is not considered a valid colour-space
 * character, so is converted to 4.
 */
int ColorSpaceCodec::csChrToInt(char tChr){
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
 * Convert base-space character to a number. Anything other than [ACGT] is not
 * considered a valid base-space character, so is converted to 4.
 */
int ColorSpaceCodec::bsChrToInt(char tChr){
	switch(toupper(tChr)){
	case 'A':
		return(0);
		break;
	case 'C':
		return(1);
		break;
	case 'G':
		return(2);
		break;
	case 'T':
		return(3);
		break;
	}
	return(4);
}



/*
 * decode color-space read into base-space, assuming it has a starting base.
 */
string ColorSpaceCodec::decodeCStoBS(string csInput, bool reverseComplement){
	if(csInput.length() == 0){
		return "";
	}
	string output("");
	output.reserve(csInput.length());
	// transfer first base directly to output
	if(reverseComplement){
		output += complementNucleotide(toupper(csInput.at(0)));
	} else {
		output += toupper(csInput.at(0));
	}
	for(unsigned int pos = 1; pos < csInput.length(); pos++){
		// get next bases, first (base) from output, second (colour) from input
		int mapX = bsChrToInt(output.at(pos-1));
		int mapY = csChrToInt(csInput.at(pos));
		if((mapX >= 4) || (mapY >= 4)){
			output += 'N';
		} else {
			// convert dimer into colour space
			if(mapX % 2 == 0){
				output += bsBases[(mapX + mapY) % 4];
			} else {
				// add 4 to avoid modulus becoming negative
				output += bsBases[((mapX - mapY) + 4) % 4];
			}
		}
	}
	if(reverseComplement){
		// reverse string in-place
		std::reverse(output.begin(), output.end());
	}
	return(output);
}

/*
 * decode color-space read into base-space, assuming it has a starting base.
 */
string ColorSpaceCodec::decodeCStoBS(string csInput){
	return(decodeCStoBS(csInput, false));
}


/*
 * Encode base-space sequence to colour space, with a starting base.
 */
string ColorSpaceCodec::encodeBStoCS(string bsInput){
	if(bsInput.length() == 0){
		return "";
	}
	// convert to A/C/G/T/N
	string output("");
	output.reserve(bsInput.length());
	// transfer first base directly to output
	output += bsInput.at(0);
	for(int pos = 0; pos < (bsInput.length()-1); pos++){
		// get next colours, both (bases) from input
		int mapX = bsChrToInt(bsInput.at(pos));
		int mapY = bsChrToInt(bsInput.at(pos+1));
		if((mapX >= 4) || (mapY >= 4)){
			output += 'N';
		} else {
			// convert dimer from colour space
			if(mapX % 2 == 0){
				output += csColours[(mapX + mapY) % 4];
			} else {
				// add 4 to avoid modulus becoming negative
				output += csColours[((mapX - mapY) + 4) % 4];
			}
		}
	}
	return(output);
}

bool ColorSpaceCodec::check(){
	ColorSpaceCodec cd;
	string csSeq1("T3.020000223213003122002213101232303020002301033000");
	string cConv1F("TANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");
	string csSeq2("T32002333220000303130320033020032123301032223033002");
	string cConv2F("TAGGGATATCTTTTTAATGCCGAAATAAGGGCTGATAACCGAGATTATTTC");
	string cConv2R("GAAATAATCTCGGTTATCAGCCCTTATTTCGGCATTAAAAAGATATCCCTA");
	cout << "1: checking colour-space decode (junk characters)... ";
	if(cd.decodeCStoBS(csSeq1).compare(cConv1F) != 0){
		cout << "Decode failed:\n"
				<< csSeq1 << endl << "->\n"
				<< cd.decodeCStoBS(csSeq1) << endl
				<< cConv1F << " expected\n";
	} else { cout << "SUCCESS!\n"; }
	cout << "2: checking colour-space decode (fully informative sequence)... ";
	if(cd.decodeCStoBS(csSeq2).compare(cConv2F) != 0){
		cout << "Decode failed:\n"
				<< csSeq2 << endl << "->\n"
				<< cd.decodeCStoBS(csSeq2) << endl
				<< cConv2F << " expected\n";
	} else { cout << "SUCCESS!\n"; }
	cout << "3: checking colour-space decode (inverse function actions)... ";
	if(cd.encodeBStoCS(cd.decodeCStoBS(csSeq2)).compare(csSeq2) != 0){
		cout << "Encode+Decode failed... encoding is not the inverse of decoding "
				<< "for fully informative input:\n"
				<< csSeq2 << endl << "->\n"
				<< cd.decodeCStoBS(csSeq2) << "->\n"
				<< cd.encodeBStoCS(cd.decodeCStoBS(csSeq2)) << endl
				<< csSeq2 << " expected\n";
	} else { cout << "SUCCESS!\n"; }
	cout << "4: checking colour-space decode (reverse decode)... ";
	if(cd.decodeCStoBS(csSeq2,true).compare(cConv2R) != 0){
		cout << "Reverse complement decode (CSRC) failed:\n"
				<< csSeq2 << endl << "->\n"
				<< cd.decodeCStoBS(csSeq2,true) << endl
				<< cConv2R << " expected\n";
	} else { cout << "SUCCESS!\n"; }
	return ((cd.decodeCStoBS(csSeq1).compare(cConv1F) == 0) &&
			(cd.decodeCStoBS(csSeq2).compare(cConv2F) == 0) &&
			(cd.encodeBStoCS(cd.decodeCStoBS(csSeq2)).compare(csSeq2) == 0) &&
			(cd.decodeCStoBS(csSeq2,true).compare(cConv2R) == 0)
			);
}
