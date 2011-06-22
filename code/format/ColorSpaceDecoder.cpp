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
const char ColorSpaceDecoder::bsBases[5] = {'A','C','G','T','N'};
const char ColorSpaceDecoder::csColours[5] = {'0','1','2','3','4'};

ColorSpaceDecoder::ColorSpaceDecoder(){
}

/*
 * Convert colour-space character to a number. Anything other than [ACGT]
 * (i.e. double encoded) or [0123] is not considered a valid colour-space
 * character, so is converted to 4.
 */
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
 * Convert base-space character to a number. Anything other than [ACGT] is not
 * considered a valid base-space character, so is converted to 4.
 */
int ColorSpaceDecoder::bsChrToInt(char tChr){
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
string ColorSpaceDecoder::decode(string x){
	if(x.length() == 0){
		return "";
	}
	string output("");
	// transfer first base directly to output
	output += toupper(x.at(0));
	for(unsigned int pos = 1; pos < x.length(); pos++){
		// get next bases, first (base) from output, second (colour) from input
		int mapX = bsChrToInt(output.at(pos-1));
		int mapY = csChrToInt(x.at(pos));
		if((mapX >= 4) || (mapY >= 4)){
			output += 'N';
		} else {
			// convert dimer into colour space
			// add 4 to avoid modulus becoming negative
			output += bsBases[((mapX - mapY) + 4) % 4];
		}
	}
	return(output);
}

/*
 * decode color-space read, assuming an sequence has a starting base,
 * and a reverse-complement starting base. The assumed input format is as follows:
 * [ACGT][0-3]+[ACGT]
 * i.e. <starting base><colour-space encoding><reverse complement starting base>
 */
string ColorSpaceDecoder::decodeCSRCtoBS(string csInputRC, bool reverseComplement){
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
 * decode color-space read (forward direction), assuming an sequence has a starting base,
 * and a reverse-complement starting base. The assumed input format is as follows:
 * [ACGT][0-3]+[ACGT]
 * i.e. <starting base><colour-space encoding><reverse complement starting base>
 */
string ColorSpaceDecoder::decodeCSRCtoBS(string csInputRC){
	return(decodeCSRCtoBS(csInputRC, false));
}

/*
 * Encode base-space sequence to colour space, with a starting base.
 */
string ColorSpaceDecoder::encode(string x){
	if(x.length() == 0){
		return "";
	}
	// convert to A/C/G/T/N
	string output("");
	output.reserve(x.length());
	// transfer first base directly to output
	output += x.at(0);
	for(int pos = 0; pos < (x.length()-1); pos++){
		// get next colours, both (bases) from input
		int mapX = bsChrToInt(x.at(pos));
		int mapY = bsChrToInt(x.at(pos+1));
		if((mapX >= 4) || (mapY >= 4)){
			output += 'N';
		} else {
			// convert dimer from colour space
			// add 4 to avoid modulus becoming negative
			output += csColours[((mapX - mapY) + 4) % 4];
		}
	}
	return(output);
}

/*
 * Convert base-pair sequence to colour space, with a starting base and ending reverse-complement base.
 */
string ColorSpaceDecoder::encodeBStoCSRC(string bsInput){
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
	string csSeq3("T32002333220000303130320033020032123301032223033002G");
	string cConv3F("TAGGGACGTCTTTTTAACACCGAAACGGAAACTGACGGCCGAGACCGTTTC");
	string cConv3R("GAAACGGTCTCGGCCGTCAGTTTCCGTTTCGGTGTTAAAAAGACGTCCCTA");
	cout << "1: checking basic decode (junk characters)... ";
	if(cd.decode(csSeq1).compare(cConv1) != 0){
		cout << "Decode failed:"
				<< csSeq1 << endl << "->\n"
				<< cd.decode(csSeq1) << endl
				<< cConv1 << " expected";
	} else { cout << "SUCCESS!\n"; }
	cout << "2: checking basic decode (fully informative sequence)... ";
	if(cd.decode(csSeq2).compare(cConv2) != 0){
		cout << "Decode failed:"
				<< csSeq2 << endl << "->\n"
				<< cd.decode(csSeq2) << endl
				<< cConv2 << " expected";
	} else { cout << "SUCCESS!\n"; }
	cout << "3: checking basic decode (inverse function actions)... ";
	if(cd.encode(cd.decode(csSeq2)).compare(csSeq2) != 0){
		cout << "Encode+Decode failed... encoding is not the inverse of decoding "
				<< "for fully informative input:"
				<< csSeq2 << endl << "->\n"
				<< cd.decode(csSeq2) << "->\n"
				<< cd.encode(cd.decode(csSeq2)) << endl
				<< csSeq2 << " expected";
	} else { cout << "SUCCESS!\n"; }
	cout << "4: checking CSRC decode (forward decode)... ";
	if(cd.decodeCSRCtoBS(csSeq3).compare(cConv3F) != 0){
		cout << "Forward decode (CSRC) failed:"
				<< csSeq3 << endl << "->\n"
				<< cd.decodeCSRCtoBS(csSeq3) << endl
				<< cConv3F << " expected";
	} else { cout << "SUCCESS!\n"; }
	cout << "5: checking CSRC decode (reverse decode)... ";
	if(cd.decodeCSRCtoBS(csSeq3,true).compare(cConv3R) != 0){
		cout << "Reverse complement decode (CSRC) failed:"
				<< csSeq3 << endl << "->\n"
				<< cd.decodeCSRCtoBS(csSeq3,true) << endl
				<< cConv3R << " expected";
	} else { cout << "SUCCESS!\n"; }
	cout << "6: checking CSRC decode (inverse function actions)... ";
	if(cd.encodeBStoCSRC(cd.decodeCSRCtoBS(csSeq3)).compare(csSeq3) != 0){
		cout << "Encode+Decode failed... encoding is not the inverse of decoding "
				<< "for fully informative CSRC input:"
				<< csSeq3 << endl << "->\n"
				<< cd.decodeCSRCtoBS(csSeq3) << "->\n"
				<< cd.encodeBStoCSRC(cd.decodeCSRCtoBS(csSeq3)) << endl
				<< csSeq3 << " expected";
	} else { cout << "SUCCESS!\n"; }
	return ((cd.decode(csSeq1).compare(cConv1) == 0) &&
			(cd.decode(csSeq2).compare(cConv2) == 0) &&
			(cd.encode(cd.decode(csSeq2)).compare(csSeq2) == 0) &&
			(cd.decodeCSRCtoBS(csSeq3,false).compare(cConv3F) == 0) &&
			(cd.decodeCSRCtoBS(csSeq3,true).compare(cConv3R) == 0) &&
			(cd.encodeBStoCSRC(cd.decodeCSRCtoBS(csSeq3,false)).compare(csSeq3) == 0)
			);
}
