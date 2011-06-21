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

#include<stdlib.h>
#include<format/ColorSpaceLoader.h>
#include<fstream>
#include<core/common_functions.h>
#include<iostream>
#include<string>
using namespace std;

#define _ENCODING_CHAR_A '0'
#define _ENCODING_CHAR_T '1'
#define _ENCODING_CHAR_C '2'
#define _ENCODING_CHAR_G '3'

int ColorSpaceLoader::open(string file){
	m_f.open(file.c_str());
	m_size=0;
	m_loaded=0;

	string bufferForLine;
	while(!m_f.eof()){
		getline(m_f,bufferForLine);
		if(bufferForLine.at(0) == '#'){
			continue;// skip csfasta comment
		}
		if(bufferForLine.at(0) == '>'){
			getline(m_f,bufferForLine);
			m_size++;
		}
	}

	m_f.close();
	m_f.open(file.c_str());
	return EXIT_SUCCESS;
}

void ColorSpaceLoader::load(int maxToLoad,ArrayOfReads*reads,MyAllocator*seqMyAllocator){
	string bufferForLine;
	int loadedSequences = 0;
	while((m_loaded < m_size) && (loadedSequences < maxToLoad)){
		getline(m_f, bufferForLine);
		if(bufferForLine.at(0) == '#'){
			continue;// skip csfasta comment
		}
		// read two lines
		if(bufferForLine.at(0) == '>'){
			getline(m_f, bufferForLine);
			string decodedLine = m_decoder.decode(bufferForLine);
			Read t;
			t.constructor(decodedLine.c_str(),seqMyAllocator,true);
			reads->push_back(&t);
			loadedSequences++;
			m_loaded++;
		}
	}
	if(m_loaded==m_size){
		m_f.close();
	}
}

int ColorSpaceLoader::getSize(){
	return m_size;
}
