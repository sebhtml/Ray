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

// Find number of sequences in file
int ColorSpaceLoader::open(string file){
	m_f.open(file.c_str());
	m_size=0;
	m_loaded=0;
	int lineMod4 = 0;
	m_ft = UNKNOWN;
	bool firstLine = true;
	string bufferForLine;
	while(!m_f.eof()){
		getline(m_f,bufferForLine);
		if(bufferForLine.length() == 0){
			// skip over empty lines (including at end of file)
			continue;
		}
		if((m_ft == UNKNOWN) && (bufferForLine.at(0) == '#')){
			continue;// skip initial comment lines
		}
		lineMod4 = ((lineMod4 + 1) % 4);
		if(firstLine){
			if(bufferForLine.at(0) == '>'){
				m_ft = FASTA;
			}
			if(bufferForLine.at(0) == '@'){
				m_ft = FASTQ;
			}
			firstLine = false;
		}
		if(m_ft == FASTA){
			if(bufferForLine.at(0) == '>'){
				getline(m_f,bufferForLine); // skip over next line -- definitely won't be ID line
				m_size++;
			}
		}
		if(m_ft == FASTQ){
			// 1: ID, 2: sequence, 3: ID, 4: quality
			if(lineMod4 == 2){
				m_size++;
			}
		}
	}
	// reset file pointer to start of file
	m_f.close();
	m_f.open(file.c_str());
	return EXIT_SUCCESS;
}

void ColorSpaceLoader::load(int maxToLoad,ArrayOfReads*reads,MyAllocator*seqMyAllocator){
	string bufferForLine;
	string sequence("");
	string id("");
	int loadedSequences = 0;
	int lineMod4 = 0;
	bool doneComments = false; // needed as # can appear at start of fastq quality string
	while(!m_f.eof() && (m_loaded < m_size) && (loadedSequences < maxToLoad)){
		getline(m_f, bufferForLine);
		if(bufferForLine.length() == 0){
			// skip over empty lines (including at end of file)
			continue;
		}
		if(!doneComments && (bufferForLine.at(0) == '#')){
			continue;// skip over comments
		} else {
			doneComments = true;
		}
		lineMod4 = ((lineMod4 + 1) % 4);
		if(m_ft == FASTA){
			// read two lines
			if(bufferForLine.at(0) == '>'){
				if(id.compare("") != 0){
					// a previous sequence has been read in
					Read t;
					t.constructor(sequence.c_str(),seqMyAllocator,true);
					reads->push_back(&t);
					loadedSequences++;
					m_loaded++;
				}
				id = bufferForLine;
				sequence.assign("");
			} else {
				sequence += m_decoder.decode(bufferForLine);
			}
		} else if(m_ft == FASTQ){
			if(lineMod4 == 2){
				string decodedLine = m_decoder.decode(bufferForLine);
				Read t;
				t.constructor(decodedLine.c_str(),seqMyAllocator,true);
				reads->push_back(&t);
				loadedSequences++;
				m_loaded++;
			}
		} else {
			/* if file type is unknown [and m_size > 0], simulate
			 * loading to break out of while loop */
			m_loaded++;
		}
	}
	if((id.compare("") != 0) && (sequence.compare("") != 0)){
		// sequence still exists in sequence buffer, so store in reads
		Read t;
		t.constructor(sequence.c_str(),seqMyAllocator,true);
		reads->push_back(&t);
		loadedSequences++;
		m_loaded++;
	}
	if(m_loaded >= m_size){
		m_f.close();
	}
}

int ColorSpaceLoader::getSize(){
	return m_size;
}
