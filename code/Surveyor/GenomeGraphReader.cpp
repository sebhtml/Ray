/*
    Copyright 2013 Sébastien Boisvert
    Copyright 2013 Université Laval
    Copyright 2013 Centre Hospitalier Universitaire de Québec

    This file is part of Ray Surveyor.

    Ray Surveyor is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    Ray Surveyor is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Ray Surveyor.  If not, see <http://www.gnu.org/licenses/>.
*/

// DONE: validate that the kmer length is the same for this file
// and the provided -k argument

#include "GenomeGraphReader.h"
#include "CoalescenceManager.h"

#include <code/Mock/constants.h>
#include <code/Mock/common_functions.h>
#include <code/KmerAcademyBuilder/Kmer.h>
#include <code/VerticesExtractor/Vertex.h>

#include <iostream>
#include <sstream>
using namespace std;

#include <string.h>

GenomeGraphReader::GenomeGraphReader() {

}

GenomeGraphReader::~GenomeGraphReader() {

}

void GenomeGraphReader::receive(Message & message) {

	int type = message.getTag();

	if(type == START_PARTY) {

		startParty(message);

	} else if(type == CoalescenceManager::PAYLOAD_RESPONSE) {

		// read the next line now !
		readLine();
	}
}

void GenomeGraphReader::startParty(Message & message) {

	char * buffer = (char*) message.getBufferBytes();

	memcpy(&m_aggregator, buffer, sizeof(int));

	m_reader.open(m_fileName.c_str());

	m_bad = false;

	if(!m_reader.isValid())
		m_bad = true;

	m_loaded = 0;

	printName();
	cout << "[GraphReader] opens file " << m_fileName << endl;

	m_parent = message.getSourceActor();

	int source = message.getSourceActor();
	Message response;
	response.setTag(START_PARTY_OK);

	send(source, response);

	readLine();
}

void GenomeGraphReader::readLine() {

	char buffer[1024];
	buffer[0] = '\0';

	bool isCurrentlyAtTheGreatEndOfTime = m_reader.eof();

	while(!m_bad && !m_reader.eof()) {
		m_reader.getline(buffer, 1024);

		// skip comment
		if(strlen(buffer) > 0 && buffer[0] == '#')
			continue;

		break;
	}

	if(m_bad || (m_reader.eof() && isCurrentlyAtTheGreatEndOfTime)) {

		m_reader.close();

		printName();

		if(m_bad) {
			cout << "[GraphReader] Error: file " << m_fileName << " does not exist";
			cout << endl;

		} else {
			cout << "[GraphReader] finished reading file " << m_fileName;
			cout << " got " << m_loaded << " objects" << endl;
		}

		Message finishedMessage;
		finishedMessage.setTag(DONE);

		send(m_parent, finishedMessage);

		die();
	} else {

		// AGCTGTGAAACTGGTGCAAGCTACCAGAATC;36;A;C
		string sequence;
		CoverageDepth coverage;
		string parents;
		string children;

		for(int i = 0 ; i < (int) strlen(buffer) ; ++i) {
			if(buffer[i] == ';')
				buffer[i] = ' ';
		}

		istringstream stringBuffer(buffer);

		stringBuffer >> sequence;
		stringBuffer >> coverage;
		stringBuffer >> parents;
		stringBuffer >> children;

		///////////////////////////////////////////////////////////////////////
		// convert the sequence to upper case

		map<char,char> translationTable;
		translationTable['a'] = 'A';
		translationTable['t'] = 'T';
		translationTable['g'] = 'G';
		translationTable['c'] = 'C';

		for(int i = 0 ; i < (int) sequence.length() ; ++i) {

			char symbol = sequence[i];

			if(translationTable.count(symbol) > 0) {
				char correct = translationTable[symbol];

				sequence [i] = correct;
			}
		}
#if 0
		cout << "DEBUG " << sequence << " with " << coverage << endl;
#endif

		// if this is the first one, send the k-mer length too
		if(m_loaded == 0) {

			Message aMessage;
			aMessage.setTag(CoalescenceManager::SET_KMER_INFO);

			int length = sequence.length();
			aMessage.setBuffer(&length);
			aMessage.setNumberOfBytes(sizeof(length));

			send(m_aggregator, aMessage);
		}

		Kmer kmer;
		kmer.loadFromTextRepresentation(sequence.c_str());

		Vertex vertex;
		vertex.setKey(kmer);
		vertex.setCoverageValue(coverage);

		// add parents
		for(int i = 0 ; i < (int)parents.length() ; ++i) {

			string parent = sequence;
			for(int j = 0 ; j < (int) parent.length()-1 ; ++j) {
				parent[j + 1] = parent[j];
			}
			parent[0] = parents[i];

			Kmer parentKmer;
			parentKmer.loadFromTextRepresentation(parent.c_str());

			vertex.addIngoingEdge(&kmer, &parentKmer, sequence.length());
		}

		// add children
		for(int i = 0 ; i < (int)children.length() ; ++i) {

			string child = sequence;
			for(int j = 0 ; j < (int) child.length()-1 ; ++j) {
				child[j] = child[j + 1];
			}
			child[child.length() - 1] = children[i];

			Kmer childKmer;
			childKmer.loadFromTextRepresentation(child.c_str());

			vertex.addOutgoingEdge(&kmer, &childKmer, sequence.length());
		}

		char messageBuffer[100];
		int position = 0;

		position += vertex.dump(messageBuffer + position);
		memcpy(messageBuffer + position, &m_sample, sizeof(m_sample));

		position += sizeof(m_sample);

		// maybe: accumulate many objects before flushing it.
		// we can go up to MAXIMUM_MESSAGE_SIZE_IN_BYTES bytes.

		// Sending PAYLOAD to the CoalescenceManager
		Message message;
		message.setTag(CoalescenceManager::PAYLOAD);
		message.setBuffer(messageBuffer);
		message.setNumberOfBytes(position);

#if 0
		printName();
		cout << "DEBUG sending PAYLOAD to " << m_aggregator;
		cout << " with " << position << " bytes ";
		vertex.print(sequence.length(), false);
		cout << endl;
#endif

		int period = 1000000;
		if(m_loaded % period == 0 && m_loaded > 0) {
			printName();
			cout << "[GraphReader] loaded " << m_loaded << " sequences" << endl;
		}
		m_loaded ++;
		send(m_aggregator, message);
	}
}

void GenomeGraphReader::setFileName(string & fileName, int sample) {

	//int nameSpace = COLOR_NAMESPACE_SAMPLE;
	//m_sample = (uint64_t)sample + nameSpace * COLOR_NAMESPACE_MULTIPLIER;

	m_sample = sample;

	m_fileName = fileName;

#if 0
	printName();
	cout << " DEBUG setFileName " << m_fileName << endl;
#endif
}
