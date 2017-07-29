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


#include "StoreKeeper.h"
#include "CoalescenceManager.h"
#include "MatrixOwner.h"
#include "KmerMatrixOwner.h"

#include <code/VerticesExtractor/Vertex.h>
#include <RayPlatform/structures/MyHashTableIterator.h>
#include <RayPlatform/core/OperatingSystem.h>

#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <bitset>

using namespace std;

#include <string.h>

#include <assert.h>


#define INPUT_TYPE_GRAPH 0
#define INPUT_FILTERIN_GRAPH 1
#define INPUT_FILTEROUT_GRAPH 2
#define INPUT_TYPE_ASSEMBLY 3
#define INPUT_FILTERIN_ASSEMBLY 4
#define INPUT_FILTEROUT_ASSEMBLY 5


StoreKeeper::StoreKeeper() {

	m_storeDataCalls = 0;

	m_receivedObjects = 0;

	m_configured = false;
	m_kmerLength = 0;

	m_receivedPushes = 0;
}

StoreKeeper::~StoreKeeper() {

	m_receivedPushes = 0;
}

void StoreKeeper::receive(Message & message) {

	int tag = message.getTag();
	int source = message.getSourceActor();
	char * buffer = message.getBufferBytes();

	if(!m_configured)
		configureHashTable();

	if(tag == PUSH_SAMPLE_VERTEX) {

		m_receivedPushes ++;

		pushSampleVertex(message);

	} else if( tag == CoalescenceManager::DIE) {

		printName();
		cout << "[StoreKeeper] received " << m_receivedObjects << " objects in total";
		cout << " with " << m_receivedPushes << " push operations" << endl;


		// * 2 because we store pairs
		uint64_t size = m_hashTable.size() * 2;

		printName();
		cout << "[StoreKeeper] has " << size << " Kmer objects in MyHashTable instance (final)" << endl;


		printName();
		cout << "[StoreKeeper] will now die " << endl;

		die();

	} else if(tag == MERGE_GRAM_MATRIX) {

#if 0
		printName();
		cout << "DEBUG at MERGE_GRAM_MATRIX message reception ";
		cout << "[StoreKeeper] received " << m_receivedObjects << " objects in total";
		cout << " with " << m_receivedPushes << " push operations" << endl;
#endif
		computeLocalGramMatrix();

		m_mother = source;

		memcpy(&m_matrixOwner, buffer, sizeof(m_matrixOwner));

		// m_iterator1 = m_localGramMatrix.begin();

		// if(m_iterator1 != m_localGramMatrix.end()) {
		// 	m_iterator2 = m_iterator1->second.begin();
		// }

		m_iterator0 = m_localGramMatrices.begin();
		m_iterator1 = m_localGramMatrices[m_iterator0->first].begin();

		if(m_iterator1 != m_localGramMatrices[m_iterator0->first].end()) {
			m_iterator2 = m_iterator1->second.begin();
		}

		sendMatrixCell();

	} else if(tag == MatrixOwner::PUSH_PAYLOAD_OK) {
		sendMatrixCell();
	} else if(tag == MERGE_KMER_MATRIX) {

		m_mother = source;

		memcpy(&m_kmerMatrixOwner, buffer, sizeof(m_kmerMatrixOwner));

		m_hashTableIterator.constructor(&m_hashTable);

		sendKmersSamples();
	} else if (tag == KmerMatrixOwner::PUSH_KMER_SAMPLES_END) {
		// empty
	} else if(tag == KmerMatrixOwner::PUSH_KMER_SAMPLES_OK) {
		sendKmersSamples();
	} else if(tag == CoalescenceManager::SET_KMER_INFO) {

		int kmerLength = 0;
		int position = 0;

		char * buffer = (char*)message.getBufferBytes();
		memcpy(&kmerLength, buffer + position, sizeof(kmerLength));
		position += sizeof(kmerLength);
		memcpy(&m_sampleInputTypes, buffer + position, sizeof(m_sampleInputTypes));
		position += sizeof(m_sampleInputTypes);
		memcpy(&m_filterMatrices, buffer + position, sizeof(m_filterMatrices));
		position += sizeof(m_filterMatrices);

		if(m_kmerLength == 0)
			m_kmerLength = kmerLength;

		if(kmerLength != m_kmerLength) {
			cout << "[StoreKeeper] ERROR: the k-mer value is different this time !" << endl;
		}

		// the color space mode is an artefact.
		m_colorSpaceMode = false;

#if 0
		std::ostringstream ss;
		for(std::vector<int>::iterator it = m_sampleInputTypes->begin(); it != m_sampleInputTypes->end(); ++it) {
			ss << *it;
		}

		cout << "DEBUG StoreKeeper SET_KMER_INFO ";
		cout << m_kmerLength << endl;

		cout << "DEBUG StoreKeeper list of filters :";
		cout << ss.str();
		cout << endl;

#endif

	}
}

void StoreKeeper::sendMatrixCell() {

	if (m_iterator0 != m_localGramMatrices.end()) {

		if(m_iterator1 != m_localGramMatrices[m_iterator0->first].end()) {

			if(m_iterator2 != m_iterator1->second.end()) {

				SampleIdentifier sample1 = m_iterator1->first;
				SampleIdentifier sample2 = m_iterator2->first;
				LargeCount count = m_iterator2->second;
				int matrixNb = m_iterator0->first;

				Message message;
				char buffer[20];
				int offset = 0;
				memcpy(buffer + offset, &sample1, sizeof(sample1));
				offset += sizeof(sample1);
				memcpy(buffer + offset, &sample2, sizeof(sample2));
				offset += sizeof(sample2);
				memcpy(buffer + offset, &count, sizeof(count));
				offset += sizeof(count);
				memcpy(buffer + offset, &matrixNb, sizeof(matrixNb));
				offset += sizeof(matrixNb);

				message.setBuffer(buffer);
				message.setNumberOfBytes(offset);
				message.setTag(MatrixOwner::PUSH_PAYLOAD);

				send(m_matrixOwner, message);

				m_iterator2++;

				// end of the line
				if(m_iterator2 == m_iterator1->second.end()) {

					m_iterator1++;

					if(m_iterator1 != m_localGramMatrices[m_iterator0->first].end()) {
						m_iterator2 = m_iterator1->second.begin();
					} else {
						m_iterator0++;
						if (m_iterator0 != m_localGramMatrices.end()) {
							m_iterator1 = m_localGramMatrices[m_iterator0->first].begin();
							m_iterator2 = m_iterator1->second.begin();
						}
					}
				}

				return;
			}
		}

	}

	// we processed all the matrix

	// free memory.
	m_localGramMatrices.clear();

	printName();
	cout << "[StoreKeeper] local Gram Matrix clearance!" << endl;

	Message response;
	response.setTag(MatrixOwner::PUSH_PAYLOAD_END);
	send(m_matrixOwner, response);
}

void StoreKeeper::configureHashTable() {

	uint64_t buckets = 268435456;

	int bucketsPerGroup = 32 + 16 + 8 + 8;

	// \see http://docs.oracle.com/javase/7/docs/api/java/util/HashMap.html
	double loadFactorThreshold = 0.75;

	int rank = getRank();

	bool showMemoryAllocation = false;

	m_hashTable.constructor(buckets,"/apps/Ray-Surveyor/actors/StoreKeeper.txt",
		showMemoryAllocation, rank,
		bucketsPerGroup,loadFactorThreshold
		);

	m_configured = true;

}

void StoreKeeper::printColorReport() {

	printName();
	cout << "Coloring Report: " << endl;

	m_colorSet.printColors(&cout);
}

void StoreKeeper::computeLocalGramMatrix() {

#if 0
	uint64_t sum = 0;
#endif

	// compute the local Gram matrix

	int colors = m_colorSet.getTotalNumberOfVirtualColors();

#if 0
	cout << "DEBUG " << colors << " virtual colors" << endl;
	cout << "DEBUG312 m_storeDataCalls " << m_storeDataCalls << endl;
#endif

#ifdef CONFIG_ASSERT
	int withZeroReferences = 0;
#endif

	for(int i = 0 ; i < colors ; ++i) {

		VirtualKmerColorHandle virtualColor = i;

		set<PhysicalKmerColor> * samples = m_colorSet.getPhysicalColors(virtualColor);

		LargeCount hits = m_colorSet.getNumberOfReferences(virtualColor);

#ifdef CONFIG_ASSERT
		if(hits == 0) {
			withZeroReferences ++;
		}
#endif

		// TODO for "Directed Surveys",
		// add a check for colors in the virtualColor that are not in a namespace.
		// This directed survey only aims at counting colored kmers with colors
		// other than sample colors

		// since people are going to use this to check
		// for genome size, don't duplicate counts
		bool reportTwoDNAStrands = false;

#if 0
		cout << "DEBUG ***********";
		cout << "virtualColor: " << i << " ";
		cout << " samples: " << samples->size() << endl;
		cout << "  Sample list:";
#endif

#if 0
		for(set<PhysicalKmerColor>:: iterator sampleIterator = samples->begin();
				sampleIterator != samples->end() ;
				++sampleIterator) {

			PhysicalKmerColor value = *sampleIterator;

			cout << " " << value;
		}
		cout << endl;

		cout << " References: " << hits << " hash table entries ";

#ifdef CONFIG_ASSERT
		cout << " DEBUG.WithZeroReferences ---> " << withZeroReferences << endl;
#endif

		cout << endl;

#endif


		// bool filterOutKmer = checkKmerFilter(samples);
		bool filterOutKmer = false;


		// we have 2 DNA strands !!!
		if(reportTwoDNAStrands)
			hits *= 2;

		// DONE: samples could be ordered and stored in a vector
		// to stop as soon as j > i

		// Complexity: quadratic in the number of samples -> physical colors of the virtual color.
		// .. Now N*LogN in the number of samples VS quadratic
		for(set<PhysicalKmerColor>::iterator sample1 = samples->begin();
				sample1 != samples->end();
				++sample1) {

			SampleIdentifier sample1Index = *sample1;

			for(set<PhysicalKmerColor>::iterator sample2 = sample1;
				sample2 != samples->end();
				++sample2) {

				SampleIdentifier sample2Index = *sample2;

				for(map< int, vector<int> >::iterator filter = m_filterMatrices->begin();
				    filter != m_filterMatrices->end(); ++filter) {

					filterOutKmer = checkKmerFilter(samples, &filter->second);

					if(filterOutKmer) {
						continue;
					}

					m_localGramMatrices[filter->first][sample1Index][sample2Index] += hits;
				}

#if 0
				sum += hits;
#endif

			}
		}
	}

#if 0
	printName();
	cout << " DEBUG checksum " << sum << endl;

	uint64_t size = m_hashTable.size();
	cout << " DEBUG m_hashTable.size() " << size << endl;
#endif

}


void StoreKeeper::printLocalGramMatrix() {

	printName();
	cout << "[StoreKeeper] Local Gram Matrix:" << endl;
	cout << endl;

	for(map<SampleIdentifier, map<SampleIdentifier, LargeCount> >::iterator column = m_localGramMatrices[-1].begin();
			column != m_localGramMatrices[-1].end(); ++column) {

		SampleIdentifier sample = column->first;

		cout << "	" << sample;
	}

	cout << endl;

	for(map<SampleIdentifier, map<SampleIdentifier, LargeCount> >::iterator row = m_localGramMatrices[-1].begin();
			row != m_localGramMatrices[-1].end(); ++row) {

		SampleIdentifier sample1 = row->first;

		cout << sample1;
		for(map<SampleIdentifier, LargeCount>::iterator cell = row->second.begin();
				cell != row->second.end(); ++cell) {

			//SampleIdentifier sample2 = cell->first;

			LargeCount hits = cell->second;

			cout << "	" << hits;
		}

		cout << endl;
	}
}

void StoreKeeper::pushSampleVertex(Message & message) {

	char * buffer = (char*)message.getBufferBytes();
	int bytes = message.getNumberOfBytes();

	int position = 0;

	int producer = -1;
	bytes -= sizeof(producer);
	memcpy(&producer, buffer + bytes, sizeof(producer));

	while(position < bytes) {
		Vertex vertex;

		position += vertex.load(buffer + position);

		int sample = -1;
		memcpy(&sample, buffer + position, sizeof(sample));
		position += sizeof(sample);

		storeData(vertex, sample);

		m_receivedObjects ++;

		if(m_receivedObjects % 1000000 == 0) {

			printStatus();
		}
	}

	int source = message.getSourceActor();

	Message response;
	response.setTag(PUSH_SAMPLE_VERTEX_OK);
	response.setBuffer(&producer);
	response.setNumberOfBytes(sizeof(producer));

#ifdef CONFIG_ASSERT
	assert(sizeof(producer) > 0);
#endif

	send(source, response);
}

void StoreKeeper::printStatus() {

	printName();
	cout << "[StoreKeeper] received " << m_receivedObjects << " objects so far !" << endl;
}

void StoreKeeper::storeData(Vertex & vertex, int & sample) {

	m_storeDataCalls++;

	Kmer kmer = vertex.getKey();
	Kmer lowerKey;
	kmer.getLowerKey(&lowerKey, m_kmerLength, m_colorSpaceMode);

	uint64_t before = m_hashTable.size();

	ExperimentVertex * graphVertex = m_hashTable.insert(&lowerKey);

	// * 2 because we store pairs
	uint64_t size = m_hashTable.size();

	// check if we inserted something.
	// if it is the case, then assign the dummy color to it.
	if(before < size) {

		set<PhysicalKmerColor> emptySet;
		VirtualKmerColorHandle noColor = m_colorSet.findVirtualColor(&emptySet);

		m_colorSet.incrementReferences(noColor);

		graphVertex->setVirtualColor(noColor);

#ifdef CONFIG_ASSERT
		assert(noColor == NULL_VIRTUAL_COLOR);
#endif
	}

	int period = 1000000;
	if(size % period == 0 && size != m_lastSize) {

		printName();
		cout << " has " << size << " Kmer objects in MyHashTable instance" << endl;

		m_lastSize = size;
	}

#if 0
	cout << "DEBUG Growth -> " << before << " -> " << size << endl;
#endif

	// add the PhysicalKmerColor to the node.

	PhysicalKmerColor sampleColor = sample;
	VirtualKmerColorHandle oldVirtualColor = graphVertex->getVirtualColor();

	if(m_colorSet.virtualColorHasPhysicalColor(oldVirtualColor, sampleColor)) {

		// Nothing to do, we already have this color
		return;
	}

#ifdef CONFIG_ASSERT

	set<PhysicalKmerColor> * theOldSamples = m_colorSet.getPhysicalColors(oldVirtualColor);
	set<PhysicalKmerColor> oldSamples = *theOldSamples;

	assert(oldSamples.count(sampleColor) == 0);
#endif

	VirtualKmerColorHandle newVirtualColor= m_colorSet.getVirtualColorFrom(oldVirtualColor, sampleColor);

#ifdef CONFIG_ASSERT
	assert(m_colorSet.virtualColorHasPhysicalColor(newVirtualColor, sampleColor));
	set<PhysicalKmerColor>* samples = m_colorSet.getPhysicalColors(newVirtualColor);


	assert(samples->count(sampleColor) > 0);
#endif


#ifdef CONFIG_ASSERT2
	if(oldVirtualColor == newVirtualColor) {

		cout << " new sampleColor " << sampleColor << endl;
		cout << endl;

		cout << " >>> Old samples " << oldSamples.size () << endl;

		for(set<PhysicalKmerColor>::iterator i = oldSamples.begin();
				i != oldSamples.end() ; ++i) {

			cout << " " << *i;

		}

		cout << endl;

		cout << " old color " << oldVirtualColor;
		cout << " refs " << m_colorSet.getNumberOfReferences(oldVirtualColor) << endl;


		set<PhysicalKmerColor>* samples = m_colorSet.getPhysicalColors(newVirtualColor);

		cout << " >>> new samples " << samples->size () << endl;

		for(set<PhysicalKmerColor>::iterator i = samples->begin();
				i != samples->end() ; ++i) {

			cout << " " << *i;

		}

		cout << endl;

		cout << " new color " << newVirtualColor;
		cout << " refs " << m_colorSet.getNumberOfReferences(newVirtualColor) << endl;


	}

	// we can reuse the same handle if it has 0 references
	// The call to decrementReferences is done above the
	// call to getVirtualColorFrom
	//assert(oldVirtualColor != newVirtualColor);
#endif

	graphVertex->setVirtualColor(newVirtualColor);

	m_colorSet.incrementReferences(newVirtualColor);
	m_colorSet.decrementReferences(oldVirtualColor);

}


void StoreKeeper::setSampleSize(int sampleSize) {
	m_sampleSize = sampleSize;
}


void StoreKeeper::sendKmersSamples() {

	char buffer[m_sampleSize+256];
	int bytes = 0;

	ExperimentVertex * currentVertex = NULL;
	VirtualKmerColorHandle currentVirtualColor = NULL_VIRTUAL_COLOR;

	char samplesArray[m_sampleSize+1];
	samplesArray[m_sampleSize] = '\0';

	// Set all samples to 0
	for (int i=0; i<m_sampleSize; i++) {
		samplesArray[i] = '0';
	}

	if(m_hashTableIterator.hasNext()){

		currentVertex = m_hashTableIterator.next();
		Kmer kmer = currentVertex->getKey();

		bytes += kmer.dump(buffer);

		currentVirtualColor = currentVertex->getVirtualColor();
		set<PhysicalKmerColor> * samples = m_colorSet.getPhysicalColors(currentVirtualColor);

		for(set<PhysicalKmerColor>:: iterator sampleIterator = samples->begin();
			sampleIterator != samples->end(); ++sampleIterator) {
			PhysicalKmerColor value = *sampleIterator;
			samplesArray[value] = '1';
		}


		for (int i=0; i<m_sampleSize; i++) {
			buffer[bytes] = samplesArray[i];
			bytes++;
		}

	}


	Message message;
	message.setNumberOfBytes(bytes);
	message.setBuffer(buffer);

	if(m_hashTableIterator.hasNext()){
		message.setTag(KmerMatrixOwner::PUSH_KMER_SAMPLES);
	}else{
		message.setTag(KmerMatrixOwner::PUSH_KMER_SAMPLES_END);
	}

	send(m_kmerMatrixOwner, message);

}



bool StoreKeeper::checkKmerFilter (set<PhysicalKmerColor> * samples, vector<int> * sampleInputTypes) {

	bool skipKmer = false;
	vector<int> samplesFILTERIN;

	for(std::vector<int>::iterator it = sampleInputTypes->begin(); it != sampleInputTypes->end(); ++it) {

		// cout << "DEBUG [storekeeper checkfilter] sampleInputTypes : " << *it << endl;

		int sample = distance(sampleInputTypes->begin(),it);

		if(*it == INPUT_TYPE_GRAPH || *it == INPUT_TYPE_ASSEMBLY){
			continue;
		}
		else if(*it == INPUT_FILTERIN_GRAPH || *it == INPUT_FILTERIN_ASSEMBLY){
			// Store the FILTERIN for later lookup
			samplesFILTERIN.push_back(sample);
		}
		else if(*it == INPUT_FILTEROUT_GRAPH || *it == INPUT_FILTEROUT_ASSEMBLY){
			// Skip the Kmer if part of a FILTEROUT
			if (samples->count(sample) > 0){
				return true;
			}
		}
	}

	// Look if the Kmer is not part of a FILTERIN
	// Only incorporates kmers from FILTERIN samples if not already Filtered OUT
	if(!samplesFILTERIN.empty()) {
		skipKmer = true;
		for(std::vector<int>::iterator it = samplesFILTERIN.begin(); it != samplesFILTERIN.end(); ++it) {
			if (samples->count(*it) > 0){
				return false;
			}
		}
	}


	return skipKmer;

}
