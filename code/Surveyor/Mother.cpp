/*
    Copyright 2014 Maxime Déraspe, Sébastien Boisvert
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


// Written by <<sebhtml>>
// 2013-10-10

#include "Mother.h"
#include "StoreKeeper.h"
#include "GenomeGraphReader.h"
#include "GenomeAssemblyReader.h"
#include "MatrixOwner.h"
#include "KmerMatrixOwner.h"

#include <RayPlatform/cryptography/crypto.h>

#include <iostream>
using namespace std;

#define PLAN_RANK_ACTORS_PER_RANK 1
#define PLAN_MOTHER_ACTORS_PER_RANK 1
#define PLAN_GENOME_GRAPH_READER_ACTORS_PER_RANK 999999
#define INPUT_TYPE_GRAPH 0
#define INPUT_FILTERIN_GRAPH 1
#define INPUT_FILTEROUT_GRAPH 2
#define INPUT_TYPE_ASSEMBLY 3
#define INPUT_FILTERIN_ASSEMBLY 4
#define INPUT_FILTEROUT_ASSEMBLY 5

Mother::Mother() {

	m_coalescenceManager = -1;
	m_matrixOwner = -1;
	m_kmerMatrixOwner = -1;

	m_parameters = NULL;
	m_bigMother = -1;

	m_finishedMothers = 0;

	m_flushedMothers = 0;
}

Mother::~Mother() {

}

void Mother::receive(Message & message) {

	int tag = message.getTag();
	int source = message.getSourceActor();
	char * buffer = message.getBufferBytes();

	if(tag == Actor::BOOT) {
		boot(message);
	} else if (tag == Mother::HELLO) {
		hello(message);
	} else if(tag == GenomeGraphReader::START_PARTY_OK) {

		// spawn the next reader now !
		spawnReader();

	} else if(tag == GenomeGraphReader::DONE) {

		m_aliveReaders--;

		if(m_aliveReaders == 0) {
			notifyController();
		}

	} else if(tag == MERGE_GRAM_MATRIX) {

		int matrixOwner = -1;
		memcpy(&matrixOwner, buffer, sizeof(matrixOwner));

#ifdef CONFIG_ASSERT
		assert(matrixOwner >= 0);
		assert(m_storeKeepers.size() == 1);
#endif

		Message theMessage;
		theMessage.setTag(StoreKeeper::MERGE_GRAM_MATRIX);
		theMessage.setBuffer(&matrixOwner);
		theMessage.setNumberOfBytes(sizeof(matrixOwner));

		int destination = m_storeKeepers[0];

		send(destination, theMessage);

		Message response;
		response.setTag(MERGE_GRAM_MATRIX_OK);
		send(source, response);

	} else if (tag == MERGE_KMER_MATRIX) {

		int kmerMatrixOwner = -1;
		memcpy(&kmerMatrixOwner, buffer, sizeof(kmerMatrixOwner));

#ifdef CONFIG_ASSERT
		assert(kmerMatrixOwner >= 0);
		assert(m_storeKeepers.size() == 1);
#endif

		Message theMessage;
		theMessage.setTag(StoreKeeper::MERGE_KMER_MATRIX);
		theMessage.setBuffer(&kmerMatrixOwner);
		theMessage.setNumberOfBytes(sizeof(kmerMatrixOwner));

		int destination = m_storeKeepers[0];

		send(destination, theMessage);

		Message response;
		response.setTag(MERGE_KMER_MATRIX_OK);
		send(source, response);

	} else if(tag == SHUTDOWN) {

		Message response;
		response.setTag(SHUTDOWN_OK);
		send(message.getSourceActor(), response);

		stop();

	} else if(tag == StoreKeeper::MERGE_GRAM_MATRIX_OK) {

		// TODO: the bug https://github.com/sebhtml/ray/issues/216
		// is caused by the fact that this message is not
		// received .

		// does nothing

	} else if(tag == FINISH_JOB) {

		m_finishedMothers++;

		if(m_finishedMothers == getSize()) {

			// all readers have finished,
			// now tell mother to flush aggregators
			sendToFirstMother(FLUSH_AGGREGATOR, FLUSH_AGGREGATOR_RETURN);
		}

	} else if(tag == FLUSH_AGGREGATOR) {

		m_bigMother = source;

		// forward the message to the aggregator
		Message newMessage;
		newMessage.setTag(CoalescenceManager::FLUSH_BUFFERS);
		send(m_coalescenceManager, newMessage);

		Message newMessage2;
		newMessage2.setTag(FLUSH_AGGREGATOR_RETURN);
		send(source, newMessage2);

	} else if(tag == CoalescenceManager::FLUSH_BUFFERS_OK) {

		// notice bigMother that the CoalescenceManager flushed its buffers
		// .. for this rank's mother
		Message response;
		response.setTag(FLUSH_AGGREGATOR_OK);
		send(m_bigMother, response);

	} else if(tag == MatrixOwner::GRAM_MATRIX_IS_READY) {

		if(m_matricesAreReady){
			sendToFirstMother(SHUTDOWN, SHUTDOWN_OK);
		}else {
			printName();
			cout << "[BigMother] GRAM_MATRIX_IS_READY" << endl;
			m_matricesAreReady = true;
		}

	} else if(tag == KmerMatrixOwner::KMER_MATRIX_IS_READY) {

		printName();
		cout << "[BigMother] KMER_MATRIX_IS_READY" << endl;

		if(m_matricesAreReady){
			sendToFirstMother(SHUTDOWN, SHUTDOWN_OK);
		}else {
			m_matricesAreReady = true;
		}

	} else if(tag == FLUSH_AGGREGATOR_OK) {

		m_flushedMothers++;

		if(m_flushedMothers < getSize())
			return;

		spawnMatrixOwner();

	} else if(tag == m_responseTag) {

		if(m_responseTag == SHUTDOWN_OK) {
			// does nothing
		} else if(m_responseTag == MERGE_GRAM_MATRIX_OK) {

			// All mothers merged their GRAM MATRIX
			// Spawn KmerMatrixOwner to print
			if(m_motherToKill < getSize() && m_printKmerMatrix){
				spawnKmerMatrixOwner();
			}

		} else if(m_responseTag == MERGE_KMER_MATRIX_OK) {
			// does nothing
		} else if(m_responseTag == FLUSH_AGGREGATOR_RETURN) {
			// does nothing
		}

		// every mother was not informed.
		if(m_motherToKill >= getSize()) {
			sendMessageWithReply(m_motherToKill, m_forwardTag);
			m_motherToKill--;
		}
	}
}

void Mother::sendToFirstMother(int forwardTag, int responseTag) {

	m_forwardTag = forwardTag;
	m_responseTag = responseTag;

	m_motherToKill = 2 * getSize() - 1;

	sendMessageWithReply(m_motherToKill, m_forwardTag);
	m_motherToKill--;

}

void Mother::sendMessageWithReply(int & actor, int tag) {

	Message message;
	message.setTag(tag);

	if(tag == MERGE_GRAM_MATRIX) {
		message.setBuffer(&m_matrixOwner);
		message.setNumberOfBytes(sizeof(m_matrixOwner));
	} else if(tag == MERGE_KMER_MATRIX) {
		message.setBuffer(&m_kmerMatrixOwner);
		message.setNumberOfBytes(sizeof(m_kmerMatrixOwner));
	} else if(tag == FLUSH_AGGREGATOR) {
		// does nothing
	}

	send(actor, message);

}

void Mother::notifyController() {

	Message message2;
	message2.setTag(FINISH_JOB);

	// first Mother
	int controller = getSize();
	send(controller, message2);

}

void Mother::stop() {

	Message kill;
	kill.setTag(CoalescenceManager::DIE);

	send(m_coalescenceManager, kill);
	m_coalescenceManager = -1;

	send(m_storeKeepers[0], kill);
	m_storeKeepers.clear();

	if(m_matrixOwner >= 0) {
		send(m_matrixOwner, kill);
		m_matrixOwner = -1;
	}

	if(m_kmerMatrixOwner >= 0) {
		send(m_kmerMatrixOwner, kill);
		m_kmerMatrixOwner = -1;
	}

	die();

}

void Mother::hello(Message & message) {
	// dump for testing purposes
}

void Mother::boot(Message & message) {

	m_aliveReaders = 0;

	Message message2;

	int joe = 9921;

	message2.setBuffer(&joe);

	message2.setNumberOfBytes( sizeof(int) * 1 );

	message2.setTag(Mother::HELLO);

	int next = getName() + 1;

	next %= (getSize() * 2);

	if(next < getSize())
		next += getSize();

	printName();
	cout << "[RankMother] friend of [RankMother] #" << next << endl;

	send(next, message2);

	if(m_parameters->hasOption("-run-surveyor")) {
		startSurveyor();
	} else {
		die();
	}
}

void Mother::startSurveyor() {

	bool isRoot = (getName() % getSize()) == 0;

	// Set matricesAreReady to true in case user doesn't want
	// to print out kmers matrix.
	m_matricesAreReady = true;

	vector<string> * commands = m_parameters->getCommands();

	for(int i = 0 ; i < (int) commands->size() ; ++i) {

		string & element = commands->at(i);

		if (element != "-write-kmer-matrix") {
			// DONE: Check bounds for file names

			map<string,int> fastTable;

			fastTable["-read-sample-graph"] = INPUT_TYPE_GRAPH;
			fastTable["-filter-in-graph"] = INPUT_FILTERIN_GRAPH;
			fastTable["-filter-out-graph"] = INPUT_FILTEROUT_GRAPH;
			fastTable["-read-sample-assembly"] = INPUT_TYPE_ASSEMBLY;
			fastTable["-filter-in-assembly"] = INPUT_FILTERIN_ASSEMBLY;
			fastTable["-filter-out-assembly"] = INPUT_FILTEROUT_ASSEMBLY;

			// Unsupported option
			if(fastTable.count(element) == 0 || i+2 > (int) commands->size())
				continue;

			string sampleName = commands->at(++i);
			string fileName = commands->at(++i);

			m_sampleNames.push_back(sampleName);

			// DONE implement this m_assemblyFileNames + type
			m_inputFileNames.push_back(fileName);

			int type = fastTable[element];

			m_sampleInputTypes.push_back(type);

		} else {
			m_matricesAreReady = false;
			m_printKmerMatrix = true;
		}

	}

	if(isRoot) {
		printName();
		cout << "[BigMother] I am the Survey's Goddess and I am watching you!" << endl;
		printName();
		cout << "[BigMother] Total number of samples to compare: " << m_sampleNames.size() << endl;
	}


	CoalescenceManager * coalescenceManager = new CoalescenceManager();
	spawn(coalescenceManager);

	m_coalescenceManager = coalescenceManager->getName();

	// spawn the local StoreKeeper and introduce the CoalescenceManager to him

	// spawn actors for storing the graph.
	for(int i = 0 ; i < PLAN_STORE_KEEPER_ACTORS_PER_RANK; ++i) {

		StoreKeeper * actor = new StoreKeeper();
		spawn(actor);

		m_storeKeepers.push_back(actor->getName());

		actor->setSampleSize(m_sampleNames.size());

		// tell the CoalescenceManager about the local StoreKeeper
		Message dummyMessage;
		int localStore = actor->getName();

		dummyMessage.setBuffer(&localStore);
		dummyMessage.setNumberOfBytes( sizeof(int) );

		dummyMessage.setTag(CoalescenceManager::INTRODUCE_STORE);

		send(m_coalescenceManager, dummyMessage);

		// send kmerLength and sampleInputTypes to localStore
		int kmerLength = m_parameters->getWordSize();
		vector<int> * sampleInputTypes = & m_sampleInputTypes;

		char buffer[32];
		int offset = 0;
		memcpy(buffer + offset, &kmerLength, sizeof(kmerLength));
		offset += sizeof(kmerLength);
		memcpy(buffer + offset, &sampleInputTypes, sizeof(sampleInputTypes));
		offset += sizeof(sampleInputTypes);

		Message kmerMessage;
		kmerMessage.setBuffer(&kmerLength);
		kmerMessage.setNumberOfBytes(sizeof(kmerLength));
		kmerMessage.setBuffer(&buffer);
		kmerMessage.setNumberOfBytes(sizeof(buffer));
		kmerMessage.setTag(CoalescenceManager::SET_KMER_INFO);
		send(localStore, kmerMessage);
	}


	// spawn an actor for each file that this actor owns
	for(int i = 0; i < (int) m_inputFileNames.size() ; ++i) {

		int mother = getName() % getSize();
		int fileMother = i % getSize();

		if(mother == fileMother) {
			m_filesToSpawn.push_back(i);
		}
	}

	printName();
	cout << "[RankMother] readers to spawn: " << m_filesToSpawn.size() << endl;

	m_fileIterator = 0;
	spawnReader();
}

void Mother::spawnReader() {

	if(m_fileIterator < (int) m_filesToSpawn.size()) {

		int sampleIdentifier = m_filesToSpawn[m_fileIterator];

		string & fileName = m_inputFileNames[sampleIdentifier];
		int type = m_sampleInputTypes[sampleIdentifier];
		m_fileIterator++;

		if(type == INPUT_TYPE_GRAPH || type == INPUT_FILTERIN_GRAPH || type == INPUT_FILTEROUT_GRAPH) {
			GenomeGraphReader * actor = new GenomeGraphReader();

			spawn(actor);
			actor->setFileName(fileName, sampleIdentifier);

			int coalescenceManagerName = m_coalescenceManager;
			int destination = actor->getName();
			Message dummyMessage;

			dummyMessage.setBuffer(&coalescenceManagerName);
			dummyMessage.setNumberOfBytes( sizeof(int) );
			dummyMessage.setTag(GenomeGraphReader::START_PARTY);

			m_aliveReaders ++;

			send(destination, dummyMessage);


		} else if(type == INPUT_TYPE_ASSEMBLY || type == INPUT_FILTERIN_ASSEMBLY || type == INPUT_FILTEROUT_ASSEMBLY){

			GenomeAssemblyReader * actor = new GenomeAssemblyReader();
			spawn(actor);
			actor->setFileName(fileName, sampleIdentifier);
			actor->setKmerSize(m_parameters->getWordSize());

			int coalescenceManagerName = m_coalescenceManager;
			int destination = actor->getName();
			Message dummyMessage;

			dummyMessage.setBuffer(&coalescenceManagerName);
			dummyMessage.setNumberOfBytes( sizeof(int) );

			dummyMessage.setTag(GenomeAssemblyReader::START_PARTY);

			m_aliveReaders ++;

			send(destination, dummyMessage);

		}
	}

	if(m_aliveReaders == 0) {
		notifyController();
	}
}


void Mother::spawnMatrixOwner() {

	// spawn the MatrixOwner here !
	MatrixOwner * matrixOwner = new MatrixOwner();
	spawn(matrixOwner);

	m_matrixOwner = matrixOwner->getName();

	printName();
	cout << "[BigMother] Spawned [MatrixOwner] actor #" << m_matrixOwner << endl;

	// tell the StoreKeeper actors to send their stuff to the MatrixOwner actor
	// bigMother will wait for a signal from MatrixOwner

	Message greetingMessage;

	vector<string> * names = & m_sampleNames;

	char buffer[32];
	int offset = 0;
	memcpy(buffer + offset, &m_parameters, sizeof(m_parameters));
	offset += sizeof(m_parameters);
	memcpy(buffer + offset, &names, sizeof(names));
	offset += sizeof(names);

	greetingMessage.setBuffer(&buffer);
	greetingMessage.setNumberOfBytes(offset);

	greetingMessage.setTag(MatrixOwner::GREETINGS);
	send(m_matrixOwner, greetingMessage);

	sendToFirstMother(MERGE_GRAM_MATRIX, MERGE_GRAM_MATRIX_OK);
}

void Mother::spawnKmerMatrixOwner() {

	// spawn the KmerMatrixOwner here !
	KmerMatrixOwner * kmerMatrixOwner = new KmerMatrixOwner();
	spawn(kmerMatrixOwner);

	m_kmerMatrixOwner = kmerMatrixOwner->getName();

	printName();
	cout << "[BigMother] Spawned [KmerMatrixOwner] actor #" << m_kmerMatrixOwner << endl;

	// tell the StoreKeeper actors to send their stuff to the KmerMatrixOwner actor
	// bigMother will wait for a signal from KmerMatrixOwner

	Message greetingMessage;

	vector<string> * names = & m_sampleNames;

	char buffer[32];
	int offset = 0;
	memcpy(buffer + offset, &m_parameters, sizeof(m_parameters));
	offset += sizeof(m_parameters);
	memcpy(buffer + offset, &names, sizeof(names));
	offset += sizeof(names);

	greetingMessage.setBuffer(&buffer);
	greetingMessage.setNumberOfBytes(offset);

	greetingMessage.setTag(KmerMatrixOwner::GREETINGS);
	send(m_kmerMatrixOwner, greetingMessage);

	sendToFirstMother(MERGE_KMER_MATRIX, MERGE_KMER_MATRIX_OK);

}


void Mother::setParameters(Parameters * parameters) {
	m_parameters = parameters;
}
