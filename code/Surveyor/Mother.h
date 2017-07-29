/*
    Copyright 2013, 2014 Sébastien Boisvert, Maxime Déraspe
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

#ifndef MotherHeader
#define MotherHeader

#include "CoalescenceManager.h"

#include <code/Mock/Parameters.h>
#include <RayPlatform/actors/Actor.h>

#include <vector>
#include <string>
#include <iostream>
using namespace std;

/**
 *
 * Map (assuming N ranks)
 *
 * The following map is not used anymore because it is stupid.
 * ----------------------------------------------------------
 * Actors			Quantity	Role
 * Start	End
 * ----------------------------------------------------------
 * 0		N - 1		N		ComputeCore
 * N		2N - 1		N		Mother
 * 2N		102N - 1	100N		StoreKeeper
 * 102N		104N - 1	2N		GenomeGraphReader
 * ----------------------------------------------------------
 *
 * It would be nice to have a number of tokens / second that can be exchanged and also
 * the point-to-point latency for this actor implementation.
 *
 * \author Sébastien Boisvert
 */
class Mother: public Actor {

private:

	int m_matrixOwner;
	int m_kmerMatrixOwner;

	int m_flushedMothers;
	int m_finishedMothers;
	bool m_matricesAreReady;
	bool m_printKmerMatrix;

	Parameters * m_parameters;

	int m_coalescenceManager;
	int m_fileIterator;
	vector<int> m_filesToSpawn;

	vector<int> m_storeKeepers;

	vector<int> m_readers;

	vector<string> m_sampleNames;
	vector<string> m_inputFileNames;

	vector<int> m_sampleInputTypes;

	map< int, vector<int> > m_filterMatrices;

	int m_bigMother;
	int m_aliveReaders;
	int m_motherToKill;

	int m_forwardTag;
	int m_responseTag;


	// Mother methods
	void startSurveyor();
	void hello(Message & message);
	void boot(Message & message);
	void stop();
	void notifyController();
	void sendMessageWithReply(int & actor, int tag);

	/*
	 * Send tag to all mothers.
	 */
	void sendToFirstMother(int forwardTag, int responseTag);

	// Actor spawning
	void spawnReader();
	void spawnMatrixOwner();
	void spawnKmerMatrixOwner();

public:

	Mother();
	~Mother();

	void receive(Message & message);

	enum {
		FIRST_TAG = 10100,
		HELLO,
		FINISH_JOB,
		SHUTDOWN,
		SHUTDOWN_OK,
		FLUSH_AGGREGATOR,
		FLUSH_AGGREGATOR_OK,
		FLUSH_AGGREGATOR_RETURN,
		MERGE_GRAM_MATRIX,
		MERGE_GRAM_MATRIX_OK,
		MERGE_KMER_MATRIX,
		MERGE_KMER_MATRIX_OK,
		LAST_TAG,
	};

	void setParameters(Parameters * parameters);
};

#endif
