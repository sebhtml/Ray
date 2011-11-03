/*
 	Ray
    Copyright (C) 2010, 2011  Sébastien Boisvert

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

#ifndef _SequencesIndexer
#define _SequencesIndexer

#include <assembler/IndexerWorker.h>
#include <communication/Message.h>
#include <communication/BufferedData.h>
#include <communication/VirtualCommunicator.h>
#include <core/Parameters.h>
#include <core/common_functions.h>
#include <memory/MyAllocator.h>
#include <memory/RingAllocator.h>
#include <structures/SplayTree.h>
#include <structures/SplayTreeIterator.h>
#include <structures/StaticVector.h>
#include <structures/ArrayOfReads.h>
#include <map>
#include <structures/Read.h>
#include <vector>
#include <fstream>
using namespace std;

/*
 * Computes optimal read markers using workers.
 * \author Sébastien Boisvert
 */
class SequencesIndexer{

	ofstream m_readMarkerFile;

	map<int,map<int,int> >m_forwardStatistics;
	map<int,map<int,int> >m_reverseStatistics;

	/** for checkpointing */
	bool m_checkedCheckpoint;

	int m_rank;
	int m_size;
	Parameters*m_parameters;
	MyAllocator m_allocator;
	MyAllocator m_workAllocator;
	int m_pendingMessages;
	int m_completedJobs;
	int m_maximumAliveWorkers;
	int m_maximumWorkers;

	VirtualCommunicator*m_virtualCommunicator;
	SplayTree<uint64_t,char> m_activeWorkers;
	SplayTreeIterator<uint64_t,char> m_activeWorkerIterator;

	SplayTree<uint64_t,IndexerWorker> m_aliveWorkers;

	bool m_communicatorWasTriggered;
	vector<uint64_t> m_workersDone;
	vector<uint64_t> m_waitingWorkers;
	vector<uint64_t> m_activeWorkersToRestore;

	bool m_initiatedIterator;
	int m_theSequenceId;
	void updateStates();

public:

	void attachReads(
ArrayOfReads*m_myReads,
				RingAllocator*m_outboxAllocator,
				StaticVector*m_outbox,
				int*m_mode,
				int m_wordSize,
				int m_size,
				int m_rank
);

	void constructor(Parameters*parameters,RingAllocator*outboxAllocator,StaticVector*inbox,StaticVector*outbox,
	VirtualCommunicator*vc);
	void setReadiness();
	MyAllocator*getAllocator();
};

#endif
