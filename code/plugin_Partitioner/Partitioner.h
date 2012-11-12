/*
 	Ray
    Copyright (C) 2011, 2012 Sébastien Boisvert

	http://DeNovoAssembler.SourceForge.Net/

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, version 3 of the License.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You have received a copy of the GNU General Public License
    along with this program (gpl-3.0.txt).  
	see <http://www.gnu.org/licenses/>

*/

#ifndef _Partitioner_H
#define _Partitioner_H

#include <code/plugin_Mock/Parameters.h>
#include <code/plugin_SequencesLoader/Loader.h>

#include <RayPlatform/structures/StaticVector.h>
#include <RayPlatform/memory/RingAllocator.h>
#include <RayPlatform/scheduling/SwitchMan.h>
#include <RayPlatform/handlers/SlaveModeHandler.h>
#include <RayPlatform/handlers/MasterModeHandler.h>
#include <RayPlatform/core/ComputeCore.h>

#include <map>
using namespace std;

__DeclarePlugin(Partitioner);

__DeclareMasterModeAdapter(Partitioner,RAY_MASTER_MODE_COUNT_FILE_ENTRIES);
__DeclareSlaveModeAdapter(Partitioner,RAY_SLAVE_MODE_COUNT_FILE_ENTRIES);

/**
 * This class counts the number of entries in each input file in parallel 
 * \author Sébastien Boisvert
 */
class Partitioner :  public CorePlugin{

	__AddAdapter(Partitioner,RAY_MASTER_MODE_COUNT_FILE_ENTRIES);
	__AddAdapter(Partitioner,RAY_SLAVE_MODE_COUNT_FILE_ENTRIES);

	MessageTag RAY_MPI_TAG_COUNT_FILE_ENTRIES;
	MessageTag RAY_MPI_TAG_COUNT_FILE_ENTRIES_REPLY;
	MessageTag RAY_MPI_TAG_FILE_ENTRY_COUNT;
	MessageTag RAY_MPI_TAG_FILE_ENTRY_COUNT_REPLY;
	MessageTag RAY_MPI_TAG_REQUEST_FILE_ENTRY_COUNTS;
	MessageTag RAY_MPI_TAG_REQUEST_FILE_ENTRY_COUNTS_REPLY;

	MasterMode RAY_MASTER_MODE_COUNT_FILE_ENTRIES;
	MasterMode RAY_MASTER_MODE_LOAD_SEQUENCES;

	SlaveMode RAY_SLAVE_MODE_COUNT_FILE_ENTRIES;
	SlaveMode RAY_SLAVE_MODE_DO_NOTHING;



	SwitchMan*m_switchMan;

	/** the loader */
	Loader m_loader;
	/** did we initialize the master peer */
	bool m_initiatedMaster;
	/** did we initialize the slave peer */
	bool m_initiatedSlave;
	/** the number of peers that have finished sending  counts */
	int m_ranksDoneSending;
	/** the number of peers that have finished counting files */
	int m_ranksDoneCounting;
	/** did the count was sent ?*/
	bool m_sentCount;
	/** the current file to send */
	int m_currentFileToSend;
	/** the current file to count */
	int m_currentFileToCount;
	/** are we currently sending counts ?*/
	bool m_currentlySendingCounts;
	/** counts for a peer slave */
	map<int,LargeCount> m_slaveCounts;
	/** counts for the master node */
	map<int,LargeCount> m_masterCounts;

	/** --- Below are the only elements necessary to bind Partitioner to the Ray software stack */
	/** allocator for outgoing messages */
	RingAllocator*m_outboxAllocator;
	/** message inbox */
	StaticVector*m_inbox;
	/** message outbox */
	StaticVector*m_outbox;
	/** parameters */
	Parameters*m_parameters;

	/** --- */

public:
	void constructor(RingAllocator*outboxAllocator,StaticVector*inbox,StaticVector*outbox,Parameters*parameters,
	SwitchMan*switchMan);
	void call_RAY_MASTER_MODE_COUNT_FILE_ENTRIES();
	void call_RAY_SLAVE_MODE_COUNT_FILE_ENTRIES();

	void registerPlugin(ComputeCore*core);
	void resolveSymbols(ComputeCore*core);
};

#endif
