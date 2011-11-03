/*
 	Ray
    Copyright (C) 2011  Sébastien Boisvert

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

#ifndef _Amos
#define _Amos

#include <core/Parameters.h>
#include <structures/StaticVector.h>
#include <memory/RingAllocator.h>
#include <stdio.h>
#include <assembler/FusionData.h>
#include <assembler/ExtensionData.h>
#include <scaffolder/Scaffolder.h>
#include <assembler/ReadFetcher.h>
#include <stdint.h>
#include <communication/VirtualCommunicator.h>
#include <vector>
using namespace std;

/**
 * AMOS specification is available : http://sourceforge.net/apps/mediawiki/amos/index.php?title=AMOS
 * \see http://sourceforge.net/apps/mediawiki/amos/index.php?title=Message_Types
 * \author Sébastien Boisvert
 */
class Amos{
	VirtualCommunicator*m_virtualCommunicator;
	StaticVector*m_inbox;
	uint64_t m_workerId;
	ReadFetcher m_readFetcher;
	int*m_master_mode;
	int*m_slave_mode;
	Scaffolder*m_scaffolder;
	FusionData*m_fusionData;
	FILE*m_amosFile;
	Parameters*m_parameters;
	RingAllocator*m_outboxAllocator;
	StaticVector*m_outbox;
	int m_contigId;
	int m_sequence_id;
	ExtensionData*m_ed;
	int m_mode_send_vertices_sequence_id_position;
	vector<uint64_t> m_activeWorkers;
public:
	void masterMode();
	void slaveMode();
	void constructor(Parameters*parameters,RingAllocator*outboxAllocator,StaticVector*outbox,
		FusionData*fusionData,ExtensionData*extensionData,int*masterMode,int*slaveMode,Scaffolder*scaffolder,
StaticVector*inbox,VirtualCommunicator*virtualCommunicator);
};


#endif

