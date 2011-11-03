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

#ifndef _NetworkTest_H
#define _NetworkTest_H

#include <structures/StaticVector.h>
#include <core/Parameters.h>
#include <assembler/TimePrinter.h>
#include <memory/RingAllocator.h>
#include <string>
#include <map>
using namespace std;

/**
 * This class tests the network
 * Tested elements:
 *
 *     - latency to get a response for a message
 *
 *     Dependencies from the Ray software stack (mostly):
 *
 *          - message inbox
 *          - message outbox
 *          - outbox memory ring allocator
 *          - slave and master modes.
 * \author Sébastien Boisvert
 */
class NetworkTest{
	/* details of the network test */

	uint64_t m_averageLatencyInMicroSeconds;

	vector<int> m_destinations;
	vector<uint64_t> m_sentMicroseconds;
	vector<uint64_t> m_receivedMicroseconds;

	/* number of words to use for network test */
	/* a word is 8 bytes */
	/* MAXIMUM_MESSAGE_SIZE_IN_BYTES is 4000 per default so
		numberOfWords must be <= 500 */
	/* this is only for the network test */
	/* default is 500 */
	int m_numberOfWords;

	TimePrinter*m_timePrinter;
	/** the message inbox */
	StaticVector*m_inbox;
	/** the message outbox */
	StaticVector*m_outbox;
	/** parameter object */
	Parameters*m_parameters;
	/** the slave mode */
	int*m_slaveMode;
	/** the master mode, always RAY_SLAVE_MODE_DO_NOTHING for rank >0 */
	int*m_masterMode;
	/** message-passing interface rank */
	int m_rank;
	/** number of ranks */
	int m_size;
	/** number of ranks who finished the tests */
	int m_doneWithNetworkTest;
	/* initialised this ? */
	bool m_initialisedNetworkTest;
	/** outbox allocator */
	RingAllocator*m_outboxAllocator;
	/** the number of test message to send */
	int m_numberOfTestMessages;
	/** the current test message */
	int m_currentTestMessage;
	/** the current test message has been sent ? */
	bool m_sentCurrentTestMessage;

	/** the starting time */
	uint64_t m_startingTimeMicroseconds;

	uint64_t m_sumOfMicroSeconds;

	/** latencies */
	map<int,int> m_latencies;
	map<int,string> m_names;

	/* processor name */
	string*m_name;
public:
	/** initialize the NetworkTest */
	void constructor(int rank,int*masterMode,int*slaveMode,int size,StaticVector*inbox,StaticVector*outbox,Parameters*parameters,RingAllocator*outboxAllocator,
		string*name,TimePrinter*timePrinter);
	/** work method for the master mode */
	void masterWork();
	/** work method for the slave mode */
	void slaveWork();

	void writeData();
};

#endif

