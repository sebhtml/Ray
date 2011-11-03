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

#include <core/constants.h>
#include <communication/MessagesHandler.h>
#include <core/common_functions.h>
#include <fstream>
#include <assert.h>
#include <memory/malloc_types.h>
#include <iostream>
#include <sstream>
#include <string.h>
using namespace std;

//#define COMMUNICATION_IS_VERBOSE

/*
 * send messages,
 */
void MessagesHandler::sendMessages(StaticVector*outbox){
	for(int i=0;i<(int)outbox->size();i++){
		Message*aMessage=((*outbox)[i]);
		int destination=aMessage->getDestination();
		void*buffer=aMessage->getBuffer();
		int count=aMessage->getCount();
		int tag=aMessage->getTag();

		#ifdef ASSERT
		assert(destination>=0);
		if(destination>=m_size){
			cout<<"Tag="<<tag<<" Destination="<<destination<<endl;
		}
		assert(destination<m_size);
		assert(!(buffer==NULL && count>0));
		assert(count<=(int)(MAXIMUM_MESSAGE_SIZE_IN_BYTES/sizeof(uint64_t)));
		#endif

		MPI_Request request;

		//  MPI_Issend
		//      Synchronous nonblocking.
		MPI_Isend(buffer,count,m_datatype,destination,tag,MPI_COMM_WORLD,&request);

		// we don't need the request object because we use a ring and this ring is never violated. */
		// this is enforced
		MPI_Request_free(&request);

		#ifdef ASSERT
		assert(request==MPI_REQUEST_NULL);
		#endif

		/** update statistics */
		m_messageStatistics[destination*RAY_MPI_TAG_DUMMY+tag]++;
		m_sentMessages++;
	}

	outbox->clear();
}

#ifdef USE_PERSISTENT_COMMUNICATION
/*
 * receiveMessages is implemented as recommanded by Mr. George Bosilca from
the University of Tennessee (via the Open-MPI mailing list)

De: George Bosilca <bosilca@…>
Reply-to: Open MPI Developers <devel@…>
À: Open MPI Developers <devel@…>
Sujet: Re: [OMPI devel] Simple program (103 lines) makes Open-1.4.3 hang
Date: 2010-11-23 18:03:04

If you know the max size of the receives I would take a different approach.
Post few persistent receives, and manage them in a circular buffer.
Instead of doing an MPI_Iprobe, use MPI_Test on the current head of your circular buffer.
Once you use the data related to the receive, just do an MPI_Start on your request.
This approach will minimize the unexpected messages, and drain the connections faster.
Moreover, at the end it is very easy to MPI_Cancel all the receives not yet matched.

    george.
 */
void MessagesHandler::pumpMessageFromPersistentRing(){
	/* persistent communication is not enabled by default */
	int flag;
	MPI_Status status;
	MPI_Test(m_ring+m_head,&flag,&status);

	if(flag){
		// get the length of the message
		// it is not necessary the same as the one posted with MPI_Recv_init
		// that one was a lower bound
		int tag=status.MPI_TAG;
		int source=status.MPI_SOURCE;
		int count;
		MPI_Get_count(&status,m_datatype,&count);
		uint64_t*filledBuffer=(uint64_t*)m_buffers+m_head*MAXIMUM_MESSAGE_SIZE_IN_BYTES/sizeof(uint64_t);

		// copy it in a safe buffer
		// internal buffers are all of length MAXIMUM_MESSAGE_SIZE_IN_BYTES
		uint64_t*incoming=(uint64_t*)m_internalBufferAllocator.allocate(MAXIMUM_MESSAGE_SIZE_IN_BYTES);

		memcpy(incoming,filledBuffer,count*sizeof(uint64_t));

		// the request can start again
		MPI_Start(m_ring+m_head);

		// add the message in the internal inbox
		// this inbox is not the real inbox
		Message aMessage(incoming,count,m_rank,tag,source);

		addMessage(&aMessage);

		if(m_head == m_ringSize-1)
			m_head = 0;
		else
			m_head++;
	}
}
#endif


#ifdef USE_PERSISTENT_COMMUNICATION
void MessagesHandler::addMessage(Message*message){

	m_bufferedMessages ++ ;

	int source=message->getSource();

	// we allocate memory for internal storage
	// this message may not be returned immediately
	MessageNode*allocatedMessage = (MessageNode*) m_internalMessageAllocator.allocate(sizeof(MessageNode));

	// the buffer of the message is already allocated in the calling function so eveyrthing is OK
	allocatedMessage->m_message = *message; // copy data
	allocatedMessage->m_next = NULL;
	allocatedMessage->m_previous = NULL;

	// there are no message from this source
	if(m_heads[source] == NULL && m_tails[source] == NULL){
		#ifdef ASSERT
		assert(m_heads[source] == NULL && m_tails[source] == NULL);
		#endif

		m_heads[source] = allocatedMessage;
		m_tails[source] = allocatedMessage;

	// This source already have messages in the linked lists
	}else{ // we add at the head
		#ifdef ASSERT
		assert(m_heads[source] != NULL && m_tails[source] != NULL);
		#endif

		// update the previous pointer of the old head

		#ifdef ASSERT
		assert(m_heads[source]->m_previous == NULL);
		#endif

		m_heads[source]->m_previous= allocatedMessage;

		// Assign the next pointer for the new head
		allocatedMessage->m_next = m_heads[source];

		// we have a new head
		m_heads[source] = allocatedMessage;
	}
}
#endif

#define USE_ROUND_ROBIN
void MessagesHandler::receiveMessages(StaticVector*inbox,RingAllocator*inboxAllocator){
	#ifdef USE_PERSISTENT_COMMUNICATION
	// pump messages with persistent communication
	pumpMessageFromPersistentRing();

	// there is at least one message to fetch
	if(m_bufferedMessages > 0){
		while(inbox->size() == 0){
			// pick a message according to a round-robin policy
			roundRobinReception(inbox,inboxAllocator);
		}
	}
	#elif defined USE_ROUND_ROBIN

	// round-robin reception seems to avoid starvation
/** round robin with Iprobe will increase the latency because there will be a lot of calls to
 MPI_Iprobe that yield no messages at all */

	roundRobinReception(inbox,inboxAllocator);

	#else

	// receive any message
	// it is assumed that MPI is fair
	// otherwise there may be some starvation

	probeAndRead(MPI_ANY_SOURCE,MPI_ANY_TAG,inbox,inboxAllocator);

	#endif
}

void MessagesHandler::roundRobinReception(StaticVector*inbox,RingAllocator*inboxAllocator){
	#ifdef USE_PERSISTENT_COMMUNICATION
	// with persistent communication, we fetch the message from the persistent ring
	// check if we have one message for m_currentRankToTryToReceiveFrom
	if(m_tails[m_currentRankToTryToReceiveFrom] != NULL){
		// we have a message

		MessageNode*messageNode = m_tails[m_currentRankToTryToReceiveFrom];
		Message*selectedMessage=&(messageNode->m_message);

		// set the new tail
		m_tails[m_currentRankToTryToReceiveFrom] = messageNode->m_previous;

		// we just took the last one
		if(m_tails[m_currentRankToTryToReceiveFrom] == NULL){
			m_heads[m_currentRankToTryToReceiveFrom] = NULL; // set the head to NULL too
		}else{
			// it was not the last one
			// we need to update the next of the new tail
			m_tails[m_currentRankToTryToReceiveFrom]->m_next = NULL;
		}

		#ifdef ASSERT
		if(m_tails[m_currentRankToTryToReceiveFrom] != NULL)
			assert(m_tails[m_currentRankToTryToReceiveFrom]->m_next == NULL);
		#endif

		int count = selectedMessage->getCount();

		// allocate the buffer using the ring buffer for that
		uint64_t*incoming=(uint64_t*)inboxAllocator->allocate(count*sizeof(uint64_t));

		// copy the data
		memcpy(incoming,selectedMessage->getBuffer(),count*sizeof(uint64_t));

		// recycle the old buffer
		// all internal buffers are MAXIMUM_MESSAGE_SIZE_IN_BYTES bytes, not count
		m_internalBufferAllocator.free(selectedMessage->getBuffer(),MAXIMUM_MESSAGE_SIZE_IN_BYTES);

		// set the newly created buffer
		selectedMessage->setBuffer(incoming);

		// push the message in the inbox
		inbox->push_back(*selectedMessage);

		// recycle the old message
		m_internalMessageAllocator.free(messageNode,sizeof(MessageNode));

		/** update statistics */
		m_receivedMessages++;

		m_bufferedMessages --;

		#ifdef ASSERT
		assert(m_bufferedMessages >= 0);
		#endif
	}

	#else

	/* otherwise, we use MPI_Iprobe + MPI_Recv */
	/* probe and read a message */

	probeAndRead(m_currentRankToTryToReceiveFrom,MPI_ANY_TAG,inbox,inboxAllocator);

	#endif


	#ifdef COMMUNICATION_IS_VERBOSE
	if(inbox->size() > 0){
		cout<<"[RoundRobin] received a message from "<<m_currentRankToTryToReceiveFrom<<endl;
	}
	#endif

	/* advance the rank to receive from */
	if(m_currentRankToTryToReceiveFrom == m_size-1){
		/* restart the loop */
		m_currentRankToTryToReceiveFrom = 0;
	}else{
		m_currentRankToTryToReceiveFrom++;
	}

}

// this code is not utilised, but is kept for reference
void MessagesHandler::probeAndRead(int source,int tag,StaticVector*inbox,RingAllocator*inboxAllocator){
	#ifdef COMMUNICATION_IS_VERBOSE
	//cout<<"call to probeAndRead source="<<source<<""<<endl;
	#endif
	int flag;
	MPI_Status status;
	MPI_Iprobe(source,tag,MPI_COMM_WORLD,&flag,&status);

	/* read at most one message */
	if(flag){
		MPI_Datatype datatype=MPI_UNSIGNED_LONG_LONG;
		int tag=status.MPI_TAG;
		int source=status.MPI_SOURCE;
		int count=-1;
		MPI_Get_count(&status,datatype,&count);

		#ifdef ASSERT
		assert(count >= 0);
		#endif

		uint64_t*incoming=NULL;
		if(count > 0){
			incoming=(uint64_t*)inboxAllocator->allocate(count*sizeof(uint64_t));
		}

		MPI_Recv(incoming,count,datatype,source,tag,MPI_COMM_WORLD,&status);

		Message aMessage(incoming,count,m_rank,tag,source);
		inbox->push_back(aMessage);

		#ifdef ASSERT
		assert(aMessage.getDestination() == m_rank);
		#endif

		m_receivedMessages++;
	}
}

void MessagesHandler::initialiseMembers(){
	#ifdef USE_PERSISTENT_COMMUNICATION
	// the ring itself  contain requests ready to receive messages
	m_ringSize=m_size+16;

	m_ring=(MPI_Request*)__Malloc(sizeof(MPI_Request)*m_ringSize,RAY_MALLOC_TYPE_PERSISTENT_MESSAGE_RING,false);
	m_buffers=(uint8_t*)__Malloc(MAXIMUM_MESSAGE_SIZE_IN_BYTES*m_ringSize,RAY_MALLOC_TYPE_PERSISTENT_MESSAGE_BUFFERS,false);
	m_head=0;

	// post a few receives.
	for(int i=0;i<m_ringSize;i++){
		uint8_t*buffer=m_buffers+i*MAXIMUM_MESSAGE_SIZE_IN_BYTES;
		MPI_Recv_init(buffer,MAXIMUM_MESSAGE_SIZE_IN_BYTES/sizeof(uint64_t),m_datatype,
			MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,m_ring+i);
		MPI_Start(m_ring+i);
	}

	m_heads= (MessageNode**) __Malloc(sizeof(MessageNode*)*m_size,RAY_MALLOC_TYPE_COMMUNICATION_LAYER,false);
	m_tails= (MessageNode**) __Malloc(sizeof(MessageNode*)*m_size,RAY_MALLOC_TYPE_COMMUNICATION_LAYER,false);

	for(int i=0;i<m_size;i++){
		m_heads[i] = NULL;
		m_tails[i] = NULL;
	}
	#endif
}

void MessagesHandler::freeLeftovers(){
	#ifdef USE_PERSISTENT_COMMUNICATION
	for(int i=0;i<m_ringSize;i++){
		MPI_Cancel(m_ring+i);
		MPI_Request_free(m_ring+i);
	}
	__Free(m_ring,RAY_MALLOC_TYPE_PERSISTENT_MESSAGE_RING,false);
	m_ring=NULL;
	__Free(m_buffers,RAY_MALLOC_TYPE_PERSISTENT_MESSAGE_BUFFERS,false);
	m_buffers=NULL;

	__Free(m_messageStatistics,RAY_MALLOC_TYPE_MESSAGE_STATISTICS,false);
	m_messageStatistics=NULL;

	__Free(m_heads,RAY_MALLOC_TYPE_COMMUNICATION_LAYER,false);
	__Free(m_tails,RAY_MALLOC_TYPE_COMMUNICATION_LAYER,false);
	m_heads= NULL;
	m_tails= NULL;
	#endif
}

void MessagesHandler::constructor(int*argc,char***argv){
	m_sentMessages=0;
	m_receivedMessages=0;
	m_datatype=MPI_UNSIGNED_LONG_LONG;
	MPI_Init(argc,argv);
	char serverName[1000];
	int len;

	/** initialize the message passing interface stack */
	MPI_Get_processor_name(serverName,&len);
	MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&m_size);
	initialiseMembers();
	m_processorName=serverName;

	/** initialize message statistics to 0 */
	m_messageStatistics=(uint64_t*)__Malloc(RAY_MPI_TAG_DUMMY*m_size*sizeof(uint64_t),RAY_MALLOC_TYPE_MESSAGE_STATISTICS,false);
	for(int rank=0;rank<m_size;rank++){
		for(int tag=0;tag<RAY_MPI_TAG_DUMMY;tag++){
			m_messageStatistics[rank*RAY_MPI_TAG_DUMMY+tag]=0;
		}
	}

	m_currentRankToTryToReceiveFrom=m_rank;

	#ifdef USE_PERSISTENT_COMMUNICATION
	int chunkSize = 4194304; // 4 MiB
	m_internalMessageAllocator.constructor(chunkSize,RAY_MALLOC_TYPE_COMMUNICATION_LAYER,false);
	m_internalBufferAllocator.constructor(chunkSize,RAY_MALLOC_TYPE_COMMUNICATION_LAYER,false);

	m_bufferedMessages=0;
	#endif
}

void MessagesHandler::destructor(){
	MPI_Finalize();
}

string*MessagesHandler::getName(){
	return &m_processorName;
}

int MessagesHandler::getRank(){
	return m_rank;
}

int MessagesHandler::getSize(){
	return m_size;
}

void MessagesHandler::barrier(){
	MPI_Barrier(MPI_COMM_WORLD);
}

void MessagesHandler::version(int*a,int*b){
	MPI_Get_version(a,b);
}

void MessagesHandler::appendStatistics(const char*file){
	/** add an entry for RAY_MPI_TAG_GOOD_JOB_SEE_YOU_SOON_REPLY
 * 	because this message is sent after writting the current file
 */
	m_messageStatistics[MASTER_RANK*RAY_MPI_TAG_DUMMY+RAY_MPI_TAG_GOOD_JOB_SEE_YOU_SOON_REPLY]++;
	m_sentMessages++;

	ofstream fp;
	fp.open(file,ios_base::out|ios_base::app);
	for(int destination=0;destination<m_size;destination++){
		for(int tag=0;tag<RAY_MPI_TAG_DUMMY;tag++){
			uint64_t count=m_messageStatistics[destination*RAY_MPI_TAG_DUMMY+tag];
			fp<<m_rank<<"\t"<<destination<<"\t"<<MESSAGES[tag]<<"\t"<<count<<"\n";
		}
	}
	fp.close();
	cout<<"Rank "<<m_rank<<": sent "<<m_sentMessages<<" messages, received "<<m_receivedMessages<<" messages."<<endl;
}

string MessagesHandler::getMessagePassingInterfaceImplementation(){
	ostringstream implementation;

	#ifdef MPICH2
        implementation<<"MPICH2 (MPICH2) "<<MPICH2_VERSION;
	#endif

	#ifdef OMPI_MPI_H
        implementation<<"Open-MPI (OMPI_MPI_H) "<<OMPI_MAJOR_VERSION<<"."<<OMPI_MINOR_VERSION<<"."<<OMPI_RELEASE_VERSION;
	#endif

	if(implementation.str().length()==0)
		implementation<<"Unknown";
	return implementation.str();
}

