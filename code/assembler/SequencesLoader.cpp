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

#include<string.h>
#include<assembler/SequencesLoader.h>
#include <core/OperatingSystem.h>
#include<communication/Message.h>
#include<assembler/BubbleData.h>
#include<assert.h>
#include<vector>
#include<core/common_functions.h>
#include<assembler/Loader.h>
#include<structures/Read.h>
#include<iostream>
using namespace std;

void SequencesLoader::registerSequence(){
	if(m_myReads->size()%100000==0){
		uint64_t amount=m_myReads->size();
		cout<<"Rank "<<m_rank<<" has "<<amount<<" sequence reads"<<endl;

		if(m_parameters->showMemoryUsage()){
			showMemoryUsage(m_rank);
		}
	}
	#ifdef ASSERT
	assert(m_distribution_sequence_id<m_loader.size());
	#endif

	Read*theRead=m_loader.at(m_distribution_sequence_id);
	char read[RAY_MAXIMUM_READ_LENGTH];
	theRead->getSeq(read,m_parameters->getColorSpaceMode(),false);

	//cout<<"DEBUG2 Read="<<m_distribution_sequence_id<<" color="<<m_parameters->getColorSpaceMode()<<" Seq= "<<read<<endl;

	Read myRead;
	myRead.constructor(read,&(*m_persistentAllocator),true);
	m_myReads->push_back(&myRead);

	if(m_LOADER_isLeftFile){
		uint64_t leftSequenceGlobalId=m_distribution_currentSequenceId;
		uint64_t leftSequenceIdOnRank=m_myReads->size()-1;

		#ifdef ASSERT
		assert(m_loader.size()!=0);
		#endif

		uint64_t rightSequenceGlobalId=leftSequenceGlobalId+m_loader.size();

		#ifdef ASSERT
		assert(leftSequenceGlobalId<rightSequenceGlobalId);
		#endif

		int rightSequenceRank=m_parameters->getRankFromGlobalId(rightSequenceGlobalId);

		#ifdef ASSERT
		if(rightSequenceRank>=m_size){
			cout<<"m_distribution_currentSequenceId="<<m_distribution_currentSequenceId<<" m_distribution_sequence_id="<<m_distribution_sequence_id<<" LoaderSize="<<m_loader.size()<<" Rank="<<rightSequenceRank<<" Size="<<m_size<<endl;
			assert(rightSequenceRank<m_size);
		}
		#endif

		uint64_t rightSequenceIdOnRank=m_parameters->getIdFromGlobalId(rightSequenceGlobalId);

		int library=m_parameters->getLibrary(m_distribution_file_id);

		(*m_myReads)[leftSequenceIdOnRank]->setLeftType();
		(*m_myReads)[leftSequenceIdOnRank]->getPairedRead()->constructor(rightSequenceRank,rightSequenceIdOnRank,library);
	}else if(m_LOADER_isRightFile){

		#ifdef ASSERT
		assert(m_loader.size()!=0);
		#endif

		uint64_t rightSequenceGlobalId=(m_distribution_currentSequenceId);
		uint64_t rightSequenceIdOnRank=m_myReads->size()-1;
		uint64_t leftSequenceGlobalId=rightSequenceGlobalId-m_loader.size();

		int leftSequenceRank=m_parameters->getRankFromGlobalId(leftSequenceGlobalId);
		#ifdef ASSERT
		if(leftSequenceRank>=m_size){
			cout<<"Global="<<leftSequenceGlobalId<<" rank="<<leftSequenceRank<<endl;
		}
		assert(leftSequenceRank<m_size);
		#endif
		uint64_t leftSequenceIdOnRank=m_parameters->getIdFromGlobalId(leftSequenceGlobalId);
		int library=m_parameters->getLibrary(m_distribution_file_id);

		(*m_myReads)[rightSequenceIdOnRank]->setRightType();
		(*m_myReads)[rightSequenceIdOnRank]->getPairedRead()->constructor(leftSequenceRank,leftSequenceIdOnRank,library);
	// left sequence in interleaved file
	}else if(m_isInterleavedFile && ((m_distribution_sequence_id)%2)==0){
		uint64_t rightSequenceGlobalId=(m_distribution_currentSequenceId)+1;
		int rightSequenceRank=m_parameters->getRankFromGlobalId(rightSequenceGlobalId);
		uint64_t rightSequenceIdOnRank=m_parameters->getIdFromGlobalId(rightSequenceGlobalId);

		uint64_t leftSequenceIdOnRank=m_myReads->size()-1;

		int library=m_parameters->getLibrary(m_distribution_file_id);

		(*m_myReads)[leftSequenceIdOnRank]->setLeftType();
		(*m_myReads)[leftSequenceIdOnRank]->getPairedRead()->constructor(rightSequenceRank,rightSequenceIdOnRank,library);

	// only the right sequence.
	}else if(m_isInterleavedFile &&((m_distribution_sequence_id)%2)==1){
		uint64_t rightSequenceGlobalId=(m_distribution_currentSequenceId);
		uint64_t rightSequenceIdOnRank=m_myReads->size()-1;
		uint64_t leftSequenceGlobalId=rightSequenceGlobalId-1;
		int leftSequenceRank=m_parameters->getRankFromGlobalId(leftSequenceGlobalId);
		uint64_t leftSequenceIdOnRank=m_parameters->getIdFromGlobalId(leftSequenceGlobalId);
		int library=m_parameters->getLibrary(m_distribution_file_id);

		(*m_myReads)[rightSequenceIdOnRank]->setRightType();
		(*m_myReads)[rightSequenceIdOnRank]->getPairedRead()->constructor(leftSequenceRank,leftSequenceIdOnRank,library);
	}
}

bool SequencesLoader::writeSequencesToAMOSFile(int rank,int size,
	StaticVector*m_outbox,
	RingAllocator*m_outboxAllocator,
	bool*m_loadSequenceStep,BubbleData*m_bubbleData,
	time_t*m_lastTime,
	Parameters*m_parameters,int*m_master_mode,int*m_mode
){
	FILE*fp=NULL;
	if(m_parameters->useAmos()){
		fp=fopen(m_parameters->getAmosFile().c_str(),"w");
		// empty the file.
		cout<<"Rank "<<m_rank<<" is adding sequences to "<<m_parameters->getAmosFile()<<endl<<endl;
	}
	m_distribution_sequence_id=0;
	vector<string> allFiles=(*m_parameters).getAllFiles();
	m_loader.constructor(m_parameters->getMemoryPrefix().c_str(),m_parameters->showMemoryAllocations());
	for(m_distribution_file_id=0;m_distribution_file_id<(int)allFiles.size();
		m_distribution_file_id++){
		int res=m_loader.load(allFiles[(m_distribution_file_id)],false);
		if(res==EXIT_FAILURE){
			return false;
		}

		fflush(stdout);
		// write Reads in AMOS format.
		if(rank==MASTER_RANK&&m_parameters->useAmos()){
			char qlt[20000];
			for(uint64_t i=0;i<m_loader.size();i++){
				uint64_t iid=m_distribution_currentSequenceId;
				m_distribution_currentSequenceId++;
				char seq[4000];
				m_loader.at(i)->getSeq(seq,m_parameters->getColorSpaceMode(),true);
				#ifdef ASSERT
				assert(seq!=NULL);
				#endif
				strcpy(qlt,seq);
				// spec: https://sourceforge.net/apps/mediawiki/amos/index.php?title=Message_Types#Sequence_t_:_Universal_t
				for(int j=0;j<(int)strlen(qlt);j++){
					qlt[j]='D';
				}
				#if defined(RAY_64_BITS)
				fprintf(fp,"{RED\niid:%lu\neid:%lu\nseq:\n%s\n.\nqlt:\n%s\n.\n}\n",iid+1,iid+1,seq,qlt);
				#elif defined(RAY_32_BITS)
				fprintf(fp,"{RED\niid:%llu\neid:%llu\nseq:\n%s\n.\nqlt:\n%s\n.\n}\n",iid+1,iid+1,seq,qlt);
				#endif
			}
			m_loader.clear();
			m_loader.load(allFiles[(m_distribution_file_id)],false);
		}
		m_loader.reset();
	}
	m_loader.clear();
	if(m_parameters->useAmos()){
		fclose(fp);
	}
	return true;
}

bool SequencesLoader::loadSequences(int rank,int size,
	StaticVector*m_outbox,
	RingAllocator*m_outboxAllocator,
	bool*m_loadSequenceStep,BubbleData*m_bubbleData,
	time_t*m_lastTime,
	Parameters*m_parameters,int*m_master_mode,int*m_mode
){
	m_rank=rank;

	this->m_parameters=m_parameters;

	printf("Rank %i is loading sequence reads\n",m_rank);
	fflush(stdout);

	/* check if the checkpoint exists */
	if(m_parameters->hasCheckpoint("Sequences")){
		cout<<"Rank "<<m_parameters->getRank()<<" is reading checkpoint Sequences"<<endl;
		cout.flush();

		ifstream f(m_parameters->getCheckpointFile("Sequences").c_str());
		uint64_t count=0;
		f.read((char*)&count,sizeof(uint64_t));
		for(uint64_t i=0;i<count;i++){
			Read myRead;
			myRead.read(&f,m_persistentAllocator);
			m_myReads->push_back(&myRead);
		}

		Message aMessage(NULL,0,MASTER_RANK,RAY_MPI_TAG_SEQUENCES_READY,rank);
		m_outbox->push_back(aMessage);
		(*m_mode)=RAY_SLAVE_MODE_DO_NOTHING;

		cout<<"Rank "<<m_parameters->getRank()<<" loaded "<<count<<" sequences from checkpoint Sequences"<<endl;

		/* true means no error */
		return true;
	}

	// count the number of sequences in all files.
	vector<string> allFiles=(*m_parameters).getAllFiles();

	uint64_t totalNumberOfSequences=0;
	for(int i=0;i<(int)m_parameters->getNumberOfFiles();i++){
		totalNumberOfSequences+=m_parameters->getNumberOfSequences(i);
	}

	uint64_t sequencesPerRank=totalNumberOfSequences/size;
	uint64_t sequencesOnRanksBeforeThisOne=rank*sequencesPerRank;

	uint64_t startingSequenceId=sequencesOnRanksBeforeThisOne;
	uint64_t endingSequenceId=startingSequenceId+sequencesPerRank-1;

	if(rank==size-1){
		endingSequenceId=totalNumberOfSequences-1;
	}

	uint64_t sequences=endingSequenceId-startingSequenceId+1;

	cout<<"Rank "<<m_rank<<" : partition is ["<<startingSequenceId+1;
	cout<<";"<<endingSequenceId+1<<"], "<<sequences<<" sequence reads"<<endl;

	m_distribution_currentSequenceId=0;
	m_loader.constructor(m_parameters->getMemoryPrefix().c_str(),m_parameters->showMemoryAllocations());
	for(m_distribution_file_id=0;m_distribution_file_id<(int)allFiles.size();
		m_distribution_file_id++){

		/** should not load more sequences than required */
		#ifdef ASSERT
		assert(m_myReads->size()<=sequences);
		#endif

		uint64_t sequencesInFile=m_parameters->getNumberOfSequences(m_distribution_file_id);

		if(!(startingSequenceId<m_distribution_currentSequenceId+sequencesInFile)){
			m_distribution_currentSequenceId+=sequencesInFile;
			continue;// skip the file
		}

		if(m_distribution_currentSequenceId>endingSequenceId){
			break;// we are done
		}

		m_loader.load(allFiles[(m_distribution_file_id)],false);

		m_isInterleavedFile=(m_LOADER_isLeftFile)=(m_LOADER_isRightFile)=false;

		if((*m_parameters).isLeftFile((m_distribution_file_id))){
			(m_LOADER_isLeftFile)=true;
		}else if(m_parameters->isRightFile((m_distribution_file_id))){
			(m_LOADER_isRightFile)=true;
		}else if((*m_parameters).isInterleavedFile((m_distribution_file_id))){
			m_isInterleavedFile=true;
		}

		for(m_distribution_sequence_id=0;
			m_distribution_sequence_id<m_loader.size();
				m_distribution_sequence_id++){

			m_loader.at(m_distribution_sequence_id);

			if(m_distribution_currentSequenceId>=startingSequenceId){
				#ifdef ASSERT
				assert(m_distribution_currentSequenceId>=startingSequenceId);
				assert(m_distribution_currentSequenceId<=endingSequenceId);
				#endif

				registerSequence();
			}

			m_distribution_currentSequenceId++;

			if(m_distribution_currentSequenceId>endingSequenceId){
				break;
			}
		}

		m_loader.reset();
	}

	m_loader.clear();
	Message aMessage(NULL,0,MASTER_RANK,RAY_MPI_TAG_SEQUENCES_READY,rank);
	m_outbox->push_back(aMessage);
	(*m_mode)=RAY_SLAVE_MODE_DO_NOTHING;

	uint64_t amount=m_myReads->size();
	cout<<"Rank "<<m_rank<<" has "<<amount<<" sequence reads (completed)"<<endl;

	/* write the checkpoint file */
	if(m_parameters->writeCheckpoints() && !m_parameters->hasCheckpoint("Sequences")){
		/* announce the user that we are writing a checkpoint */
		cout<<"Rank "<<m_parameters->getRank()<<" is writing checkpoint Sequences"<<endl;
		cout.flush();

		ofstream f(m_parameters->getCheckpointFile("Sequences").c_str());

		uint64_t count=m_myReads->size();
		f.write((char*)&count,sizeof(uint64_t));
		for(uint64_t i=0;i<count;i++){
			m_myReads->at(i)->write(&f);
		}

		f.close();
	}

	return true;
}

void SequencesLoader::constructor(int size,MyAllocator*allocator,ArrayOfReads*reads){
	m_size=size;
	m_persistentAllocator=allocator;
	m_myReads=reads;
	m_myReads->constructor(allocator);
}
