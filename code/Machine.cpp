/*
 	Ray
    Copyright (C) 2010  Sébastien Boisvert

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


 	Funding:

Sébastien Boisvert has a scholarship from the Canadian Institutes of Health Research (Master's award: 200910MDR-215249-172830 and Doctoral award: 200902CGM-204212-172830).

*/

#include<mpi.h>
#include<EdgesExtractor.h>
#include<Machine.h>
#include<VerticesExtractor.h>
#include<sstream>
#include<Message.h>
#include<time.h>
#include<RepeatedVertexWatchdog.h>
#include<TipWatchdog.h>
#include<BubbleTool.h>
#include<assert.h>
#include<common_functions.h>
#include<iostream>
#include<fstream>
#include<CoverageDistribution.h>
#include<string.h>
#include<SplayTreeIterator.h>
#include<Read.h>
#include<Loader.h>
#include<MyAllocator.h>
#include<algorithm>
#include<unistd.h>

bool myComparator_sort(const vector<VERTEX_TYPE>&a,const vector<VERTEX_TYPE>&b){
	return a.size()>b.size();
}

using namespace std;

void Machine::showUsage(){
	cout<<endl;
	cout<<"Usage:"<<endl<<endl;
	cout<<"Supported sequences file format: "<<endl;

	cout<<".fasta"<<endl;
	#ifdef HAVE_ZLIB
	cout<<".fasta.gz"<<endl;
	#endif
	#ifdef HAVE_LIBBZ2
	cout<<".fasta.bz2"<<endl;
	#endif
	cout<<".fastq"<<endl;
	#ifdef HAVE_ZLIB
	cout<<".fastq.gz"<<endl;
	#endif
	#ifdef HAVE_LIBBZ2
	cout<<".fastq.bz2"<<endl;
	#endif
	cout<<".sff (paired reads must be extracted manually)"<<endl;

	cout<<endl;

	cout<<"Parameters:"<<endl;
	cout<<endl;

	cout<<"Single-end reads"<<endl;
    	cout<<" -s <sequencesFile>"<<endl;
	cout<<endl;
	cout<<"Paired-end reads:"<<endl;
	cout<<" -p <leftSequencesFile> <rightSequencesFile> [ <fragmentLength> <standardDeviation> ]"<<endl;
	cout<<endl;
	cout<<"Paired-end reads:"<<endl;
	cout<<" -i <interleavedFile> [ <fragmentLength> <standardDeviation> ]"<<endl;
	cout<<endl;
	cout<<"Output (default: RayOutput)"<<endl;
	cout<<" -o <outputPrefix>"<<endl;
	cout<<endl;	
	cout<<"AMOS output"<<endl;
	cout<<" -a  "<<endl;
    	cout<<endl;
	cout<<"k-mer size (default: 21)"<<endl;
	cout<<" -k <kmerSize>"<<endl;
	cout<<endl;
	cout<<"Ray writes a contigs file, a coverage distribution file, and an AMOS file (if -a is provided)."<<endl;
	cout<<"The name of these files is based on the value provided with -o."<<endl;
	cout<<endl;
	cout<<"use --help to show this help"<<endl;
	cout<<endl;
}

void Machine::sendLibraryDistances(){
	if(!m_ready){
		return;
	}
	if(m_libraryIterator==(int)m_libraryDistances.size()){

		m_bufferedData.flushAll(TAG_LIBRARY_DISTANCE,&m_outboxAllocator,&m_outbox,getRank());

		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,MASTER_RANK,TAG_ASK_LIBRARY_DISTANCES_FINISHED,getRank());
		m_outbox.push_back(aMessage);
		m_mode=MODE_DO_NOTHING;
	}else if(m_libraryIndex==m_libraryDistances[m_libraryIterator].end()){
		m_libraryIterator++;
		m_libraryIndexInitiated=false;
	}else{
		if(!m_libraryIndexInitiated){
			m_libraryIndexInitiated=true;
			m_libraryIndex=m_libraryDistances[m_libraryIterator].begin();
		}
		int library=m_libraryIterator;
		int distance=m_libraryIndex->first;
		int count=m_libraryIndex->second;
		m_bufferedData.addAt(MASTER_RANK,library);
		m_bufferedData.addAt(MASTER_RANK,distance);
		m_bufferedData.addAt(MASTER_RANK,count);
		if(m_bufferedData.flush(MASTER_RANK,3,TAG_LIBRARY_DISTANCE,&m_outboxAllocator,&m_outbox,getRank(),false)){
			m_ready=false;
		}

		m_libraryIndex++;
	}
}

/*
 * get the Directions taken by a vertex.
 *
 * m_Machine_getPaths_INITIALIZED must be set to false before any calls.
 * also, you must set m_Machine_getPaths_DONE to false;
 *
 * when done, m_Machine_getPaths_DONE is true
 * and
 * the result is in m_Machine_getPaths_result (a vector<Direction>)
 */
void Machine::getPaths(VERTEX_TYPE vertex){
	if(!m_Machine_getPaths_INITIALIZED){
		m_Machine_getPaths_INITIALIZED=true;
		m_fusionData->m_FUSION_paths_requested=false;
		m_Machine_getPaths_DONE=false;
		m_Machine_getPaths_result.clear();
		return;
	}

	if(!m_fusionData->m_FUSION_paths_requested){
		VERTEX_TYPE theVertex=vertex;
		VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
		message[0]=theVertex;
		Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(theVertex),TAG_ASK_VERTEX_PATHS_SIZE,getRank());
		m_outbox.push_back(aMessage);
		m_fusionData->m_FUSION_paths_requested=true;
		m_fusionData->m_FUSION_paths_received=false;
		m_fusionData->m_FUSION_path_id=0;
		m_fusionData->m_FUSION_path_requested=false;
		m_fusionData->m_FUSION_receivedPaths.clear();
	}else if(m_fusionData->m_FUSION_paths_received){
		if(m_fusionData->m_FUSION_path_id<m_fusionData->m_FUSION_numberOfPaths){
			if(!m_fusionData->m_FUSION_path_requested){
				VERTEX_TYPE theVertex=vertex;
				VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
				message[0]=m_fusionData->m_FUSION_path_id;
				Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(theVertex),TAG_ASK_VERTEX_PATH,getRank());
				m_outbox.push_back(aMessage);
				m_fusionData->m_FUSION_path_requested=true;
				m_fusionData->m_FUSION_path_received=false;
			}else if(m_fusionData->m_FUSION_path_received){
				m_fusionData->m_FUSION_path_id++;
				m_Machine_getPaths_result.push_back(m_fusionData->m_FUSION_receivedPath);
				m_fusionData->m_FUSION_path_requested=false;
			}
		}else{
			m_Machine_getPaths_DONE=true;
		}
	}
}

Machine::Machine(int argc,char**argv){
	m_argc=argc;
	m_argv=argv;
	m_bubbleData=new BubbleData();
	#ifdef SHOW_SENT_MESSAGES
	m_stats=new StatisticsData();
	#endif
	m_dfsData=new DepthFirstSearchData();
	m_fusionData=new FusionData();
	m_disData=new DistributionData();
	m_seedingData=new SeedingData();
	m_ed=new ExtensionData();
	m_sd=new ScaffolderData();
	m_cd=new ChooserData();


	m_master_methods[MASTER_MODE_LOAD_CONFIG]=&Machine::call_MASTER_MODE_LOAD_CONFIG;
	m_master_methods[MASTER_MODE_LOAD_SEQUENCES]=&Machine::call_MASTER_MODE_LOAD_SEQUENCES;
	m_master_methods[MASTER_MODE_TRIGGER_VERTICE_DISTRIBUTION]=&Machine::call_MASTER_MODE_TRIGGER_VERTICE_DISTRIBUTION;
	m_master_methods[MASTER_MODE_SEND_COVERAGE_VALUES]=&Machine::call_MASTER_MODE_SEND_COVERAGE_VALUES;
	m_master_methods[MASTER_MODE_TRIGGER_EDGES_DISTRIBUTION]=&Machine::call_MASTER_MODE_TRIGGER_EDGES_DISTRIBUTION;
	m_master_methods[MASTER_MODE_START_EDGES_DISTRIBUTION]=&Machine::call_MASTER_MODE_START_EDGES_DISTRIBUTION;
	m_master_methods[MASTER_MODE_DO_NOTHING]=&Machine::call_MASTER_MODE_DO_NOTHING;
	m_master_methods[MASTER_MODE_UPDATE_DISTANCES]=&Machine::call_MASTER_MODE_UPDATE_DISTANCES;
	m_master_methods[MASTER_MODE_ASK_EXTENSIONS]=&Machine::call_MASTER_MODE_ASK_EXTENSIONS;
	m_master_methods[MASTER_MODE_AMOS]=&Machine::call_MASTER_MODE_AMOS;
	m_master_methods[MASTER_MODE_PREPARE_DISTRIBUTIONS]=&Machine::call_MASTER_MODE_PREPARE_DISTRIBUTIONS;
	m_master_methods[MASTER_MODE_TRIGGER_EDGES]=&Machine::call_MASTER_MODE_TRIGGER_EDGES;
	m_master_methods[MASTER_MODE_TRIGGER_INDEXING]=&Machine::call_MASTER_MODE_TRIGGER_INDEXING;
	m_master_methods[MASTER_MODE_INDEX_SEQUENCES]=&Machine::call_MASTER_MODE_INDEX_SEQUENCES;
	m_master_methods[MASTER_MODE_PREPARE_DISTRIBUTIONS_WITH_ANSWERS]=&Machine::call_MASTER_MODE_PREPARE_DISTRIBUTIONS_WITH_ANSWERS;
	m_master_methods[MASTER_MODE_PREPARE_SEEDING]=&Machine::call_MASTER_MODE_PREPARE_SEEDING;
	m_master_methods[MASTER_MODE_TRIGGER_SEEDING]=&Machine::call_MASTER_MODE_TRIGGER_SEEDING;
	m_master_methods[MASTER_MODE_TRIGGER_DETECTION]=&Machine::call_MASTER_MODE_TRIGGER_DETECTION;
	m_master_methods[MASTER_MODE_ASK_DISTANCES]=&Machine::call_MASTER_MODE_ASK_DISTANCES;
	m_master_methods[MASTER_MODE_START_UPDATING_DISTANCES]=&Machine::call_MASTER_MODE_START_UPDATING_DISTANCES;
	m_master_methods[MASTER_MODE_TRIGGER_EXTENSIONS]=&Machine::call_MASTER_MODE_TRIGGER_EXTENSIONS;
	m_master_methods[MASTER_MODE_TRIGGER_FUSIONS]=&Machine::call_MASTER_MODE_TRIGGER_FUSIONS;
	m_master_methods[MASTER_MODE_TRIGGER_FIRST_FUSIONS]=&Machine::call_MASTER_MODE_TRIGGER_FIRST_FUSIONS;
	m_master_methods[MASTER_MODE_START_FUSION_CYCLE]=&Machine::call_MASTER_MODE_START_FUSION_CYCLE;

	m_slave_methods[MODE_START_SEEDING]=&Machine::call_MODE_START_SEEDING;
	m_slave_methods[MODE_DO_NOTHING]=&Machine::call_MODE_DO_NOTHING;
	m_slave_methods[MODE_SEND_EXTENSION_DATA]=&Machine::call_MODE_SEND_EXTENSION_DATA;
	m_slave_methods[MODE_ASSEMBLE_WAVES]=&Machine::call_MODE_ASSEMBLE_WAVES;
	m_slave_methods[MODE_FUSION]=&Machine::call_MODE_FUSION;
	m_slave_methods[MODE_PERFORM_CALIBRATION]=&Machine::call_MODE_PERFORM_CALIBRATION;
	m_slave_methods[MODE_INDEX_SEQUENCES]=&Machine::call_MODE_INDEX_SEQUENCES;
	m_slave_methods[MODE_FINISH_FUSIONS]=&Machine::call_MODE_FINISH_FUSIONS;
	m_slave_methods[MODE_DISTRIBUTE_FUSIONS]=&Machine::call_MODE_DISTRIBUTE_FUSIONS;
	m_slave_methods[MODE_AUTOMATIC_DISTANCE_DETECTION]=&Machine::call_MODE_AUTOMATIC_DISTANCE_DETECTION;
	m_slave_methods[MODE_SEND_LIBRARY_DISTANCES]=&Machine::call_MODE_SEND_LIBRARY_DISTANCES;
	m_slave_methods[MODE_EXTRACT_VERTICES]=&Machine::call_MODE_EXTRACT_VERTICES;
	m_slave_methods[MODE_SEND_DISTRIBUTION]=&Machine::call_MODE_SEND_DISTRIBUTION;
	m_slave_methods[MODE_PROCESS_INGOING_EDGES]=&Machine::call_MODE_PROCESS_INGOING_EDGES;
	m_slave_methods[MODE_PROCESS_OUTGOING_EDGES]=&Machine::call_MODE_PROCESS_OUTGOING_EDGES;
	m_slave_methods[MODE_EXTENSION]=&Machine::call_MODE_EXTENSION;


}


void Machine::start(){
	m_ready=true;
	m_maxCoverage=0;
	m_maxCoverage--;// underflow.
	int numberOfTrees=_FOREST_SIZE;

	#ifdef SHOW_PROGRESS
	#endif
	srand(m_lastTime);
	m_fusionData->m_fusionStarted=false;
	m_ed->m_EXTENSION_numberOfRanksDone=0;
	m_colorSpaceMode=false;
	m_messageSentForEdgesDistribution=false;
	m_numberOfRanksDoneSeeding=0;
	m_numberOfMachinesReadyForEdgesDistribution=0;
	m_ed->m_mode_EXTENSION=false;
	m_aborted=false;
	m_readyToSeed=0;
	m_wordSize=-1;
	m_reverseComplementVertex=false;
	m_last_value=0;
	m_mode_send_ingoing_edges=false;
	m_mode_send_vertices=false;
	m_mode_sendDistribution=false;
	m_mode_send_outgoing_edges=false;
	m_mode_send_vertices_sequence_id=0;
	m_mode_send_vertices_sequence_id_position=0;
	m_numberOfMachinesDoneSendingVertices=0;
	m_numberOfMachinesDoneSendingEdges=0;
	m_numberOfMachinesReadyToSendDistribution=0;
	m_numberOfMachinesDoneSendingCoverage=0;
	m_machineRank=0;
	m_messageSentForVerticesDistribution=false;
	m_sequence_ready_machines=0;
	m_isFinalFusion=false;

	m_outboxAllocator.constructor(MAX_ALLOCATED_MESSAGES_IN_OUTBOX,MPI_BTL_SM_EAGER_LIMIT);
	m_inboxAllocator.constructor(MAX_ALLOCATED_MESSAGES_IN_INBOX,MPI_BTL_SM_EAGER_LIMIT);
	m_persistentAllocator.constructor(PERSISTENT_ALLOCATOR_CHUNK_SIZE);
	m_directionsAllocator.constructor(PERSISTENT_ALLOCATOR_CHUNK_SIZE);

	m_mode=MODE_DO_NOTHING;
	m_master_mode=MASTER_MODE_DO_NOTHING;
	m_mode_AttachSequences=false;
	m_startEdgeDistribution=false;

	m_ranksDoneAttachingReads=0;

	char serverName[1000];
	int len;
	MPI_Init(&m_argc,&m_argv);
	MPI_Get_processor_name(serverName,&len);

	MPI_Comm_rank(MPI_COMM_WORLD,&m_rank);
	MPI_Comm_size(MPI_COMM_WORLD,&m_size);
	m_sl.constructor(m_size);

	assert(getSize()<=MAX_NUMBER_OF_MPI_PROCESSES);
	
	m_fusionData->constructor(getSize(),MPI_BTL_SM_EAGER_LIMIT);

	if(isMaster()){
		cout<<"Bienvenue !"<<endl;
		cout<<endl;

		cout<<"Rank 0: Ray "<<RAY_VERSION<<endl;

		#ifdef MPICH2
                cout<<"Rank 0: compiled with MPICH2 "<<MPICH2_VERSION<<endl;
		#endif

		#ifdef OMPI_MPI_H
                cout<<"Rank 0: compiled with Open-MPI "<<OMPI_MAJOR_VERSION<<"."<<OMPI_MINOR_VERSION<<"."<<OMPI_RELEASE_VERSION<<endl;
		#endif

		#ifdef HAVE_ZLIB
		cout<<"Rank 0: compiled with GZIP"<<endl;
		#endif

		#ifdef HAVE_LIBBZ2
		cout<<"Rank 0: compiled with BZIP2"<<endl;
		#endif
		cout<<endl;
		m_timePrinter.printElapsedTime("Beginning of computation");
		cout<<endl;
	}
	m_disData->constructor(getSize(),MPI_BTL_SM_EAGER_LIMIT,&m_persistentAllocator);
	m_bufferedData.constructor(getSize(),MPI_BTL_SM_EAGER_LIMIT);
	m_library.constructor(getRank(),&m_outbox,&m_outboxAllocator,&m_bufferedData,&m_sequence_id,&m_sequence_idInFile,
		m_ed,&m_readsPositions,getSize(),&m_timePrinter,&m_mode,&m_master_mode,
	&m_parameters,&m_fileId,m_seedingData,&m_libraryDistances);


	m_subgraph.constructor(numberOfTrees,&m_persistentAllocator);
	
	m_edgesExtractor.m_disData=m_disData;
	m_edgesExtractor.getRank=getRank();
	m_edgesExtractor.getSize=getSize();
	m_edgesExtractor.m_outboxAllocator=&m_outboxAllocator;
	m_edgesExtractor.m_outbox=&m_outbox;
	m_edgesExtractor.m_mode_send_outgoing_edges=&m_mode_send_outgoing_edges;
	m_edgesExtractor.m_mode_send_ingoing_edges=&m_mode_send_ingoing_edges;
	m_edgesExtractor.m_colorSpaceMode=m_colorSpaceMode;
	m_edgesExtractor.m_myReads=&m_myReads;
	m_edgesExtractor.m_mode=&m_mode;

	if(isMaster()){
		cout<<endl<<"**************************************************"<<endl;
    		cout<<"This program comes with ABSOLUTELY NO WARRANTY."<<endl;
    		cout<<"This is free software, and you are welcome to redistribute it"<<endl;
    		cout<<"under certain conditions; see \"COPYING\" for details."<<endl;
		cout<<"**************************************************"<<endl;
		cout<<endl;
		cout<<"Ray Copyright (C) 2010  Sébastien Boisvert, Jacques Corbeil, François Laviolette"<<endl;
		cout<<"Centre de recherche en infectiologie de l'Université Laval"<<endl;
		cout<<"Project funded by the Canadian Institutes of Health Research (Doctoral award 200902CGM-204212-172830 to S.B.)"<<endl;
 		cout<<"http://denovoassembler.sf.net/"<<endl<<endl;

		cout<<"Reference to cite: "<<endl<<endl;
		cout<<"Sébastien Boisvert, François Laviolette & Jacques Corbeil."<<endl;
		cout<<"Ray: simultaneous assembly of reads from a mix of high-throughput sequencing technologies."<<endl;
		cout<<"Journal of Computational Biology (Mary Ann Liebert, Inc. publishers, New York, U.S.A.)."<<endl;
		cout<<"November 2010, Volume 17, Issue 11, Pages 1519-1533."<<endl;
		cout<<"doi:10.1089/cmb.2009.0238"<<endl;
		cout<<"http://dx.doi.org/doi:10.1089/cmb.2009.0238"<<endl;
		cout<<endl;

		cout<<"Rank "<<getRank()<<" welcomes you to the MPI_COMM_WORLD"<<endl;
	}

	cout<<"Rank "<<getRank()<<" is running as UNIX process "<<getpid()<<" on "<<serverName<<endl;
	m_alive=true;
	m_welcomeStep=true;
	m_loadSequenceStep=false;
	m_totalLetters=0;

	MPI_Barrier(MPI_COMM_WORLD);

	m_mp.constructor(
&m_messagesHandler,
&m_library,&m_ready,
&m_verticesExtractor,
&m_edgesExtractor,
&m_sl,
			m_ed,
			&m_numberOfRanksDoneDetectingDistances,
			&m_numberOfRanksDoneSendingDistances,
			&m_parameters,
			&m_libraryIterator,
			&m_libraryIndexInitiated,
			&m_subgraph,
			&m_outboxAllocator,
			getRank(),
			&m_ed->m_EXTENSION_receivedReads,
			&m_numberOfMachinesDoneSendingEdges,
			m_fusionData,
			&m_ed->m_EXTENSION_contigs,
			&m_wordSize,
			&m_minimumCoverage,
			&m_seedCoverage,
			&m_peakCoverage,
			&m_myReads,
			&m_ed->m_EXTENSION_currentRankIsDone,
			&m_FINISH_newFusions,
		getSize(),
	&m_inboxAllocator,
	&m_persistentAllocator,
	&m_identifiers,
	&m_mode_sendDistribution,
	&m_alive,
	&(m_seedingData->m_SEEDING_receivedIngoingEdges),
	&(m_seedingData->m_SEEDING_receivedKey),
	&(m_seedingData->m_SEEDING_i),
	&m_colorSpaceMode,
	&m_FINISH_fusionOccured,
	&m_Machine_getPaths_INITIALIZED,
	&m_mode,
	&m_allPaths,
	&m_ed->m_EXTENSION_VertexAssembled_received,
	&m_ed->m_EXTENSION_numberOfRanksDone,
	&m_ed->m_EXTENSION_currentPosition,
	&m_last_value,
	&m_ed->m_EXTENSION_identifiers,
	&m_ranksDoneAttachingReads,
	&(m_seedingData->m_SEEDING_edgesReceived),
	&m_ed->m_EXTENSION_pairedRead,
	&m_ed->m_mode_EXTENSION,
	&(m_seedingData->m_SEEDING_receivedOutgoingEdges),
	&m_DISTRIBUTE_n,
	&(m_seedingData->m_SEEDING_nodes),
	&m_ed->m_EXTENSION_hasPairedReadReceived,
	&m_numberOfRanksDoneSeeding,
	&(m_seedingData->m_SEEDING_vertexKeyAndCoverageReceived),
	&(m_seedingData->m_SEEDING_receivedVertexCoverage),
	&m_ed->m_EXTENSION_readLength_received,
	&m_Machine_getPaths_DONE,
	&m_CLEAR_n,
	&m_FINISH_vertex_received,
	&m_ed->m_EXTENSION_initiated,
	&m_readyToSeed,
	&(m_seedingData->m_SEEDING_NodeInitiated),
	&m_FINISH_n,
	&m_nextReductionOccured,
	&m_ed->m_EXTENSION_hasPairedReadAnswer,
	&m_directionsAllocator,
	&m_FINISH_pathLengths,
	&m_ed->m_EXTENSION_pairedSequenceReceived,
	&m_ed->m_EXTENSION_receivedLength,
	&m_mode_send_coverage_iterator,
	&m_coverageDistribution,
	&m_FINISH_received_vertex,
	&m_ed->m_EXTENSION_read_vertex_received,
	&m_sequence_ready_machines,
	&(m_seedingData->m_SEEDING_InedgesReceived),
	&m_ed->m_EXTENSION_vertexIsAssembledResult,
	&(m_seedingData->m_SEEDING_vertexCoverageReceived),
	&m_ed->m_EXTENSION_receivedReadVertex,
	&m_numberOfMachinesReadyForEdgesDistribution,
	&m_numberOfMachinesReadyToSendDistribution,
	&m_mode_send_outgoing_edges,
	&m_mode_send_vertices_sequence_id,
	&m_mode_send_vertices,
	&m_numberOfMachinesDoneSendingVertices,
	&m_numberOfMachinesDoneSendingCoverage,
	&m_ed->m_EXTENSION_reads_received,
				&m_outbox,
	&m_sd->m_allIdentifiers,&m_oa,
	&m_numberOfRanksWithCoverageData,&m_seedExtender,
	&m_master_mode,&m_isFinalFusion);

	m_messagesHandler.constructor(getRank(),getSize());
	if(m_argc==1 or ((string)m_argv[1])=="--help"){
		if(isMaster()){
			m_aborted=true;
			showUsage();
		}
	}else{
		if(isMaster()){
			cout<<"Rank "<<getRank()<<": I am the master among "<<getSize()<<" ranks in the MPI_COMM_WORLD."<<endl;
			m_master_mode=MASTER_MODE_LOAD_CONFIG;
		}
		run();
	}

	if(isMaster() && !m_aborted){
		m_timePrinter.printElapsedTime("Collection of fusions");
		m_timePrinter.printDurations();

		cout<<endl;
		string file=m_parameters.getReceivedMessagesFile();
		const char*tmp=file.c_str();
		m_messagesHandler.writeStats(tmp);

		cout<<"Rank "<<getRank()<<" wrote "<<m_parameters.getCoverageDistributionFile()<<""<<endl;
		m_parameters.printFinalMessage();
		cout<<"Rank "<<getRank()<<" wrote "<<m_parameters.getOutputFile()<<endl;
		cout<<"Rank "<<getRank()<<" wrote "<<m_parameters.getReceivedMessagesFile()<<endl;

		cout<<endl;
		cout<<"Au revoir !"<<endl;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();
}

/*
 * this is the function that runs a lots
 *
 * it
 * 	1) receives messages
 * 	3) process message. The function that deals with a message is selected with the message's tag
 * 	4) process data, this depends on the master-mode and slave-mode states.
 * 	5) send messages
 */
void Machine::run(){
	while(isAlive()){
		receiveMessages(); 
		processMessages();
		processData();
		sendMessages();
	}
}

int Machine::getRank(){
	return m_rank;
}

/*
 * finish hyper fusions now!
 */
void Machine::finishFusions(){
	if(m_seedingData->m_SEEDING_i==(int)m_ed->m_EXTENSION_contigs.size()){
		VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
		message[0]=m_FINISH_fusionOccured;
		cout<<"Rank "<<getRank()<<" is finishing fusions "<<m_ed->m_EXTENSION_contigs.size()<<"/"<<m_ed->m_EXTENSION_contigs.size()<<" (completed)"<<endl;

	
		/*
		char number[10];
		sprintf(number,"%d",m_rank);
		string theNumber=number;
		string file="Rank_"+theNumber+".fasta";
		ofstream f(file.c_str());

		for(int i=0;i<(int)m_FINISH_newFusions.size();i++){
			string contig=convertToString(&(m_FINISH_newFusions[i]),m_wordSize);
			f<<">contig-"<<i<<" "<<contig.length()<<" nucleotides"<<endl<<addLineBreaks(contig);
		}
		f.close();
		*/

		Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,MASTER_RANK,TAG_FINISH_FUSIONS_FINISHED,getRank());
		m_outbox.push_back(aMessage);
		m_mode=MODE_DO_NOTHING;
		return;
	}
	int overlapMinimumLength=1000;
	if((int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()<overlapMinimumLength){
		#ifdef SHOW_PROGRESS
		#endif
		m_seedingData->m_SEEDING_i++;
		m_FINISH_vertex_requested=false;
		m_ed->m_EXTENSION_currentPosition=0;
		m_fusionData->m_FUSION_pathLengthRequested=false;
		m_Machine_getPaths_INITIALIZED=false;
		m_Machine_getPaths_DONE=false;
		m_checkedValidity=false;
		return;
	}
	// check if the path begins with someone else.
	
	int currentId=m_ed->m_EXTENSION_identifiers[m_seedingData->m_SEEDING_i];
	// don't do it if it is removed.

	// start threading the extension
	// as the algorithm advance on it, it stores the path positions.
	// when it reaches a choice, it will use the available path as basis.
	
	// we have the extension in m_ed->m_EXTENSION_contigs[m_SEEDING_i]
	// we get the paths with getPaths
	bool done=false;
	if(m_ed->m_EXTENSION_currentPosition<(int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()){
		if(!m_Machine_getPaths_DONE){
			getPaths(m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][m_ed->m_EXTENSION_currentPosition]);
		}else{
			vector<Direction> a;
			for(int i=0;i<(int)m_Machine_getPaths_result.size();i++){
				if(m_Machine_getPaths_result[i].getWave()!=currentId){
					a.push_back(m_Machine_getPaths_result[i]);
				}
			}
			m_FINISH_pathsForPosition.push_back(a);
			if(m_ed->m_EXTENSION_currentPosition==0){
				if(m_seedingData->m_SEEDING_i%10==0){
					cout<<"Rank "<<getRank()<<" is finishing fusions "<<m_seedingData->m_SEEDING_i+1<<"/"<<m_ed->m_EXTENSION_contigs.size()<<endl;

				}
				vector<VERTEX_TYPE> a;
				m_FINISH_newFusions.push_back(a);
				m_FINISH_vertex_requested=false;
				m_fusionData->m_FUSION_eliminated.insert(currentId);
				m_fusionData->m_FUSION_pathLengthRequested=false;
				m_checkedValidity=false;
			}
			VERTEX_TYPE vertex=m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][m_ed->m_EXTENSION_currentPosition];
			m_FINISH_newFusions[m_FINISH_newFusions.size()-1].push_back(vertex);
			m_ed->m_EXTENSION_currentPosition++;
			m_Machine_getPaths_DONE=false;
			m_Machine_getPaths_INITIALIZED=false;
		}
	}else if(!m_checkedValidity){
		// basically, directions1 contains the paths at a particular vertex in the path
		// directions2 contains the paths at another vertex in the path
		// both vertices are distanced by overlapMinimumLength, or so
		// basically, here we say we have a hit if and only if
		// there is a pair x,y with x in directions1 ad y in directions2
		// with the property that the difference of progressions are exactly overlapMinimumLength (progressions
		// are simply positions of these vertices on another path.)
		// 
		int hits=0;

		done=true;
		vector<Direction> directions1=m_FINISH_pathsForPosition[m_FINISH_pathsForPosition.size()-1];
		vector<Direction> directions2=m_FINISH_pathsForPosition[m_FINISH_pathsForPosition.size()-overlapMinimumLength];

		map<int,vector<int> > indexOnDirection2;
		
		// index the index for each wave
		for(int j=0;j<(int)directions2.size();j++){
			int waveId=directions2[j].getWave();
			if(indexOnDirection2.count(waveId)==0){
				vector<int> emptyVector;
				indexOnDirection2[waveId]=emptyVector;
			}
			indexOnDirection2[waveId].push_back(j);
		}

		// find all hits
		//
		for(int i=0;i<(int)directions1.size();i++){
			int wave1=directions1[i].getWave();
			if(indexOnDirection2.count(wave1)==0){
				continue;
			}
			vector<int> searchResults=indexOnDirection2[wave1];
			int progression1=directions1[i].getProgression();
			for(int j=0;j<(int)searchResults.size();j++){
				int index2=searchResults[j];
				int otherProgression=directions2[index2].getProgression();
				if(progression1-otherProgression+1==overlapMinimumLength){
					// this is 
					done=false;
					hits++;
					m_selectedPath=wave1;
					m_selectedPosition=progression1;
				}
			}
		}

		indexOnDirection2.clear();

		/**
 *		if there is more than one hit, they must be repeated regions. (?)
 *
 */
		if(hits>1){// we don't support that right now.
			done=true;
		}
		m_checkedValidity=true;

	}else{
		// check if it is there for at least overlapMinimumLength
		int pathId=m_selectedPath;
		int progression=m_selectedPosition;

		// only one path, just go where it goes...
		// except if it has the same number of vertices and
		// the same start and end.
		if(m_FINISH_pathLengths.count(pathId)==0){
			if(!m_fusionData->m_FUSION_pathLengthRequested){
				int rankId=pathId%MAX_NUMBER_OF_MPI_PROCESSES;
				VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(sizeof(VERTEX_TYPE));
				message[0]=pathId;
				Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,rankId,TAG_GET_PATH_LENGTH,getRank());
				m_outbox.push_back(aMessage);
				m_fusionData->m_FUSION_pathLengthRequested=true;
				m_fusionData->m_FUSION_pathLengthReceived=false;
			}else if(m_fusionData->m_FUSION_pathLengthReceived){
				m_FINISH_pathLengths[pathId]=m_fusionData->m_FUSION_receivedLength;
			}
		}else if(m_FINISH_pathLengths[pathId]!=(int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()){// avoid fusion of same length.
			int nextPosition=progression+1;
			if(nextPosition<m_FINISH_pathLengths[pathId]){
				// get the vertex
				// get its paths,
				// and continue...
				if(!m_FINISH_vertex_requested){
					int rankId=pathId%MAX_NUMBER_OF_MPI_PROCESSES;
					VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(sizeof(VERTEX_TYPE)*2);
					message[0]=pathId;
					message[1]=nextPosition;
					Message aMessage(message,2,MPI_UNSIGNED_LONG_LONG,rankId,TAG_GET_PATH_VERTEX,getRank());
					m_outbox.push_back(aMessage);
					m_FINISH_vertex_requested=true;
					m_FINISH_vertex_received=false;
				}else if(m_FINISH_vertex_received){
					if(!m_Machine_getPaths_DONE){
						getPaths(m_FINISH_received_vertex);
					}else{
						m_FINISH_pathsForPosition.push_back(m_Machine_getPaths_result);
						m_FINISH_newFusions[m_FINISH_newFusions.size()-1].push_back(m_FINISH_received_vertex);
						m_FINISH_vertex_requested=false;
						m_Machine_getPaths_INITIALIZED=false;
						m_Machine_getPaths_DONE=false;
						m_selectedPosition++;
						m_FINISH_fusionOccured=true;
					}
				}
			}else{
				#ifdef SHOW_FUSION
				cout<<"Ray says: extension-"<<m_ed->m_EXTENSION_identifiers[m_seedingData->m_SEEDING_i]<<" ("<<m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()<<" vertices) and extension-"<<pathId<<" ("<<m_FINISH_pathLengths[pathId]<<" vertices) make a fusion, result: "<<m_FINISH_newFusions[m_FINISH_newFusions.size()-1].size()<<" vertices."<<endl;
				#endif

				done=true;
			}
		}else{
			done=true;
		}
	}
	if(done){
		// there is nothing we can do.
		m_seedingData->m_SEEDING_i++;
		m_FINISH_vertex_requested=false;
		m_ed->m_EXTENSION_currentPosition=0;
		m_fusionData->m_FUSION_pathLengthRequested=false;
		m_Machine_getPaths_INITIALIZED=false;
		m_Machine_getPaths_DONE=false;
		m_checkedValidity=false;
	}

}

void Machine::makeFusions(){
	// fusion.
	// find a path that matches directly.
	// if a path is 100% included in another, but the other is longer, keep the longest.
	// if a path is 100% identical to another one, keep the one with the lowest ID
	// if a path is 100% identical to another one, but is reverse-complement, keep the one with the lowest ID
	
	int END_LENGTH=100;
	// avoid duplication of contigs.
	if(m_seedingData->m_SEEDING_i<(int)m_ed->m_EXTENSION_contigs.size()){
		if((int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()<=END_LENGTH){
			END_LENGTH=m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()-1;
		}
	}
	if(m_seedingData->m_SEEDING_i==(int)m_ed->m_EXTENSION_contigs.size()){





		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,MASTER_RANK,TAG_FUSION_DONE,getRank());
		m_outbox.push_back(aMessage);
		m_mode=MODE_DO_NOTHING;
		#ifdef SHOW_PROGRESS
		int seedIndex=m_seedingData->m_SEEDING_i-1;
		if(m_ed->m_EXTENSION_contigs.size()==0){
			seedIndex++;
		}
		cout<<"Rank "<<getRank()<<" is computing fusions "<<m_ed->m_EXTENSION_contigs.size()<<"/"<<m_ed->m_EXTENSION_contigs.size()<<" (completed)"<<endl;
		#endif
		#ifdef ASSERT
		//cout<<"Rank "<<getRank()<<" eliminated: "<<m_fusionData->m_FUSION_eliminated.size()<<endl;
		#endif
		return;
	}else if((int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()<=END_LENGTH){
		#ifdef SHOW_PROGRESS
		cout<<"No fusion for me. "<<m_seedingData->m_SEEDING_i<<" "<<m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()<<" "<<m_ed->m_EXTENSION_identifiers[m_seedingData->m_SEEDING_i]<<endl;
		#endif
		m_fusionData->m_FUSION_direct_fusionDone=false;
		m_fusionData->m_FUSION_first_done=false;
		m_fusionData->m_FUSION_paths_requested=false;
		m_seedingData->m_SEEDING_i++;
		return;
	}else if(!m_fusionData->m_FUSION_direct_fusionDone){
		int currentId=m_ed->m_EXTENSION_identifiers[m_seedingData->m_SEEDING_i];
		if(!m_fusionData->m_FUSION_first_done){
			if(!m_fusionData->m_FUSION_paths_requested){
				#ifdef SHOW_PROGRESS
				if(m_seedingData->m_SEEDING_i%10==0){
					cout<<"Rank "<<getRank()<<" is computing fusions "<<m_seedingData->m_SEEDING_i+1<<"/"<<m_ed->m_EXTENSION_contigs.size()<<endl;
				}
				#endif
				// get the paths going on the first vertex
				#ifdef ASSERT
				assert((int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()>END_LENGTH);
				#endif
				VERTEX_TYPE theVertex=m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][END_LENGTH];
				VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
				message[0]=theVertex;
				Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(theVertex),TAG_ASK_VERTEX_PATHS_SIZE,getRank());
				m_outbox.push_back(aMessage);
				m_fusionData->m_FUSION_paths_requested=true;
				m_fusionData->m_FUSION_paths_received=false;
				m_fusionData->m_FUSION_path_id=0;
				m_fusionData->m_FUSION_path_requested=false;
			}else if(m_fusionData->m_FUSION_paths_received){
				if(m_fusionData->m_FUSION_path_id<m_fusionData->m_FUSION_numberOfPaths){
					if(!m_fusionData->m_FUSION_path_requested){
						VERTEX_TYPE theVertex=m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][END_LENGTH];
						VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
						message[0]=m_fusionData->m_FUSION_path_id;
						Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(theVertex),TAG_ASK_VERTEX_PATH,getRank());
						m_outbox.push_back(aMessage);
						m_fusionData->m_FUSION_path_requested=true;
						m_fusionData->m_FUSION_path_received=false;
					}else if(m_fusionData->m_FUSION_path_received){
						m_fusionData->m_FUSION_path_id++;
						m_fusionData->m_FUSION_receivedPaths.push_back(m_fusionData->m_FUSION_receivedPath);
						m_fusionData->m_FUSION_path_requested=false;
					}
				}else{
					m_fusionData->m_FUSION_first_done=true;
					m_fusionData->m_FUSION_paths_requested=false;
					m_fusionData->m_FUSION_last_done=false;
					m_fusionData->m_FUSION_firstPaths=m_fusionData->m_FUSION_receivedPaths;
					#ifdef ASSERT
					assert(m_fusionData->m_FUSION_numberOfPaths==(int)m_fusionData->m_FUSION_firstPaths.size());
					#endif
				}
			}
		}else if(!m_fusionData->m_FUSION_last_done){
			// get the paths going on the last vertex.

			if(!m_fusionData->m_FUSION_paths_requested){
				// get the paths going on the lastvertex<
				#ifdef ASSERT
				assert((int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()>=END_LENGTH);
				#endif
				VERTEX_TYPE theVertex=m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()-END_LENGTH];
				VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
				message[0]=theVertex;
				Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(theVertex),TAG_ASK_VERTEX_PATHS_SIZE,getRank());
				m_outbox.push_back(aMessage);
				m_fusionData->m_FUSION_paths_requested=true;
				m_fusionData->m_FUSION_paths_received=false;
				m_fusionData->m_FUSION_path_id=0;
				m_fusionData->m_FUSION_path_requested=false;
			}else if(m_fusionData->m_FUSION_paths_received){
				if(m_fusionData->m_FUSION_path_id<m_fusionData->m_FUSION_numberOfPaths){
					if(!m_fusionData->m_FUSION_path_requested){
						VERTEX_TYPE theVertex=m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()-END_LENGTH];
						VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
						message[0]=m_fusionData->m_FUSION_path_id;
						Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(theVertex),TAG_ASK_VERTEX_PATH,getRank());
						m_outbox.push_back(aMessage);
						m_fusionData->m_FUSION_path_requested=true;
						m_fusionData->m_FUSION_path_received=false;
					}else if(m_fusionData->m_FUSION_path_received){
						m_fusionData->m_FUSION_path_id++;
						m_fusionData->m_FUSION_receivedPaths.push_back(m_fusionData->m_FUSION_receivedPath);
						m_fusionData->m_FUSION_path_requested=false;
					}
				}else{
					m_fusionData->m_FUSION_last_done=true;
					m_fusionData->m_FUSION_paths_requested=false;
					m_fusionData->m_FUSION_lastPaths=m_fusionData->m_FUSION_receivedPaths;
					m_fusionData->m_FUSION_matches_done=false;
					m_fusionData->m_FUSION_matches.clear();

					#ifdef ASSERT
					assert(m_fusionData->m_FUSION_numberOfPaths==(int)m_fusionData->m_FUSION_lastPaths.size());
					#endif
				}
			}


		}else if(!m_fusionData->m_FUSION_matches_done){
			m_fusionData->m_FUSION_matches_done=true;
			map<int,int> index;
			map<int,vector<int> > starts;
			map<int,vector<int> > ends;


			// extract those that are on both starting and ending vertices.
			for(int i=0;i<(int)m_fusionData->m_FUSION_firstPaths.size();i++){
				index[m_fusionData->m_FUSION_firstPaths[i].getWave()]++;
				int pathId=m_fusionData->m_FUSION_firstPaths[i].getWave();
				int progression=m_fusionData->m_FUSION_firstPaths[i].getProgression();
				starts[pathId].push_back(progression);
			}

			vector<int> matches;

			for(int i=0;i<(int)m_fusionData->m_FUSION_lastPaths.size();i++){
				index[m_fusionData->m_FUSION_lastPaths[i].getWave()]++;
				
				int pathId=m_fusionData->m_FUSION_lastPaths[i].getWave();
				int progression=m_fusionData->m_FUSION_lastPaths[i].getProgression();
				ends[pathId].push_back(progression);
			}
			

			
			for(map<int,int>::iterator i=index.begin();i!=index.end();++i){
				int otherPathId=i->first;
				if(i->second>=2 and otherPathId != currentId){
					// try to find a match with the current size.
					for(int k=0;k<(int)starts[otherPathId].size();k++){
						bool found=false;
						for(int p=0;p<(int)ends[otherPathId].size();p++){
							int observedLength=ends[otherPathId][p]-starts[otherPathId][k]+1;
							int expectedLength=m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()-2*END_LENGTH+1;
							//cout<<observedLength<<" versus "<<expectedLength<<endl;
							if(observedLength==expectedLength){
								m_fusionData->m_FUSION_matches.push_back(otherPathId);
								found=true;
								break;
							}
						}
						if(found)
							break;
					}
				}
			}
			if(m_fusionData->m_FUSION_matches.size()==0){ // no match, go next.
				m_fusionData->m_FUSION_direct_fusionDone=true;
				m_fusionData->m_FUSION_reverse_fusionDone=false;
				m_fusionData->m_FUSION_first_done=false;
				m_fusionData->m_FUSION_paths_requested=false;
			}
			m_fusionData->m_FUSION_matches_length_done=false;
			m_fusionData->m_FUSION_match_index=0;
			m_fusionData->m_FUSION_pathLengthRequested=false;
		}else if(!m_fusionData->m_FUSION_matches_length_done){
			int currentId=m_ed->m_EXTENSION_identifiers[m_seedingData->m_SEEDING_i];
			if(m_fusionData->m_FUSION_match_index==(int)m_fusionData->m_FUSION_matches.size()){// tested all matches, and nothing was found.
				m_fusionData->m_FUSION_matches_length_done=true;
			}else if(!m_fusionData->m_FUSION_pathLengthRequested){
				int uniquePathId=m_fusionData->m_FUSION_matches[m_fusionData->m_FUSION_match_index];
				int rankId=uniquePathId%MAX_NUMBER_OF_MPI_PROCESSES;
				VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(sizeof(VERTEX_TYPE));
				message[0]=uniquePathId;
				Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,rankId,TAG_GET_PATH_LENGTH,getRank());
				m_outbox.push_back(aMessage);
				m_fusionData->m_FUSION_pathLengthRequested=true;
				m_fusionData->m_FUSION_pathLengthReceived=false;
			}else if(m_fusionData->m_FUSION_pathLengthReceived){
				if(m_fusionData->m_FUSION_receivedLength==0){
				}else if(m_fusionData->m_FUSION_matches[m_fusionData->m_FUSION_match_index]<currentId and m_fusionData->m_FUSION_receivedLength == (int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()){
					m_fusionData->m_FUSION_eliminated.insert(currentId);
					m_fusionData->m_FUSION_direct_fusionDone=false;
					m_fusionData->m_FUSION_first_done=false;
					m_fusionData->m_FUSION_paths_requested=false;
					m_seedingData->m_SEEDING_i++;
				}else if(m_fusionData->m_FUSION_receivedLength>(int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size() ){
					m_fusionData->m_FUSION_eliminated.insert(currentId);
					m_fusionData->m_FUSION_direct_fusionDone=false;
					m_fusionData->m_FUSION_first_done=false;
					m_fusionData->m_FUSION_paths_requested=false;
					m_seedingData->m_SEEDING_i++;
				}
				m_fusionData->m_FUSION_match_index++;
				m_fusionData->m_FUSION_pathLengthRequested=false;
			}
		}else if(m_fusionData->m_FUSION_matches_length_done){ // no candidate found for fusion.
			m_fusionData->m_FUSION_direct_fusionDone=true;
			m_fusionData->m_FUSION_reverse_fusionDone=false;
			m_fusionData->m_FUSION_first_done=false;
			m_fusionData->m_FUSION_paths_requested=false;
		}
	}else if(!m_fusionData->m_FUSION_reverse_fusionDone){
		int currentId=m_ed->m_EXTENSION_identifiers[m_seedingData->m_SEEDING_i];
		if(!m_fusionData->m_FUSION_first_done){
			if(!m_fusionData->m_FUSION_paths_requested){
				// get the paths going on the first vertex
				VERTEX_TYPE theVertex=complementVertex(m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()-END_LENGTH],m_wordSize,m_colorSpaceMode);
				VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
				message[0]=theVertex;
				Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(theVertex),TAG_ASK_VERTEX_PATHS_SIZE,getRank());
				m_outbox.push_back(aMessage);
				m_fusionData->m_FUSION_paths_requested=true;
				m_fusionData->m_FUSION_paths_received=false;
				m_fusionData->m_FUSION_path_id=0;
				m_fusionData->m_FUSION_path_requested=false;
			}else if(m_fusionData->m_FUSION_paths_received){
				if(m_fusionData->m_FUSION_path_id<m_fusionData->m_FUSION_numberOfPaths){
					if(!m_fusionData->m_FUSION_path_requested){
						VERTEX_TYPE theVertex=complementVertex(m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()-END_LENGTH],m_wordSize,m_colorSpaceMode);
						VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
						message[0]=m_fusionData->m_FUSION_path_id;
						Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(theVertex),TAG_ASK_VERTEX_PATH,getRank());
						m_outbox.push_back(aMessage);
						m_fusionData->m_FUSION_path_requested=true;
						m_fusionData->m_FUSION_path_received=false;
					}else if(m_fusionData->m_FUSION_path_received){
						m_fusionData->m_FUSION_path_id++;
						m_fusionData->m_FUSION_receivedPaths.push_back(m_fusionData->m_FUSION_receivedPath);
						m_fusionData->m_FUSION_path_requested=false;
					}
				}else{
					m_fusionData->m_FUSION_first_done=true;
					m_fusionData->m_FUSION_paths_requested=false;
					m_fusionData->m_FUSION_last_done=false;
					m_fusionData->m_FUSION_firstPaths=m_fusionData->m_FUSION_receivedPaths;
					#ifdef ASSERT
					assert(m_fusionData->m_FUSION_numberOfPaths==(int)m_fusionData->m_FUSION_firstPaths.size());
					#endif
				}
			}
		}else if(!m_fusionData->m_FUSION_last_done){
			// get the paths going on the last vertex.

			if(!m_fusionData->m_FUSION_paths_requested){
				// get the paths going on the first vertex
				VERTEX_TYPE theVertex=complementVertex(m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][END_LENGTH],m_wordSize,m_colorSpaceMode);
				VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
				message[0]=theVertex;
				Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(theVertex),TAG_ASK_VERTEX_PATHS_SIZE,getRank());
				m_outbox.push_back(aMessage);
				m_fusionData->m_FUSION_paths_requested=true;
				m_fusionData->m_FUSION_paths_received=false;
				m_fusionData->m_FUSION_path_id=0;
				m_fusionData->m_FUSION_path_requested=false;
			}else if(m_fusionData->m_FUSION_paths_received){
				if(m_fusionData->m_FUSION_path_id<m_fusionData->m_FUSION_numberOfPaths){
					if(!m_fusionData->m_FUSION_path_requested){
						VERTEX_TYPE theVertex=complementVertex(m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][END_LENGTH],m_wordSize,m_colorSpaceMode);
						VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
						message[0]=m_fusionData->m_FUSION_path_id;
						Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(theVertex),TAG_ASK_VERTEX_PATH,getRank());
						m_outbox.push_back(aMessage);
						m_fusionData->m_FUSION_path_requested=true;
						m_fusionData->m_FUSION_path_received=false;
					}else if(m_fusionData->m_FUSION_path_received){
						m_fusionData->m_FUSION_path_id++;
						m_fusionData->m_FUSION_receivedPaths.push_back(m_fusionData->m_FUSION_receivedPath);
						m_fusionData->m_FUSION_path_requested=false;
					}
				}else{
					m_fusionData->m_FUSION_last_done=true;
					m_fusionData->m_FUSION_paths_requested=false;
					m_fusionData->m_FUSION_lastPaths=m_fusionData->m_FUSION_receivedPaths;
					m_fusionData->m_FUSION_matches_done=false;
					m_fusionData->m_FUSION_matches.clear();

					#ifdef ASSERT
					assert(m_fusionData->m_FUSION_numberOfPaths==(int)m_fusionData->m_FUSION_lastPaths.size());
					#endif
				}
			}



		}else if(!m_fusionData->m_FUSION_matches_done){
			m_fusionData->m_FUSION_matches_done=true;
			map<int,int> index;
			map<int,vector<int> > starts;
			map<int,vector<int> > ends;
			for(int i=0;i<(int)m_fusionData->m_FUSION_firstPaths.size();i++){
				index[m_fusionData->m_FUSION_firstPaths[i].getWave()]++;
				int pathId=m_fusionData->m_FUSION_firstPaths[i].getWave();
				int progression=m_fusionData->m_FUSION_firstPaths[i].getProgression();
				starts[pathId].push_back(progression);
			}
			for(int i=0;i<(int)m_fusionData->m_FUSION_lastPaths.size();i++){
				index[m_fusionData->m_FUSION_lastPaths[i].getWave()]++;
				
				int pathId=m_fusionData->m_FUSION_lastPaths[i].getWave();
				int progression=m_fusionData->m_FUSION_lastPaths[i].getProgression();
				ends[pathId].push_back(progression);
			}
			vector<int> matches;
			for(map<int,int>::iterator i=index.begin();i!=index.end();++i){
				int otherPathId=i->first;
				if(i->second>=2 and i->first != currentId){
					// try to find a match with the current size.
					for(int k=0;k<(int)starts[otherPathId].size();k++){
						bool found=false;
						for(int p=0;p<(int)ends[otherPathId].size();p++){
							int observedLength=ends[otherPathId][p]-starts[otherPathId][k]+1;
							int expectedLength=m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()-2*END_LENGTH+1;
							//cout<<observedLength<<" versus "<<expectedLength<<endl;
							if(observedLength==expectedLength){
								m_fusionData->m_FUSION_matches.push_back(otherPathId);
								found=true;
								break;
							}
						}
						if(found)
							break;
					}
				}
			}
			if(m_fusionData->m_FUSION_matches.size()==0){ // no match, go next.
				m_fusionData->m_FUSION_direct_fusionDone=false;
				m_fusionData->m_FUSION_first_done=false;
				m_fusionData->m_FUSION_paths_requested=false;
				m_seedingData->m_SEEDING_i++;
			}
			m_fusionData->m_FUSION_matches_length_done=false;
			m_fusionData->m_FUSION_match_index=0;
			m_fusionData->m_FUSION_pathLengthRequested=false;
		}else if(!m_fusionData->m_FUSION_matches_length_done){
			int currentId=m_ed->m_EXTENSION_identifiers[m_seedingData->m_SEEDING_i];
			if(m_fusionData->m_FUSION_match_index==(int)m_fusionData->m_FUSION_matches.size()){
				m_fusionData->m_FUSION_matches_length_done=true;
			}else if(!m_fusionData->m_FUSION_pathLengthRequested){
				int uniquePathId=m_fusionData->m_FUSION_matches[m_fusionData->m_FUSION_match_index];
				int rankId=uniquePathId%MAX_NUMBER_OF_MPI_PROCESSES;
				VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(sizeof(VERTEX_TYPE));
				message[0]=uniquePathId;
				Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,rankId,TAG_GET_PATH_LENGTH,getRank());
				m_outbox.push_back(aMessage);
				m_fusionData->m_FUSION_pathLengthRequested=true;
				m_fusionData->m_FUSION_pathLengthReceived=false;
			}else if(m_fusionData->m_FUSION_pathLengthReceived){
				if(m_fusionData->m_FUSION_receivedLength==0){
				}else if(m_fusionData->m_FUSION_matches[m_fusionData->m_FUSION_match_index]<currentId and m_fusionData->m_FUSION_receivedLength == (int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()){
					m_fusionData->m_FUSION_eliminated.insert(currentId);
					m_fusionData->m_FUSION_direct_fusionDone=false;
					m_fusionData->m_FUSION_first_done=false;
					m_fusionData->m_FUSION_paths_requested=false;
					m_seedingData->m_SEEDING_i++;
				}else if(m_fusionData->m_FUSION_receivedLength>(int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()){
					m_fusionData->m_FUSION_eliminated.insert(currentId);
					m_fusionData->m_FUSION_direct_fusionDone=false;
					m_fusionData->m_FUSION_first_done=false;
					m_fusionData->m_FUSION_paths_requested=false;
					m_seedingData->m_SEEDING_i++;
				}
				m_fusionData->m_FUSION_match_index++;
				m_fusionData->m_FUSION_pathLengthRequested=false;
			}
		}else if(m_fusionData->m_FUSION_matches_length_done){ // no candidate found for fusion.
			m_fusionData->m_FUSION_direct_fusionDone=false;
			m_fusionData->m_FUSION_first_done=false;
			m_fusionData->m_FUSION_paths_requested=false;
			m_seedingData->m_SEEDING_i++;
		}
	}
}

void Machine::processMessages(){
	for(int i=0;i<(int)m_inbox.size();i++){
		m_mp.processMessage((m_inbox[i]));
	}
	m_inbox.clear();
}

void Machine::sendMessages(){
	m_messagesHandler.sendMessages(&m_outbox,getRank());
}

void Machine::receiveMessages(){
	m_messagesHandler.receiveMessages(&m_inbox,&m_inboxAllocator,getRank());
}



void Machine::call_MASTER_MODE_LOAD_CONFIG(){
	if(m_argc==2){
		ifstream f(m_argv[1]);
		if(!f){
			cout<<"Rank "<<getRank()<<" invalid input file."<<endl;
			m_alive=false;
			m_aborted=true;
			f.close();
			killRanks();
			return;
		}
	}



	m_parameters.load(m_argc,m_argv);
	if(m_parameters.getError()){
		m_master_mode=MASTER_MODE_DO_NOTHING;
		m_aborted=true;
		killRanks();
		return;
	}
	if(m_parameters.useAmos()){
		// empty the file.
		cout<<"Preparing AMOS file "<<m_parameters.getAmosFile()<<endl;
		m_bubbleData->m_amos=fopen(m_parameters.getAmosFile().c_str(),"w");
	}


	VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
	message[0]=m_parameters.getWordSize();
	VERTEX_TYPE*message2=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
	message2[0]=m_parameters.getColorSpaceMode();
	for(int i=0;i<getSize();i++){
		Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,i,TAG_SET_WORD_SIZE,getRank());
		m_outbox.push_back(aMessage);
		
		Message aMessage2(message2,1,MPI_UNSIGNED_LONG_LONG,i,TAG_SET_COLOR_MODE,getRank());
		m_outbox.push_back(aMessage2);
	}
	m_master_mode=MASTER_MODE_LOAD_SEQUENCES;
}

void Machine::call_MASTER_MODE_LOAD_SEQUENCES(){
	bool res=m_sl.loadSequences(getRank(),getSize(),
	&m_outbox,
	m_disData,&m_outboxAllocator,
	&m_loadSequenceStep,
	m_bubbleData,
	&m_lastTime,
	&m_parameters,&m_master_mode
);
	if(!res){
		m_aborted=true;
		killRanks();
		m_mode=MODE_DO_NOTHING;
		m_master_mode=MASTER_MODE_DO_NOTHING;
	}
}

void Machine::call_MASTER_MODE_TRIGGER_VERTICE_DISTRIBUTION(){
	#ifdef SHOW_PROGRESS
	m_timePrinter.printElapsedTime("Distribution of sequence reads");
	cout<<endl;
	cout<<"Rank "<<getRank()<<": starting vertices distribution."<<endl;
	#else
	cout<<"\r"<<"Counting vertices"<<endl;
	#endif
	for(int i=0;i<getSize();i++){
		Message aMessage(NULL, 0, MPI_UNSIGNED_LONG_LONG,i,TAG_START_VERTICES_DISTRIBUTION,getRank());
		m_outbox.push_back(aMessage);
	}
	m_messageSentForVerticesDistribution=true;
	m_master_mode=MASTER_MODE_DO_NOTHING;
}

void Machine::call_MASTER_MODE_TRIGGER_EDGES_DISTRIBUTION(){
	#ifndef SHOW_PROGRESS
	cout<<"\r"<<"Connecting vertices"<<endl;
	#endif
	for(int i=0;i<getSize();i++){
		Message aMessage(NULL, 0, MPI_UNSIGNED_LONG_LONG,i,TAG_START_EDGES_DISTRIBUTION_ASK,getRank());
		m_outbox.push_back(aMessage);
	}
	m_startEdgeDistribution=false;
}

void Machine::call_MASTER_MODE_START_EDGES_DISTRIBUTION(){
	m_timePrinter.printElapsedTime("Calculation of coverage distribution");
	cout<<endl;
	cout<<"Rank 0 tells its friends to proceed with the distribution of edges."<<endl;
	m_numberOfMachinesReadyForEdgesDistribution=-1;
	for(int i=0;i<getSize();i++){
		Message aMessage(NULL, 0, MPI_UNSIGNED_LONG_LONG,i,TAG_START_EDGES_DISTRIBUTION,getRank());
		m_outbox.push_back(aMessage);
	}
	m_messageSentForEdgesDistribution=true;
	m_master_mode=MASTER_MODE_DO_NOTHING;
}

void Machine::call_MASTER_MODE_SEND_COVERAGE_VALUES(){
	m_numberOfMachinesDoneSendingCoverage=-1;
	string file=m_parameters.getCoverageDistributionFile();
	CoverageDistribution distribution(&m_coverageDistribution,&file);
	m_minimumCoverage=distribution.getMinimumCoverage();
	m_peakCoverage=distribution.getPeakCoverage();
	m_seedCoverage=(m_minimumCoverage+m_peakCoverage)/2;

	m_coverageDistribution.clear();


	if(m_minimumCoverage > m_peakCoverage or m_peakCoverage==m_maxCoverage){
		killRanks();
		cout<<"Error: no enrichment observed."<<endl;
		return;
	}

	// see these values to everyone.
	
	VERTEX_TYPE*buffer=(VERTEX_TYPE*)m_outboxAllocator.allocate(3*sizeof(VERTEX_TYPE));
	buffer[0]=m_minimumCoverage;
	buffer[1]=m_seedCoverage;
	buffer[2]=m_peakCoverage;
	m_numberOfRanksWithCoverageData=0;
	for(int i=0;i<getSize();i++){
		Message aMessage(buffer,3,MPI_UNSIGNED_LONG_LONG,i,TAG_SEND_COVERAGE_VALUES,getRank());
		m_outbox.push_back(aMessage);
	}
	m_master_mode=MASTER_MODE_DO_NOTHING;
}

void Machine::call_MASTER_MODE_DO_NOTHING(){
}

void Machine::call_MODE_DO_NOTHING(){
}

void Machine::call_MODE_EXTRACT_VERTICES(){
	m_verticesExtractor.process(		&m_mode_send_vertices_sequence_id,
			&m_myReads,
			&m_reverseComplementVertex,
			&m_mode_send_vertices_sequence_id_position,
			getRank(),
			&m_outbox,
			&m_mode_send_vertices,
			m_wordSize,
			m_disData,
			getSize(),
			&m_outboxAllocator,
			m_colorSpaceMode,&m_mode
		);
}

void Machine::call_MASTER_MODE_TRIGGER_EDGES(){
	m_numberOfRanksWithCoverageData=-1;
	m_numberOfMachinesReadyForEdgesDistribution=getSize();
	m_master_mode=MASTER_MODE_START_EDGES_DISTRIBUTION;
}

void Machine::call_MASTER_MODE_TRIGGER_INDEXING(){
	m_numberOfMachinesDoneSendingEdges=-9;
	m_master_mode=MASTER_MODE_DO_NOTHING;
	m_timePrinter.printElapsedTime("Distribution of edges");
	cout<<endl;
	for(int i=0;i<getSize();i++){
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_START_INDEXING_SEQUENCES,getRank());
		m_outbox.push_back(aMessage);
	}
}

void Machine::call_MASTER_MODE_PREPARE_DISTRIBUTIONS(){
	m_numberOfMachinesDoneSendingVertices=-1;
	for(int i=0;i<getSize();i++){
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG, i, TAG_PREPARE_COVERAGE_DISTRIBUTION_QUESTION,getRank());
		m_outbox.push_back(aMessage);
	}
	m_master_mode=MASTER_MODE_DO_NOTHING;
}

void Machine::call_MASTER_MODE_PREPARE_DISTRIBUTIONS_WITH_ANSWERS(){
	m_numberOfMachinesReadyToSendDistribution=-1;
	m_timePrinter.printElapsedTime("Distribution of vertices");
	cout<<endl;
	cout<<"Rank 0 computes the coverage distribution."<<endl;


	for(int i=0;i<getSize();i++){
		Message aMessage(NULL, 0, MPI_UNSIGNED_LONG_LONG, i, TAG_PREPARE_COVERAGE_DISTRIBUTION,getRank());
		m_outbox.push_back(aMessage);
	}
	m_master_mode=MASTER_MODE_DO_NOTHING;
}

void Machine::call_MASTER_MODE_PREPARE_SEEDING(){
	m_ranksDoneAttachingReads=-1;
	m_readyToSeed=getSize();
	m_master_mode=MASTER_MODE_TRIGGER_SEEDING;
}

void Machine::call_MODE_ASSEMBLE_WAVES(){
	// take each seed, and extend it in both direction using previously obtained information.
	if(m_seedingData->m_SEEDING_i==(int)m_seedingData->m_SEEDING_seeds.size()){
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,MASTER_RANK,TAG_ASSEMBLE_WAVES_DONE,getRank());
		m_outbox.push_back(aMessage);
	}else{
	}
}

void Machine::call_MODE_PERFORM_CALIBRATION(){
	int rank=rand()%getSize();
	VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
	Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,rank,TAG_CALIBRATION_MESSAGE,getRank());
	m_outbox.push_back(aMessage);
}

void Machine::call_MODE_FINISH_FUSIONS(){
	finishFusions();
}

void Machine::call_MODE_DISTRIBUTE_FUSIONS(){
	m_fusionData->distribute(m_seedingData,m_ed,m_rank,&m_outboxAllocator,&m_outbox,getSize(),&m_mode);
}

void Machine::call_MODE_SEND_DISTRIBUTION(){
	if(m_distributionOfCoverage.size()==0){
		for(int i=0;i<m_subgraph.getNumberOfTrees();i++){
			SplayTreeIterator<VERTEX_TYPE,Vertex> iterator(m_subgraph.getTree(i));
			while(iterator.hasNext()){
				int coverage=iterator.next()->getValue()->getCoverage();
				m_distributionOfCoverage[coverage]++;
			}
		}
	}

	int*data=(int*)m_outboxAllocator.allocate(sizeof(int)*2*m_maxCoverage);
	int j=0;
	data[j++]=m_distributionOfCoverage.size();
	for(map<int,VERTEX_TYPE>::iterator i=m_distributionOfCoverage.begin();i!=m_distributionOfCoverage.end();i++){
		int coverage=i->first;
		VERTEX_TYPE count=i->second;
		data[j++]=coverage;
		data[j++]=count;
	}
	Message aMessage(data,MPI_BTL_SM_EAGER_LIMIT/sizeof(VERTEX_TYPE), MPI_UNSIGNED_LONG_LONG, MASTER_RANK, TAG_COVERAGE_DATA,getRank());
	m_outbox.push_back(aMessage);

	m_distributionOfCoverage.clear();
	m_mode=MODE_DO_NOTHING;
}

void Machine::call_MODE_PROCESS_OUTGOING_EDGES(){
	m_edgesExtractor.m_wordSize=m_wordSize;
	m_edgesExtractor.processOutgoingEdges();
}

void Machine::call_MODE_PROCESS_INGOING_EDGES(){
	m_edgesExtractor.processIngoingEdges();
}

void Machine::call_MASTER_MODE_TRIGGER_SEEDING(){
	m_timePrinter.printElapsedTime("Indexing of sequence reads");
	cout<<endl;
	cout<<"Rank 0 tells other ranks to calculate their seeds."<<endl;
	m_readyToSeed=-1;
	m_numberOfRanksDoneSeeding=0;
	// tell everyone to seed now.
	for(int i=0;i<getSize();i++){
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_START_SEEDING,getRank());
		m_outbox.push_back(aMessage);
	}
	m_master_mode=MASTER_MODE_DO_NOTHING;
}

void Machine::call_MODE_START_SEEDING(){
	// assign a first vertex
	if(!m_seedingData->m_SEEDING_NodeInitiated){
		if(m_seedingData->m_SEEDING_i==(int)m_subgraph.size()-1){

			m_mode=MODE_DO_NOTHING;
			cout<<"Rank "<<getRank()<<" is creating seeds "<<m_seedingData->m_SEEDING_i+1<<"/"<<m_subgraph.size()<<" (completed)"<<endl;
			Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,MASTER_RANK,TAG_SEEDING_IS_OVER,getRank());
			m_seedingData->m_SEEDING_nodes.clear();
			m_outbox.push_back(aMessage);
		}else{
			if(m_seedingData->m_SEEDING_i % 100000 ==0){
				cout<<"Rank "<<getRank()<<" is creating seeds "<<m_seedingData->m_SEEDING_i+1<<"/"<<m_subgraph.size()<<endl;
			}
			m_seedingData->m_SEEDING_currentVertex=m_seedingData->m_SEEDING_nodes[m_seedingData->m_SEEDING_i];
			m_seedingData->m_SEEDING_first=m_seedingData->m_SEEDING_currentVertex;
			m_seedingData->m_SEEDING_testInitiated=false;
			m_seedingData->m_SEEDING_1_1_test_done=false;
			m_seedingData->m_SEEDING_i++;
			m_seedingData->m_SEEDING_NodeInitiated=true;
			m_seedingData->m_SEEDING_firstVertexTestDone=false;
		}
	// check that this node has 1 ingoing edge and 1 outgoing edge.
	}else if(!m_seedingData->m_SEEDING_firstVertexTestDone){
		if(!m_seedingData->m_SEEDING_1_1_test_done){
			do_1_1_test();
		}else{
			if(!m_seedingData->m_SEEDING_1_1_test_result){
				m_seedingData->m_SEEDING_NodeInitiated=false;// abort
			}else{
				m_seedingData->m_SEEDING_firstVertexParentTestDone=false;
				m_seedingData->m_SEEDING_firstVertexTestDone=true;
				m_seedingData->m_SEEDING_currentVertex=m_seedingData->m_SEEDING_currentParentVertex;
				m_seedingData->m_SEEDING_testInitiated=false;
				m_seedingData->m_SEEDING_1_1_test_done=false;
			}
		}

	// check that the parent does not have 1 ingoing edge and 1 outgoing edge
	}else if(!m_seedingData->m_SEEDING_firstVertexParentTestDone){
		if(!m_seedingData->m_SEEDING_1_1_test_done){
			do_1_1_test();
		}else{
			if(m_seedingData->m_SEEDING_1_1_test_result){
				m_seedingData->m_SEEDING_NodeInitiated=false;//abort
			}else{
				m_seedingData->m_SEEDING_firstVertexParentTestDone=true;
				m_seedingData->m_SEEDING_vertices.clear();
				m_seedingData->m_SEEDING_seed.clear();
				// restore original starter.
				m_seedingData->m_SEEDING_currentVertex=m_seedingData->m_SEEDING_first;
				m_seedingData->m_SEEDING_testInitiated=false;
				m_seedingData->m_SEEDING_1_1_test_done=false;
			}
		}

	// check if currentVertex has 1 ingoing edge and 1 outgoing edge, if yes, add it
	}else{
		// attempt to add m_SEEDING_currentVertex
		if(!m_seedingData->m_SEEDING_1_1_test_done){
			do_1_1_test();
		}else{
			if(m_seedingData->m_SEEDING_vertices.count(m_seedingData->m_SEEDING_currentVertex)>0){
				m_seedingData->m_SEEDING_1_1_test_result=false;
			}
			if(!m_seedingData->m_SEEDING_1_1_test_result){
				m_seedingData->m_SEEDING_NodeInitiated=false;
				int nucleotides=m_seedingData->m_SEEDING_seed.size()+m_wordSize-1;
				// only consider the long ones.
				if(nucleotides>=m_parameters.getMinimumContigLength()){
					m_seedingData->m_SEEDING_seeds.push_back(m_seedingData->m_SEEDING_seed);
				}
			}else{
				m_seedingData->m_SEEDING_seed.push_back(m_seedingData->m_SEEDING_currentVertex);
				m_seedingData->m_SEEDING_vertices.insert(m_seedingData->m_SEEDING_currentVertex);
				m_seedingData->m_SEEDING_currentVertex=m_seedingData->m_SEEDING_currentChildVertex;
				m_seedingData->m_SEEDING_testInitiated=false;
				m_seedingData->m_SEEDING_1_1_test_done=false;
			}
		}
	}

}

void Machine::call_MASTER_MODE_TRIGGER_DETECTION(){
	m_timePrinter.printElapsedTime("Computation of seeds");
	cout<<endl;
	cout<<"Rank 0 asks others to approximate library sizes."<<endl;
	m_numberOfRanksDoneSeeding=-1;
	for(int i=0;i<getSize();i++){
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_AUTOMATIC_DISTANCE_DETECTION,getRank());
		m_outbox.push_back(aMessage);
	}
	m_numberOfRanksDoneDetectingDistances=0;
	m_master_mode=MASTER_MODE_DO_NOTHING;
}

void Machine::call_MASTER_MODE_ASK_DISTANCES(){
	m_numberOfRanksDoneDetectingDistances=-1;
	m_numberOfRanksDoneSendingDistances=0;
	for(int i=0;i<getSize();i++){
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_ASK_LIBRARY_DISTANCES,getRank());
		m_outbox.push_back(aMessage);
	}
	m_master_mode=MASTER_MODE_DO_NOTHING;
}

void Machine::call_MASTER_MODE_START_UPDATING_DISTANCES(){
	m_numberOfRanksDoneSendingDistances=-1;
	m_parameters.computeAverageDistances();
	m_mode=MODE_DO_NOTHING;
	m_master_mode=MASTER_MODE_UPDATE_DISTANCES;
	m_fileId=0;
	m_sequence_idInFile=0;
	m_sequence_id=0;
}

void Machine::call_MASTER_MODE_INDEX_SEQUENCES(){
}

void Machine::call_MODE_INDEX_SEQUENCES(){
	m_si.attachReads(&m_myReads,&m_outboxAllocator,&m_outbox,&m_mode,m_wordSize,
	&m_bufferedData,m_size,m_rank,m_colorSpaceMode);
}

void Machine::call_MASTER_MODE_TRIGGER_EXTENSIONS(){
	for(int i=0;i<getSize();i++){
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_ASK_EXTENSION,getRank());
		m_outbox.push_back(aMessage);
	}
	m_master_mode=MASTER_MODE_DO_NOTHING;
}

void Machine::call_MODE_SEND_EXTENSION_DATA(){
	if(!m_ready){
		return;
	}
	if(m_seedingData->m_SEEDING_i==(int)m_ed->m_EXTENSION_contigs.size()){
		m_mode=MODE_DO_NOTHING;
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,MASTER_RANK,TAG_EXTENSION_DATA_END,getRank());
		m_outbox.push_back(aMessage);
	}else{
		if(m_fusionData->m_FUSION_eliminated.count(m_ed->m_EXTENSION_identifiers[m_seedingData->m_SEEDING_i])>0){ // skip merged paths.
			m_seedingData->m_SEEDING_i++;
			m_ed->m_EXTENSION_currentPosition=0;
		}else{
			if(m_ed->m_EXTENSION_currentPosition==0){
				VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(sizeof(VERTEX_TYPE)*1);
				int theId=m_ed->m_EXTENSION_identifiers[m_seedingData->m_SEEDING_i];
				message[0]=theId;
				Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,MASTER_RANK,TAG_EXTENSION_START,getRank());
				m_outbox.push_back(aMessage);
			}
			VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(MPI_BTL_SM_EAGER_LIMIT);

			int count=0;
			for(int i=0;i<(int)(MPI_BTL_SM_EAGER_LIMIT/sizeof(VERTEX_TYPE));i++){
				if(m_ed->m_EXTENSION_currentPosition==(int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()){
					break;
				}
				message[i+0]=m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i][m_ed->m_EXTENSION_currentPosition];
				m_ed->m_EXTENSION_currentPosition++;
				count++;
			}
			
			Message aMessage(message,count,MPI_UNSIGNED_LONG_LONG,MASTER_RANK,TAG_EXTENSION_DATA,getRank());
			m_outbox.push_back(aMessage);
			m_ready=false;
			if(m_ed->m_EXTENSION_currentPosition==(int)m_ed->m_EXTENSION_contigs[m_seedingData->m_SEEDING_i].size()){
				Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,MASTER_RANK,TAG_EXTENSION_END,getRank());
				m_outbox.push_back(aMessage);
				m_seedingData->m_SEEDING_i++;
				m_ed->m_EXTENSION_currentPosition=0;
			}
		}
	}
}

void Machine::call_MODE_FUSION(){
	makeFusions();
}

void Machine::call_MODE_AUTOMATIC_DISTANCE_DETECTION(){
	m_library.detectDistances();
}

void Machine::call_MODE_SEND_LIBRARY_DISTANCES(){
	sendLibraryDistances();
}

void Machine::call_MASTER_MODE_UPDATE_DISTANCES(){
	m_library.updateDistances();
}

void Machine::call_MASTER_MODE_TRIGGER_FUSIONS(){
	m_timePrinter.printElapsedTime("Extension of seeds");
	cout<<endl;

	// ask one at once to do the fusion
	// because otherwise it may lead to hanging of the program for unknown reasons
	m_ed->m_EXTENSION_numberOfRanksDone=-1;
	m_fusionData->m_FUSION_numberOfRanksDone=0;
	for(int i=0;i<(int)getSize();i++){// start fusion.
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_START_FUSION,getRank());
		m_outbox.push_back(aMessage);
	}
	m_fusionData->m_fusionStarted=true;
	m_master_mode=MASTER_MODE_DO_NOTHING;
}

void Machine::call_MASTER_MODE_TRIGGER_FIRST_FUSIONS(){

	m_reductionOccured=true;
	m_master_mode=MASTER_MODE_START_FUSION_CYCLE;
	m_cycleStarted=false;
	m_cycleNumber=0;
}

void Machine::call_MASTER_MODE_START_FUSION_CYCLE(){
	// the finishing is
	//
	//  * a clear cycle
	//  * a distribute cycle
	//  * a finish cycle
	//  * a clear cycle
	//  * a distribute cycle
	//  * a fusion cycle


	if(!m_cycleStarted){
		#ifdef SHOW_PROGRESS
		#endif
		m_nextReductionOccured=false;
		m_cycleStarted=true;
		m_isFinalFusion=false;
		for(int i=0;i<getSize();i++){
			Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_CLEAR_DIRECTIONS,getRank());
			m_outbox.push_back(aMessage);
		}
		//cout<<"Cycle "<<m_cycleNumber<<" sending 1) TAG_CLEAR_DIRECTIONS"<<endl;
		m_currentCycleStep=1;
		m_CLEAR_n=0;
	}else if(m_CLEAR_n==getSize() and !m_isFinalFusion and m_currentCycleStep==1){
		#ifdef SHOW_PROGRESS
		#endif
		m_currentCycleStep++;
		m_CLEAR_n=-1;

		for(int i=0;i<getSize();i++){
			Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_DISTRIBUTE_FUSIONS,getRank());
			m_outbox.push_back(aMessage);
		}
		m_DISTRIBUTE_n=0;
		//cout<<"Cycle "<<m_cycleNumber<<" sending 2) TAG_DISTRIBUTE_FUSIONS"<<endl;
	}else if(m_DISTRIBUTE_n==getSize() and !m_isFinalFusion and m_currentCycleStep==2){
		#ifdef SHOW_PROGRESS
		#endif
		m_currentCycleStep++;
		m_DISTRIBUTE_n=-1;
		m_isFinalFusion=true;
		for(int i=0;i<getSize();i++){
			Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_FINISH_FUSIONS,getRank());
			m_outbox.push_back(aMessage);
		}
		//cout<<"Cycle "<<m_cycleNumber<<" sending 3) TAG_FINISH_FUSIONS"<<endl;
		m_FINISH_n=0;
	}else if(m_FINISH_n==getSize() and m_isFinalFusion and m_currentCycleStep==3){
		#ifdef SHOW_PROGRESS
		#endif
		m_currentCycleStep++;
		for(int i=0;i<getSize();i++){
			Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_CLEAR_DIRECTIONS,getRank());
			m_outbox.push_back(aMessage);
		}
		//cout<<"Cycle "<<m_cycleNumber<<" sending 4) TAG_CLEAR_DIRECTIONS"<<endl;
		m_FINISH_n=-1;
		m_CLEAR_n=0;
	}else if(m_CLEAR_n==getSize() and m_isFinalFusion && m_currentCycleStep==4){
		m_CLEAR_n=-1;
		m_currentCycleStep++;
		#ifdef SHOW_PROGRESS
		#endif

		for(int i=0;i<getSize();i++){
			Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_DISTRIBUTE_FUSIONS,getRank());
			m_outbox.push_back(aMessage);
		}
		m_DISTRIBUTE_n=0;
		//cout<<"Cycle "<<m_cycleNumber<<" sending 5) TAG_DISTRIBUTE_FUSIONS"<<endl;

	}else if(m_DISTRIBUTE_n==getSize() and m_isFinalFusion && m_currentCycleStep==5){
		m_currentCycleStep++;
		#ifdef SHOW_PROGRESS
		cout<<"Rank 0 tells others to compute fusions."<<endl;

		#endif
		m_fusionData->m_FUSION_numberOfRanksDone=0;
		m_DISTRIBUTE_n=-1;
		for(int i=0;i<(int)getSize();i++){// start fusion.
			Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_START_FUSION,getRank());
			m_outbox.push_back(aMessage);
		}
		
	}else if(m_fusionData->m_FUSION_numberOfRanksDone==getSize() && m_isFinalFusion && m_currentCycleStep==6){
		
		m_reductionOccured=m_nextReductionOccured;
		m_fusionData->m_FUSION_numberOfRanksDone=-1;
		if(!m_reductionOccured or m_cycleNumber ==5){ 
			m_timePrinter.printElapsedTime("Computation of fusions");
			cout<<endl;
			cout<<"Rank 0 is "<<"collecting fusions"<<endl;
			m_master_mode=MASTER_MODE_ASK_EXTENSIONS;

			m_sd->m_computedTopology=false;

			m_sd->m_pathId=0;
			m_sd->m_visitedVertices.clear();
			while(!m_sd->m_verticesToVisit.empty())
				m_sd->m_verticesToVisit.pop();
			while(!m_sd->m_depthsToVisit.empty())
				m_sd->m_depthsToVisit.pop();
			m_sd->m_processedLastVertex=false;
			m_ed->m_EXTENSION_currentRankIsSet=false;
			m_ed->m_EXTENSION_rank=-1;
		}else{
			// we continue now!
			m_cycleStarted=false;
			m_cycleNumber++;
		}
	}
}

void Machine::call_MASTER_MODE_ASK_EXTENSIONS(){
	#ifndef SHOW_PROGRESS
	time_t tmp=time(NULL);
	if(tmp>m_lastTime){
		m_lastTime=tmp;
		showProgress(m_lastTime);
	}
	#endif

	// ask ranks to send their extensions.
	if(!m_ed->m_EXTENSION_currentRankIsSet){
		m_ed->m_EXTENSION_currentRankIsSet=true;
		m_ed->m_EXTENSION_currentRankIsStarted=false;
		m_ed->m_EXTENSION_rank++;
	}
	if(m_ed->m_EXTENSION_rank==getSize()){
		#ifdef SHOW_SCAFFOLDER
		#endif
		int minimumLength=500;
		if(!m_sd->m_computedTopology){ // in development.
			// for each contig path, take the last vertex, and search for other contig paths 
			// reachable from it.
			if(false and m_sd->m_pathId<(int)m_allPaths.size()){
				int currentPathId=m_identifiers[m_sd->m_pathId];
				if(!m_sd->m_processedLastVertex){
					if((int)m_allPaths[m_sd->m_pathId].size()<minimumLength){
						m_sd->m_pathId++;
						return;
					}
				
					#ifdef SHOW_SCAFFOLDER
					//cout<<"push last vertex."<<endl;
					#endif
					m_sd->m_processedLastVertex=true;
					int theLength=m_allPaths[m_sd->m_pathId].size();
					VERTEX_TYPE lastVertex=m_allPaths[m_sd->m_pathId][theLength-1];
					#ifdef SHOW_SCAFFOLDER
					cout<<"contig-"<<currentPathId<<" Last="<<idToWord(lastVertex,m_wordSize)<<" "<<theLength<<" vertices"<<endl;
	
					#endif
					m_sd->m_verticesToVisit.push(lastVertex);
					m_sd->m_depthsToVisit.push(0);
					m_seedingData->m_SEEDING_edgesRequested=false;
					m_sd->m_visitedVertices.clear();
				}else if(!m_sd->m_verticesToVisit.empty()){
					VERTEX_TYPE theVertex=m_sd->m_verticesToVisit.top();
					int theDepth=m_sd->m_depthsToVisit.top();
					if(!m_seedingData->m_SEEDING_edgesRequested){
						#ifdef SHOW_SCAFFOLDER
						//cout<<"Asking for arcs. "<<theVertex<<endl;
						#endif
						m_seedingData->m_SEEDING_edgesReceived=false;
						m_seedingData->m_SEEDING_edgesRequested=true;
						VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
						message[0]=(VERTEX_TYPE)theVertex;
						Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(message[0]),TAG_REQUEST_VERTEX_OUTGOING_EDGES,getRank());
						m_outbox.push_back(aMessage);
						m_Machine_getPaths_DONE=false;
						m_Machine_getPaths_INITIALIZED=false;
						m_Machine_getPaths_result.clear();
						m_sd->m_visitedVertices.insert(theVertex);
					}else if(m_seedingData->m_SEEDING_edgesReceived){
						if(!m_Machine_getPaths_DONE){
							getPaths(theVertex);
						}else{
							vector<Direction> nextPaths;
							#ifdef SHOW_SCAFFOLDER
							if(nextPaths.size()>0){
								cout<<"We have "<<nextPaths.size()<<" paths with "<<idToWord(theVertex,m_wordSize)<<endl;	
							}
							#endif
							for(int i=0;i<(int)m_Machine_getPaths_result.size();i++){
								int pathId=m_Machine_getPaths_result[i].getWave();
								if(pathId==currentPathId)
									continue;
								// this one is discarded.
								if(m_sd->m_allIdentifiers.count(pathId)==0){
									continue;
								}
								// not at the front.
								if(m_Machine_getPaths_result[i].getProgression()>0)
									continue;
								int index=m_sd->m_allIdentifiers[pathId];

								// too small to be relevant.
								if((int)m_allPaths[index].size()<minimumLength)
									continue;
								
								#ifdef SHOW_SCAFFOLDER
								#endif
								nextPaths.push_back(m_Machine_getPaths_result[i]);
							}

							m_sd->m_verticesToVisit.pop();
							m_sd->m_depthsToVisit.pop();
							m_seedingData->m_SEEDING_edgesRequested=false;

							if(nextPaths.size()>0){// we found a path
								for(int i=0;i<(int)nextPaths.size();i++){
									cout<<"contig-"<<m_identifiers[m_sd->m_pathId]<<" -> "<<"contig-"<<nextPaths[i].getWave()<<" ("<<theDepth<<","<<nextPaths[i].getProgression()<<") via "<<idToWord(theVertex,m_wordSize)<<endl;
									#ifdef ASSERT
									assert(m_sd->m_allIdentifiers.count(nextPaths[i].getWave())>0);
									#endif
								}
							}else{// continue the visit.
								for(int i=0;i<(int)m_seedingData->m_SEEDING_receivedOutgoingEdges.size();i++){
									VERTEX_TYPE newVertex=m_seedingData->m_SEEDING_receivedOutgoingEdges[i];
									if(m_sd->m_visitedVertices.count(newVertex)>0)
										continue;
									int d=theDepth+1;
									if(d>3000)
										continue;
									m_sd->m_verticesToVisit.push(newVertex);
									m_sd->m_depthsToVisit.push(d);
								}
							}

						}
					}
				}else{
					m_sd->m_processedLastVertex=false;
					m_sd->m_pathId++;
					#ifdef SHOW_SCAFFOLDER
					//cout<<"Processing next."<<endl;
					#endif
				}
			}else{
				m_sd->m_computedTopology=true;
			}
			return;
		}

		m_master_mode=MASTER_MODE_DO_NOTHING;

		int totalLength=0;
		
		#ifdef ASSERT
		assert(m_allPaths.size()==m_identifiers.size());
		#endif
		ofstream f(m_parameters.getOutputFile().c_str());
		for(int i=0;i<(int)m_allPaths.size();i++){
			string contig=convertToString(&(m_allPaths[i]),m_wordSize);
			#ifdef ASSERT
			assert(i<(int)m_identifiers.size());
			#endif
			int id=m_identifiers[i];
			#ifdef ASSERT
			int theRank=id%MAX_NUMBER_OF_MPI_PROCESSES;
			assert(theRank<getSize());
			#endif
			f<<">contig-"<<id<<" "<<contig.length()<<" nucleotides"<<endl<<addLineBreaks(contig);
			totalLength+=contig.length();
		}
		f.close();
		#ifdef SHOW_PROGRESS
		#else
		cout<<"\r"<<"              "<<endl<<"Writing "<<m_parameters.getOutputFile()<<endl;
		#endif
		cout<<endl<<"Rank 0: "<<m_allPaths.size()<<" contigs/"<<totalLength<<" nucleotides"<<endl;
		if(m_parameters.useAmos()){
			m_master_mode=MASTER_MODE_AMOS;
			m_seedingData->m_SEEDING_i=0;
			m_mode_send_vertices_sequence_id_position=0;
			m_ed->m_EXTENSION_reads_requested=false;
			cout<<"\rCompleting "<<m_parameters.getAmosFile()<<endl;
		}else{// we are done.
			killRanks();
		}
		
	}else if(!m_ed->m_EXTENSION_currentRankIsStarted){
		m_ed->m_EXTENSION_currentRankIsStarted=true;
		#ifdef SHOW_PROGRESS
		cout<<"Rank "<<getRank()<<" asks "<<m_ed->m_EXTENSION_rank<<" its fusions"<<endl;
		#endif
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,m_ed->m_EXTENSION_rank,TAG_ASK_EXTENSION_DATA,getRank());
		m_outbox.push_back(aMessage);
		m_ed->m_EXTENSION_currentRankIsDone=false;
	}else if(m_ed->m_EXTENSION_currentRankIsDone){
		m_ed->m_EXTENSION_currentRankIsSet=false;
	}
}

void Machine::call_MASTER_MODE_AMOS(){
	// in development.
	/*
	* use m_allPaths and m_identifiers
	*
	* iterators: m_SEEDING_i: for the current contig
	*            m_mode_send_vertices_sequence_id_position: for the current position in the current contig.
	*/
	if(m_seedingData->m_SEEDING_i==(int)m_allPaths.size()){// all contigs are processed
		killRanks();
		m_master_mode=MASTER_MODE_DO_NOTHING;
		fclose(m_bubbleData->m_amos);
	}else if(m_mode_send_vertices_sequence_id_position==(int)m_allPaths[m_seedingData->m_SEEDING_i].size()){// iterate over the next one
		m_seedingData->m_SEEDING_i++;
		m_mode_send_vertices_sequence_id_position=0;
		m_ed->m_EXTENSION_reads_requested=false;
		
		FILE*fp=m_bubbleData->m_amos;
		fprintf(fp,"}\n");
	}else{
		if(!m_ed->m_EXTENSION_reads_requested){
			if(m_mode_send_vertices_sequence_id_position==0){
				FILE*fp=m_bubbleData->m_amos;
				string seq=convertToString(&(m_allPaths[m_seedingData->m_SEEDING_i]),m_wordSize);
				char*qlt=(char*)__Malloc(seq.length()+1);
				strcpy(qlt,seq.c_str());
				for(int i=0;i<(int)strlen(qlt);i++)
					qlt[i]='D';
				fprintf(fp,"{CTG\niid:%i\neid:contig-%i\ncom:\nRay\n.\nseq:\n%s\n.\nqlt:\n%s\n.\n",
					m_seedingData->m_SEEDING_i+1,
					m_identifiers[m_seedingData->m_SEEDING_i],
					seq.c_str(),
					qlt
					);
				__Free(qlt);
			}

			m_ed->m_EXTENSION_reads_requested=true;
			m_ed->m_EXTENSION_reads_received=false;
			VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
			message[0]=m_allPaths[m_seedingData->m_SEEDING_i][m_mode_send_vertices_sequence_id_position];
			Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(message[0]),TAG_REQUEST_READS,getRank());
			m_outbox.push_back(aMessage);

			// iterator on reads
			m_fusionData->m_FUSION_path_id=0;
			m_ed->m_EXTENSION_readLength_requested=false;
		}else if(m_ed->m_EXTENSION_reads_received){
			if(m_fusionData->m_FUSION_path_id<(int)m_ed->m_EXTENSION_receivedReads.size()){
				int readRank=m_ed->m_EXTENSION_receivedReads[m_fusionData->m_FUSION_path_id].getRank();
				char strand=m_ed->m_EXTENSION_receivedReads[m_fusionData->m_FUSION_path_id].getStrand();
				int idOnRank=m_ed->m_EXTENSION_receivedReads[m_fusionData->m_FUSION_path_id].getReadIndex();
				if(!m_ed->m_EXTENSION_readLength_requested){
					m_ed->m_EXTENSION_readLength_requested=true;
					m_ed->m_EXTENSION_readLength_received=false;
					VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
					message[0]=idOnRank;
					Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,readRank,TAG_ASK_READ_LENGTH,getRank());
					m_outbox.push_back(aMessage);
				}else if(m_ed->m_EXTENSION_readLength_received){
					int readLength=m_ed->m_EXTENSION_receivedLength;
					int globalIdentifier=idOnRank*getSize()+readRank;
					FILE*fp=m_bubbleData->m_amos;
					int start=0;
					int theEnd=readLength-1;
					int offset=m_mode_send_vertices_sequence_id_position;
					if(strand=='R'){
						int t=start;
						start=theEnd;
						theEnd=t;
						offset++;
					}
					fprintf(fp,"{TLE\nsrc:%i\noff:%i\nclr:%i,%i\n}\n",globalIdentifier+1,offset,
						start,theEnd);
		
					// increment to get the next read.
					m_fusionData->m_FUSION_path_id++;
					m_ed->m_EXTENSION_readLength_requested=false;
				}
			}else{
				// continue.
				m_mode_send_vertices_sequence_id_position++;
				m_ed->m_EXTENSION_reads_requested=false;
			}
		}

	}
}

void Machine::call_MODE_EXTENSION(){
	int maxCoverage=m_maxCoverage;
	m_seedExtender.extendSeeds(&(m_seedingData->m_SEEDING_seeds),m_ed,getRank(),&m_outbox,&(m_seedingData->m_SEEDING_currentVertex),
	m_fusionData,&m_outboxAllocator,&(m_seedingData->m_SEEDING_edgesRequested),&(m_seedingData->m_SEEDING_outgoingEdgeIndex),
	&m_last_value,&(m_seedingData->m_SEEDING_vertexCoverageRequested),m_wordSize,&m_colorSpaceMode,getSize(),&(m_seedingData->m_SEEDING_vertexCoverageReceived),
	&(m_seedingData->m_SEEDING_receivedVertexCoverage),&m_repeatedLength,&maxCoverage,&(m_seedingData->m_SEEDING_receivedOutgoingEdges),&m_c,
	m_cd,m_bubbleData,m_dfsData,
m_minimumCoverage,&m_oa,&(m_seedingData->m_SEEDING_edgesReceived),&m_mode);
}

void Machine::call_MASTER_MODE_ASSEMBLE_WAVES(){
	// ask ranks to send their extensions.
	if(!m_ed->m_EXTENSION_currentRankIsSet){
		m_ed->m_EXTENSION_currentRankIsSet=true;
		m_ed->m_EXTENSION_currentRankIsStarted=false;
		m_ed->m_EXTENSION_rank++;
	}
	if(m_ed->m_EXTENSION_rank==getSize()){
		m_master_mode=MASTER_MODE_DO_NOTHING;
		cout<<"Rank "<<getRank()<<" contigs computed."<<endl;
		killRanks();
	}else if(!m_ed->m_EXTENSION_currentRankIsStarted){
		m_ed->m_EXTENSION_currentRankIsStarted=true;
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,m_ed->m_EXTENSION_rank,TAG_ASSEMBLE_WAVES,getRank());
		m_outbox.push_back(aMessage);
		m_ed->m_EXTENSION_currentRankIsDone=false;
	}else if(m_ed->m_EXTENSION_currentRankIsDone){
		m_ed->m_EXTENSION_currentRankIsSet=false;
	}
}

void Machine::processData(){
	MachineMethod masterMethod=m_master_methods[m_master_mode];
	(this->*masterMethod)();
	MachineMethod slaveMethod=m_slave_methods[m_mode];
	(this->*slaveMethod)();
}

/*
 * check if (m_SEEDING_currentRank,m_SEEDING_currentPointer) has
 * 1 ingoing edge and 1 outgoing edge
 *
 * before entering the first call, m_SEEDING_testInitiated and m_SEEDING_1_1_test_done must be false
 *
 * outputs:
 *
 *  m_SEEDING_1_1_test_done
 *  m_SEEDING_currentChildVertex
 *  m_SEEDING_currentChildRank
 *  m_SEEDING_currentChildPointer
 *  m_SEEDING_currentParentRank
 *  m_SEEDING_currentParentPointer
 *
 *
 *  internals:
 *
 *  m_SEEDING_InedgesRequested
 *  m_SEEDING_InedgesReceived
 *  m_SEEDING_Inedge
 *  m_SEEDING_edgesRequested
 *  m_SEEDING_edgesReceived
 */
void Machine::do_1_1_test(){
	if(m_seedingData->m_SEEDING_1_1_test_done){
		return;
	}else if(!m_seedingData->m_SEEDING_testInitiated){
		m_seedingData->m_SEEDING_testInitiated=true;
		m_seedingData->m_SEEDING_ingoingEdgesDone=false;
		m_seedingData->m_SEEDING_InedgesRequested=false;
	}else if(!m_seedingData->m_SEEDING_ingoingEdgesDone){
		if(!m_seedingData->m_SEEDING_InedgesRequested){
			VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
			message[0]=(VERTEX_TYPE)m_seedingData->m_SEEDING_currentVertex;
			Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(m_seedingData->m_SEEDING_currentVertex),TAG_REQUEST_VERTEX_INGOING_EDGES,getRank());
			m_outbox.push_back(aMessage);
			m_seedingData->m_SEEDING_numberOfIngoingEdges=0;
			m_seedingData->m_SEEDING_numberOfIngoingEdgesWithSeedCoverage=0;
			m_seedingData->m_SEEDING_vertexCoverageRequested=false;
			m_seedingData->m_SEEDING_InedgesReceived=false;
			m_seedingData->m_SEEDING_InedgesRequested=true;
			m_seedingData->m_SEEDING_ingoingEdgeIndex=0;
		}else if(m_seedingData->m_SEEDING_InedgesReceived){
			if(m_seedingData->m_SEEDING_ingoingEdgeIndex<(int)m_seedingData->m_SEEDING_receivedIngoingEdges.size()){
				if(!m_seedingData->m_SEEDING_vertexCoverageRequested){
					VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
					message[0]=(VERTEX_TYPE)m_seedingData->m_SEEDING_receivedIngoingEdges[m_seedingData->m_SEEDING_ingoingEdgeIndex];
					Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(message[0]),TAG_REQUEST_VERTEX_COVERAGE,getRank());
					m_outbox.push_back(aMessage);
					m_seedingData->m_SEEDING_vertexCoverageRequested=true;
					m_seedingData->m_SEEDING_vertexCoverageReceived=false;
					m_seedingData->m_SEEDING_receivedVertexCoverage=-1;
				}else if(m_seedingData->m_SEEDING_vertexCoverageReceived){
					if(m_seedingData->m_SEEDING_receivedIngoingEdges.size()==1){//there is only one anyway
						m_seedingData->m_SEEDING_currentParentVertex=m_seedingData->m_SEEDING_receivedIngoingEdges[m_seedingData->m_SEEDING_ingoingEdgeIndex];
					}
					if(m_seedingData->m_SEEDING_receivedVertexCoverage>=m_seedCoverage){
						m_seedingData->m_SEEDING_currentParentVertex=m_seedingData->m_SEEDING_receivedIngoingEdges[m_seedingData->m_SEEDING_ingoingEdgeIndex];
						m_seedingData->m_SEEDING_numberOfIngoingEdgesWithSeedCoverage++;
					}
					m_seedingData->m_SEEDING_ingoingEdgeIndex++;
					m_seedingData->m_SEEDING_numberOfIngoingEdges++;
					m_seedingData->m_SEEDING_vertexCoverageRequested=false;
				}
			}else{// done analyzing ingoing edges.
				m_seedingData->m_SEEDING_ingoingEdgesDone=true;
				m_seedingData->m_SEEDING_outgoingEdgesDone=false;
				m_seedingData->m_SEEDING_edgesRequested=false;
			}
		}
	}else if(!m_seedingData->m_SEEDING_outgoingEdgesDone){
		if(!m_seedingData->m_SEEDING_edgesRequested){
			VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
			message[0]=(VERTEX_TYPE)m_seedingData->m_SEEDING_currentVertex;
			Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(m_seedingData->m_SEEDING_currentVertex),TAG_REQUEST_VERTEX_OUTGOING_EDGES,getRank());
			m_outbox.push_back(aMessage);
			m_seedingData->m_SEEDING_edgesRequested=true;
			m_seedingData->m_SEEDING_numberOfOutgoingEdges=0;
			m_seedingData->m_SEEDING_numberOfOutgoingEdgesWithSeedCoverage=0;
			m_seedingData->m_SEEDING_vertexCoverageRequested=false;
			m_seedingData->m_SEEDING_edgesReceived=false;
			m_seedingData->m_SEEDING_outgoingEdgeIndex=0;
		}else if(m_seedingData->m_SEEDING_edgesReceived){
			if(m_seedingData->m_SEEDING_outgoingEdgeIndex<(int)m_seedingData->m_SEEDING_receivedOutgoingEdges.size()){
				// TODO: don't check the coverage if there is only one
				if(!m_seedingData->m_SEEDING_vertexCoverageRequested){
					VERTEX_TYPE*message=(VERTEX_TYPE*)m_outboxAllocator.allocate(1*sizeof(VERTEX_TYPE));
					message[0]=(VERTEX_TYPE)m_seedingData->m_SEEDING_receivedOutgoingEdges[m_seedingData->m_SEEDING_outgoingEdgeIndex];
					Message aMessage(message,1,MPI_UNSIGNED_LONG_LONG,vertexRank(message[0]),TAG_REQUEST_VERTEX_COVERAGE,getRank());
					m_outbox.push_back(aMessage);
					m_seedingData->m_SEEDING_vertexCoverageRequested=true;
					m_seedingData->m_SEEDING_vertexCoverageReceived=false;
					m_seedingData->m_SEEDING_receivedVertexCoverage=-1;
				}else if(m_seedingData->m_SEEDING_vertexCoverageReceived){
					if(m_seedingData->m_SEEDING_receivedOutgoingEdges.size()==1){//there is only one anyway
						m_seedingData->m_SEEDING_currentChildVertex=m_seedingData->m_SEEDING_receivedOutgoingEdges[m_seedingData->m_SEEDING_outgoingEdgeIndex];
					}
					if(m_seedingData->m_SEEDING_receivedVertexCoverage>=m_seedCoverage){
						m_seedingData->m_SEEDING_currentChildVertex=m_seedingData->m_SEEDING_receivedOutgoingEdges[m_seedingData->m_SEEDING_outgoingEdgeIndex];
						m_seedingData->m_SEEDING_numberOfOutgoingEdgesWithSeedCoverage++;
					}
					m_seedingData->m_SEEDING_outgoingEdgeIndex++;
					m_seedingData->m_SEEDING_numberOfOutgoingEdges++;
					m_seedingData->m_SEEDING_vertexCoverageRequested=false;
				}
			}else{// done analyzing ingoing edges.
				m_seedingData->m_SEEDING_outgoingEdgesDone=true;
			}
		}


	}else{
		m_seedingData->m_SEEDING_1_1_test_done=true;
		m_seedingData->m_SEEDING_1_1_test_result=(m_seedingData->m_SEEDING_numberOfIngoingEdgesWithSeedCoverage==1)and
			(m_seedingData->m_SEEDING_numberOfOutgoingEdgesWithSeedCoverage==1);
	}
}

void Machine::killRanks(){
	for(int i=getSize()-1;i>=0;i--){
		Message aMessage(NULL,0,MPI_UNSIGNED_LONG_LONG,i,TAG_GOOD_JOB_SEE_YOU_SOON,getRank());
		m_outbox.push_back(aMessage);
	}
}

bool Machine::isMaster(){
	return getRank()==MASTER_RANK;
}





int Machine::getSize(){
	return m_size;
}



bool Machine::isAlive(){
	return m_alive;
}

void Machine::printStatus(){
	cout<<"********"<<endl;
	cout<<"Rank: "<<getRank()<<endl;
	cout<<"Reads: "<<m_myReads.size()<<endl;
	cout<<"Inbox: "<<m_inbox.size()<<endl;
	cout<<"Outbox: "<<m_outbox.size()<<endl;
}



Machine::~Machine(){
	// do nothing.
	delete m_dfsData;
	delete m_bubbleData;
	m_dfsData=NULL;
	m_bubbleData=NULL;
}


int Machine::vertexRank(VERTEX_TYPE a){
	return hash_VERTEX_TYPE(a)%(getSize());
}


