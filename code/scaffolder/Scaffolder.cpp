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

#include <iostream>
#include <assembler/ReadFetcher.h>
#include <scaffolder/Scaffolder.h>
#include <communication/Message.h>
#include <algorithm> /* for sort */
#include <vector>
#include <fstream>
#include <scaffolder/ScaffoldingVertex.h>
#include <scaffolder/ScaffoldingEdge.h>
#include <scaffolder/ScaffoldingAlgorithm.h>
#include <sstream>
#include <core/statistics.h>
#include <assert.h>
#include <math.h> /* for sqrt */
using namespace std;

void Scaffolder::addMasterLink(SummarizedLink*a){
	m_masterLinks.push_back(*a);
}

void Scaffolder::addMasterContig(uint64_t name,int length){
	m_masterContigs.push_back(name);
	m_masterLengths.push_back(length);
}

void Scaffolder::solve(){

/*
 *  The value of the minimum number of raw links is completely arbitrary.
 *  However, it reduces significantly the number of scaffolds.
 *  Therefore, its usage is warranted.
 */
	int minimumNumberOfRawLinks=3;

	map<uint64_t,map<char,map<uint64_t,map<char,vector<int> > > > > keys;
	for(int i=0;i<(int)m_masterLinks.size();i++){
		uint64_t leftContig=m_masterLinks[i].getLeftContig();
		char leftStrand=m_masterLinks[i].getLeftStrand();
		uint64_t rightContig=m_masterLinks[i].getRightContig();
		char rightStrand=m_masterLinks[i].getRightStrand();
		int average=m_masterLinks[i].getAverage();
		int number=m_masterLinks[i].getCount();
		int standardDeviation=m_masterLinks[i].getStandardDeviation();

		if(number<minimumNumberOfRawLinks)
			continue;

		keys[leftContig][leftStrand][rightContig][rightStrand].push_back(number);
		keys[leftContig][leftStrand][rightContig][rightStrand].push_back(average);
		keys[leftContig][leftStrand][rightContig][rightStrand].push_back(standardDeviation);
	}


	ostringstream linkFile;
	linkFile<<m_parameters->getPrefix()<<"ScaffoldLinks.txt";
	ofstream f(linkFile.str().c_str());

	/* Prototype */
	vector<ScaffoldingVertex> scaffoldingVertices;
	vector<ScaffoldingEdge> scaffoldingEdges;

	for(map<uint64_t,map<char,map<uint64_t,map<char,vector<int> > > > >::iterator i=
		keys.begin();i!=keys.end();i++){
		uint64_t leftContig=i->first;
		for(map<char,map<uint64_t,map<char,vector<int> > > >::iterator j=i->second.begin();
			j!=i->second.end();j++){
			char leftStrand=j->first;
			for(map<uint64_t,map<char,vector<int> > >::iterator k=j->second.begin();
				k!=j->second.end();k++){
				uint64_t rightContig=k->first;
				for(map<char,vector<int> >::iterator l=k->second.begin();
					l!=k->second.end();l++){
					char rightStrand=l->first;
					int sum=0;
					int n=0;
					int pos=0;

					vector<int> countValues;
					vector<int> averageValues;
					vector<int> standardDeviationValues;

					for(vector<int>::iterator m=l->second.begin();m!=l->second.end();m++){
						/* +0 is average, +1 is the count */
						if(pos%3==1){
							int average=*m;
							sum+=average;
							n++;
							averageValues.push_back(average);
						}else if(pos%3==0){
							int count=*m;
							countValues.push_back(count);
						}else if(pos%3==2){
							int standardDeviation=*m;
							standardDeviationValues.push_back(standardDeviation);
						}
						pos++;
					}

					/* we want contig A to reach contig B and
						we want contig B to reach contig A */
					if(n==2){
						int average=sum/n;
						f<<"contig-"<<leftContig<<"\t"<<leftStrand<<"\tcontig-"<<rightContig<<"\t"<<rightStrand<<"\t"<<average;

						#ifdef ASSERT
						assert(averageValues.size() == countValues.size());
						assert(standardDeviationValues.size() == countValues.size());
						#endif

						f<<"	"<<averageValues.size();
						for(int summarizedLinkIterator=0;
							summarizedLinkIterator<(int)averageValues.size();
							summarizedLinkIterator++){
							f<<"\t"<<summarizedLinkIterator;
							f<<"	"<<countValues[summarizedLinkIterator];
							f<<"	"<<averageValues[summarizedLinkIterator];
							f<<"	"<<standardDeviationValues[summarizedLinkIterator];
						}
						f<<endl;


						ScaffoldingEdge scaffoldingEdge(leftContig,leftStrand,rightContig,rightStrand,average,averageValues[0],countValues[0],standardDeviationValues[0],
averageValues[1],countValues[1],standardDeviationValues[1]);

						scaffoldingEdges.push_back(scaffoldingEdge);
					}
				}
			}
		}
	}
	f.close();

	// write contig list
	ostringstream contigList;
	contigList<<m_parameters->getPrefix()<<"ContigLengths.txt";
	ofstream f2(contigList.str().c_str());
	for(int i=0;i<(int)m_masterContigs.size();i++){
		int length=m_masterLengths[i]+m_parameters->getWordSize()-1;
		f2<<"contig-"<<m_masterContigs[i]<<"\t"<<length<<endl;
		m_allContigLengths.push_back(length);
		ScaffoldingVertex scaffoldingVertex(m_masterContigs[i],length);
		scaffoldingVertices.push_back(scaffoldingVertex);
	}
	f2.close();

	/* run the greedy solver */
	ScaffoldingAlgorithm solver;
	solver.setVertices(&scaffoldingVertices);
	solver.setEdges(&scaffoldingEdges);
	solver.solve(&m_scaffoldContigs,&m_scaffoldStrands,&m_scaffoldGaps);

	for(int i=0;i<(int)m_masterLengths.size();i++){
		m_contigLengths[m_masterContigs[i]]=m_masterLengths[i];
	}

	// write scaffold list
	ostringstream scaffoldList;
	scaffoldList<<m_parameters->getPrefix()<<"ScaffoldComponents.txt";
	ostringstream scaffoldLengths;
	scaffoldLengths<<m_parameters->getPrefix()<<"ScaffoldLengths.txt";
	ofstream f3(scaffoldLengths.str().c_str());
	ofstream f4(scaffoldList.str().c_str());
	for(int i=0;i<(int)m_scaffoldContigs.size();i++){
		int scaffoldName=i;
		int length=0;
		for(int j=0;j<(int)m_scaffoldContigs[i].size();j++){
			uint64_t contigName=m_scaffoldContigs[i][j];
			char contigStrand=m_scaffoldStrands[i][j];
			int theLength=m_contigLengths[contigName]+m_parameters->getWordSize()-1;
			f4<<"scaffold-"<<scaffoldName<<"\t"<<"contig-"<<contigName<<"\t"<<contigStrand<<"\t"<<theLength<<endl;
			length+=theLength;
			if(j!=(int)m_scaffoldContigs[i].size()-1){
				int theLength=m_scaffoldGaps[i][j];
				f4<<"scaffold-"<<scaffoldName<<"\tgap\t-\t"<<theLength<<endl;
				length+=theLength;
			}
		}
		f3<<"scaffold-"<<scaffoldName<<"\t"<<length<<endl;
		f4<<endl;
		m_allScaffoldLengths.push_back(length);
	}
	f2.close();
	f3.close();
}

void Scaffolder::printFinalMessage(){
	ostringstream outputStat;
	outputStat<<m_parameters->getPrefix()<<"OutputNumbers.txt";
	ofstream f5(outputStat.str().c_str());

	printInStream(&cout);
	printInStream(&f5);

	f5.close();
}

void Scaffolder::computeStatistics(vector<int>*lengths,int minimumLength,ostream*outputStream){
	vector<int> accepted;
	uint64_t totalLength=0;
	for(int i=0;i<(int)lengths->size();i++){
		int contigLength=lengths->at(i);
		if(contigLength>=minimumLength){
			accepted.push_back(contigLength);
			totalLength+=contigLength;
		}
	}
	if(accepted.size()==0)
		return;

	int average=totalLength/accepted.size();

	(*outputStream)<<" Number: "<<accepted.size()<<endl;
	(*outputStream)<<" Total length: "<<totalLength<<endl;
	(*outputStream)<<" Average: "<<average<<endl;

	sort(accepted.begin(), accepted.end());

	uint64_t sumOfLengths=0;
	int i=0;
	sumOfLengths += accepted[i];
	while(sumOfLengths < totalLength/2){
		i++;
		if(i<(int)accepted.size())
			sumOfLengths += accepted[i];
		else
			break;
	}

	#ifdef ASSERT
	assert(i<(int)accepted.size());
	#endif
	int n50=accepted[i];

	(*outputStream)<<" N50: "<<n50<<endl;

	int median=0;
	if(accepted.size()%2==0){
		median=(accepted[accepted.size()/2] + accepted[accepted.size()/2 -1]) / 2;
	}else{
		median=accepted[accepted.size()/2];
	}
	(*outputStream)<<" Median: "<<median<<endl;

	(*outputStream)<<" Largest: "<<accepted[accepted.size()-1]<<endl;
}

void Scaffolder::printInStream(ostream*outputStream){
	(*outputStream)<<"Contigs >= "<<m_parameters->getMinimumContigLength()<<" nt"<<endl;
	computeStatistics(&m_allContigLengths,m_parameters->getMinimumContigLength(),outputStream);
	(*outputStream)<<"Contigs >= "<<m_parameters->getLargeContigThreshold()<<" nt"<<endl;
	computeStatistics(&m_allContigLengths,m_parameters->getLargeContigThreshold(),outputStream);
	(*outputStream)<<"Scaffolds >= "<<m_parameters->getMinimumContigLength()<<" nt"<<endl;
	computeStatistics(&m_allScaffoldLengths,m_parameters->getMinimumContigLength(),outputStream);
	(*outputStream)<<"Scaffolds >= "<<m_parameters->getLargeContigThreshold()<<" nt"<<endl;
	computeStatistics(&m_allScaffoldLengths,m_parameters->getLargeContigThreshold(),outputStream);
}

void Scaffolder::constructor(StaticVector*outbox,StaticVector*inbox,RingAllocator*outboxAllocator,Parameters*parameters,
	int*slaveMode,VirtualCommunicator*vc){
	m_slave_mode=slaveMode;
	m_virtualCommunicator=vc;
	m_outbox=outbox;
	m_inbox=inbox;
	m_outboxAllocator=outboxAllocator;
	m_parameters=parameters;
	m_initialised=false;
	m_workerId=0;
}

void Scaffolder::run(){
	if(!m_initialised){
		m_initialised=true;
		m_ready=true;
		m_contigId=0;
		m_positionOnContig=0;
		m_forwardDone=false;
		m_coverageRequested=false;
	}

	m_virtualCommunicator->forceFlush();
	m_virtualCommunicator->processInbox(&m_activeWorkers);
	m_activeWorkers.clear();

	if(m_contigId<(int)(*m_contigs).size()){
		processContig();
	}else{
		if(m_parameters->hasOption("-debug-scaffolder")){
			cout<<" sending MASTER_RANK,RAY_MPI_TAG_I_FINISHED_SCAFFOLDING"<<endl;
		}
		Message aMessage(NULL,0,MASTER_RANK,RAY_MPI_TAG_I_FINISHED_SCAFFOLDING,
			m_parameters->getRank());
		m_outbox->push_back(aMessage);
		(*m_slave_mode)=RAY_SLAVE_MODE_DO_NOTHING;
	}
}

void Scaffolder::setContigPaths(vector<uint64_t>*names,vector<vector<Kmer> >*paths){
	m_contigNames=names;
	m_contigs=paths;
}

void Scaffolder::processContig(){
	if(m_positionOnContig<(int)(*m_contigs)[m_contigId].size()){
		processContigPosition();
	}else if(!m_summaryPerformed){
		performSummary();
	}else if(!m_summarySent){
		sendSummary();
	}else if(!m_sentContigMeta){
		sendContigInfo();
	}else{
		m_contigId++;
		m_positionOnContig=0;
	}
}

void Scaffolder::sendContigInfo(){
	if(!m_sentContigInfo){
		uint64_t*message=(uint64_t*)m_outboxAllocator->allocate(MAXIMUM_MESSAGE_SIZE_IN_BYTES);
		message[0]=(*m_contigNames)[m_contigId];
		message[1]=(*m_contigs)[m_contigId].size();
		Message aMessage(message,2,
			MASTER_RANK,RAY_MPI_TAG_CONTIG_INFO,m_parameters->getRank());
		m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
		m_sentContigInfo=true;
	}else if(m_virtualCommunicator->isMessageProcessed(m_workerId)){
		vector<uint64_t> response;
		m_virtualCommunicator->getMessageResponseElements(m_workerId,&response);
		m_sentContigMeta=true;
	}
}

void Scaffolder::sendSummary(){
	if(m_summaryIterator<(int)m_summary.size()){
		if(!m_entrySent){
			uint64_t*message=(uint64_t*)m_outboxAllocator->allocate(MAXIMUM_MESSAGE_SIZE_IN_BYTES);

			int i=0;
			m_summary[m_summaryIterator].pack(message,&i);

			Message aMessage(message,i,
				MASTER_RANK,RAY_MPI_TAG_SCAFFOLDING_LINKS,m_parameters->getRank());
			m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
			m_entrySent=true;
		}else if(m_virtualCommunicator->isMessageProcessed(m_workerId)){
			vector<uint64_t> response;
			m_virtualCommunicator->getMessageResponseElements(m_workerId,&response);
			m_summaryIterator++;
			m_entrySent=false;
		}
	}else{
		m_summarySent=true;
		m_sentContigMeta=false;
		m_sentContigInfo=false;
	}
}

void Scaffolder::performSummary(){

	uint64_t sum=0;

	#ifdef ASSERT
	assert(m_contigId < (int)m_contigs->size());
	assert(m_vertexCoverageValues.size() == (*m_contigs)[m_contigId].size());
	#endif

	map<int,int> distribution;
	int n=0;
	for(int i=0;i<(int)m_vertexCoverageValues.size();i++){
		int coverageValue=m_vertexCoverageValues[i];
		distribution[coverageValue]++;

		if(coverageValue < m_parameters->getRepeatCoverage()){
			sum+=coverageValue;
			n++;
		}
	}

	int mean=sum;

	if(n>0){
		mean /= n;
	}

	uint64_t sumOfSquares=0;
	for(int i=0;i<(int)m_vertexCoverageValues.size();i++){
		int coverageValue=m_vertexCoverageValues[i];
		int diff=coverageValue-mean;
		if(coverageValue < m_parameters->getRepeatCoverage()){
			sumOfSquares+= diff*diff;
		}
	}

	if(n > 0){
		sumOfSquares /= n;
	}

	int standardDeviation=(int)sqrt(sumOfSquares);

	cout<<"contig: "<<(*m_contigNames)[m_contigId]<<" vertices: "<<m_vertexCoverageValues.size()<<" averageCoverage: "<<mean<<" standardDeviation: "<<standardDeviation<<" peakCoverage: "<<m_parameters->getPeakCoverage()<<" repeatCoverage: "<<m_parameters->getRepeatCoverage()<<endl;

	#ifdef SCAFFOLDER_SHOW_DISTRIBUTION
	cout<<"Distribution "<<endl;
	for(map<int,int>::iterator i=distribution.begin();i!=distribution.end();i++){
		cout<<" "<<i->first<<"	"<<i->second<<endl;
	}
	#endif

	/* write coverage values to a file if requested */
	if(m_parameters->hasOption("-write-contig-paths")){
		ostringstream fileName;
		fileName<<m_parameters->getPrefix()<<"Rank"<<m_parameters->getRank()<<".RayContigPaths.txt";
		ofstream fp;
		if(m_contigId==0){
			fp.open(fileName.str().c_str());
		}else{
			fp.open(fileName.str().c_str(),ios_base::out|ios_base::app);
		}

		#ifdef ASSERT
		assert(m_contigId < (int)m_contigNames->size());
		#endif

		uint64_t contigName=(*m_contigNames)[m_contigId];
		int vertices=m_vertexCoverageValues.size();
		fp<<"contig-"<<contigName<<endl;
		fp<<vertices<<" vertices"<<endl;
		fp<<"#Index	Vertex	Coverage"<<endl;
		for(int i=0;i<vertices;i++){
			Kmer kmer=(*m_contigs)[m_contigId][i];
			int coverage=m_vertexCoverageValues[i];

			fp<<i<<"	"<<kmer.idToWord(m_parameters->getWordSize(),m_parameters->getColorSpaceMode())<<"	"<<coverage<<endl;
		}

		fp.close();
	}

	m_summary.clear();
	m_summaryIterator=0;
	for(map<uint64_t,map<char,map<uint64_t,map<char,vector<ScaffoldingLink> > > > >::iterator i=
		m_scaffoldingSummary.begin();i!=m_scaffoldingSummary.end();i++){
		uint64_t leftContig=i->first;
		for(map<char,map<uint64_t,map<char,vector<ScaffoldingLink> > > >::iterator j=i->second.begin();
			j!=i->second.end();j++){
			char leftStrand=j->first;
			for(map<uint64_t,map<char,vector<ScaffoldingLink> > >::iterator k=j->second.begin();
				k!=j->second.end();k++){
				uint64_t rightContig=k->first;
				for(map<char,vector<ScaffoldingLink> >::iterator l=k->second.begin();
					l!=k->second.end();l++){
					char rightStrand=l->first;

					vector<int> veryRawDistances;

					for(vector<ScaffoldingLink>::iterator m=l->second.begin();m!=l->second.end();m++){
						int distance=(*m).getDistance();
						int coverage1=(*m).getCoverage1();
						int coverage2=(*m).getCoverage2();

						int numberOfStandardDeviations=1;

						/* only pick up things that are not repeated */
						if((mean-numberOfStandardDeviations*standardDeviation) <= coverage1 && coverage1 <= (mean+numberOfStandardDeviations*standardDeviation)
						  && (mean-numberOfStandardDeviations*standardDeviation) <= coverage2 && coverage2 <= (mean+numberOfStandardDeviations*standardDeviation)){
							veryRawDistances.push_back(distance);
						}
					}

					int count=veryRawDistances.size();

					/* no links are valid */
					if(count == 0){
						continue;
					}

					int averageValue=getAverage(&veryRawDistances);
					int standardDeviationValue=getStandardDeviation(&veryRawDistances);
					/* this summary information will be sent to MASTER later */
					SummarizedLink entry(leftContig,leftStrand,rightContig,rightStrand,averageValue,count,standardDeviationValue);
					m_summary.push_back(entry);
				}
			}
		}
	}
	m_summaryPerformed=true;
	m_summarySent=false;
	m_entrySent=false;
}

void Scaffolder::processContigPosition(){
	#ifdef ASSERT
	assert(m_contigId<(int)(*m_contigs).size());
	assert(m_positionOnContig<(int)(*m_contigs)[m_contigId].size());
	#endif

	Kmer vertex=(*m_contigs)[m_contigId][m_positionOnContig];
	#ifdef ASSERT
	assert(m_parameters!=NULL);
	#endif
	if(!m_forwardDone){
		processVertex(&vertex);
	}else if(!m_reverseDone){
		// get the coverage
		// if < maxCoverage
		// 	get read markers
		// 	for each read marker
		// 		if it is paired
		// 			get its pair
		// 				get the vertex for the opposite strand of the first read
		// 				get the coverage of this vertex
		// 				if < maxCoverage
		// 					get the paths that goes on them
		// 					print the linking information

		m_reverseDone=true;
	}else{
		m_positionOnContig++;
		m_forwardDone=false;
		m_coverageRequested=false;
	}
}

void Scaffolder::processVertex(Kmer*vertex){
	// get the coverage
	// if < maxCoverage
	// 	get read markers
	// 	for each read marker
	// 		if it is paired
	// 			get its pair
	// 				get the vertex for the opposite strand of the first read
	// 				get the coverage of this vertex
	// 				if < maxCoverage
	// 					get the paths that goes on them
	// 					print the linking information
	if(!m_coverageRequested){
		uint64_t*buffer=(uint64_t*)m_outboxAllocator->allocate(1*sizeof(Kmer));
		int bufferPosition=0;
		vertex->pack(buffer,&bufferPosition);
		Message aMessage(buffer,m_virtualCommunicator->getElementsPerQuery(RAY_MPI_TAG_GET_VERTEX_EDGES_COMPACT),
			m_parameters->_vertexRank(vertex),RAY_MPI_TAG_GET_VERTEX_EDGES_COMPACT,m_parameters->getRank());
		m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
		m_coverageRequested=true;
		m_coverageReceived=false;
		if(m_positionOnContig==0){
			m_scaffoldingSummary.clear();
			m_summaryPerformed=false;

			m_vertexCoverageValues.clear();
		}
		if(m_positionOnContig==(int)(*m_contigs)[m_contigId].size()-1){
			printf("Rank %i: gathering scaffold links [%i/%i] [%i/%i] (completed)\n",m_parameters->getRank(),
				m_contigId+1,(int)(*m_contigs).size(),
				m_positionOnContig+1,(int)(*m_contigs)[m_contigId].size());
		}else if(m_positionOnContig%10000==0){
			printf("Rank %i: gathering scaffold links [%i/%i] [%i/%i]\n",m_parameters->getRank(),
				m_contigId+1,(int)(*m_contigs).size(),
				m_positionOnContig+1,(int)(*m_contigs)[m_contigId].size());
		}
	}else if(!m_coverageReceived
		&&m_virtualCommunicator->isMessageProcessed(m_workerId)){
		vector<uint64_t> elements;
		m_virtualCommunicator->getMessageResponseElements(m_workerId,&elements);

		#ifdef ASSERT
		assert(elements.size() > 0);
		#endif

		uint8_t edges=elements[0];
		int coverage=elements[1];

		m_receivedCoverage=coverage;
		m_coverageReceived=true;
		m_initialisedFetcher=false;

		/* receive the coverage value at this position */
		m_vertexCoverageValues.push_back(m_receivedCoverage);

		/* here, make sure that the vertex has exactly 1 parent and 1 child */
		int parents=vertex->_getIngoingEdges(edges,m_parameters->getWordSize()).size();
		int children=vertex->_getOutgoingEdges(edges,m_parameters->getWordSize()).size();

		#ifdef SHOW_EDGES
		cout<<"/ "<<coverage<<" "<<parents<<" "<<children<<endl;
		#endif

		bool invalidVertex=(!(parents == 1 && children == 1));

		// don't judge a vertex by its parents and children
		// TODO: dump invalid parents and children by using their coverage (see SeedWorker.cpp)
		invalidVertex=false;

		if(invalidVertex){
			m_forwardDone=true;
			m_reverseDone=false;
		}

	}else if(m_coverageReceived){
		if(m_receivedCoverage<m_parameters->getRepeatCoverage()){
			if(!m_initialisedFetcher){
				m_readFetcher.constructor(vertex,m_outboxAllocator,m_inbox,
				m_outbox,m_parameters,m_virtualCommunicator,m_workerId);
				m_readAnnotationId=0;
				m_initialisedFetcher=true;
				m_hasPairRequested=false;
			}else if(!m_readFetcher.isDone()){
				m_readFetcher.work();
			}else{
				processAnnotations();
			}
		}else{
			m_forwardDone=true;
			m_reverseDone=false;
		}
	}
}

void Scaffolder::processAnnotations(){
	if(m_readAnnotationId<(int)m_readFetcher.getResult()->size()){
		processAnnotation();
	}else{
		m_forwardDone=true;
		m_reverseDone=false;
	}
}

void Scaffolder::processAnnotation(){
	// if is paired
	// 	get the forward and the reverse markers
	// 	get the coverage of the forward vertex
	// 	if < maxCoverage
	//	 	get the Direction
	//	 	if only 1 Direction
	//	 		if contig is not self
	//	 			get its length
	//	 			print link information
	//
	// 	get the coverage of the reverse vertex
	// 	if < maxCoverage
	//	 	get the Direction
	//	 	if only 1 Direction
	//	 		if contig is not self
	//	 			get its length
	//	 			print link information
	//

	ReadAnnotation*a=&(m_readFetcher.getResult()->at(m_readAnnotationId));
	int rank=a->getRank();
	int sequenceId=a->getReadIndex();
	char strand=a->getStrand();
	int positionOnStrand=a->getPositionOnStrand();
	if(!m_hasPairRequested){
		uint64_t*buffer=(uint64_t*)m_outboxAllocator->allocate(1*sizeof(Kmer));
		buffer[0]=sequenceId;
		Message aMessage(buffer,1,rank,RAY_MPI_TAG_HAS_PAIRED_READ,m_parameters->getRank());
		m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
		m_hasPairRequested=true;
		m_hasPairReceived=false;
	}else if(!m_hasPairReceived
	&&m_virtualCommunicator->isMessageProcessed(m_workerId)){
		vector<uint64_t> response;
		m_virtualCommunicator->getMessageResponseElements(m_workerId,&response);
		m_hasPair=response[0];
		m_hasPairReceived=true;
		m_pairRequested=false;
	}else if(!m_hasPairReceived){
		return;
	}else if(!m_hasPair){
		m_readAnnotationId++;
		m_hasPairRequested=false;
	}else if(!m_pairRequested){
		uint64_t*buffer=(uint64_t*)m_outboxAllocator->allocate(1*sizeof(Kmer));
		buffer[0]=sequenceId;
		Message aMessage(buffer,1,
		rank,RAY_MPI_TAG_GET_READ_MATE,m_parameters->getRank());
		m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
		m_pairRequested=true;
		m_pairReceived=false;
	}else if(!m_pairReceived
	&&m_virtualCommunicator->isMessageProcessed(m_workerId)){
		vector<uint64_t> response;
		m_virtualCommunicator->getMessageResponseElements(m_workerId,&response);
		m_readLength=response[0];
		m_pairedReadRank=response[1];
		m_pairedReadIndex=response[2];
		m_pairedReadLibrary=response[3];
		m_pairReceived=true;
		m_markersRequested=false;
	}else if(!m_pairReceived){
		return;
	}else if(!m_markersRequested){
		uint64_t*buffer=(uint64_t*)m_outboxAllocator->allocate(1*sizeof(Kmer));
		buffer[0]=m_pairedReadIndex;
		Message aMessage(buffer,1,
		m_pairedReadRank,RAY_MPI_TAG_GET_READ_MARKERS,m_parameters->getRank());
		m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
		m_markersRequested=true;
		m_markersReceived=false;
	}else if(!m_markersReceived
	&&m_virtualCommunicator->isMessageProcessed(m_workerId)){
		vector<uint64_t> response;
		m_virtualCommunicator->getMessageResponseElements(m_workerId,&response);
		int bufferPosition=0;
		m_pairedReadLength=response[bufferPosition++];
		m_pairedForwardMarker.unpack(&response,&bufferPosition);
		m_pairedReverseMarker.unpack(&response,&bufferPosition);
		m_pairedForwardOffset=response[bufferPosition++];
		m_pairedReverseOffset=response[bufferPosition++];
		m_markersReceived=true;
		m_forwardDirectionsRequested=false;
	}else if(!m_markersReceived){
		return;
/***
 *
 * Forward Directions
 *
 ***/
	}else if(!m_forwardDirectionsRequested){
		// skip unrelated marker
		if(m_pairedForwardOffset>m_pairedReadLength-m_parameters->getWordSize()){
			m_forwardDirectionsRequested=true;
			m_forwardDirectionsReceived=true;
			m_forwardDirectionLengthRequested=true;
			m_forwardDirectionLengthReceived=true;
		}
		uint64_t*buffer=(uint64_t*)m_outboxAllocator->allocate(1*sizeof(Kmer));
		int bufferPosition=0;
		m_pairedForwardMarker.pack(buffer,&bufferPosition);

		int elementsPerQuery=m_virtualCommunicator->getElementsPerQuery(RAY_MPI_TAG_GET_COVERAGE_AND_DIRECTION);

		Message aMessage(buffer,elementsPerQuery,
			m_parameters->_vertexRank(&m_pairedForwardMarker),
			RAY_MPI_TAG_GET_COVERAGE_AND_DIRECTION,m_parameters->getRank());
		m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
		m_forwardDirectionsRequested=true;
		m_forwardDirectionsReceived=false;
	}else if(!m_forwardDirectionsReceived
	&&m_virtualCommunicator->isMessageProcessed(m_workerId)){
		vector<uint64_t> response;
		m_virtualCommunicator->getMessageResponseElements(m_workerId,&response);
		m_pairedForwardMarkerCoverage=response[0];

		m_pairedForwardHasDirection=response[1];
		m_pairedForwardDirectionName=response[2];
		m_pairedForwardDirectionPosition=response[3];

		uint8_t edges=response[4];

		m_forwardDirectionsReceived=true;
		m_reverseDirectionsRequested=false;
		m_forwardDirectionLengthRequested=false;

		/* here, make sure that the vertex has exactly 1 parent and 1 child */
		int parents=m_pairedForwardMarker._getIngoingEdges(edges,m_parameters->getWordSize()).size();
		int children=m_pairedForwardMarker._getOutgoingEdges(edges,m_parameters->getWordSize()).size();

		#ifdef SHOW_EDGES
		cout<<"/ "<<m_pairedForwardMarkerCoverage<<" "<<parents<<" "<<children<<endl;
		#endif

		bool invalidVertex=(!(parents == 1 && children == 1));

		// don't judge a vertex by its parents and children
		// TODO: dump invalid parents and children by using their coverage (see SeedWorker.cpp)
		invalidVertex=false;

		/* the hit is invalid */
		if((*m_contigNames)[m_contigId]==m_pairedForwardDirectionName
		||!(m_pairedForwardMarkerCoverage<m_parameters->getRepeatCoverage())
		|| !m_pairedForwardHasDirection || invalidVertex){
			m_forwardDirectionLengthRequested=true;
			m_forwardDirectionLengthReceived=true;
		}
	}else if(!m_forwardDirectionsReceived){
		return;
	}else if(!m_forwardDirectionLengthRequested){
		uint64_t*buffer=(uint64_t*)m_outboxAllocator->allocate(1*sizeof(Kmer));
		int rankId=getRankFromPathUniqueId(m_pairedForwardDirectionName);
		buffer[0]=m_pairedForwardDirectionName;
		Message aMessage(buffer,1,
		rankId,
		RAY_MPI_TAG_GET_PATH_LENGTH,m_parameters->getRank());
		m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
		m_forwardDirectionLengthRequested=true;
		m_forwardDirectionLengthReceived=false;

	}else if(!m_forwardDirectionLengthReceived
	&&m_virtualCommunicator->isMessageProcessed(m_workerId)){
		vector<uint64_t> response;
		m_virtualCommunicator->getMessageResponseElements(m_workerId,&response);
		m_pairedForwardDirectionLength=response[0];
		m_forwardDirectionLengthReceived=true;

		int range=m_parameters->getLibraryMaxAverageLength(m_pairedReadLibrary)+3*m_parameters->getLibraryMaxStandardDeviation(m_pairedReadLibrary);

		if(m_pairedForwardDirectionLength<range
		||(int)(*m_contigs)[m_contigId].size()<range
		|| 2*m_receivedCoverage<m_pairedForwardMarkerCoverage
		|| 2*m_pairedForwardMarkerCoverage<m_receivedCoverage ){
			return;
		}

		if(m_parameters->hasOption("-debug-scaffolder")){
			cout<<endl;
			cout<<"AverageDistance: "<<m_parameters->getLibraryMaxAverageLength(m_pairedReadLibrary)<<endl;
			cout<<"StandardDeviation: "<<m_parameters->getLibraryMaxStandardDeviation(m_pairedReadLibrary)<<endl;
			cout<<"Path1: "<<(*m_contigNames)[m_contigId]<<endl;
			cout<<" Length: "<<(*m_contigs)[m_contigId].size()<<endl;
			cout<<" Position: "<<m_positionOnContig<<endl;
			cout<<" Coverage: "<<m_receivedCoverage<<endl;
			cout<<" PathStrand: F"<<endl;
			cout<<" ReadStrand: "<<strand<<endl;
			cout<<" ReadLength: "<<m_readLength<<endl;
			cout<<" PositionInRead: "<<positionOnStrand<<endl;
			cout<<"Path2: "<<m_pairedForwardDirectionName<<endl;
			cout<<" Length: "<<m_pairedForwardDirectionLength<<endl;
			cout<<" Position: "<<m_pairedForwardDirectionPosition<<endl;
			cout<<" Coverage: "<<m_pairedForwardMarkerCoverage<<endl;
			cout<<" PathStrand: F"<<endl;
			cout<<" ReadStrand: F"<<endl;
			cout<<" ReadLength: "<<m_pairedReadLength<<endl;
			cout<<" PositionInRead: "<<m_pairedForwardOffset<<endl;
		}

		bool path1IsLeft=false;
		bool path1IsRight=false;
		bool path2IsLeft=false;
		bool path2IsRight=false;
		if(m_positionOnContig<range)
			path1IsLeft=true;
		if(m_positionOnContig>(int)(*m_contigs)[m_contigId].size()-range)
			path1IsRight=true;
		if(m_pairedForwardDirectionPosition<range)
			path2IsLeft=true;
		if(m_pairedForwardDirectionPosition>m_pairedForwardDirectionLength-range)
			path2IsRight=true;

		/* don't try to solve the problem if both paths allow both sides */
		if((path1IsLeft&&path1IsRight)&&(path2IsLeft&&path2IsRight))
			return;
/*
Case 6. (allowed)

                    ---->
                                                           ---->
------------------------>              ------------------------>
*/

		if(path1IsRight&&path2IsRight&&strand=='F'){
			int distanceIn1=(*m_contigs)[m_contigId].size()-m_positionOnContig+positionOnStrand;
			int distanceIn2=m_pairedForwardDirectionLength-m_pairedForwardDirectionPosition+m_pairedForwardOffset;
			int distance=range-distanceIn1-distanceIn2;
			if(distance>0){
				ScaffoldingLink hit;
				hit.constructor(distance,m_receivedCoverage,m_pairedForwardMarkerCoverage);
				m_scaffoldingSummary[(*m_contigNames)[m_contigId]]['F'][m_pairedForwardDirectionName]['R'].push_back(hit);
				if(m_parameters->hasOption("-debug-scaffolder")){
					cout<<"LINK06 "<<(*m_contigNames)[m_contigId]<<",F,"<<m_pairedForwardDirectionName<<",R,"<<distance<<endl;
				}
			}
/*
Case 1. (allowed)

---->
                                       ---->
------------------------>              ------------------------>
*/
		}else if(path1IsLeft&&path2IsLeft&&strand=='F'){
			int distanceIn1=m_positionOnContig+m_readLength-positionOnStrand;
			int distanceIn2=m_pairedForwardDirectionPosition+m_pairedReadLength-m_pairedForwardOffset;
			int distance=range-distanceIn1-distanceIn2;
			if(distance>0){
				if(m_parameters->hasOption("-debug-scaffolder")){
					cout<<"LINK01 "<<(*m_contigNames)[m_contigId]<<",R,"<<m_pairedForwardDirectionName<<",F,"<<distance<<endl;
				}
				ScaffoldingLink hit;
				hit.constructor(distance,m_receivedCoverage,m_pairedForwardMarkerCoverage);
				m_scaffoldingSummary[(*m_contigNames)[m_contigId]]['R'][m_pairedForwardDirectionName]['F'].push_back(hit);
			}
/*
Case 10. (allowed)

<----
                                                           ---->
------------------------>              ------------------------>

                   ---->              <----
<-----------------------              <-------------------------
*/
		}else if(path1IsLeft&&path2IsRight&&strand=='R'){
			int distanceIn1=m_positionOnContig+positionOnStrand;
			int distanceIn2=m_pairedForwardDirectionLength-m_pairedForwardDirectionPosition+m_pairedForwardOffset;
			int distance=range-distanceIn1-distanceIn2;
			if(distance>0){
				if(m_parameters->hasOption("-debug-scaffolder")){
					cout<<"LINK10 "<<(*m_contigNames)[m_contigId]<<",R,"<<m_pairedForwardDirectionName<<",R,"<<distance<<endl;
				}
				ScaffoldingLink hit;
				hit.constructor(distance,m_receivedCoverage,m_pairedForwardMarkerCoverage);
				m_scaffoldingSummary[(*m_contigNames)[m_contigId]]['R'][m_pairedForwardDirectionName]['R'].push_back(hit);
			}

/*
Case 13. (allowed)

                    <----
                                       ---->
------------------------>              ------------------------>
*/
		}else if(path1IsRight&&path2IsLeft&&strand=='R'){
			int distanceIn1=(*m_contigs)[m_contigId].size()-m_positionOnContig-positionOnStrand+m_readLength;
			int distanceIn2=m_pairedForwardDirectionPosition+m_pairedReadLength-m_pairedForwardOffset;
			int distance=range-distanceIn1-distanceIn2;
			if(distance>0){
				if(m_parameters->hasOption("-debug-scaffolder")){
					cout<<"LINK13 "<<(*m_contigNames)[m_contigId]<<",F,"<<m_pairedForwardDirectionName<<",F,"<<distance<<endl;
				}

				ScaffoldingLink hit;
				hit.constructor(distance,m_receivedCoverage,m_pairedForwardMarkerCoverage);
				m_scaffoldingSummary[(*m_contigNames)[m_contigId]]['F'][m_pairedForwardDirectionName]['F'].push_back(hit);
			}
		}

	}else if(!m_forwardDirectionLengthReceived){
		return;

/***
 *
 * Reverse Directions
 *
 ***/
	}else if(!m_reverseDirectionsRequested){
		if(m_pairedReverseOffset>m_pairedReadLength-m_parameters->getWordSize()){
			m_reverseDirectionsRequested=true;
			m_reverseDirectionsReceived=true;
			m_reverseDirectionLengthRequested=true;
			m_reverseDirectionLengthReceived=true;
		}

		uint64_t*buffer=(uint64_t*)m_outboxAllocator->allocate(1*sizeof(Kmer));
		int bufferPosition=0;
		m_pairedReverseMarker.pack(buffer,&bufferPosition);

		int elementsPerQuery=m_virtualCommunicator->getElementsPerQuery(RAY_MPI_TAG_GET_COVERAGE_AND_DIRECTION);

		Message aMessage(buffer,elementsPerQuery,
			m_parameters->_vertexRank(&m_pairedReverseMarker),
			RAY_MPI_TAG_GET_COVERAGE_AND_DIRECTION,m_parameters->getRank());
		m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
		m_reverseDirectionsRequested=true;
		m_reverseDirectionsReceived=false;
	}else if(!m_reverseDirectionsReceived
	&&m_virtualCommunicator->isMessageProcessed(m_workerId)){
		vector<uint64_t> response;
		m_virtualCommunicator->getMessageResponseElements(m_workerId,&response);
		m_pairedReverseMarkerCoverage=response[0];
		m_pairedReverseHasDirection=response[1];
		m_pairedReverseDirectionName=response[2];
		m_pairedReverseDirectionPosition=response[3];

		uint8_t edges=response[4];

		m_reverseDirectionsReceived=true;
		m_reverseDirectionLengthRequested=false;

		/* here, make sure that the vertex has exactly 1 parent and 1 child */
		int parents=m_pairedReverseMarker._getIngoingEdges(edges,m_parameters->getWordSize()).size();
		int children=m_pairedReverseMarker._getOutgoingEdges(edges,m_parameters->getWordSize()).size();

		#ifdef SHOW_EDGES
		cout<<"/ "<<m_pairedReverseMarkerCoverage<<" "<<parents<<" "<<children<<endl;
		#endif

		bool invalidVertex=(!(parents == 1 && children == 1));

		// don't judge a vertex by its parents and children
		// TODO: dump invalid parents and children by using their coverage (see SeedWorker.cpp)
		invalidVertex=false;

		/* the hit is invalid */
		if((*m_contigNames)[m_contigId]==m_pairedReverseDirectionName
		||!(m_pairedReverseMarkerCoverage<m_parameters->getRepeatCoverage())
		|| !m_pairedReverseHasDirection || invalidVertex){
			m_reverseDirectionLengthRequested=true;
			m_reverseDirectionLengthReceived=true;
		}
	}else if(!m_reverseDirectionsReceived){
		return;
	}else if(!m_reverseDirectionLengthRequested){
		uint64_t*buffer=(uint64_t*)m_outboxAllocator->allocate(1*sizeof(Kmer));
		int rankId=getRankFromPathUniqueId(m_pairedReverseDirectionName);
		buffer[0]=m_pairedReverseDirectionName;
		Message aMessage(buffer,1,
		rankId,RAY_MPI_TAG_GET_PATH_LENGTH,m_parameters->getRank());
		m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
		m_reverseDirectionLengthRequested=true;
		m_reverseDirectionLengthReceived=false;
	}else if(!m_reverseDirectionLengthReceived
	&&m_virtualCommunicator->isMessageProcessed(m_workerId)){
		vector<uint64_t> response;
		m_virtualCommunicator->getMessageResponseElements(m_workerId,&response);
		m_pairedReverseDirectionLength=response[0];
		m_reverseDirectionLengthReceived=true;

		int range=m_parameters->getLibraryMaxAverageLength(m_pairedReadLibrary)+3*m_parameters->getLibraryMaxStandardDeviation(m_pairedReadLibrary);

		if(m_pairedReverseDirectionLength<range
		||(int)(*m_contigs)[m_contigId].size()<range
		|| 2*m_receivedCoverage<m_pairedReverseMarkerCoverage
		|| 2*m_pairedReverseMarkerCoverage<m_receivedCoverage){
			return;
		}

		if(m_parameters->hasOption("-debug-scaffolder")){
			cout<<endl;
			cout<<"AverageDistance: "<<m_parameters->getLibraryMaxAverageLength(m_pairedReadLibrary)<<endl;
			cout<<"StandardDeviation: "<<m_parameters->getLibraryMaxStandardDeviation(m_pairedReadLibrary)<<endl;
			cout<<"Path1: "<<(*m_contigNames)[m_contigId]<<endl;
			cout<<" Length: "<<(*m_contigs)[m_contigId].size()<<endl;
			cout<<" Position: "<<m_positionOnContig<<endl;
			cout<<" Coverage: "<<m_receivedCoverage<<endl;
			cout<<" PathStrand: F"<<endl;
			cout<<" ReadStrand: "<<strand<<endl;
			cout<<" ReadLength: "<<m_readLength<<endl;
			cout<<" PositionInRead: "<<positionOnStrand<<endl;
			cout<<"Path2: "<<m_pairedReverseDirectionName<<endl;
			cout<<" Length: "<<m_pairedReverseDirectionLength<<endl;
			cout<<" Position: "<<m_pairedReverseDirectionPosition<<endl;
			cout<<" Coverage: "<<m_pairedReverseMarkerCoverage<<endl;
			cout<<" PathStrand: F"<<endl;
			cout<<" ReadStrand: R"<<endl;
			cout<<" ReadLength: "<<m_pairedReadLength<<endl;
			cout<<" PositionInRead: "<<m_pairedReverseOffset<<endl;
		}

		bool path1IsLeft=false;
		bool path1IsRight=false;
		bool path2IsLeft=false;
		bool path2IsRight=false;
		if(m_positionOnContig<range)
			path1IsLeft=true;
		if(m_positionOnContig>(int)(*m_contigs)[m_contigId].size()-range)
			path1IsRight=true;
		if(m_pairedReverseDirectionPosition<range)
			path2IsLeft=true;
		if(m_pairedReverseDirectionPosition>m_pairedReverseDirectionLength-range)
			path2IsRight=true;

		/* don't try to solve the problem if both paths allow both sides */
		if((path1IsLeft&&path1IsRight)&&(path2IsLeft&&path2IsRight))
			return;

/*
Case 4. (allowed)

---->
                                                           <----
------------------------>              ------------------------>
*/

		if(path1IsLeft&&path2IsRight&&strand=='F'){
			int distanceIn1=m_positionOnContig+m_readLength-positionOnStrand;
			int distanceIn2=m_pairedReverseDirectionLength-m_pairedReverseDirectionPosition-m_pairedReverseOffset+m_pairedReadLength;
			int distance=range-distanceIn1-distanceIn2;
			if(distance>0){
				if(m_parameters->hasOption("-debug-scaffolder")){
					cout<<"LINK04 "<<(*m_contigNames)[m_contigId]<<",R,"<<m_pairedReverseDirectionName<<",R,"<<distance<<endl;
				}
				ScaffoldingLink hit;
				hit.constructor(distance,m_receivedCoverage,m_pairedReverseMarkerCoverage);
				m_scaffoldingSummary[(*m_contigNames)[m_contigId]]['R'][m_pairedReverseDirectionName]['R'].push_back(hit);
			}


/*
Case 7. (allowed)

                    ---->
                                       <----
------------------------>              ------------------------>
*/
		}else if(path1IsRight&&path2IsLeft&&strand=='F'){
			int distanceIn1=(*m_contigs)[m_contigId].size()-m_positionOnContig+positionOnStrand;
			int distanceIn2=m_pairedReverseDirectionPosition+m_pairedReverseOffset;
			int distance=range-distanceIn1-distanceIn2;
			if(distance>0){
				if(m_parameters->hasOption("-debug-scaffolder")){
					cout<<"LINK07 "<<(*m_contigNames)[m_contigId]<<",F,"<<m_pairedReverseDirectionName<<",F,"<<distance<<endl;
				}
				ScaffoldingLink hit;
				hit.constructor(distance,m_receivedCoverage,m_pairedReverseMarkerCoverage);
				m_scaffoldingSummary[(*m_contigNames)[m_contigId]]['F'][m_pairedReverseDirectionName]['F'].push_back(hit);
			}


/*
Case 11. (allowed)

<----
                                       <----
------------------------>              ------------------------>
*/
		}else if(path1IsLeft&&path2IsLeft&&strand=='R'){
			int distanceIn1=m_positionOnContig+positionOnStrand;
			int distanceIn2=m_pairedReverseDirectionPosition+m_pairedReverseOffset;
			int distance=range-distanceIn1-distanceIn2;
			if(distance>0){
				if(m_parameters->hasOption("-debug-scaffolder")){
					cout<<"LINK11 "<<(*m_contigNames)[m_contigId]<<",R,"<<m_pairedReverseDirectionName<<",F,"<<distance<<endl;
				}
				ScaffoldingLink hit;
				hit.constructor(distance,m_receivedCoverage,m_pairedReverseMarkerCoverage);
				m_scaffoldingSummary[(*m_contigNames)[m_contigId]]['R'][m_pairedReverseDirectionName]['F'].push_back(hit);
			}

/*
Case 16. (allowed)

                    <----
                                                           <----
------------------------>              ------------------------>
*/
		}else if(path1IsRight&&path2IsRight&&strand=='R'){
			int distanceIn1=(*m_contigs)[m_contigId].size()-m_positionOnContig-positionOnStrand+m_readLength;
			int distanceIn2=m_pairedReverseDirectionLength-m_pairedReverseDirectionPosition-m_pairedReverseOffset+m_pairedReadLength;
			int distance=range-distanceIn1-distanceIn2;
			if(distance>0){
				if(m_parameters->hasOption("-debug-scaffolder")){
					cout<<"LINK16 "<<(*m_contigNames)[m_contigId]<<",F,"<<m_pairedReverseDirectionName<<",R,"<<distance<<endl;
				}
				ScaffoldingLink hit;
				hit.constructor(distance,m_receivedCoverage,m_pairedReverseMarkerCoverage);
				m_scaffoldingSummary[(*m_contigNames)[m_contigId]]['F'][m_pairedReverseDirectionName]['R'].push_back(hit);
			}
		}
	}else if(!m_reverseDirectionLengthReceived){
		return;

	}else if(m_reverseDirectionLengthReceived){
		m_readAnnotationId++;
		m_hasPairRequested=false;
	}
}

void Scaffolder::getContigSequence(uint64_t id){
	if(!m_hasContigSequence_Initialised){
		m_hasContigSequence_Initialised=true;
		m_rankIdForContig=getRankFromPathUniqueId(id);
		m_theLength=m_contigLengths[id];
		m_position=0;
		m_contigPath.clear();
		m_requestedContigChunk=false;
	}

	if(m_position<m_theLength){
		if(!m_requestedContigChunk){
			m_requestedContigChunk=true;
			uint64_t*message=(uint64_t*)m_outboxAllocator->allocate(MAXIMUM_MESSAGE_SIZE_IN_BYTES);
			message[0]=id;
			message[1]=m_position;
			Message aMessage(message,2,
				m_rankIdForContig,RAY_MPI_TAG_GET_CONTIG_CHUNK,m_parameters->getRank());
			m_virtualCommunicator->pushMessage(m_workerId,&aMessage);
		}else if(m_virtualCommunicator->isMessageProcessed(m_workerId)){
			vector<uint64_t> data;
			m_virtualCommunicator->getMessageResponseElements(m_workerId,&data);
			/* the position in the message buffer */
			int pos=0;
			/* the first element is the number of Kmer */
			int count=data[pos++];
			int kmerIterator=0;
			while(kmerIterator<count){
				Kmer a;
				a.unpack(&data,&pos);
				m_contigPath.push_back(a);
				kmerIterator++;
			}
			m_position+=count;
			m_requestedContigChunk=false;
		}
	}else{
		/* we should receive a correct number of vertices */
		#ifdef ASSERT
		assert((int)m_contigPath.size()==m_theLength);
		#endif

		m_contigSequence=convertToString(&m_contigPath,m_parameters->getWordSize(),m_parameters->getColorSpaceMode());
		m_hasContigSequence=true;
	}
}

void Scaffolder::writeScaffolds(){
	if(!m_initialised){
		m_initialised=true;
		m_scaffoldId=0;
		m_contigId=0;
		/* actually it is a position on the scaffold */
		m_positionOnScaffold=0;
		m_hasContigSequence=false;
		m_hasContigSequence_Initialised=false;
		string file=m_parameters->getScaffoldFile();
		m_fp=fopen(file.c_str(),"w");
	}

	m_virtualCommunicator->forceFlush();
	m_virtualCommunicator->processInbox(&m_activeWorkers);
	m_activeWorkers.clear();

	if(m_scaffoldId<(int)m_scaffoldContigs.size()){
		if(m_contigId<(int)m_scaffoldContigs[m_scaffoldId].size()){
			uint64_t contigNumber=m_scaffoldContigs[m_scaffoldId][m_contigId];
			if(!m_hasContigSequence){
				getContigSequence(contigNumber);
			}else{ /* at this point, m_contigSequence is filled. */
				if(m_contigId==0){
					fprintf(m_fp,">scaffold-%i\n",m_scaffoldId);
					m_positionOnScaffold=0;
				}

				int contigPosition=0;
				char strand=m_scaffoldStrands[m_scaffoldId][m_contigId];
				if(strand=='R'){
					m_contigSequence=reverseComplement(&m_contigSequence);
				}

				int length=m_contigSequence.length();

				#ifdef ASSERT
				int theLength=m_contigLengths[contigNumber]+m_parameters->getWordSize()-1;
				assert(length==theLength);
				#endif

				int columns=m_parameters->getColumns();
				ostringstream outputBuffer;
				while(contigPosition<length){
					char nucleotide=m_contigSequence[contigPosition];
					outputBuffer<<nucleotide;
					contigPosition++;
					m_positionOnScaffold++;
					if(m_positionOnScaffold%columns==0){
						outputBuffer<<"\n";
					}
				}

				fprintf(m_fp,"%s",outputBuffer.str().c_str());

				if(m_contigId<(int)m_scaffoldContigs[m_scaffoldId].size()-1){
					int gapSize=m_scaffoldGaps[m_scaffoldId][m_contigId];
					int i=0;
					int columns=m_parameters->getColumns();
					ostringstream outputBuffer2;
					while(i<gapSize){
						outputBuffer2<<"N";
						i++;
						m_positionOnScaffold++;
						if(m_positionOnScaffold%columns==0){
							outputBuffer2<<"\n";
						}
					}
					fprintf(m_fp,"%s",outputBuffer2.str().c_str());
				}
				m_contigId++;
				m_hasContigSequence=false;
				m_hasContigSequence_Initialised=false;
			}
		}else{
			fprintf(m_fp,"\n");
			m_scaffoldId++;
			m_contigId=0;
			m_positionOnScaffold=0;
			m_hasContigSequence=false;
			m_hasContigSequence_Initialised=false;
		}
	}else{
		fclose(m_fp);
		m_parameters->setMasterMode(RAY_MASTER_MODE_KILL_RANKS);
	}
}


