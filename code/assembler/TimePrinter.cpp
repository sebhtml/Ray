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

#include<stdio.h>
#include<assembler/TimePrinter.h>
#include<iostream>
#include <sstream>
using namespace std;

void TimePrinter::printElapsedTime(string description){

	m_endingTime=time(NULL);
	int differenceWithLast=m_endingTime-m_lastTime;
	m_lastTime=m_endingTime;
	struct tm * timeinfo;
	timeinfo=localtime(&m_endingTime);

	m_descriptions.push_back(description);
	m_durations.push_back(differenceWithLast);

	printElapsedTimeInStream(&cout,description,timeinfo,differenceWithLast);

	if(m_fileSet){
		printElapsedTimeInStream(&m_file,description,timeinfo,differenceWithLast);
	}
}

void TimePrinter::printElapsedTimeInStream(ostream*stream,string description,struct tm*timeinfo,
int differenceWithLast){
	(*stream)<<endl;
	(*stream)<<"***"<<endl;
	(*stream)<<"Step: "<<description<<endl;
	(*stream)<<"Date: "<<asctime(timeinfo);
	(*stream)<<"Elapsed time: ";
	printDifference(differenceWithLast,stream);
	int totalSeconds=m_endingTime-m_startingTime;
	(*stream)<<endl;
	(*stream)<<"Since beginning: ";
	printDifference(totalSeconds,stream);
	(*stream)<<endl;
	(*stream)<<"***"<<endl;
	(*stream)<<endl;
	fflush(stdout);
}

void TimePrinter::printDifferenceFromStart(int rank){
	time_t endingTime=time(NULL);
	if(endingTime==m_last){
		return;
	}
	if((endingTime-m_startingTime)%3600!=0){
		return;
	}
	m_last=endingTime;
	int differenceWithLast=endingTime-m_startingTime;
	cout<<"Rank "<<rank<<": I am still running... ";
	printDifference(differenceWithLast,&cout);
	cout<<endl;
}

void TimePrinter::constructor(){
	m_startingTime=m_lastTime=m_endingTime=time(NULL);
	m_descriptions.clear();
	m_durations.clear();

	m_fileSet=false;
}

void TimePrinter::setFile(string prefix){
	ostringstream fileName;
	fileName<<prefix<<"ElapsedTime.txt";

	m_file.open(fileName.str().c_str(),ios_base::out);

	m_fileSet=true;
}

void TimePrinter::printDifference(int difference,ostream*stream){
	int minutes=difference/60;
	int seconds=difference%60;
	int hours=minutes/60;
	minutes=minutes%60;
	int days=hours/24;
	hours=hours%24;

	bool printed=false;

	if(days>0){
		(*stream)<<days<<" days";
		printed=true;
	}
	if(hours>0){
		if(printed){
			(*stream)<<", ";
		}
		printed=true;
		(*stream)<<hours<<" hours";
	}
	if(minutes>0){
		if(printed){
			(*stream)<<", ";
		}
		printed=true;
		(*stream)<<minutes<<" minutes";
	}

	if(printed){
		(*stream)<<", ";
	}
	(*stream)<<seconds<<" seconds";
}

void TimePrinter::printDurations(){

	m_endingTime=time(NULL);
	struct tm * timeinfo;
	timeinfo=localtime(&m_endingTime);
	m_descriptions.push_back("Total");
	m_durations.push_back(m_endingTime-m_startingTime);

	printDurationsInStream(&cout,timeinfo);

	if(m_fileSet){
		printDurationsInStream(&m_file,timeinfo);
		m_file.close();
		cout<<"Closing the file"<<endl;
	}
}

void TimePrinter::printDurationsInStream(ostream*stream,struct tm*timeinfo){
	(*stream)<<"\nElapsed time for each step, "<<asctime(timeinfo)<<endl;
	for(int i=0;i<(int)m_descriptions.size();i++){
		string text=m_descriptions[i];
		int seconds=m_durations[i];
		(*stream)<<" "<<text<<": ";
		printDifference(seconds,stream);
		(*stream)<<endl;
	}
}
