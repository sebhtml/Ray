/*
    Copyright 2014 Maxime Déraspe
    Copyright 2013 Sébastien Boisvert
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

#include "MatrixOwner.h"
#include "CoalescenceManager.h" // for DIE

#include <RayPlatform/core/OperatingSystem.h>

#include <fstream>
#include <string>
#include <sstream>
using namespace std;

#include <math.h>
#include <numeric>

#define INPUT_TYPE_GRAPH 0
#define INPUT_FILTERIN_GRAPH 1
#define INPUT_FILTEROUT_GRAPH 2
#define INPUT_TYPE_ASSEMBLY 3
#define INPUT_FILTERIN_ASSEMBLY 4
#define INPUT_FILTEROUT_ASSEMBLY 5

MatrixOwner::MatrixOwner() {

	m_completedStoreActors = 0;

	m_receivedPayloads = 0;

}

MatrixOwner::~MatrixOwner() {

}


void MatrixOwner::receive(Message & message) {

	int tag = message.getTag();
	char * buffer = message.getBufferBytes();
	int source = message.getSourceActor();

	if( tag == CoalescenceManager::DIE) {

		die();

	} else if(tag == GREETINGS) {

		int offset = 0;
		memcpy(&m_parameters, buffer + offset, sizeof(m_parameters));
		offset += sizeof(m_parameters);
		memcpy(&m_sampleNames, buffer + offset, sizeof(m_sampleNames));
		offset += sizeof(m_sampleNames);
		memcpy(&m_filterMatrices, buffer + offset, sizeof(m_filterMatrices));
		offset += sizeof(m_filterMatrices);

		// Build a m_sampleByFilter to only print samples and filters that belong to a filtered matrix
		for (map< int, vector<int> >::iterator it = m_filterMatrices->begin(); it!=m_filterMatrices->end(); ++it) {
			// fill sampleByFilter with 1 (sample present in that matrix)
			for (int z=0;z<it->second.size();z++) {
				m_sampleByFilter[it->first].push_back(1);
			}
		}

		for (map< int, vector<int> >::iterator it = m_filterMatrices->begin(); it!=m_filterMatrices->end(); ++it) {

			// skip -1 complete no filter matrix
			if (it->first == -1)
				continue;

			// look for filter samples to hide in other matrices
			vector<int> samples_to_hide_in_other_matrix;

			for(vector<int>::iterator sampleIt = it->second.begin(); sampleIt != it->second.end(); ++sampleIt) {
				if (*sampleIt != INPUT_TYPE_GRAPH && *sampleIt != INPUT_TYPE_ASSEMBLY) {
					samples_to_hide_in_other_matrix.push_back(distance(it->second.begin(), sampleIt));
				}
			}

			for (map< int, vector<int> >::iterator it2 = m_sampleByFilter.begin(); it2!=m_sampleByFilter.end(); ++it2) {
				if(it2->first == -1 || it2->first == it->first)
					continue;

				for(vector<int>::iterator sampleIt2 = samples_to_hide_in_other_matrix.begin(); sampleIt2 != samples_to_hide_in_other_matrix.end(); ++sampleIt2) {
					it2->second.at(*sampleIt2) = 0;
				}
			}

		}


#ifdef CONFIG_ASSERT
		assert(m_parameters != NULL);
		assert(m_sampleNames != NULL);
#endif
		m_mother = source;

	} else if(tag == PUSH_PAYLOAD) {

		SampleIdentifier sample1 = -1;
		SampleIdentifier sample2 = -1;
		LargeCount count = 0;
		int matrixNb;

		int offset = 0;

		memcpy(&sample1, buffer + offset, sizeof(sample1));
		offset += sizeof(sample1);
		memcpy(&sample2, buffer + offset, sizeof(sample2));
		offset += sizeof(sample2);
		memcpy(&count, buffer + offset, sizeof(count));
		offset += sizeof(count);
		memcpy(&matrixNb, buffer + offset, sizeof(matrixNb));
		offset += sizeof(matrixNb);


#ifdef CONFIG_ASSERT
		assert(sample1 >= 0);
		assert(sample2 >= 0);
		assert(count >= 0);
#endif

		m_receivedPayloads ++;

		m_gramMatrices[matrixNb][sample1][sample2] += count;

		if (sample1 != sample2) {
			m_gramMatrices[matrixNb][sample2][sample1] += count;
		}

		Message response;
		response.setTag(PUSH_PAYLOAD_OK);
		send(source, response);
        } else if(tag == PUSH_PAYLOAD_END) {

		m_completedStoreActors++;

		if(m_completedStoreActors == getSize()) {

			printName();
			cout << "[MatrixOwner] received " << m_receivedPayloads << " payloads" << endl;

			// create directory for Surveyor
			ostringstream matrixFile;
			matrixFile << m_parameters->getPrefix() << "/Surveyor/";
			string surveyorDirectory = matrixFile.str();

			// avoid a race condition
			// we don't know which actor will win the race of life.
			if(!fileExists(surveyorDirectory.c_str())) {
				createDirectory(surveyorDirectory.c_str());
			}


			// print all SimilarityMatrices
			for (map<int,localGramMatrix>::iterator it=m_gramMatrices.begin(); it!=m_gramMatrices.end(); ++it) {

				// matrixFile << "SimilarityMatrix.tsv";
				string matrixNb = static_cast<ostringstream*>( &(ostringstream() << it->first) )->str();

				string similarityMatrix = "";

				if (it->first != -1) {
				        similarityMatrix = (matrixFile.str() + "SimilarityMatrix.filter-" +  matrixNb + ".tsv");
				} else {
					similarityMatrix = (matrixFile.str() + "SimilarityMatrix.global.tsv");
				}

				ofstream similarityFile;

				similarityFile.open(similarityMatrix.c_str());
				printLocalGramMatrix(similarityFile, it->second, m_sampleByFilter[it->first]);
				similarityFile.close();

				printName();
				cout << "[MatrixOwner] printed the Similarity Matrix: ";
				cout << similarityMatrix << endl;
			}


			// normalize the similarity matrix

			normalizeMatrix();

			ostringstream matrixFileForNormalized;
			matrixFileForNormalized << m_parameters->getPrefix() << "/Surveyor/";
			matrixFileForNormalized << "SimilarityMatrix.global.normalized.tsv";

			string normalizedMatrix = matrixFileForNormalized.str();
			ofstream normalizedFile;
			normalizedFile.open(normalizedMatrix.c_str());
			printLocalGramMatrix(normalizedFile, m_normalizedSimilarityMatrix, m_sampleByFilter[-1]);
			normalizedFile.close();

			printName();
			cout << "[MatrixOwner] printed the normalized Similarity Matrix: ";
			cout << normalizedMatrix << endl;


			// create DistanceMatrix

			computeDistanceMatrix();

			ostringstream matrixFileForDistances;
			matrixFileForDistances << m_parameters->getPrefix() << "/Surveyor/";
			matrixFileForDistances << "DistanceMatrix.global.euclidean_raw.tsv";

			string distanceMatrix = matrixFileForDistances.str();
			ofstream distanceFile;
			distanceFile.open(distanceMatrix.c_str());
			printLocalGramMatrix(distanceFile, m_kernelDistanceMatrix, m_sampleByFilter[-1]);
			distanceFile.close();

			printName();
			cout << "[MatrixOwner] printed the Distance Matrix: ";
			cout << distanceMatrix << endl;

			ostringstream matrixFileForNormDistances;
			matrixFileForNormDistances << m_parameters->getPrefix() << "/Surveyor/";
			matrixFileForNormDistances << "DistanceMatrix.global.euclidean_normalized.tsv";

			string normDistanceMatrix = matrixFileForNormDistances.str();
			ofstream normDistanceFile;
			normDistanceFile.open(normDistanceMatrix.c_str());
			printLocalGramMatrix(normDistanceFile, m_normalizedDistanceMatrix, m_sampleByFilter[-1]);
			normDistanceFile.close();

			printName();
			cout << "[MatrixOwner] printed the normalized Distance Matrix: ";
			cout << normDistanceMatrix << endl;


			// tell Mother that the matrix is ready now.
                        Message coolMessage;
                        coolMessage.setTag(GRAM_MATRIX_IS_READY);
                        send(m_mother, coolMessage);


			// clear matrices
		        m_gramMatrices.clear();
			m_kernelDistanceMatrix.clear();
		}
	}
}


// TODO: save time by only computing the lower triangle.
void MatrixOwner::computeDistanceMatrix() {

	// raw matrix
	for(map<SampleIdentifier, map<SampleIdentifier, LargeCount> >::iterator row = m_gramMatrices[-1].begin();
			row != m_gramMatrices[-1].end(); ++row) {

		SampleIdentifier sample1 = row->first;

		for(map<SampleIdentifier, LargeCount>::iterator cell = row->second.begin();
				cell != row->second.end(); ++cell) {

			SampleIdentifier sample2 = cell->first;

			// This is not Euclidean distance .. d(x, x') = sqrt( k(x,x) + k(x', x') - 2 k (x, x'))
			LargeCount distance = 0;
			distance += m_gramMatrices[-1][sample1][sample1];
			distance += m_gramMatrices[-1][sample2][sample2];
			distance -= 2 * m_gramMatrices[-1][sample1][sample2];

			distance = (LargeCount) sqrt((double)distance);

			m_kernelDistanceMatrix[sample1][sample2] = distance;
			m_kernelDistanceMatrix[sample2][sample1] = distance;

		}

	}


	// normalized matrix
	// for(map<SampleIdentifier, map<SampleIdentifier, double> >::iterator row = m_normalizedSimilarityMatrix.begin();
	// 		row != m_normalizedSimilarityMatrix.end(); ++row) {

	// 	SampleIdentifier sample1 = row->first;

	// 	for(map<SampleIdentifier, double>::iterator cell = row->second.begin();
	// 			cell != row->second.end(); ++cell) {

	// 		SampleIdentifier sample2 = cell->first;


	// 		// This is not Euclidean distance .. d(x, x') = sqrt( k(x,x) + k(x', x') - 2 k (x, x'))
	// 		double distance = 0;
	// 		distance += m_normalizedSimilarityMatrix[sample1][sample1];
	// 		distance += m_normalizedSimilarityMatrix[sample2][sample2];
	// 		distance -= 2 * m_normalizedSimilarityMatrix[sample1][sample2];

	// 		distance = (double) sqrt((double)distance);

	// 		m_normalizedDistanceMatrix[sample1][sample2] = distance;
	// 	}

	// }


	for(map<SampleIdentifier, map<SampleIdentifier, double> >::iterator row = m_normalizedSimilarityMatrix.begin();
	    row != m_normalizedSimilarityMatrix.end(); ++row) {

		SampleIdentifier sample1 = row->first;
		double v1[row->second.size()];

		int x = 0;
		for(map<SampleIdentifier, double>::iterator cell = row->second.begin(); cell != row->second.end(); ++cell) {
			v1[x] = double(cell->second);
			x += 1;
		}

		// 		SampleIdentifier sample2 = cell->first;

		for(map<SampleIdentifier, map<SampleIdentifier, double> >::iterator row2 = row;
		    row2 != m_normalizedSimilarityMatrix.end(); ++row2) {

			SampleIdentifier sample2 = row2->first;
			double v2[row->second.size()];

			int y = 0;
			for(map<SampleIdentifier, double>::iterator cell2 = row2->second.begin(); cell2 != row2->second.end(); ++cell2) {
				v2[y] = double(cell2->second);
				y += 1;
			}

			double inner_product = 0;
			for(unsigned int i = 0; i != row->second.size(); i++) {
				double diff = (v1[i]-v2[i]);
				inner_product += (diff*diff);
			}

			m_normalizedDistanceMatrix[sample1][sample2] = sqrt(inner_product);
			m_normalizedDistanceMatrix[sample2][sample1] = sqrt(inner_product);
		}

	}

}


void MatrixOwner::normalizeMatrix() {


	for(map<SampleIdentifier, map<SampleIdentifier, LargeCount> >::iterator row = m_gramMatrices[-1].begin();
			row != m_gramMatrices[-1].end(); ++row) {

		SampleIdentifier sample1 = row->first;

		for(map<SampleIdentifier, LargeCount>::iterator cell = row->second.begin();
				cell != row->second.end(); ++cell) {

			SampleIdentifier sample2 = cell->first;

			double count = (double)(double(m_gramMatrices[-1][sample1][sample2]+1)/sqrt(double(m_gramMatrices[-1][sample1][sample1]+1)*double(m_gramMatrices[-1][sample2][sample2]+1)));

			m_normalizedSimilarityMatrix[sample1][sample2] = count;
		}

	}

}



void MatrixOwner::printLocalGramMatrix(ostream & stream, map<SampleIdentifier, map<SampleIdentifier, LargeCount> > & matrix, vector<int> & samplesToInclude) {

	int numberOfSamples = m_sampleNames->size();

	for(int i = 0 ; i < numberOfSamples ; ++i) {

		if (samplesToInclude[i] == 0)
			continue;

		string & sampleName1 = m_sampleNames->at(i);

		stream << "	" << sampleName1;
	}

	stream << endl;


	for(int i = 0 ; i < numberOfSamples ; ++i) {

		if (samplesToInclude[i] == 0)
			continue;

		string & sampleName1 = m_sampleNames->at(i);

		stream << sampleName1;

		for(int j = 0 ; j < numberOfSamples ; ++j) {

			if (samplesToInclude[j] == 0)
				continue;

			//string & sampleName2 = m_sampleNames->at(j);

			LargeCount hits = 0;

			if(matrix.count(i) > 0 && matrix[i].count(j) > 0) {

				hits = matrix[i][j];
			}

			stream << "	" << hits;
		}

		stream << endl;
	}
}



void MatrixOwner::printLocalGramMatrix(ostream & stream, map<SampleIdentifier, map<SampleIdentifier, double> > & matrix, vector<int> & samplesToInclude) {

	int numberOfSamples = m_sampleNames->size();

	for(int i = 0 ; i < numberOfSamples ; ++i) {

		if (samplesToInclude[i] == 0)
			continue;

		string & sampleName1 = m_sampleNames->at(i);

		stream << "	" << sampleName1;
	}

	stream << endl;


	for(int i = 0 ; i < numberOfSamples ; ++i) {

		if (samplesToInclude[i] == 0)
			continue;

		string & sampleName1 = m_sampleNames->at(i);

		stream << sampleName1;

		for(int j = 0 ; j < numberOfSamples ; ++j) {

			if (samplesToInclude[j] == 0)
				continue;

			//string & sampleName2 = m_sampleNames->at(j);

			double hits = 0;

			if(matrix.count(i) > 0 && matrix[i].count(j) > 0) {

				hits = matrix[i][j];
			}

			stream << "	" << hits;
		}

		stream << endl;
	}
}




/**
 * Write it in RaySurveyorResults/SurveyorMatrix.tsv
 * Also write a distance matrix too !
 */
void MatrixOwner::printLocalGramMatrixWithHash(ostream & stream, map<SampleIdentifier, map<SampleIdentifier, LargeCount> > & matrix) {

	for(map<SampleIdentifier, map<SampleIdentifier, LargeCount> >::iterator column = matrix.begin();
			column != matrix.end(); ++column) {

		SampleIdentifier sample = column->first;

		stream << "	" << m_sampleNames->at(sample);
	}

	stream << endl;

	for(map<SampleIdentifier, map<SampleIdentifier, LargeCount> >::iterator row = matrix.begin();
			row != matrix.end(); ++row) {

		SampleIdentifier sample1 = row->first;

		stream << m_sampleNames->at(sample1);

		for(map<SampleIdentifier, LargeCount>::iterator cell = row->second.begin();
				cell != row->second.end(); ++cell) {

			//SampleIdentifier sample2 = cell->first;

			LargeCount hits = cell->second;

			stream << "	" << hits;
		}

		stream << endl;
	}
}
