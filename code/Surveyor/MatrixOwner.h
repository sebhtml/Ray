/*
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

#ifndef MatrixOwnerHeader
#define MatrixOwnerHeader

#include <code/Mock/constants.h>
#include <code/Mock/Parameters.h>

#include <RayPlatform/actors/Actor.h>

#include <map>
#include <iostream>
#include <sstream>
using namespace std;

class MatrixOwner : public Actor {
private:

	int m_receivedPayloads;

	Parameters * m_parameters;
	vector<string> * m_sampleNames;

	// Distance matrix
	map<SampleIdentifier, map<SampleIdentifier, LargeCount> > m_kernelDistanceMatrix;

	// Normalized matrices
	map<SampleIdentifier, map<SampleIdentifier, double> > m_normalizedSimilarityMatrix;
	map<SampleIdentifier, map<SampleIdentifier, double> > m_normalizedDistanceMatrix;

	typedef map<SampleIdentifier, map<SampleIdentifier, LargeCount> > localGramMatrix;
	map <int, localGramMatrix > m_gramMatrices;

	map< int, vector<int> > * m_filterMatrices;
	map <int, vector<int> > m_sampleByFilter;

	int m_mother;
	int m_completedStoreActors;

	void printLocalGramMatrix(ostream & stream, map<SampleIdentifier, map<SampleIdentifier, LargeCount> > & matrix, vector<int> & samplesToInclude);
	void printLocalGramMatrix(ostream & stream, map<SampleIdentifier, map<SampleIdentifier, double> > & matrix, vector<int> & samplesToInclude);

	void printLocalGramMatrixWithHash(ostream & stream, map<SampleIdentifier, map<SampleIdentifier, LargeCount> > & matrix);

	void computeDistanceMatrix();
	void normalizeMatrix();

public:

	MatrixOwner();
	~MatrixOwner();

	void receive(Message & message);

	enum {
		FIRST_TAG = 10300,
		GREETINGS,
		PUSH_PAYLOAD,
		PUSH_PAYLOAD_OK,
		PUSH_PAYLOAD_END,
		GRAM_MATRIX_IS_READY,
		LAST_TAG
	};

};

#endif

