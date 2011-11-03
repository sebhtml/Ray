/*
Ray
Copyright (C)  2011  Sébastien Boisvert

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

#ifndef _constants
#define _constants

#ifndef RAY_VERSION
#define RAY_VERSION "Unknown"
#endif

/*
 * Define the maximum k-mer length when
 * the compiler/make does not.
 */
#ifndef MAXKMERLENGTH
#define MAXKMERLENGTH 32
#endif

#include <stdlib.h> /* for __WORDSIZE hopefully */
#include <stdint.h>

/* exit codes */

/* EXIT_SUCCESS 0 (defined in stdlib.h) */
#define EXIT_NEEDS_ARGUMENTS 5
#define EXIT_NO_MORE_MEMORY 42

/*
 * Include those libraries for Microsoft Visual C++
 */
#ifdef _MSC_VER
#include <xiosbase>
#include <stdexcept>
/* http://msdn.microsoft.com/en-us/library/b0084kay%28VS.80%29.aspx */
#define __func__ __FUNCTION__
#endif

#ifdef FORCE_PACKING
/*
 * With gcc, one can pack data structures.
 */
	#ifdef __GNUC__
		#define ATTRIBUTE_PACKED  __attribute__ ((packed))
/*
 * For Microsoft Visual C++
 */
	#elif defined(_MSC_VER)
		#define ATTRIBUTE_PACKED /* sorry, not available yet */
	#else
		#define ATTRIBUTE_PACKED
	#endif
#else
	#define ATTRIBUTE_PACKED
#endif

#define DUMMY_LIBRARY 40000

#define RAY_NUCLEOTIDE_A 0 /* ~00 == 11 */
#define RAY_NUCLEOTIDE_C 1 /* ~01 == 10 */
#define RAY_NUCLEOTIDE_G 2 /* ~10 == 01 */
#define RAY_NUCLEOTIDE_T 3 /* ~11 == 00 */

#define DOUBLE_ENCODING_A_COLOR '0'
#define DOUBLE_ENCODING_C_COLOR '1'
#define DOUBLE_ENCODING_G_COLOR '2'
#define DOUBLE_ENCODING_T_COLOR '3'


/* the maximum of processes is utilized to construct unique hyperfusions IDs */
#define MAX_NUMBER_OF_MPI_PROCESSES 1000000
#define INVALID_RANK MAX_NUMBER_OF_MPI_PROCESSES

#define MAX_ALLOCATED_OUTPUT_BUFFERS 17

/* maximum value for a uint16_t */
#define RAY_MAXIMUM_READ_LENGTH 65535

#define MAX_VERTICES_TO_VISIT 500
#define TIP_LIMIT 40

/*
 Open-MPI eager threshold is 4k (4096), and this include Open-MPI's metadata.
 tests show that 4096-100 bytes are sent eagerly, too.
 divide that by eight and you get the number of 64-bit integers
 allowed in a single eager communication

 * "4096 is rendezvous. For eager, try 4000 or lower. "
 *  --Eugene Loh  (Oracle)
 *  http://www.open-mpi.org/community/lists/devel/2010/11/8700.php
 *
 */

#define MAXIMUM_MESSAGE_SIZE_IN_BYTES 4000

#define MASTER_RANK 0

/*
 * this is the type used to store coverage values
 */
#define COVERAGE_TYPE uint16_t

/** 32-bit or 64-bit system */

#if defined(__WORDSIZE)
/** use __WORDSIZE */
#define NUMBER_OF_BITS __WORDSIZE

/** assume 64 bits */
/* you may get some compilation warnings about printf and fprintf */
#else
#define NUMBER_OF_BITS 64
#endif

/* 64-bit system */
#if NUMBER_OF_BITS == 64
#define RAY_64_BITS

/* 32-bit system */
#elif NUMBER_OF_BITS == 32
#define RAY_32_BITS

/* assume a 64-bit system */
#else
#define RAY_64_BITS
#endif

#endif
