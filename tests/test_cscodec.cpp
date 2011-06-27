/*
 * test_cscodec.cpp
 *
 *  Created on: 22/06/2011
 *      Author: David Eccles (gringer) <david.eccles@mpi-muenster.mpg.de>
 *
 *  Checks to make sure colour space encoding / decoding works correctly
 */

#include<format/ColorSpaceCodec.h>

int main(){
	ColorSpaceCodec cd;
	return(cd.check()?0:1);
}
