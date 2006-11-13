#include "Mutator.h"
#include <vector>
#include <cmath>
using namespace std;

gnSeqI uniformSample( gnSeqI min, gnSeqI max );
gnSeqI uniformSample( gnSeqI min, gnSeqI max ){
	gnSeqI sample = rand() % (max - min);
	return sample + min;
}

double exponentialSample( double theta ){
	double z_samp = (double)rand() / (double)RAND_MAX;
	z_samp = log( z_samp );
	return z_samp / theta;
}

gnSeqI poissonSample( double p ){
	gnSeqI count = 0;
	double lambda = 1;
	double sum = 0;
	for( ;; count++ ){
		sum += exponentialSample( lambda );
		if( sum >= p )
			break;
	}
	return count;
}

void IndelInserter::getLocation( gnSeqI& source_start, gnSeqI& source_end, gnSeqI& dest, gnSeqI dest_len ) {
	gnSeqI len = poissonSample( size );
	source_start = uniformSample( 0, donor.alignedSeqsSize() - len - 1 );
	source_end = source_start + len;
	dest = uniformSample( 0, dest_len );
}

void SmallHTInserter::getLocation( gnSeqI& source_start, gnSeqI& source_end, gnSeqI& dest, gnSeqI dest_len ) {
	gnSeqI len = exponentialSample( size );
	source_start = uniformSample( 0, donor.alignedSeqsSize() - len - 1 );
	source_end = source_start + len;
	dest = uniformSample( 0, dest_len );
}

void LargeHTInserter::getLocation( gnSeqI& source_start, gnSeqI& source_end, gnSeqI& dest, gnSeqI dest_len ) {
	gnSeqI len = uniformSample( min_size, max_size );
	source_start = uniformSample( 0, donor.alignedSeqsSize() - len - 1 );
	source_end = source_start + len;
	dest = uniformSample( 0, dest_len );
}

void IndelDeleter::getLocation( gnSeqI& start, gnSeqI& end, gnSeqI dest_len ) {
	gnSeqI len = poissonSample( size );
	start = uniformSample( 0, dest_len - len - 1 );
	end = start + len;
}

void SmallHTDeleter::getLocation( gnSeqI& start, gnSeqI& end, gnSeqI dest_len ) {
	gnSeqI len = exponentialSample( size );
	start = uniformSample( 0, dest_len - len - 1 );
	end = start + len;
}

void LargeHTDeleter::getLocation( gnSeqI& start, gnSeqI& end, gnSeqI dest_len ) {
	gnSeqI len = uniformSample( min_size, max_size );
	start = uniformSample( 0, dest_len - len - 1 );
	end = start + len;
}

void Inverter::getLocation( gnSeqI& start, gnSeqI& end, gnSeqI dest_len ) {
	gnSeqI len = exponentialSample( size );
	start = uniformSample( 0, dest_len - len - 1 );
	end = start + len;
}

void Inserter::mutate( node_id_t nodeI, const PhyloTree& tree, std::vector< gnAlignedSequences >& evolved_alignment ) {

}

void Deleter::mutate( node_id_t nodeI, const PhyloTree& tree, std::vector< gnAlignedSequences >& evolved_alignment ) {

}

void Inverter::mutate( node_id_t nodeI, const PhyloTree& tree, std::vector< gnAlignedSequences >& evolved_alignment ) {

}
