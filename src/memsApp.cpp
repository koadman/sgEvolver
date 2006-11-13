#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <cstdio>
#include <stack>
#include "Matrix.h"
#include "NumericMatrix.h"
#include "getopt.h"
#include "gn/gnSequence.h"
#include "gnAlignedSequences.h"
#include "Mutator.h"

/*
typedef struct mutation_event_s {
	char* name;
	double frequency;
	Mutator& mutator;
} mutation_event_t;
*/

class MutationEvent {
public:
	MutationEvent( const char* mut_name, double mut_frequency, Mutator& mut_mutator );
	MutationEvent( const MutationEvent& me );
	MutationEvent& operator=( const MutationEvent& me );
	string name;
	double frequency;
	Mutator& mutator;
private:
	MutationEvent();
};

MutationEvent::MutationEvent( const char* mut_name, double mut_frequency, Mutator& mut_mutator ) :
name( mut_name ),
frequency( mut_frequency ),
mutator( mut_mutator )
{
}

MutationEvent::MutationEvent( const MutationEvent& me ) :
name( me.name ),
frequency( me.frequency ),
mutator( me.mutator )
{
}

MutationEvent& MutationEvent::operator=( const MutationEvent& me )
{
	name = me.name;
	frequency = me.frequency;
}


class AlignmentEvolver {
public:
	AlignmentEvolver( vector< MutationEvent >& mutation_events );
	void evolve( const PhyloTree& tree, const gnAlignedSequences& destination, vector< gnAlignedSequences >& evolved_alignment );
protected:
	vector< MutationEvent >& mut_events;
};

AlignmentEvolver::AlignmentEvolver( vector< MutationEvent >& mutation_events ) :
mut_events( mutation_events )
{
}

void AlignmentEvolver::evolve( const PhyloTree& tree, const gnAlignedSequences& destination, vector< gnAlignedSequences >& evolved_alignment )
{

	vector< gnAlignedSequences >& ev = evolved_alignment;
	ev.clear();
	ev.push_back( destination );
	
	// scale the mutation frequencies to the tree length,
	// assuming the branch lengths given are expected number of
	// substitutions per site.
	
	double tree_length = tree.getHeight();
	
	// scale mutation rates by tree length
	int mutI;
	for( mutI = 0; mutI < mut_events.size(); mutI++ ){
		mut_events[ mutI ].frequency *= tree_length;
	}

	// get the total mutation rate so we can sample mutation event types from a distribution
	double total_mutation_rate = 0;
	for( mutI = 0; mutI < mut_events.size(); mutI++ )
		total_mutation_rate += mut_events[ mutI ].frequency;
	
	double total_branch_len = 0;
	vector< double > branch_lengths;
	node_id_t nodeI = 0;
	for( nodeI = 0; nodeI < tree.size(); nodeI++ ){
		branch_lengths.push_back( tree[ nodeI ].distance );
		total_branch_len += tree[ nodeI ].distance;
	}

	// calculate the total number of events
	gnSeqI total_events = total_mutation_rate * total_branch_len * destination.alignedSeqsSize();
	
	// normalize branch lengths to a distribution
	for( nodeI = 0; nodeI < tree.size(); nodeI++ )
		branch_lengths[ nodeI ] /= total_branch_len;
	
	// sample branches for events
	vector< vector< uint > > branch_events( branch_lengths.size() );
	gnSeqI eventI;
	double randy, sum;
	uint typeI;
//	srand( time(NULL) );
	srand( 1000 );
	for( eventI = 0; eventI < total_events; eventI++ ){
		randy = (double)rand() / (double)RAND_MAX;
		sum = 0;
		for( nodeI = 0; nodeI < branch_lengths.size(); nodeI++ ){
			sum += branch_lengths[ nodeI ];
			if( sum >= randy )
				break;
		}
		// watch out for floating point imprecision
		if( nodeI == branch_lengths.size() )
			nodeI--;
		
		// so the event will occur on branch 'nodeI'.
		// now sample the event type
		randy = (double)rand() / (double)RAND_MAX;
		sum = 0;
		for( typeI = 0; typeI < mut_events.size(); typeI++ ){
			sum += mut_events[ typeI ].frequency;
			if( sum >= randy )
				break;
		}
		// watch out for floating point imprecision
		if( typeI == mut_events.size() )
			typeI--;

		branch_events[ nodeI ].push_back( typeI );
	}
	
	// now perform mutations in tree-branching order
	stack< node_id_t > node_stack;
	stack< node_id_t > child_stack;
	node_stack.push( tree.root );

	while( node_stack.size() > 0 ){
		if( node_stack.size() != child_stack.size() ){
			// perform mutations on this node
			for( eventI = 0; eventI < branch_events[ node_stack.top() ].size(); eventI++ ){
				mut_events[ branch_events[ node_stack.top() ][ eventI ] ].mutator.mutate( 
						node_stack.top(), tree, evolved_alignment );
			}
			
			// push a child if possible
			if( tree[ node_stack.top() ].children.size() > 0 ){
				node_stack.push( tree[ node_stack.top() ].children[ 0 ] );
				child_stack.push( 0 );
			}else	// done with this leaf node
				node_stack.pop();
			continue;
		}
		if( child_stack.top() != tree[ node_stack.top() ].children.size() ){
			// push a child
			node_stack.push( tree[ node_stack.top() ].children[ child_stack.top() ] );
			child_stack.top()++;
			continue;
		}
		// no more children under this parent node
		node_stack.pop();
		child_stack.pop();
	}
}

void print_usage( const char* pname ){
	cerr << "Usage:" << endl;
	cerr << pname << "[options] <tree file> <alignment file> <evolved output file>" << endl;
	cerr << "Options:" << endl;
	cerr << "\t    --indel-freq=<number> Frequency of indel events" << endl;
	cerr << "\t    --small-ht-freq=<number> Frequency of small horizontal transfer and deletion events" << endl;
	cerr << "\t    --large-ht-freq=<number> Frequency of large horizontal transfer and deletion events" << endl;
	cerr << "\t    --inversion-freq=<number> Frequency of inversion events" << endl;
	cerr << "\t    --indel-size=<number> Poisson parameter for size of indels [3]" << endl;
	cerr << "\t    --small-ht-size=<number> Average length of small horizontal transfer events (sampled from an exponential distribution) [200]" << endl;
	cerr << "\t    --large-ht-min=<number> Minimum length of large horizontal transfer events (sampled uniformly between min and max) [10000]" << endl;
	cerr << "\t    --large-ht-max=<number> Maximum length of large horizontal transfer events (sampled uniformly between min and max) [60000]" << endl;
	cerr << "\t    --inversion-size=<number> Average size of inversion events, samples will be taken from an exponential distribution [50000]" << endl;
	cerr << endl;
}

int viceroy( int argc, char* argv[] );
#define NELEMS(a) ( sizeof( a ) / sizeof( *a ) )
#if defined(__MWERKS__) && defined(__GNDEBUG__)
int main( int argc, char* argv[] ){
#else
int viceroy( int argc, char* argv[] ){
#endif
	char* m_argv[] = {
		"sgEvolver",
		"--indel-freq=.05",
		"--small-ht-freq=.0005",
		"--large-ht-freq=.00001",
		"--inversion-freq=.00001",
		"example.tree",
		"example.dat",
		"evolved.dat",
	};

	int m_argc = NELEMS( m_argv );
#if defined(__MWERKS__) && defined(__GNDEBUG__)
	viceroy( m_argc, m_argv );
#endif
}

/**
 * 
 */
#if defined(__MWERKS__) && defined(__GNDEBUG__)
int viceroy( int argc, char* argv[] ){
#else
int main( int argc, char* argv[] ){
#endif
try{
	if( argc <= 0 ){
		print_usage( "sgEvolver" );
		return -1;
	}
	if( argc == 1 ){
		print_usage( argv[0] );
		return -1;
	}
	
	
	double indel_frequency;
	double small_ht_frequency;
	double large_ht_frequency;
	double inversion_frequency;
	gnSeqI indel_size = 3;
	gnSeqI small_ht_size = 200;
	gnSeqI large_ht_min = 10000;
	gnSeqI large_ht_max = 60000;
	gnSeqI inversion_size = 50000;
	string tree_filename;
	string alignment_filename;
	string output_filename;
	bool print_version = false;
	char* tmp;

	// parse command line with gnu getopt
	int opt;
	int long_opt;
	int ac = argc;
	char** av = argv;
	
	const char* short_args= "";
	struct option long_opts[] = {
		{"indel-freq", required_argument, &long_opt, 1 },
		{"small-ht-freq", required_argument, &long_opt, 2 },
		{"large-ht-freq", required_argument, &long_opt, 3 },
		{"inversion-freq", required_argument, &long_opt, 4 },
		{"indel-size", required_argument, &long_opt, 5 },
		{"small-ht-size", required_argument, &long_opt, 6 },
		{"large-ht-min", required_argument, &long_opt, 7 },
		{"large-ht-max", required_argument, &long_opt, 8 },
		{"inversion-size", required_argument, &long_opt, 9 },
		{"version", no_argument, NULL, 10 },
		{0, 0, 0, 0}	// for correct termination of option list
						// getopt_long can segfault without this
	};

	int indexptr;
	while( (opt = getopt_long( ac, av, short_args, long_opts, &indexptr )) != EOF ){
		switch( opt ){
			case 0:
				switch( long_opt ){
					case 1:
						indel_frequency = strtold( optarg, &tmp );
						break;
					case 2:
						small_ht_frequency = strtold( optarg, &tmp );
						break;
					case 3:
						large_ht_frequency = strtold( optarg, &tmp );
						break;
					case 4:
						inversion_frequency = strtold( optarg, &tmp );
						break;
					case 5:
						indel_size = strtoll( optarg, &tmp, 10 );
						break;
					case 6:
						small_ht_size = strtoll( optarg, &tmp, 10 );
						break;
					case 7:
						large_ht_min = strtoll( optarg, &tmp, 10 );
						break;
					case 8:
						large_ht_max = strtoll( optarg, &tmp, 10 );
						break;
					case 9:
						inversion_size = strtoll( optarg, &tmp, 10 );
						break;
					case 10:
						print_version = true;
						break;
					default:
						print_usage( argv[0] );
						return -1;
						break;
				}
		}
	}
	
	// print the version if the user requested it
	if( print_version ){
		cerr << "sgEvolver " << " build date " << __DATE__ << " at " << __TIME__ << endl;
	}

	if( optind + 2 >= argc ){
		cerr << "You must specify a tree file name, and alignment file name, and an output file name\n";
		return -1;
	}
	tree_filename = av[ optind ];
	alignment_filename = av[ optind + 1];
	output_filename = av[ optind + 2];

	// open the tree file
	ifstream tree_file( tree_filename.c_str() );
	if( !tree_file.is_open() ){
		cerr << "Unable to open " << tree_filename << endl;
		return -1;
	}

	// open the alignment file
	ifstream alignment_file( alignment_filename.c_str() );
	if( !alignment_file.is_open() ){
		cerr << "Unable to open " << alignment_filename << endl;
		return -1;
	}

	// parse the tree file
	PhyloTree tree;
	tree.readTree( tree_file );
	ofstream tree_out( "tree_debug.out" );
	tree.writeTree( tree_out );
	
	// parse the alignment file
	gnAlignedSequences destination;
	gnAlignedSequences donor;
	string line;
	getline( alignment_file, line );
	if( line != "#NEXUS" ){
		cerr << "Error: Alignment file is not in SeqGen NEXUS format.\n";
		return -1;
	}
	destination.constructFromRelaxedNexus( alignment_file );
	donor.constructFromRelaxedNexus( alignment_file );


	// initialize the mutators
	IndelInserter ii( donor, indel_size );
	IndelDeleter id( donor, indel_size );
	SmallHTInserter shti( donor, small_ht_size );
	SmallHTDeleter shtd( donor, small_ht_size );
	LargeHTInserter lhti( donor, large_ht_min, large_ht_max );
	LargeHTDeleter lhtd( donor, large_ht_min, large_ht_max );
	Inverter inv( donor, inversion_size );
/*	mutation_event_t mut_events[] = {
		"insert", indel_frequency, ii,
		"deletion", indel_frequency, id, 
		"small_ht_insert", small_ht_frequency, shti, 
		"small_ht_deletion", small_ht_frequency, shtd,
		"large_ht_insert", large_ht_frequency, lhti,
		"large_ht_deletion", large_ht_frequency, lhtd,
		"inversion", inversion_frequency, inv,
	};
*/
	vector< MutationEvent > mutation_events;
	mutation_events.push_back( MutationEvent( "indel_insert", indel_frequency, ii ) );
	mutation_events.push_back( MutationEvent( "indel_deletion", indel_frequency, id ) );
	mutation_events.push_back( MutationEvent( "small_ht_insert", small_ht_frequency, shti ) );
	mutation_events.push_back( MutationEvent( "small_ht_deletion", small_ht_frequency, shtd ) );
	mutation_events.push_back( MutationEvent( "large_ht_insert", large_ht_frequency, lhti ) );
	mutation_events.push_back( MutationEvent( "large_ht_deletion", large_ht_frequency, lhtd ) );
	mutation_events.push_back( MutationEvent( "inversion", inversion_frequency, inv ) );
	
	// start evolving the sequence alignment
	AlignmentEvolver ae( mutation_events );
	
	vector< gnAlignedSequences > evolved_alignment;
	ae.evolve( tree, destination, evolved_alignment );
	
	
}catch( gnException& gne ) {
	cerr << "Unhandled gnException: " << gne << endl;
}catch( exception& e ) {
	cerr << "Unhandled exception: " << e.what() << endl;
}catch(...){
	cerr << "Unknown exception occurred.\n";
}

	return 0;
}
