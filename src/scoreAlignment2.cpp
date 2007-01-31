/*******************************************************************************
 * $Id: scoreAlignment.cpp,v 1.14 2004/02/28 00:01:31 darling Exp $
 * This file is copyright 2002-2004 Aaron Darling.  All rights reserved.
 * Please see the file called COPYING for licensing, copying, and modification
 * rights.  Redistribution of this file, in whole or in part is prohibited
 * without express permission.
 ******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include "libGenome/gnFilter.h"
#include "libMems/IntervalList.h"
#include "libMems/MemSubsets.h"
#include "libMems/MatchList.h"
#include "libMems/GappedAlignment.h"
#include "libMems/Matrix.h"
#include "libMems/MatchProjectionAdapter.h"
#include "libMems/Aligner.h"

using namespace std;
using namespace genome;
using namespace mems;

// basic data structures

//
// baaad:  copied functions from ProgressiveAligner.h
//

template< typename PairType >
class LabelSort 
{
public:
	LabelSort( uint seqI ) : ssc( seqI ) {};
	bool operator()( const PairType& pt1, const PairType& pt2 )
	{
		return ssc( pt1.first, pt2.first );
	}
private:
	LabelSort();
	mems::SSC<mems::AbstractMatch> ssc;
};

template<class MatchVector>
void IdentifyBreakpoints( MatchVector& mlist, std::vector<gnSeqI>& breakpoints )
{
	if( mlist.size() == 0 )
		return;
	breakpoints = std::vector<gnSeqI>(1, mlist.size()-1);

	mems::SSC<mems::AbstractMatch> ssc(0);
	std::sort( mlist.begin(), mlist.end(), ssc );
	typedef typename MatchVector::value_type value_type;
	typedef std::pair< value_type, size_t > LabelPairType;
	std::vector< LabelPairType > label_list;
	typename MatchVector::iterator cur = mlist.begin();
	typename MatchVector::iterator end = mlist.end();
	size_t i = 0;
	for( ;cur != end; ++cur )
	{
		label_list.push_back( std::make_pair( *cur, i ) );
		++i;
	}

	uint seq_count = mlist[0]->SeqCount();
	// check for breakpoints in each sequence
	for( uint seqI = 1; seqI < seq_count; seqI++ )
	{
		LabelSort< LabelPairType > ls(seqI); 
		std::sort( label_list.begin(), label_list.end(), ls );

		typename std::vector< LabelPairType >::const_iterator prev = label_list.begin();
		typename std::vector< std::pair< typename MatchVector::value_type, size_t > >::const_iterator iter = label_list.begin();
		typename std::vector< std::pair< typename MatchVector::value_type, size_t > >::const_iterator lab_end = label_list.end();

		bool prev_orient = (*prev).first->Orientation(seqI) == (*prev).first->Orientation(0);
		if( !prev_orient )	// if we start in a different orientation than the ref seq there's a bp here
			breakpoints.push_back(prev->second);

		for( ++iter; iter != lab_end; ++iter )
		{
			bool cur_orient = (*iter).first->Orientation(seqI) == (*iter).first->Orientation(0);
			if( prev_orient == cur_orient &&
				( ( prev_orient && (*prev).second + 1 == (*iter).second) ||
				  ( !prev_orient && (*prev).second - 1 == (*iter).second) 
				)
			  )
			{
				prev_orient = cur_orient;
				++prev;
				continue;	// no breakpoint here
			}

			// always add the last match in a new block (scanning from left to right in seq 0)
			if( prev_orient )
				breakpoints.push_back( prev->second );
			if( !cur_orient )
				breakpoints.push_back( iter->second );

			prev_orient = cur_orient;
			++prev;
		}
		if( prev_orient )
			breakpoints.push_back( prev->second );
	}
	std::sort( breakpoints.begin(), breakpoints.end() );
	std::vector<gnSeqI>::iterator uni = std::unique( breakpoints.begin(), breakpoints.end() );
	breakpoints.erase( uni, breakpoints.end() );
}


template< class MatchVector >
void ComputeLCBs_v2( const MatchVector& meml, const std::vector<gnSeqI>& breakpoints, std::vector< MatchVector >& lcb_list )
{
	// there must be at least one end of a block defined
	if( breakpoints.size() < 1 )
		return;
		
	lcb_list.clear();
	
	// organize the LCBs into different MatchVector instances
	std::vector<gnSeqI>::const_iterator break_iter = breakpoints.begin();
	uint prev_break = 0;	// prev_break is the first match in the current block
	MatchVector lcb;
	for( ; break_iter != breakpoints.end(); ++break_iter ){
		// add the new MatchList to the set if it made the cut
		lcb_list.push_back( lcb );
		lcb_list.back().insert( lcb_list.back().end(), meml.begin() + prev_break, meml.begin() + *break_iter + 1 );
		prev_break = *break_iter + 1;
	}
}


template <class MatchVector>
void computeLCBAdjacencies_v3( const std::vector< MatchVector >& lcb_list, std::vector< double >& weights, std::vector< mems::LCB >& adjacencies )
{
	adjacencies.clear(); // start with no LCB adjacencies
	if( lcb_list.size() == 0 )
		return;	// there aren't any LCBs so there aren't any adjacencies!

	uint seq_count = lcb_list.front().front()->SeqCount();
	uint seqI;
	uint lcbI;
	for( lcbI = 0; lcbI < lcb_list.size(); ++lcbI ){
		mems::LCB lcb;
		std::vector<gnSeqI> left_end;
		std::vector<gnSeqI> length;
		std::vector<bool> orientation;
		FindBoundaries( lcb_list[lcbI], left_end, length, orientation );

		lcb.left_adjacency = std::vector<uint>( left_end.size(), -1 );
		lcb.right_adjacency = std::vector<uint>( left_end.size(), -1 );
		lcb.left_end = std::vector<int64>( left_end.size(), 0 );
		lcb.right_end = std::vector<int64>( left_end.size(), 0 );

		for( seqI = 0; seqI < seq_count; seqI++ ){
			// support "ragged edges" on the ends of LCBs
			if( left_end[seqI] == mems::NO_MATCH )
				continue;
			lcb.left_end[seqI] = left_end[seqI];
			lcb.right_end[seqI] = left_end[seqI] + length[seqI];
			if( !orientation[seqI] )
			{
				lcb.left_end[seqI] = -lcb.left_end[seqI];
				lcb.right_end[seqI] = -lcb.right_end[seqI];
			}
		}
		lcb.lcb_id = adjacencies.size();
		lcb.weight = weights[ lcbI ];
		lcb.to_be_deleted = false;
		adjacencies.push_back( lcb );
	}

	for( seqI = 0; seqI < seq_count; seqI++ ){
		mems::LCBLeftComparator llc( seqI );
		std::sort( adjacencies.begin(), adjacencies.end(), llc );
		for( lcbI = 1; lcbI + 1 < lcb_list.size(); lcbI++ ){
			adjacencies[ lcbI ].left_adjacency[ seqI ] = adjacencies[ lcbI - 1 ].lcb_id;
			adjacencies[ lcbI ].right_adjacency[ seqI ] = adjacencies[ lcbI + 1 ].lcb_id;
		}
		if( lcbI == lcb_list.size() )
			lcbI--;	// need to decrement when there is only a single LCB

		// set first and last lcb adjacencies to -1
		adjacencies[ 0 ].left_adjacency[ seqI ] = (uint)-1;
		adjacencies[ lcbI ].right_adjacency[ seqI ] = (uint)-1;
		if( lcbI > 0 ){
			adjacencies[ 0 ].right_adjacency[ seqI ] = adjacencies[ 1 ].lcb_id;
			adjacencies[ lcbI ].left_adjacency[ seqI ] = adjacencies[ lcbI - 1 ].lcb_id;
		}
	}
	mems::LCBIDComparator lic;
	std::sort( adjacencies.begin(), adjacencies.end(), lic );

}

//
// end baaad
//

/** store a pair of aligned positions and the characters */
typedef struct aligned_coords_s {
	int64 pos1;
	int64 pos2;
	char char1;
	char char2;
} aligned_coords_t;


class AlignedCoordSeqIComparator {
public:
	boolean operator()(const aligned_coords_t& a, const aligned_coords_t& b) const{
		if( absolut(a.pos1) == absolut(b.pos1) )
			return absolut(a.pos2) < absolut(b.pos2);
		return absolut(a.pos1) < absolut(b.pos1);
	}
};


void constructCoordList( uint seqI, uint seqJ, IntervalList& iv_list, vector< aligned_coords_t >& coord_list, vector< gnSequence* >& seq_table ){

	//
	// pre-allocate the vector
	//
	gnSeqI ij_vec_size = 0;
	for( int ivI = 0; ivI < iv_list.size(); ivI++ ){
		ij_vec_size += iv_list[ivI].AlignmentLength();
	}
	coord_list = vector< aligned_coords_t >( ij_vec_size );

	//
	// fill in the vector with all aligned pairs
	//
	gnSeqI vecI = 0;	// current place in vector
	for( int ivI = 0; ivI < iv_list.size(); ivI++ ){
		GappedAlignment* aln;
		aln = dynamic_cast< GappedAlignment* >( iv_list[ ivI ].GetMatches()[0] );
		if( aln == NULL ){
			throw "Error:  expecting interval to contain a single GappedAlignment";
		}
		int64 pos1 = aln->Start( seqI );
		int64 pos2 = aln->Start( seqJ );
		// if rev. comp then we're starting at the other (right) side
		if( pos1 < 0 )
			pos1 -= aln->Length( seqI ) - 1;
		if( pos2 < 0 )
			pos2 -= aln->Length( seqJ ) - 1;

		const std::vector< std::string >& align_matrix = GetAlignment( *aln, seq_table );
		for( gnSeqI colI = 0; colI < aln->Length(); colI++ ){
			aligned_coords_t act;
			act.char1 = align_matrix[ seqI ][ colI ];
			act.char2 = align_matrix[ seqJ ][ colI ];
			act.pos1 = act.char1 == '-' ? 0 : pos1;
			act.pos2 = act.char2 == '-' ? 0 : pos2;
			coord_list[ vecI++ ] = act;
			if( act.char1 != '-' )
				pos1++;
			if( act.char2 != '-' )
				pos2++;
		}
	}

	//
	// sort the vector on aligned position
	//
	AlignedCoordSeqIComparator acsc;
	sort( coord_list.begin(), coord_list.end(), acsc );
}


const gnFilter* comp_filter = gnFilter::DNAComplementFilter();

void sanityCheck( string& aln_name, uint seqI, uint seqJ, aligned_coords_t& act, vector< string >& evolved_seqs ){
	if( act.pos1 != 0 ){
		char seqI_char = evolved_seqs[ seqI ][ absolut( act.pos1 ) - 1 ];
		if( act.pos1 < 0 )
			seqI_char = comp_filter->Filter( seqI_char );
		if( act.char1 != seqI_char ){
			cerr << "Error: " << aln_name << " character " << act.char1 << " instead of " << seqI_char << " at position " << act.pos1 << endl;
		}
	}

	if( act.pos2 != 0 ){
		char seqJ_char = evolved_seqs[ seqJ ][ absolut( act.pos2 ) - 1 ];
		if( act.pos2 < 0 )
			seqJ_char = comp_filter->Filter( seqJ_char );

		if( act.char2 != seqJ_char ){
			cerr << "Error: " << aln_name << " character " << act.char2 << " instead of " << seqJ_char << " at position " << act.pos2 << endl;
		}
	}
}

bool warn_missing = true;
int warn_count = 0;
const int WARN_MAX = 1000;

void compareAlignments( IntervalList& correct, IntervalList& calculated, vector< string >& evolved_seqs, vector< gnSequence* >& seq_table ){
	
	string cor_name = "correct";
	string calc_name = "calculated";
	gnSeqI sp_truepos = 0;
	gnSeqI sp_falsepos = 0;
	gnSeqI sp_trueneg = 0;
	gnSeqI sp_falseneg = 0;
	gnSeqI sp_extraneg = 0;
	gnSeqI wrong_strand = 0;

	// calculate total possible score
	gnSeqI sp_possible = 0;

	uint seqI = 0;
	uint seqJ = 0;

	// now for every entry in the correct alignment, look for a corresponding
	// entry in the calculated alignment.  Do sanity checking while we're at it
	for( seqI = 0; seqI < seq_table.size(); seqI++ ){
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ ){
			//
			// the correct coordinate matrix is expected to have one entry for
			// each character in seqI, plus an arbitrary number of "zero" entries
			// for other sequences (seqJ > seqI) that couldn't align to anything
			// in seqI
			//
			vector< aligned_coords_t > cor;
			constructCoordList( seqI, seqJ, correct, cor, seq_table );

			//
			// the calculated alignment is expected to account for every residue
			// at least once, either aligning it to another residue or a gap
			//
			vector< aligned_coords_t > calc;
			constructCoordList( seqI, seqJ, calculated, calc, seq_table );

			gnSeqI corI = 0;
			gnSeqI calcI = 0;

			// skip any gaps aligned to gaps
			while( corI < cor.size() && cor[ corI ].pos1 == 0 && cor[ corI ].pos2 == 0 )
				corI++;
			while( calcI < calc.size() && calc[ calcI ].pos1 == 0 && calc[ calcI ].pos2 == 0 )
				calcI++;

			//
			// get thru the unaligned region of seqI (where pos1 is zero)
			//
			while( corI < cor.size() && calcI < calc.size() && 
				cor[ corI ].pos1 == 0 && calc[ calcI ].pos1 == 0 )
			{
				// sanity check the "correct" aligned characters...
				sanityCheck( cor_name, seqI, seqJ, cor[ corI ], evolved_seqs );

				while( calcI < calc.size() && calc[ calcI ].pos1 == 0 &&
					absolut(calc[ calcI ].pos2) < absolut(cor[ corI ].pos2) )
				{
					// sanity check the aligned characters...
					sanityCheck( calc_name, seqI, seqJ, calc[ calcI ], evolved_seqs );
					// these are all false positives, catch them later
					calcI++;
				}

				uint match_count = 0;
				while( calcI < calc.size() && calc[ calcI ].pos1 == 0 &&
					absolut(calc[ calcI ].pos2) == absolut(cor[ corI ].pos2) ){
					// sanity check the aligned characters...
					sanityCheck( calc_name, seqI, seqJ, calc[ calcI ], evolved_seqs );
					match_count++;
					calcI++;
				}
				sp_trueneg += match_count > 0 ? 1 : 0;
				sp_extraneg += match_count > 0 ? match_count - 1 : 0;
				// catch false positives later
				corI++;
			}

			// scan thru any remaining false negatives
			// that were in the calculated alignment
			while( calcI < calc.size() && calc[ calcI ].pos1 == 0 ){
				// sanity check the aligned characters...
				sanityCheck( calc_name, seqI, seqJ, calc[ calcI ], evolved_seqs );
				calcI++;
				sp_falseneg++;
			}
			while( corI < cor.size() && cor[ corI ].pos1 == 0 ){
				sanityCheck( cor_name, seqI, seqJ, cor[ corI ], evolved_seqs );
				// this is a false positive...
				corI++;
			}
			
			//
			// now look at the rest of seqI
			//
			for( ; corI < cor.size(); corI++ ){
				// sanity check the "correct" aligned characters
				sanityCheck( cor_name, seqI, seqJ, cor[ corI ], evolved_seqs );

				// if this is an aligned pair it gets counted towards
				// the total possible points in the sum of pairs score
				if( cor[ corI ].pos2 != 0 )
					sp_possible++;

				// make sure calc is up to where cor is...
				while( calcI < calc.size() && (absolut(calc[ calcI ].pos1) < absolut(cor[ corI ].pos1)) ){
					// sanity check the aligned characters
					sanityCheck( calc_name, seqI, seqJ, calc[ calcI ], evolved_seqs );

					// cor skipped more than one char and is missing something...
					cerr << "correct alignment is missing characters\n";
					calcI++;
				}

				// make sure cor is up to where calc is...
				if( calcI == calc.size() || (absolut(calc[ calcI ].pos1) > absolut(cor[ corI ].pos1)) ){
					// calc is missing something!
					// assume it's a gap
					if( warn_missing )
					{
						cerr << "Warning, calculated alignment missing prediction for seqI " << seqI << " position " << absolut(cor[ corI ].pos1) << ".  Assuming gap prediction\n";
						warn_count++;
						if( warn_count > WARN_MAX )
						{
							cerr << "Too many warnings... no longer reporting\n";
							warn_missing = false;
						}
					}
					if( cor[ corI ].pos2 == 0 )
						sp_trueneg++;
					else
						sp_falseneg++;
					continue;
				}

				// if all is well, we should have found the corresponding entries in
				// cor and calc
				int cur_count = 0;
				// examing all predictions for this residue in the calculated alignment...
				while( calcI < calc.size() && (absolut(calc[ calcI ].pos1) == absolut(cor[ corI ].pos1)) ){
					// sanity check the aligned characters...
					sanityCheck( calc_name, seqI, seqJ, calc[ calcI ], evolved_seqs );
					
					// check the strands that were aligned
					bool cor_parity = ((cor[ corI ].pos1 > 0 && cor[ corI ].pos2 > 0) ||
										(cor[ corI ].pos1 < 0 && cor[ corI ].pos2 < 0)) ?
										true : false;

					bool calc_parity = ((calc[ calcI ].pos1 > 0 && calc[ calcI ].pos2 > 0) ||
										(calc[ calcI ].pos1 < 0 && calc[ calcI ].pos2 < 0)) ?
										true : false;

					// if they match then increment the sp_score
					if( absolut(calc[ calcI ].pos2) == absolut(cor[ corI ].pos2) ){
						if( calc[ calcI ].pos2 == 0 )
							sp_trueneg++;	// correct alignment to a gap
						else if( cor_parity == calc_parity )
							sp_truepos++;	// correct positions on correct strand!!
						else
						{
							sp_falsepos++;	// correct positions, wrong strand
							wrong_strand++;
						}
					}else{
						if( calc[ calcI ].pos2 == 0 )
							sp_falseneg++;	// incorrect alignment to a gap
						else
							sp_falsepos++;	// incorrectly aligned

					}
					calcI++;
				}
			}
		}
	}
	
	
	double sp_accuracy = (double)sp_truepos / (double)sp_possible;
	cout << "Sum of pairs accuracy: " << sp_accuracy << endl;
	cout << "Sum of pairs error rate: " << (double)sp_falsepos / (double)sp_truepos << endl;
	cout << "The error rate gives the number of incorrect orthology predictions per correct prediction." << endl;
	cout << "Sum of pairs positive predictive value: " << (double)sp_truepos / (double)(sp_truepos + sp_falsepos) << endl;
	cout << "trueneg: " << sp_trueneg << "  falseneg: " << sp_falseneg << "  extraneg: " << sp_extraneg << endl;
	cout << "wrong_strand: " << wrong_strand << endl;
}








template< typename FwdIt >
bool findFirstCorrectlyAlignedPair( 
	const FwdIt& cor_first,
	const FwdIt& cor_last,
	const FwdIt& calc_first,
	const FwdIt& calc_last,
	FwdIt& cor_iter,
	FwdIt& calc_iter,
	int dir
	)
{
	while( dir == -1 || dir == 1 && calc_iter != calc_last && cor_iter != cor_last )
	{
		if( genome::absolut(calc_iter->pos1) < genome::absolut(cor_iter->pos1) )
		{
			if( dir == 1 )
				++calc_iter;
			else if( cor_iter == cor_last )
				return false;
			else
				--cor_iter;
			continue;
		}
		if( genome::absolut(calc_iter->pos1) > genome::absolut(cor_iter->pos1) )
		{
			if( dir == 1 )
				++cor_iter;
			else if( calc_iter == calc_last )
				return false;
			else
				--calc_iter;
			continue;
		}
		if( genome::absolut(calc_iter->pos2) != genome::absolut(cor_iter->pos2) ||
			calc_iter->pos2 == 0)
		{
			if( calc_iter == calc_last || cor_iter == cor_last )
				return false;
			calc_iter += dir;
			cor_iter += dir;
			continue;
		}
		bool cor_a = cor_iter->pos1 < 0;
		bool cor_b = cor_iter->pos2 < 0;
		bool cor_parity = cor_a == cor_b;
		bool calc_a = calc_iter->pos1 < 0;
		bool calc_b = calc_iter->pos2 < 0;
		bool calc_parity = calc_a == calc_b;
		if( cor_parity != calc_parity )
		{
			// calc has the opposite strand aligned
			++calc_iter;
			++cor_iter;
			continue;
		}
		break;	// these are aligned correctly
	}
	if( dir == 1 && (cor_iter == cor_last || calc_iter == calc_last) )
	{
		calc_iter = calc_last;
		return false;
	}
}

size_t findLcb( vector< LCB >& calc_pair_adj, int64 pos )
{
	size_t calc_ivI = 0;
	for( ; calc_ivI < calc_pair_adj.size(); ++calc_ivI )
	{
		if( calc_pair_adj[calc_ivI].left_end[0] == NO_MATCH )
			continue;
		if( calc_pair_adj[calc_ivI].left_end[0] <= genome::absolut( pos ) &&
			genome::absolut( pos ) < calc_pair_adj[calc_ivI].right_end[0] )
			break;
	}
	return calc_ivI;
}


void markAllCorrectlyAlignedLcbs( 
	const vector< aligned_coords_t >::iterator& cor_first,
	const vector< aligned_coords_t >::iterator& cor_last,
	const vector< aligned_coords_t >::iterator& calc_first,
	const vector< aligned_coords_t >::iterator& calc_last,
	vector< LCB >& calc_pair_adj,
	boost::dynamic_bitset<>& found_lcbs
	)
{

	vector< aligned_coords_t >::iterator cor_iter = cor_first;
	vector< aligned_coords_t >::iterator calc_iter = calc_first;
	while( calc_iter != calc_last && cor_iter != cor_last )
	{
		if( genome::absolut(calc_iter->pos1) < genome::absolut(cor_iter->pos1) )
		{
			++calc_iter;
			continue;
		}
		if( genome::absolut(calc_iter->pos1) > genome::absolut(cor_iter->pos1) )
		{
			++cor_iter;
			continue;
		}
		if( genome::absolut(calc_iter->pos2) != genome::absolut(cor_iter->pos2)  ||
			calc_iter->pos2 == 0)
		{
			++calc_iter;
			++cor_iter;
			continue;
		}
		bool cor_a = cor_iter->pos1 < 0;
		bool cor_b = cor_iter->pos2 < 0;
		bool cor_parity = cor_a == cor_b;
		bool calc_a = calc_iter->pos1 < 0;
		bool calc_b = calc_iter->pos2 < 0;
		bool calc_parity = calc_a == calc_b;
		if( cor_parity != calc_parity )
		{
			// calc has the opposite strand aligned
			++calc_iter;
			++cor_iter;
			continue;
		}
		// these are aligned correctly

		size_t calc_ivI = findLcb( calc_pair_adj, calc_iter->pos1 );
		found_lcbs[calc_ivI] = true;
		// move calc_iter to the next LCB
		aligned_coords_t next_act;
		next_act.pos1 = calc_pair_adj[calc_ivI].right_end[0] + 1;
		AlignedCoordSeqIComparator acsi;
		calc_iter = std::lower_bound( calc_first, calc_last, next_act, acsi );
	}
}


void computeLCBaccuracy( IntervalList& correct, IntervalList& calculated, vector< string >& evolved_seqs, vector< gnSequence* >& seq_table )
{	
	string cor_name = "correct";
	string calc_name = "calculated";
	double bp_sp_truepos = 0;
	double bp_sp_falsepos = 0;
	double bp_sp_falseneg = 0;
	uint seqI = 0;
	uint seqJ = 0;


	// find all pairwise LCBs
	Matrix< vector< LCB > > cor_pairwise_adjs(seq_table.size(),seq_table.size());
	Matrix< vector< LCB > > calc_pairwise_adjs(seq_table.size(),seq_table.size());
	for( seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ )
		{
			vector< size_t > projection( 2 );
			projection[0] = seqI;
			projection[1] = seqJ;

			// construct pairwise Interval projections
			vector< AbstractMatch* > cor_mpa_list;
			vector< AbstractMatch* > calc_mpa_list;
			for( size_t corI = 0; corI < correct.size(); corI++ )
			{
				if( correct[corI].LeftEnd(seqI) == NO_MATCH || correct[corI].LeftEnd(seqJ) == NO_MATCH )
					continue;
				MatchProjectionAdapter mpa_tmp( &correct[corI], projection );
				cor_mpa_list.push_back( mpa_tmp.Copy() );
				if( cor_mpa_list.back()->Orientation(0) == AbstractMatch::reverse )
					cor_mpa_list.back()->Invert();
			}
			for( size_t calcI = 0; calcI < calculated.size(); calcI++ )
			{
				if( calculated[calcI].LeftEnd(seqI) == NO_MATCH || calculated[calcI].LeftEnd(seqJ) == NO_MATCH )
					continue;
				MatchProjectionAdapter mpa_tmp( &calculated[calcI], projection );
				calc_mpa_list.push_back( mpa_tmp.Copy() );
				if( calc_mpa_list.back()->Orientation(0) == AbstractMatch::reverse )
					calc_mpa_list.back()->Invert();
			}
			vector< vector< AbstractMatch* > > LCB_list;
			vector< gnSeqI > breakpoints;
			IdentifyBreakpoints( cor_mpa_list, breakpoints );
			ComputeLCBs_v2( cor_mpa_list, breakpoints, LCB_list );
			vector< double > lcb_scores( LCB_list.size(), 0 );
			computeLCBAdjacencies_v3( LCB_list, lcb_scores, cor_pairwise_adjs(seqI,seqJ) );

			breakpoints.clear();
			LCB_list.clear();
			IdentifyBreakpoints( calc_mpa_list, breakpoints );
			ComputeLCBs_v2( calc_mpa_list, breakpoints, LCB_list );
			lcb_scores = vector< double >( LCB_list.size(), 0 );
			computeLCBAdjacencies_v3( LCB_list, lcb_scores, calc_pairwise_adjs(seqI,seqJ) );
		}
	}

	// calculate total possible score
	double bp_sp_possible = 0;
	double bp_dist_sum = 0;
	double bp_dist_count = 0;
	vector< double > bp_dist;

	Matrix< vector<double> > left_dists( seq_table.size(), seq_table.size() );
	Matrix< vector<double> > right_dists( seq_table.size(), seq_table.size() );

	Matrix< boost::dynamic_bitset<> > found_lcbs( seq_table.size(), seq_table.size() );
	for( seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ )
		{
			boost::dynamic_bitset<> tmp_bs ( calc_pairwise_adjs(seqI,seqJ).size(), false );	// tracks whether calculated LCBs were found in the correct alignment
			found_lcbs(seqI,seqJ) = tmp_bs;
			vector<double> asdf( calc_pairwise_adjs(seqI,seqJ).size(), (std::numeric_limits<double>::max)() );
			left_dists(seqI, seqJ) = asdf;
			left_dists(seqJ, seqI) = asdf;
			right_dists(seqI, seqJ) = asdf;
			right_dists(seqJ, seqI) = asdf;
		}
	}

	// now for every entry in the correct alignment, look for a corresponding
	// entry in the calculated alignment.  Do sanity checking while we're at it
	for( seqI = 0; seqI < seq_table.size(); seqI++ ){
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ ){
			//
			// the correct coordinate matrix is expected to have one entry for
			// each character in seqI, plus an arbitrary number of "zero" entries
			// for other sequences (seqJ > seqI) that couldn't align to anything
			// in seqI
			//
			vector< aligned_coords_t > cor;
			constructCoordList( seqI, seqJ, correct, cor, seq_table );

			//
			// the calculated alignment is expected to account for every residue
			// at least once, either aligning it to another residue or a gap
			//
			vector< aligned_coords_t > calc;
			constructCoordList( seqI, seqJ, calculated, calc, seq_table );

			gnSeqI corI = 0;
			gnSeqI calcI = 0;

			vector< LCB >& cor_pair_adj = cor_pairwise_adjs(seqI,seqJ);
			vector< LCB >& calc_pair_adj = calc_pairwise_adjs(seqI,seqJ);
			for( size_t cor_ivI = 0; cor_ivI < cor_pair_adj.size(); cor_ivI++ )
			{
				bp_sp_possible++;

				int64 lend_seqI = cor_pair_adj[cor_ivI].left_end[0];
				int64 rend_seqI = cor_pair_adj[cor_ivI].right_end[0];

				aligned_coords_t act_left;
				aligned_coords_t act_right;
				act_left.pos1 = lend_seqI;
				act_left.pos2 = (std::numeric_limits<int64>::min)();
				act_right.pos1 = rend_seqI;
				act_right.pos2 = (std::numeric_limits<int64>::max)();
				AlignedCoordSeqIComparator acsi;
				vector< aligned_coords_t >::iterator cor_first;
				vector< aligned_coords_t >::iterator cor_last;
				cor_first = std::lower_bound( cor.begin(), cor.end(), act_left, acsi );
				cor_last = std::upper_bound( cor.begin(), cor.end(), act_right, acsi );

				vector< aligned_coords_t >::iterator calc_first;
				vector< aligned_coords_t >::iterator calc_last;
				calc_first = std::lower_bound( calc.begin(), calc.end(), act_left, acsi );
				calc_last = std::upper_bound( calc.begin(), calc.end(), act_right, acsi );

				// find the left-most LCB that has a correctly aligned character
				vector< aligned_coords_t >::iterator calc_iter = calc_first;
				vector< aligned_coords_t >::iterator cor_iter = cor_first;
				findFirstCorrectlyAlignedPair( cor_first, cor_last, calc_first, calc_last, cor_iter, calc_iter, 1 );
				if( calc_iter == calc_last )
				{
					bp_sp_falseneg++;	// didn't align anything correctly in this region
					continue;
				}
				bp_sp_truepos++;	// found this LCB!

				markAllCorrectlyAlignedLcbs( cor_first, cor_last, calc_first, calc_last, calc_pair_adj, found_lcbs(seqI,seqJ) );

				// determine which LCB in the calculated ivs we're inside of
				size_t calc_ivI = findLcb( calc_pair_adj, calc_iter->pos1 );
				bp_dist.push_back( (double)lend_seqI - (double)calc_pair_adj[calc_ivI].left_end[0] );
				bp_dist_sum += bp_dist.back();
				bp_dist_count++;
				if( genome::absolut( bp_dist.back() ) < genome::absolut( left_dists(seqI,seqJ)[calc_ivI] ) )
					left_dists(seqI,seqJ)[calc_ivI] = bp_dist.back();

				bool parity = (cor_pair_adj[cor_ivI].left_end[1] < 0) == (calc_pair_adj[calc_ivI].left_end[1] < 0);
				if( parity )
					bp_dist.push_back( (double)cor_pair_adj[cor_ivI].left_end[1] - (double)calc_pair_adj[calc_ivI].left_end[1] );
				else
				{
					bp_dist.push_back( (double)calc_pair_adj[calc_ivI].right_end[1] + (double)cor_pair_adj[cor_ivI].right_end[1] );
					if( calc_pair_adj[calc_ivI].left_end[1] < 0 )
						bp_dist.back() *= -1;
				}
				bp_dist_sum += bp_dist.back();
				bp_dist_count++;

				if( genome::absolut( bp_dist.back() ) < genome::absolut( left_dists(seqJ,seqI)[calc_ivI] ) )
					left_dists(seqJ,seqI)[calc_ivI] = bp_dist.back();

				// now do the same for the other side
				// determine how far away the end of the LCB prediction is
				vector< aligned_coords_t >::iterator rcor_first;
				vector< aligned_coords_t >::iterator rcor_last;
				rcor_first = cor_last - 1;
				rcor_last = cor_first;

				vector< aligned_coords_t >::iterator rcalc_first;
				vector< aligned_coords_t >::iterator rcalc_last;
				rcalc_first = calc_last - 1;
				rcalc_last = calc_first;

				// find the left-most LCB that has a correctly aligned character
				vector< aligned_coords_t >::iterator rcalc_iter = rcalc_first;
				vector< aligned_coords_t >::iterator rcor_iter = rcor_first;
				bool found = findFirstCorrectlyAlignedPair( rcor_first, rcor_last, rcalc_first, rcalc_last, rcor_iter, rcalc_iter, -1 );

				// determine which LCB in the calculated ivs we're inside of
				calc_ivI = findLcb( calc_pair_adj, rcalc_iter->pos1 );
				bp_dist.push_back( (double)calc_pair_adj[calc_ivI].right_end[0] - (double)cor_pair_adj[cor_ivI].right_end[0] );
				bp_dist_sum += bp_dist.back();
				bp_dist_count++;
				if( genome::absolut( bp_dist.back() ) < genome::absolut( right_dists(seqI,seqJ)[calc_ivI] ) )
					right_dists(seqI,seqJ)[calc_ivI] = bp_dist.back();

				parity = (cor_pair_adj[cor_ivI].left_end[1] < 0) == (calc_pair_adj[calc_ivI].left_end[1] < 0);
				if( !parity )
				{
					bp_dist.push_back( (double)cor_pair_adj[cor_ivI].left_end[1] + (double)calc_pair_adj[calc_ivI].left_end[1] );
					if( calc_pair_adj[calc_ivI].left_end[1] > 0 )
						bp_dist.back() *= -1;
				}else
					bp_dist.push_back( (double)calc_pair_adj[calc_ivI].right_end[1] - (double)cor_pair_adj[cor_ivI].right_end[1] );
				bp_dist_sum += bp_dist.back();
				bp_dist_count++;

				if( genome::absolut( bp_dist.back() ) < genome::absolut( right_dists(seqJ,seqI)[calc_ivI] ) )
					right_dists(seqJ,seqI)[calc_ivI] = bp_dist.back();
			}
		}
	}

	for( seqI = 0; seqI < seq_table.size(); seqI++ ){
		for( seqJ = seqI+1; seqJ < seq_table.size(); seqJ++ ){
			// compute falsepos as the number of lcbs that weren't found in the correct aln
			bp_sp_falsepos += found_lcbs(seqI, seqJ).size() - found_lcbs(seqI, seqJ).count();
		}
	}

	double bp_sp_sensitivity = bp_sp_truepos / bp_sp_possible;
	double bp_sp_ppv = bp_sp_truepos / (bp_sp_truepos + bp_sp_falsepos);
	double bp_sp_mean_dist = bp_dist_sum / bp_dist_count;
	double bp_sp_stddev_dist = 0;
	for( size_t i = 0; i < bp_dist.size(); ++i )
		bp_sp_stddev_dist += (bp_dist[i] - bp_sp_mean_dist) * (bp_dist[i] - bp_sp_mean_dist);
	bp_sp_stddev_dist /= bp_dist.size()-1;
	bp_sp_stddev_dist = sqrt( bp_sp_stddev_dist );

	double bp_sp_mean_2 = 0;
	double bp_sp_count_2 = 0;
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		for( uint seqJ = 0; seqJ < seq_table.size(); seqJ++ )
		{
			for( size_t lcbI = 0; lcbI < left_dists(seqI,seqJ).size(); ++lcbI )
			{
				if( left_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
				{
					bp_sp_mean_2 += left_dists(seqI,seqJ)[lcbI];
					bp_sp_count_2++;
				}
				if( right_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
				{
					bp_sp_mean_2 += right_dists(seqI,seqJ)[lcbI];
					bp_sp_count_2++;
				}
			}
		}
	}
	bp_sp_mean_2 /= bp_sp_count_2;

	double bp_sp_stddev_2 = 0;
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		for( uint seqJ = 0; seqJ < seq_table.size(); seqJ++ )
		{
			for( size_t lcbI = 0; lcbI < left_dists(seqI,seqJ).size(); ++lcbI )
			{
				if( left_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
					bp_sp_stddev_2 += (left_dists(seqI,seqJ)[lcbI] - bp_sp_mean_2) * (left_dists(seqI,seqJ)[lcbI] - bp_sp_mean_2);
				if( right_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
					bp_sp_stddev_2 += (right_dists(seqI,seqJ)[lcbI] - bp_sp_mean_2) * (right_dists(seqI,seqJ)[lcbI] - bp_sp_mean_2);
			}
		}
	}
	bp_sp_stddev_2 /= (bp_sp_count_2 - 1);
	bp_sp_stddev_2 = sqrt(bp_sp_stddev_2);

	cout << "Sum of pairs LCB sensitivity: " << bp_sp_sensitivity << std::endl;
	cout << "Sum of pairs LCB positive predictive value: " << bp_sp_ppv << std::endl;
//	cout << "Mean distance between predicted breakpoint and true breakpoint: " << bp_sp_mean_dist << std::endl;
//	cout << "Standard deviation: " << bp_sp_stddev_dist << std::endl;

//	cout << "\nSecond bp accuracy method:\n";
	cout << "Mean distance between predicted breakpoint and true breakpoint: " << bp_sp_mean_2 << std::endl;
	cout << "Standard deviation: " << bp_sp_stddev_2 << std::endl;

	// calculate quartile statistics
	vector< gnSeqI > all_dists;
	for( uint seqI = 0; seqI < seq_table.size(); seqI++ )
	{
		for( uint seqJ = 0; seqJ < seq_table.size(); seqJ++ )
		{
			for( size_t lcbI = 0; lcbI < left_dists(seqI,seqJ).size(); ++lcbI )
			{
				if( left_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
					all_dists.push_back( genome::absolut( left_dists(seqI,seqJ)[lcbI] ));
				if( right_dists(seqI,seqJ)[lcbI] != (std::numeric_limits<double>::max)() )
					all_dists.push_back( genome::absolut( right_dists(seqI,seqJ)[lcbI] ) );
			}
		}
	}
	std::sort( all_dists.begin(), all_dists.end() );
	if( all_dists.size() > 0 )
	{
		cout << "Absolute distance, min: " << all_dists.front() << endl;
		cout << "Absolute distance, first quartile: " << all_dists[ (size_t)((float)all_dists.size() * 0.25) ] << endl;
		cout << "Absolute distance, second quartile: " << all_dists[ (size_t)((float)all_dists.size() * 0.5) ] << endl;
		cout << "Absolute distance, third quartile: " << all_dists[ (size_t)((float)all_dists.size() * 0.75) ] << endl;
		cout << "Absolute distance, max: " << all_dists.back() << endl;
	}
}










/**
 * program to score alignments
 * reads in a "correct" alignment and a calculated alignment
 * scores the calculated alignment based on the correct one
 */
int main( int argc, char* argv[] ){
	
	if( argc < 3 ){
		cout << "scoreAlignment <correct alignment> <calculated alignment> <evolved sequence file> [--disable-lcb-scoring]\n";
		cout << "Use --disable-lcb-scoring when scoring an alignment that may align a given nucleotide more than once.\n";
		return -1;
	}

	boolean score_lcbs = true;
	boolean debug_mismatches = false;	/**< turns on code to debug mismatches in evolved and aligned base pairs */
	string correct_fname = argv[ 1 ];
	string calculated_fname = argv[ 2 ];
	string evolved_fname;
	if( argc > 3 ){
		debug_mismatches = true;
		evolved_fname = argv[ 3 ];
	}

	ifstream correct_in;
	correct_in.open( correct_fname.c_str() );
	if( !correct_in.is_open() ){
		cerr << "Error opening " << correct_fname << endl;
		return -1;
	}
	ifstream calculated_in;
	calculated_in.open( calculated_fname.c_str() );
	if( !calculated_in.is_open() ){
		cerr << "Error opening " << calculated_fname << endl;
		return -1;
	}
	if( argc > 4 )
	{
		string lcb_score_arg = argv[4];
		if( lcb_score_arg == "--disable-lcb-scoring" )
			score_lcbs = false;
	}
//try{
	IntervalList correct_ivs;
	IntervalList calculated_ivs;
	correct_ivs.ReadStandardAlignment( correct_in );
	correct_in.close();
	calculated_ivs.ReadStandardAlignment( calculated_in );
	calculated_in.close();

	if( calculated_ivs.size() == 0 )
	{
		cerr << "WARNING!  No alignment found in input file.  Assuming gap predictions!\n";
		warn_missing = false;
	}

	gnSequence empty_seq;
	vector< gnSequence* > seq_table( correct_ivs[0].SeqCount(), &empty_seq );
	uint seq_count = seq_table.size();
	const gnFilter* comp_filter = gnFilter::DNAComplementFilter();
	
	gnSequence evolved_gnseqs;
	vector< string > evolved_seqs( seq_count );
	if( debug_mismatches ){
		evolved_gnseqs.LoadSource( evolved_fname );
		for( uint i = 0; i < seq_count; i++ ){
			evolved_seqs[ i ] = evolved_gnseqs.contig( i ).ToString();
		}
	}
	
	compareAlignments( correct_ivs, calculated_ivs, evolved_seqs, seq_table );

	if( score_lcbs )
		computeLCBaccuracy( correct_ivs, calculated_ivs, evolved_seqs, seq_table );
/*	
}catch( gnException& gne ){
	cerr << gne << endl;
}catch( exception& e ){
	cerr << e.what() << endl;
}catch( char const* c ){
	cerr << c << endl;
}catch(...){
	cerr << "Unhandled exception" << endl;
}
*/
}


