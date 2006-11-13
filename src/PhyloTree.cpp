#include "PhyloTree.h"
#include <sstream>
#include <stack>
using namespace std;

typedef unsigned uint;

PhyloTree::PhyloTree() : vector< TreeNode >() {
	weight = 0;
	root = 0;
}

void PhyloTree::clear()
{
	vector< TreeNode >::clear();
	weight = 0;
	root = 0;
}

/**
 *  read in a phylogenetic tree in the Newick file format.
 *  Currently this only works on rooted trees...
 */
void PhyloTree::readTree( istream& tree_file ){
	string line;
	clear();
	if( !getline( tree_file, line ) )
		return;

	stringstream line_str( line );

	// look for a weight
	string::size_type open_bracket_pos = line.find( "[" );
	string::size_type bracket_pos = line.find( "]" );
	if( open_bracket_pos != string::npos && bracket_pos != string::npos && 
		open_bracket_pos < bracket_pos && bracket_pos < line.find( "(" ) ){
		// read in a weight
		getline( line_str, line, '[' );
		getline( line_str, line, ']' );
		stringstream weight_str( line );
		weight_str >> weight;
	}

	stack< node_id_t > parent_stack;
	while( getline( line_str, line, '(' ) ){
		// if this is a lone open parens then simply create a new
		// parent node and push it on the parent stack
		if( line.size() == 0 ){
			TreeNode new_parent;
			if( parent_stack.size() > 0 ){
				new_parent.parents.push_back( parent_stack.top() );
				(*this)[ parent_stack.top() ].children.push_back( (*this).size() );
			}
			parent_stack.push( (*this).size() );
			push_back( new_parent );
			continue;
		}
		
		// this is not a lone parens,
		// process the data between the parenthesis
		stringstream data_str( line );
		while( getline( data_str, line, ')' ) ){
			// this is a group of leaf nodes
			stringstream leaf_str( line );
			// read in each leaf in this grouping
			// leaf nodes are in name:distance format and are separated
			// by commas
			while( getline( leaf_str, line, ',' ) ){
				stringstream leaf_data_str( line );
				getline( leaf_data_str, line, ':' );
				TreeNode new_leaf;
				new_leaf.name = line;
				leaf_data_str >> new_leaf.distance;
				// link the parent and child
				new_leaf.parents.push_back( parent_stack.top() );
				(*this)[ parent_stack.top() ].children.push_back( (*this).size() );
				push_back( new_leaf );
			}
			
			// now read the distance for the current parent and 
			// pop it off the stack
			getline( data_str, line, ',' );
			if( line.find( ":" ) == string::npos ){
				// default the distance to 0
				(*this)[ parent_stack.top() ].distance = 0;
			}else{
				stringstream parent_str( line );
				getline( parent_str, line, ':' );
				(*this)[ parent_stack.top() ].name = line;
				parent_str >> (*this)[ parent_stack.top() ].distance;
			}
			parent_stack.pop();
		}
	}
}

void PhyloTree::writeTree( ostream& os ) const{
	stack< node_id_t > node_stack;
	stack< uint > child_stack;
	node_stack.push( root );
	child_stack.push( 0 );

	if( (*this).weight != 0 )
		os << "[" << weight << "]";
	os << "(";

	while( node_stack.size() > 0 ) {
		if( (*this)[ node_stack.top() ].children.size() != 0 ){
			// this is a parent node
			// if we have scanned all its children then pop it
			if( child_stack.top() == (*this)[ node_stack.top() ].children.size() ){
				os << ")";
				if( node_stack.size() > 1 )
					os << ":" << (*this)[ node_stack.top() ].distance;
				node_stack.pop();
				child_stack.pop();
				continue;
			}
			// try to recurse to its children
			// if the child is a parent as well spit out a paren
			node_id_t child = (*this)[ node_stack.top() ].children[ child_stack.top() ];
			node_stack.push( child );
			child_stack.top()++;
			// print a comma to separate multiple children
			if( child_stack.top() > 1 )
				os << ",";
			if( (*this)[ child ].children.size() > 0 ){
				child_stack.push( 0 );
				os << "(";
			}
			continue;
		}
		
		// this is a leaf node
		os << (*this)[ node_stack.top() ].name << ":" << (*this)[ node_stack.top() ].distance;
		
		// pop the child
		node_stack.pop();
	}
	os << ";" << endl;
}


double PhyloTree::getHeight() const
{
	return getHeight( root );
}
double PhyloTree::getHeight( node_id_t nodeI ) const
{
	if( (*this)[ nodeI ].children.size() == 0 )
		return (*this)[ nodeI ].distance;
	return (*this)[ nodeI ].distance + getHeight( (*this)[ nodeI ].children[ 0 ] );
}
