#ifndef __PhyloTree_h__
#define __PhyloTree_h__

#include <vector>
#include <string>
#include <iostream>

typedef int node_id_t;
struct TreeNode {
	std::string name;
	double distance;
	std::vector< node_id_t > parents;	/**< if parents.size() == 0 this is a root node */
	std::vector< node_id_t > children;	/**< if children.size() == 0 this is a leaf node */
};

class PhyloTree : public std::vector< TreeNode > {
public:
	PhyloTree();
	double weight;
	node_id_t root;
	void clear();
	/**
	 * Reads a tree in Newick format.  WARNING:  only reads rooted trees correctly
	 */
	void readTree( std::istream& tree_file );
	/**
	 * Writes a tree in Newick format
	 */
	void writeTree( std::ostream& os ) const;
	/**
	 * Determines the height of the tree along the path from the root to the left-most leaf node
	 */
	double getHeight() const;
	/**
	 * Determines the height of the tree along the path from nodeI to its left-most descendant leaf node
	 */
	double getHeight( node_id_t nodeI ) const;
protected:
};



#endif // __PhyloTree_h__
