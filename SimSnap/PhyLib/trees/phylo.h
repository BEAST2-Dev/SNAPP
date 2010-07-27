/**
 * @file phylo.h
 * @brief Header file for the generic phylogeny class phylo
 * @author David Bryant
 * @version 0.1
 */

#ifndef PHYLO_H
#define PHYLO_H

#include <functional>
#include <algorithm>
#include <iterator>
#include <climits>
//#include "memman.h"
#include "../utilities/phylibException.h"

#ifdef __ITERATOR_NEEDED
namespace std {

	// According to ANSI C++ this should be defined but it is missing
	// from current GNU C++. So we add it here, if __ITERATOR_NEEDED is defined.
	// (MSVC++ 6.0 defines it, but incorrectly.)

	template <class _Category, class _Tp, class _Distance = ptrdiff_t,
	class _Pointer = _Tp*, class _Reference = _Tp&>
	struct iterator {
		typedef _Category iterator_category;
		typedef _Tp value_type;
		typedef _Distance difference_type;
		typedef _Pointer pointer;
		typedef _Reference reference;
	};
};
#endif

namespace Phylib {

	using namespace std;

	/**
	 * @class phylo 
	 *
	 * An STL inspired class for phylogenies. 
	 * 
	 * The phylogeny is seen as a container to hold information at the nodes, and the type of information stored is specified using 
	 * the template. For example, a phylo<int> would store a single integer for every node. This makes it easy to write tree based
	 * algorithms: decide on what information is to be stored at a node and include this in the class. For example, an algorithm
	 * for Fitch parsimony might store:
	 * class {
	 * 		BitSet states;
	 * 		int taxon id;
	 * }
	 * 
	 * at the nodes.
	 * 
	 * 
	 * Unlike lists and vectors and so on, the structure of the container (i.e. the underlying tree) is important. We therefore add convertors
	 * to go from two different types of phylo: to copy from phylo<type1> to phylo<type2>. If you are implementing several algorithms in 
	 * one program you can  either define a different phylo type for each algorithm, or put the data from all the algorithms into one class.
	 * The advantage of the former is simplicty and memory reduction. The advantage of the latter is that you can reduce the number of times
	 * phylogenies are copied.
	 * 
	 * We note that the ordering of children is maintained, even though it (technically) doesn't affect the meaning of the phylogeny. Thus two
	 * phylogenies are equal if and only if the trees are the same, their contents are the same, and the ordering of the children of 
	 * each node is the same. See the tree utilities section for routines to compare trees that do not respect child ordering (a much more
	 * computationally difficult issue).
	 * 
	 * The storage of the tree is based on a data structure designed (I think) by Dave Swofford, and used to great effect in PAUP. There is an
	 * 'immortal' header node that can never be removed. It usually has exactly one child, the root of the tree.[CHECK]
	 * A node in the tree has three links:
	 * 	one to its parent (or to the header if this is the root)
	 *   one to its leftmost child (or to null if this node has no children)
	 *   one to its next sibling (or to null if this is the rightmost child)
	 * I have programmed using both this data structure and one where children are stored as an array, and find that Dave's structure
	 * is easier to work with, after a bit of a steep learning curve. In future, I'll implement a binary tree only version that just
	 * stores the left and right child.
	 * 
	 * The basic phylo class provides two iterators:
	 *  iterator: a pre-order traversal of the tree (parents before children)
	 *  post_iterator: a post-order traversal of the tree
	 * For both, the time to go from one node to the next is linear in the length of the path connecting the nodes, which is (amortised) constant
	 * time over the whole traversal. If you need faster than this, you may want to precompute the traversal and store the sequence of nodes.
	 * There are also new versions of begin()  and end() allowing you to easily traverse a given subtree. 
	 * 
	 * The manipulation of trees differs substantially from that of lists (and vectors etc.) In this class, there are only 
	 *
	 * I based the implementation of the phylo on the SGI list code, as adapted by D.Musser for a tutorial on lists. As such, I include the
	 * HP copyright header. Once I figure out how, I'll place this code under the GNU lesser license.
	 * 
	 *
	 * Copyright (c) 1994
	 * Hewlett-Packard Company
	 *
	 * Permission to use, copy, modify, distribute and sell this software
	 * and its documentation for any purpose is hereby granted without fee,
	 * provided that the above copyright notice appear in all copies and
	 * that both that copyright notice and this permission notice appear
	 * in supporting documentation.  Hewlett-Packard Company makes no
	 * representations about the suitability of this software for any
	 * purpose.  It is provided "as is" without express or implied warranty.
	 *
	 * Adapted by D.Musser for a tutorial on lists.
	 *
	 * @version 0.1 27/03/08 Initial implementation by DJB
	 */
	template<typename T> class phylo {
protected:

		class phylo_node {
	public:
			phylo_node* par;
			phylo_node* left;
			phylo_node* right;
			T data;
			phylo_node() :
				par(0), left(0), right(0), data() {
			}
		};

public:
		//Types:
		typedef T value_type;
		typedef T& reference;
		typedef const T& const_reference;
		typedef size_t size_type;
		typedef ptrdiff_t difference_type;

private:
		typedef T* pointer;
		typedef phylo_node* link_type;
		typedef const T* const_pointer;
		//1typedef memory_manager<phylo_node> memory_manager_type;

protected:
		link_type _header;
		link_type first_leaf;
		//1static int number_of_phylos;
		//1static memory_manager< phylo_node> manager;

public:

		class const_iterator; // incomplete definition, for use in friend decl.

		/**
		 *@class iterator
		 * Pointer type for phylo<T> 
		 */
		class iterator : public std::iterator<std::forward_iterator_tag, value_type, void, reference> {
			friend class phylo<T>;
			
			friend class const_iterator;
	protected:
			link_type current;
	protected:
			iterator(link_type x) :
				current(x) {
			}
	public:

			/**
			 * Default constructor. Iterator is NULL.
			 * */
			iterator() {
				current = 0;
			}

			/**
			 * Test equality.
			 * @param x iterator to compare
			 * @return true if this iterator is pointing to the same node as the given iterator.
			 */
			bool operator==(const iterator& x) const {
				return current == x.current;
			}

			/**
			 * Test inequality
			 * @param x Iterator to compare to.
			 * @return true if this iterator not pointing to the same node as the given iterator.
			 * 
			 */
			bool operator!=(const iterator& x) const {
				return current != x.current;
			}

			/**
			 * Arrow operator. 
			 * @return reference to data pointed to by this iterator. 
			 */
			inline pointer operator->() {

				return &(current->data);
			}

			/**
			 * Access  the data at the node the iterator is pointing to
			 * @return reference to data
			 */
			typename iterator::reference operator*() const {
				return current->data;
			}

			/**
			 * Test if null.
			 * @return true if this iterator is pointing to null.
			 */
			bool null() {
				return current == 0;
			}

			/**
			 * Make iterator null.
			 */
			void set_null() {
				current = 0;
			}

			/**
			 * test if leaf
			 * @return true if iterator points to a leaf.
			 */
			bool leaf() const {
				return left().null();
			}

			/**
			 * test if header
			 * @return true if iterator points to the header
			 */
			bool header() const {
				return par().null();
			}

			
			
			/**
			 * test if root
			 * @return true if this is a child of the header node.
			 */
			bool root() const {
				return !header() && (par().header());
			}
			
			
			
			
			/**
			 * Leftmost child.
			 * @return iterator pointing to leftmost child, or null iterator if this is a leaf.
			 */
			iterator left() const {
				return iterator(current->left);
			}

			/**
			 * Leftmost child.
			 * @return iterator pointing to leftmost child, or null iterator if this is a leaf.
			 */
			inline iterator lchild() const {
				return left();
			}

			/**
			 * Sibling to right
			 *@return iterator pointing to sibling to the right of this one, or null if this is the rightmost child of its parent.
			 */
			iterator right() const {
				return iterator(current->right);
			}

			/**
			 * Sibling to right
			 *@return iterator pointing to sibling to the right of this one, or null if this is the rightmost child of its parent.
			 */
			inline iterator rsib() const {
				return right();
			}

			/**
			 * Parent node
			 * @return iterator pointing to parent node, or null iterator if this is the header
			 */
			iterator par() const {
				return iterator(current->par);
			}

			iterator& operator++() { // Pre-increment
				iterator tmp = this->next_post();
				current = tmp.current;
				return *this;
			}

			iterator operator++(int) { // Post-increment
				iterator tmp = *this;
				operator++();
				return tmp;
			}

			iterator next_pre() const {
				if (this->header()) {
					throw PhylibException("Trying to apply next_pre to the header node");
				}

				if (!this->leaf()) {
					return this->left();
				}

				iterator p = *this;
				while (!p.header() && p.right().null()) {
					p=p.par();
				}
				if (p.header()) {
					return iterator(0);
				} else {
					return p.right();
				}
			}

			iterator next_post() {

				if (this->header()) {
					throw PhylibException("Trying to apply next_post to the header node");
				}

				if (this->right().null()) {
					if (this->par().header())
						return iterator(0);
					else
						return this->par();
				}
				iterator p = this->right();
				while (!p.leaf()) {
					p = p.left();
				}
				return p;
			}

			int numChildren()  {
				int n=0;
				for(iterator p = left();!p.null();p=p.right())
					n++;
				return n;
			}
			
		};

		/**
		 *@class const_iterator
		 * Constant iterator for phylo<T> implementing pre-order traversals.
		 * the time to go from one node to the next is linear in the length of the path connecting the nodes, which is (amortised) constant
		 * time over the whole traversal. If you need faster than this, you may want to precompute the traversal and store the sequence of nodes.
		 */
		class const_iterator : public std::iterator<forward_iterator_tag,value_type, pointer, const_reference> {
			friend class phylo<T>;
	public:
			link_type current;
	protected:
			const_iterator(link_type x) :
				current(x) {
			}

	public:
			/**
			 * Default constructor.
			 * */
			const_iterator() {
			}

			/**
			 * Convert iterator to const_iterator
			 * @param x iterator
			 */
			const_iterator(const iterator& x) :
				current(x.current) {
			}

			/**
			 * Test equality.
			 * @param x const_iterator to compare
			 * @return true if this const_iterator is pointing to the same node as the given const_iterator.
			 */
			bool operator==(const const_iterator& x) const {
				return current == x.current;
			}

			/**
			 * Test inequality
			 * @param x Iterator to compare to.
			 * @return true if this const_iterator not pointing to the same node as the given const_iterator.
			 * 
			 */
			bool operator!=(const const_iterator& x) const {
				return current != x.current;
			}

			/**
			 * Arrow operator. 
			 * @return reference to data pointed to by this iterator. 
			 */
			inline const_pointer operator->() {
				
				return &(current->data);
			}

			
			/**
			 * Access  the data at the node the const_iterator is pointing to
			 * @return reference to data
			 */
			typename const_iterator::reference operator*() const {
				return current->data;
			}

			/**
			 * Test if null.
			 * @return true if this const_iterator is pointing to null.
			 */
			bool null() const {
				return current == 0;
			}

			/**
			 * Make const_iterator null.
			 */
			void set_null() {
				current = 0;
			}

			/**
			 * test if leaf
			 * @return true if iterator points to a leaf.
			 */
			bool leaf() const {
				return left().null();
			}

			/**
			 * test if header
			 * @return true if iterator points to the header
			 */
			bool header() const {
				return par().null();
			}

			/**
			 * test if root
			 * @return true if this is a child of the header node.
			 */
			bool root() const {
				return !header() && (par().header());
			}
			/**
			 * Leftmost child.
			 * @return const_iterator pointing to leftmost child, or null const_iterator if this is a leaf.
			 */
			const_iterator left() const {
				return const_iterator(current->left);
			}

			/**
			 * Leftmost child.
			 * @return const_iterator pointing to leftmost child, or null const_iterator if this is a leaf.
			 */
			inline const_iterator lchild() const {
				return left();
			}

			/**
			 * Sibling to right
			 *@return const_iterator pointing to sibling to the right of this one, or null if this is the rightmost child of its parent.
			 */
			const_iterator right() const {
				return const_iterator(current->right);
			}

			/**
			 * Sibling to right
			 *@return const_iterator pointing to sibling to the right of this one, or null if this is the rightmost child of its parent.
			 */
			inline const_iterator rsib() const {
				return right();
			}

			/**
			 * Parent node
			 * @return const_iterator pointing to parent node, or null const_iterator if this is the header
			 */
			const_iterator par() const {
				return const_iterator(current->par);
			}

			const_iterator& operator++() { // Pre-increment
				const_iterator tmp = this->next_post();
				current = tmp.current;
				return *this;
			}

			const_iterator operator++(int) { // Post-increment
				const_iterator tmp = *this;
				operator++();
				return tmp;
			}

			const_iterator next_pre() const {
				if (this->header()) {
					throw PhylibException("Trying to apply next_pre to the header node");
				}

				if (!this->leaf()) {
					return this->left();
				}

				const_iterator p = *this;
				while (!p.header() && p.right().null()) {
					p=p.par();
				}
				if (p.header()) {
					return const_iterator(0);
				} else {
					return p.right();
				}
			}

			const_iterator next_post() {
				if (this->header()) {
					throw PhylibException("Trying to apply next_post to the header node");
				}

				if (this->right().null()) {
					if (this->par().header())
						return const_iterator(0);
					else
						return this->par();
				}
				const_iterator p = this->right();
				while (!p.leaf()) {
					p = p.left();
				}
				return p;
			}
			
			int numChildren() const {
				int n=0;
				for(const_iterator p = left();!p.null();p=p.right())
					n++;
				return n;
			}
			

		};

		//Constructors

		/**
		 * Construct an empty phylogeny
		 */
		phylo() {
			//1++number_of_phylos;
			//1_header = manager.get_node();
			_header = new phylo_node;
			_header->left = _header->right = _header->par = 0;
			//1new (&(_header->data)) T();

			//cout<<"initialise header = "<<_header<<endl;
		}

		/**
		 * Copy constructor
		 * Create a nw phylogeny that is a copy of a given phylogeny.
		 * @param x phylo to copy
		 */
		phylo(const phylo<T>& x) {
			//1++number_of_phylos;
			//1_header = manager.get_node();
			//1new (&(_header->data)) T();
			_header = new phylo_node; 
			
			_header->left = _header->right = _header->par = 0;
			operator=(x);
			//cout<<"initialise header = "<<_header<<endl;
		}

		/** 
		 * Destructor
		 * If phylo is non-empty, data data fields are destructed
		 * */
		~phylo() {
			erase(_header);
			//cout<<"removing "<<_header<<endl;
			//1(_header->data).~T();
			//1manager.put_node(_header);
			//1if (--number_of_phylos == 0)
			//1	manager.deallocate_buffers();
			delete _header;
		}
		/**
		 * assignment.
		 * makes this a copy of x
		 * @param x phylo to copy
		 * @return copy of x
		 */
		phylo<T>& operator=(const phylo<T>& x) {
			if (this!=&x) {
				erase(_header);
				copy_children(x._header, _header);
			}
			return *this;
		}

		/**
		 * swap.
		 * swaps data between this phylo and x. Note that iterators currently pointing to this this phylo
		 * will now be pointing to x and vice versa. Operation takes constant time.
		 * @param x phylo to swap with
		 */
		void swap(phylo<T>& x) {
			link_type tmp = x._header;
			x._header = _header;
			_header = tmp;
		}

		/**
		 * Returns the root of the tree. This is the first node of pre-order traversal. Returns
		 * null iterator if the tree is empty. Constant time.
		 * @return iterator pointing to the root (the first root if there are multiple)
		 */
		iterator root() const {return iterator(_header->left);}

		/**
		 * returns the leftmost leaf in the tree. This is the first node of a post-order traversal.
		 *  Returns null iterator if tree is empty. Time linear in the length of the path from the
		 * root to this leaf.
		 * @return iterator pointing to leftmost leaf in the tree.
		 */
		iterator leftmost_leaf() {
			iterator p = iterator(_header->left);
			while(!p.null() && !p.leaf()) {
				p = p.left();
			}
			return p;
		}

		/**
		 * Returns iterator pointing to the header node. 
		 */
		iterator header() const {
			return iterator(_header);
		}

		/**
		 * test if empty.
		 * @return true if the phylo has no nodes (apart from the header)
		 */
		bool empty() const {return header().left().null();}

		/**
		 * compute number of nodes.
		 * this takes linear time in the size of the tree.
		 * @return number of nodes in the phylo
		 */
		size_type size() const {

			size_type n=0;
			for(const_iterator p = root();!p.null();p++)
			n++;
			return n;
		}

		/**
		 * maximum possible size given memory limits
		 * @return maximum size for a phylo.
		 */
		size_type max_size() const {
			return max(size_type(1), size_type(UINT_MAX/sizeof(T)));
		}

		/**
		 * insert new child.
		 * Creates new leftmost child of node. Assigns data with empty constructor
		 * @param position parent of new node
		 * @throw PhylibException if position is null
		 * @return iterator to newly inserted node
		 */
		iterator insert_child(iterator position) {
			if (position.null())
			throw PhylibException("Invalid insertion: position iterator is null");
			//1link_type tmp = manager.get_node();
			//1new (&(tmp->data)) T();
			link_type tmp = new phylo_node;
			tmp->par = position.current;
			tmp->right = position.current->left;
			tmp->left = 0;
			position.current->left = tmp;
			return iterator(tmp);
		}

		/**
		 * insert new child.
		 * Creates new leftmost child of node and assigns x to it.
		 * @param position parent of new node
		 * @param x data to assign to new node.
		 * @throw PhylibException if position is null
		 * @return iterator to newly inserted node
		 */
		iterator insert_child(iterator position, const T& x) {
			if (position.null())
			throw Phylib::PhylibException("Invalid insertion: position iterator is null");
			//link_type tmp = manager.get_node();
			//new (&(tmp->data)) T(x);
			link_type tmp = new phylo_node;
			tmp->data = x;
			
			
			tmp->par = position.current;
			tmp->right = position.current->left;
			tmp->left = 0;
			position.current->left = tmp;
			return iterator(tmp); //Could also return position.left()
		}

		/**
		 * insert new sibling.
		 * creates a new sibling of the node position and assigns data with null constructor.
		 * throws exception if position in the header node
		 * @param position pointing to position in the tree
		 * @throw Phylib::PhylibException if position equals header or is null
		 * @return iterator to newly inserted node
		 */
		iterator insert_sibling(iterator position) {
			if (position.header())
			throw Phylib::PhylibException("Invalid insertion: attempting to create sibling of header node in phylo.");
			if (position.null())
			throw Phylib::PhylibException("Invalid insertion: position iterator is null");
			//1link_type tmp = manager.get_node();
			//1new (&(tmp->data)) T();
			link_type tmp = new phylo_node;
			
			tmp->par = position.current->par;
			tmp->right = position.current->right;
			tmp->left = 0;
			(position.current)->right = tmp;
			return iterator(tmp);
		}

		/**
		 * insert new sibling.
		 * creates a new sibling of the node position and assigns the data x to it.
		 * throws exception if position in the header node
		 * @param position pointing to position in the tree
		 * @param x data to assign to the new node
		 * @throw Phylib::PhylibException if position equals header or is null
		 * @return iterator to newly inserted node
		 */
		iterator insert_sibling(iterator position, const T&x) {

			if (position.header())
			throw Phylib::PhylibException("Invalid insertion: attempting to create sibling of header node in phylo.");
			if (position.null())
			throw Phylib::PhylibException("Invalid insertion: position iterator is null");
			//1link_type tmp = manager.get_node();
			//1new (&(tmp->data)) T(x);
			link_type tmp = new phylo_node;
			tmp->data = x;
			tmp->par = position.current->par;
			tmp->right = position.current->right;
			tmp->left = 0;
			(position.current)->right = tmp;
			return iterator(tmp);
		}

		/**
		 * insert subtree as child.
		 * Creates new leftmost child of node copies the subtree rooted at x to it.
		 * @param position parent of new node
		 * @param x root of subtree to be copied from
		 * @throw Phylib::PhylibException if node or x are null
		 * @return iterator pointing to copied subtree
		 */
		iterator insert_child(iterator position, const_iterator x) {
			if (x.null())
			throw Phylib::PhylibException("Invalid insertion: subtree iterator is null");
			iterator tmp = insert_child(position,x.current->data);
			copy_children(x,tmp);
			return tmp;
		}

		/**
		 * insert subtree as new sibling.
		 * creates a new sibling of the node position and assigns the data x to it.
		 * throws exception if position in the header node
		 * @param position iterator pointing to position in the tree
		 * @param x root of the subtree to be copied
		 * @throw Phylib::PhylibException if position equals header or if either iterator is null
		 */
		iterator insert_sibling(iterator& position, const_iterator x) {
			if (x.null())
			throw Phylib::PhylibException("Invalid insertion: subtree iterator is null");
			iterator tmp = insert_sibling(position,x.current->data);
			copy_children(x,tmp);
			return position.right();
		}

		/**
		 Contracts the edge above this node. 
		 @param position node below edge to be contracted.
		 @returns iterator pointing to parent of the node.
		 **/
		iterator contract(iterator position) {
			//first copy children.
			while (!position.left().null()) 
				this->graft_sibling(position,*this,position.left());
			iterator result = position.par();
			erase(position);
			return result;
		}
		 
		
		/**
		 * Clear all nodes from the tree.
		 */
		void clear() {
			erase(header());
		}

		/**
		 * erase subtree.
		 * deletes all nodes below position. If position is not the header, then this node is deleted as well. 
		 * any data stored at this nodes will be destructed,
		 * Any iterators pointing to nodes in the subtree will become invalidated.
		 * @param position root of subtree to be deleted.
		 * @throw Phylib::PhylibException if position is null
		 */
		void erase(iterator position) {

			if (position.null()) {
				throw Phylib::PhylibException("Illegal erase: position is null");
			}
			while(!position.left().null()) {
				erase(position.left()); //Erase the children
			}
			if (!position.header()) {

				if (position == position.par().left()) {
					(position.par().current)->left = (position.current)->right; //Deleting leftmost child of parent
				} else {
					iterator leftsib = position.par().left(); //Find sibling to the left of position.
					while(leftsib.right()!=position) {
						leftsib = leftsib.right();
					}
					(leftsib.current)->right = (position.current)->right;
				}
				//1(position.current)->data.~T();
				//1manager.put_node(position.current);
				delete position.current;
			}
		}

		/**
		 * Inserts the subtree rooted at node from phylogeny x as the first child of position. If x is a forest, these become new children of position.
		 * Returns an iterator pointing to the first tree inserted (i.e. the leftmost child of position). Operation takes
		 * constant time as elements are transferred rather than copied.
		 * @param position iterator in destination tree
		 * @param x source tree. 
		 * @param node root of subtree that is removed from x and inserted below position
		 * @return iterator pointer to first tree inserted.
		 */
		iterator graft_child(iterator position, phylo<T>& x, iterator node) {
			if (node.null())
			throw PhylibException("Error: trying to graft null source node");

			detach(node);
			(node.current)->right = (position.current)->left;
			(node.current)->par = (position.current);
			(position.current)->left = (node.current);
			return position.left();
		}

		/**
		 * Inserts the subtree rooted at node from phylogeny x as the right sibling of position. 
		 * Returns an iterator pointing to the first tree inserted (i.e. the new right sibling of position). Operation takes
		 * constant time as elements are transferred rather than copied.
		 * @param position iterator in destination tree
		 * @param x source tree. Left empty.
		 * @return iterator pointer to first tree inserted.
		 */

		iterator graft_sibling(iterator position, phylo<T>& x, iterator node) {
			if (node.null())
				throw PhylibException("Error: trying to graft null source node");
			if (position.par().null())
				throw PhylibException("Error: trying to graft node as sibling of header node");

			detach(node);
			(node.current)->right = (position.current)->right;
			(node.current)->par = (position.current)->par;
			(position.current)->right = (node.current);
			return position.right();
		}

		/**
		 * Inserts the phylogeny x as the first child of position. If x is a forest, these become new children of position.
		 * Returns an iterator pointing to the first tree inserted (i.e. the leftmost child of position). Operation takes\
		 * constant time as elements are transferred rather than copied.
		 * @param position iterator in destination tree
		 * @param x source tree. Left empty.
		 * @return iterator pointer to first tree inserted.
		 */

		iterator graft_child(iterator position, phylo<T>& x) {

			typename phylo<T>::iterator from = x.header().left().right();
			typename phylo<T>::iterator to_sib = graft_child(position,x,x.header().left());

			while(!from.null()) {
				iterator from_next = from.right();
				to_sib = graft_sibling(to_sib,x,from);
				from = from_next;
			}
			return position.left();
		}

		/**
		 * Inserts phylogeny x as the right sibling of position. If x is a forest, then they become consecutive siblings.
		 * Returns an iterator pointing to the first tree inserted (i.e. the new right sibling of position). Operation takes
		 * constant time as elements are transferred rather than copied.
		 * @param position iterator in destination tree
		 * @param x source tree. Left empty.
		 * @return iterator pointer to first tree inserted.
		 */

		iterator graft_sibling(iterator position, phylo<T>& x) { //removes the subtree at node from x and attaches it to this
			
			if (x.empty())
				throw PhylibException("Trying to graft empty tree");
			if (position.par().null())
				throw PhylibException("Error: trying to graft tree as sibling of header node");

			
			typename phylo<T>::iterator from = x.header().left();
			typename phylo<T>::iterator to_sib = position;

			while(!from.null()) {

				iterator from_next = from.right();
				to_sib = graft_sibling(to_sib,x,from);
				from = from_next;

			}
			return position.right();
		}

	protected:

		/**
		 * copies children of from to children of to.
		 * @param from iterator pointing to node above children that are copied from
		 * @param to iterator pointing to node to which the children are copied (right children of current siblings).
		 */
		void copy_children(const_iterator from, iterator to) {
			iterator to_child = to.left();
			for (const_iterator from_child = from.left();!from_child.null();from_child = from_child.right()) {
				//Create a new node, with the data from from_child
				if (to_child.null()) {
					to_child = insert_child(to,(*from_child)); //Insert a single node that has a copy of the data at from_child
				}
				else {
					to_child = insert_sibling(to_child,(*from_child));
				}

				copy_children(from_child,to_child);
			}
		}

		/**
		 * detach a node.
		 * fix up links pointing to a node so that it is effectively removed from the tree.
		 * @param  position node in tree.
		 */
		void detach(iterator position) {
			if (position == position.par().left()) {
				(position.par().current)->left = (position.current)->right; //Deleting leftmost child of parent
			} else {
				iterator leftsib = position.par().left(); //Find sibling to the left of position.
				while(leftsib.right()!=position) {
					leftsib = leftsib.right();
				}
				(leftsib.current)->right = (position.current)->right;
			}
			(position.current)->par = 0;

		}

		bool test_if_equal(iterator x, iterator y) {
			if (!x.header() && (*x != *y))
			return false;
			iterator q = y.left();
			for(iterator p=x.left();!p.null();p=p.right()) {
				if (q.null() || !test_if_equal(p,q)) {
					return false;
				}
				q = q.right();
			}
			if (!q.null()) {
				return false;
			}
			return true;
		}

	};

	template<typename T> inline bool operator==(const phylo<T>& x, const phylo<T>& y) {
		return test_if_equal(x._header,y._header);
	}

	//template<typename T> inline bool operator<(const list<T>&x, const list<T>& y);


	//template <typename T>
	//1int phylo<T>::number_of_phylos = 0;

	//1template <typename T>
	//1memory_manager< typename phylo<T>::phylo_node> phylo<T>::manager;

	/**
	 * Copies one phylogeny to another. The destination phylogeny is first cleared, and the copy
	 * is non-recursive (so suitable for large trees). This method assumes that type T2 can
	 * be instantiated from type T1; a compilation error will result if this is not the case.
	 * To copy just the topology, and not the content at each node, use copy_topology.
	 * @param from phylogeny to copy from (source)
	 * @param to phylogeny to copy to (destination). Any initial contents of this phylogeny are cleared.
	 */
	template <typename T1, typename T2> void copy(const phylo<T1>& from, phylo<T2>& to) {
		to.clear();
		if (from.empty())
		return;
		//Non-recursive copy algorithm. Copies tree while followin pre-order traversal of from tree.
		typename phylo<T1>::const_iterator v_from = from.header().left();
		typename phylo<T2>::iterator v_to = to.insert_child(to.header(),*v_from);
		while(!v_from.header()) {
			if (!v_from.leaf()) {
				v_from = v_from.left();
				v_to = to.insert_child(v_to,(*v_from));
				
				//if (v_from.leaf())
					//cerr<<"id from = "<<v_from->id<<" \t id to = "<<v_to->id<<endl;
				
			} else {
				while(!v_from.header()&&v_from.right().null()) {
					v_from=v_from.par();
					v_to = v_to.par();
				}
				if (!v_from.header()) {
					v_from = v_from.right();
					v_to = to.insert_sibling(v_to,(*v_from));
					
				//	if (v_from.leaf())
						//cerr<<"id from = "<<v_from->id<<" \t id to = "<<v_to->id<<endl;

					
				}
			}
		}
	}

	/**
	 * Copies the topology of one phylogeny to another. The destination phylogeny is first cleared, and the copy
	 * is non-recursive (so suitable for large trees). This method does not copy the data stored
	 * at each node, for this, use the copy function. As such, it can be used for phylogenies of any type.
	 * @param from phylogeny to copy from (source)
	 * @param to phylogeny to copy to (destination). Any initial contents of this phylogeny are cleared.
	 */
	template <typename T1, typename T2> void copy_topology(const phylo<T1>& from, phylo<T2>& to) {
		to.clear();
		if (from.empty())
		return;
		//Non-recursive copy algorithm. Copies tree while followin pre-order traversal of from tree.
		typename phylo<T1>::const_iterator v_from = from.header().left();
		typename phylo<T2>::iterator v_to = to.insert_child(to.header());
		while(!v_from.header()) {
			if (!v_from.leaf()) {
				v_from = v_from.left();
				v_to = to.insert_child(v_to);
			} else {
				while(!v_from.header()&&v_from.right().null()) {
					v_from=v_from.par();
					v_to = v_to.par();
				}
				if (!v_from.header()) {
					v_from = v_from.right();
					v_to = to.insert_sibling(v_to);
				}
			}
		}
	}
	
	
	
	
	
	
	
	
	
}

#endif


