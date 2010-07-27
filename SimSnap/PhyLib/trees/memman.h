/*
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
 */

/*
 * Adapted for tutorial purposes from the 10/31/95 HP STL release:
 *   Removed all dependence on allocators
 *   Moved nested iterator classes to a separate file (listiter.h)
 *   Moved special list memory management code to a separate class
 *     defined in file memman.h (this file)
 *   Rearranged the order of some of the member function definitions
 *
 * D. Musser, 3/14/96, revised 9/15/2000
 */

//#define NOMANAGEMENT   ...bug if you uncomment this... data field can be constructed twice.

#include <new>
using namespace std;

template <class T> inline T* block_allocate(ptrdiff_t size, T*) {
	return (T*)(::operator new((size_t)(size * sizeof(T))));
}

template <typename phylo_node> class memory_manager {

	typedef phylo_node* link_type;

	struct node_buffer {
		node_buffer* next_buffer;
		phylo_node* buffer;
	};

	typedef node_buffer* buffer_pointer;

	buffer_pointer buffer_list;
	link_type free_list;
	link_type next_avail;
	link_type last;
	ptrdiff_t buffer_size() {
		return 1024;
	}

	void add_new_buffer() {
#ifndef NOMANAGEMENT
		buffer_pointer tmp = new node_buffer;
		tmp->buffer = block_allocate(buffer_size(), (phylo_node*)0);
		
	//	tmp->buffer = new phylo_node[buffer_size()];
		
		
		tmp->next_buffer = buffer_list;
		buffer_list = tmp;
		next_avail = buffer_list->buffer;
		last = next_avail + buffer_size();
#endif
	}

public:
	memory_manager() {
		//cout<<"Constructing memory manager"<<endl;
	}
	void deallocate_buffers() {
#ifndef NOMANAGMENT
		while (buffer_list) {
			buffer_pointer tmp = buffer_list;
			buffer_list = (buffer_pointer)(buffer_list->next_buffer);
			
			
			delete  (tmp->buffer);
			delete tmp;
		}
			
		
		free_list = 0;
		next_avail = 0;
		last = 0;
#endif
	}

	link_type get_node() {
#ifndef NOMANAGEMENT
		link_type tmp = free_list;
		if (free_list != NULL) {
			free_list = (link_type)(free_list->left);
			//cout<<"returned node "<<tmp<<endl;
			return tmp;
		}
		if (next_avail == last)
			add_new_buffer();
		//cout<<"returned node2 "<<next_avail<<endl;
		return next_avail++;

		
#else
		return new phylo_node;
#endif
		
	}

	void put_node(link_type p) {
#ifndef NOMANAGEMENT

		//cout<<"put node "<<p<<endl;
		p->left = free_list;
		free_list = p;
#else
		delete p;
#endif
	}
};


