/*
 *  sorting.h
 *
 *  Created by David Bryant on 9/06/09.
 *  Copyright 2009 University of Auckland. All rights reserved.
 *
 */

#ifndef SORTING_H_INCLUDE
#define SORTING_H_INCLUDE  

#include"../global/stdIncludes.h"

namespace Phylib {
	
	template<typename A, typename B> bool compareFirst(const pair<A,B>& x, const pair<A,B>& y) {
		return (x.first < y.first);
	}
	
	template<typename A, typename B> bool compareSecond(const pair<A,B>& x, const pair<A,B>& y) {
		return (x.second < y.second);
	}
	
	
	/**
	 Sort the elements of data and keys according to the values in keys. Both vectors must be the
	 same length.
	 **/
	template<typename DATA, typename KEY> void sort_vectors(vector<DATA>& data, vector<KEY>& keys) {
		assert(data.size()==keys.size());
		int n = data.size();
		vector< pair< DATA,  KEY> > both(data.size());
		for(int i=0;i<n;i++) 
			both[i] = pair<DATA,KEY>(data[i],keys[i]);
		std::sort(both.begin(),both.end(),compareSecond<DATA,KEY>);
		for(int i=0;i<n;i++) {
			data[i] = both[i].first;
			keys[i] = both[i].second;
		}
	}
	
		
}


#endif


