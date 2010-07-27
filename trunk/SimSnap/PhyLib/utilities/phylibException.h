/*
 *  phyloException.h
 *  PhyloLib
 *
 *  Created by David Bryant on 4/03/08.
 *  Copyright 2008 David Bryant
 *
 * Implements an exception designed specifically for debugging algorithms on
 * trees
 */

#ifndef PHYLOEXCEPTION_H_INCLUDE
#define PHYLOEXCEPTION_H_INCLUDE

using namespace std;

namespace Phylib {
	class PhylibException {
		public:
		PhylibException(const string& msg) :
			_msg(msg) {
		}
		~PhylibException() {
		}
		string getMessage() const {
			return (_msg);
		}
		private:
		string _msg;
	};
}

#endif


