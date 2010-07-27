/**
 * @file phylo_io.h
 * @brief Header file generic input/output of phylos
 * @author David Bryant
 * @version 0.1
 */



#ifndef PHYLO_IO_INCLUDE
#define PHYLO_IO_INCLUDE

#include"../global/stdIncludes.h"
#include"../utilities/phylibException.h"
#include"phylo.h"


namespace Phylib {
	class basic_newick {
	public:
		int id;
		double length;
		string meta_data;
		
		//TODO: Check which of these constructors we actually use/need.
		basic_newick() :
		id(-1), length(0.0), meta_data("") {
		}
		basic_newick(int val) :
		id(val), length(0.0), meta_data("") {
		}
		basic_newick(int val, double val2) :
		id(val), length(val2) , meta_data("") {
		}
		basic_newick(int val, double val2, const string& val3):
		id(val), length(val2) , meta_data(val3) {
		}
		
		void copy(const basic_newick& data) {
			if (this!=&data) {
				id = data.id;
				length = data.length;
				meta_data = data.meta_data;
			}
		}
		
		basic_newick& operator=(const basic_newick& data) {
			copy(data);
			return *this;
		}
	};
	
	/**
	 * Prints out Newick format for a subtree in a phylogeny for which the node inherits the basic_newick class. Edge lengths are printed
	 * unless print_lengths is set to false. Only non-negative id values will be printed.
	 * @param os output stream
	 * @param v root of subtree
	 * @param print_lengths true (default) if branch lengths are to be printed.
	 * @param print_meta Print metadata in brackets [] after each node with metadata.
	 */
	template<typename T> void print_newick(ostream& os, typename phylo<T>::const_iterator v, bool print_lengths = true, bool print_meta = true) {
		if (!v.leaf()) {
			//Recurse
			os<<"(";
			for (typename phylo<T>::const_iterator p = v.left(); !p.null(); p =p.right()) {
				print_newick<T>(os, p, print_lengths);
				if (!p.right().null())
					os<<",";
			}
			os<<")";
		}
		if ((*v).id>=0)
			os<<(*v).id;
		if (print_meta && !(*v).meta_data.empty())
			os<<"["<<(*v).meta_data<<"]";
		if (print_lengths && !v.root()) {
			os<<":"<<(*v).length;
		}
	}
	
	/**
	 * Prints out Newick format for any phylogeny for which the node inherits the basic_newick class. Edge lengths are printed
	 * unless print_lengths is set to false. Only non-negative id values will be printed.
	 * @param os output stream
	 * @param tree phylogeny. Node type must inherit basic_newick (or equivalent)
	 * @param print_lengths true (default) if branch lengths are to be printed.
	 * @param print_meta Print metadata in brackets [] after each node with metadata.
	 */
	template<typename T> void print_newick(ostream & os, const phylo<T>& tree, bool print_lengths = true,bool print_meta = true) {
		if (tree.empty())
			return;
		for (typename phylo<T>::const_iterator p = tree.header().left(); !p.null(); p=p.right()) {
			print_newick<T>(os, p, print_lengths,print_meta);
			os<<";";
			if (!p.right().null())
				os<<endl;
		}
	}
	
	/**
	 * Prints out Newick format for a subtree in a phylogeny for which the node inherits the basic_newick class. Edge lengths are printed
	 * unless print_lengths is set to false. Taxon names for non-negative ids will be printed
	 * @param os output stream
	 * @param v root of subtree
	 * @param taxa_names vector of taxa names
	 * @param print_lengths true (default) if branch lengths are to be printed.
	 * @param print_meta Print metadata in brackets [] after each node with metadata.
	 */
	template<typename T> void print_newick(ostream& os, typename phylo<T>::const_iterator v, const vector<string>& taxa_names, bool print_lengths = true, bool print_meta = true) {
		if (!v.leaf()) {
			//Recurse
			os<<"(";
			for (typename phylo<T>::const_iterator p = v.left(); !p.null(); p =p.right()) {
				print_newick<T>(os, p, taxa_names, print_lengths);
				if (!p.right().null())
					os<<",";
			}
			os<<")";
		}
		if ((*v).id>=0)
			os<<taxa_names[(*v).id];
		if (print_meta && !(*v).meta_data.empty())
			os<<"["<<(*v).meta_data<<"]";
		if (print_lengths && !v.root()) {
			os<<":"<<(*v).length;
		}
	}
	
	/**
	 * Prints out Newick format for any phylogeny for which the node inherits the basic_newick class. Edge lengths are printed
	 * unless print_lengths is set to false. Taxon names for non-negative id values will be printed.
	 * @param os output stream
	 * @param tree phylogeny. Node type must inherit basic_newick (or equivalent)
	 * @param print_lengths true (default) if branch lengths are to be printed.
	 * @param print_meta Print metadata in brackets [] after each node with metadata.
	 */
	template<typename T> void print_newick(ostream & os, const phylo<T>& tree, const vector<string>& taxa_names, bool print_lengths = true, bool print_meta = true) {
		if (tree.empty())
			return;
		for (typename phylo<T>::const_iterator p = tree.header().left(); !p.null(); p=p.right()) {
			print_newick<T>(os, p, taxa_names, print_lengths,print_meta);
			os<<endl;
		}
	}
	
	static const int io_no_token = 0;
	static const int io_tree_left=1;
	static const int io_tree_right=2;
	static const int io_tree_comma=3;
	static const int io_tree_name=4;
	static const int io_tree_eof=0;
	static const int io_tree_colon=5;
	static const int io_tree_semicolon=6;
	static const int io_tree_comment = 7;
	
	template<typename T> int peek_next_tree_token(istream& is) {
		
		is >> ws; //extract whitespace
		
		int ch = is.peek(); //FIX THIS.
		
		if (ch==EOF || ch==10)
			return io_tree_eof;
		
		if (ch=='(')
			return io_tree_left;
		if (ch==')')
			return io_tree_right;
		if (ch==',')
			return io_tree_comma;
		if (ch==':')
			return io_tree_colon;
		if (ch==';')
			return io_tree_semicolon;
		if (ch=='[')
			return io_tree_comment;

		
		return io_tree_name;
	}
	
	template<typename T> int next_tree_token(istream& is, string& next_word) {
		
		int val = peek_next_tree_token<T>(is);
		char ch = is.get();
		next_word.erase();
		
		if (val==io_tree_comment) { //TODO: handle nested []
			is>>ch; //skip the '['
			while(is&&(ch!=0)&&(ch!=']')) {
				next_word+=ch;
				is>>ch;
			}
			return val;
		}
		
		if (val!= io_tree_name)
			return val;
		
		while (is&&(ch!=0)&&(!isspace(ch) ) && (ch!='(') && (ch!=')') && (ch!=',') && (ch!=':') && (ch!=';')&&(ch!='[')) {
			next_word+=ch;
			is>>ch;
		}
		is.putback(ch);
		
		return io_tree_name;
	}
	
	/**
	 *parse_newick	
	 *
	 *Reads in a tree or subtree in newick format. A subtree is defined by:
	 * tree = subtree;
	 * subtree =  leaf | node
	 * leaf = taxa_name[:length]
	 * node = (subtree,....,subtree)[taxa_name][:length]
	 * Except: if we are at the toplevel then we cannot have a length and can only have a taxa name 
	 *@param is	input stream
	 *@param taxa_names  uses these names to set the id field for leaves. If a new taxon is encountered, it is added to the list.
	 *@param defaultLength	branch length used if none is specified. Default = 0.0.
	 */
	//TODO: Add support for internal labels when reading trees.
	template<typename T> void read_newick(istream& is, phylo<T>& tree, vector<string>& taxa_names, double defaultLength, bool toplevel) { //reads in tree in phylip format.
		
		string next_word;
		int next_tok;
		next_tok=peek_next_tree_token<T>(is);
		tree.clear();
		
		bool at_leaf = true;
		
		//Read in a subtree if there is one.
		if (next_tok == io_tree_left) {
			next_tree_token<T>(is, next_word); //Swallow '('
			at_leaf = false;
			typename phylo<T>::iterator new_node = tree.insert_child(tree.header()); //Parent node.
			typename phylo<T>::iterator new_child;
			new_child.set_null();
			while (true) { //Loop through subtrees
				phylo<T> subtree;
				read_newick<T>(is, subtree, taxa_names, defaultLength, false);
				if (new_child.null()) {
					new_child = tree.graft_child(new_node, subtree);
				} else {
					new_child = tree.graft_sibling(new_child, subtree);
				}
				next_tok = next_tree_token<T>(is, next_word);
				if (next_tok==io_tree_right) {
					if (new_node.left().null() )
						throw new PhylibException("Error: Error in Newick syntax");
					else
						break;
				} else if (next_tok!=io_tree_comma) {
					throw new PhylibException("Error: Expected ',' or ')' ");
				}
			}
			next_tok = peek_next_tree_token<T>(is);
		}
		
		//Check for taxa
		if (next_tok == io_tree_name) {
			next_tree_token<T>(is, next_word);
			unsigned int id;
			for (id=0; (id<taxa_names.size())&&(next_word!=taxa_names[id]); id++)
				;
			if (id==taxa_names.size()) {
				taxa_names.push_back(next_word); //Can't find the taxa name. Make a new one
				cerr<<"Found taxon '"<<next_word<<"'"<<endl;
			}
			if (at_leaf) {
				typename phylo<T>::iterator tmp = tree.insert_child(tree.header());
				tmp->id = id;
			} else
				tree.header().left()->id = id;
			next_tok = peek_next_tree_token<T>(is);
		} else {
			if (at_leaf) {
				throw new PhylibException("Error: Expected taxon name");
			}
		}
		//Check for metadata
		if (next_tok == io_tree_comment) {
			next_tree_token<T>(is, next_word);
			(*(tree.header().left())).meta_data=next_word;
			next_tok = peek_next_tree_token<T>(is);
		}
			
		//Check for branch length
		if (next_tok == io_tree_colon) {
			next_tree_token<T>(is, next_word);
			double weight;
			is>>weight;
			(*(tree.header().left())).length=weight;
			next_tok = peek_next_tree_token<T>(is);
		} else {
			(*(tree.header().left())).length=defaultLength; //if not, give default weight
		}
		
		//Check for semi-colon
		if (next_tok == io_tree_semicolon) {
			next_tree_token<T>(is, next_word);
			if (!toplevel)
				throw new PhylibException("Error: Read ';' in middle of tree description");
			return;
		}
		
		return;
	}
	
	template<typename T> void read_newick(string s, phylo<T>& tree, vector<string>& taxa_names, double default_branchlength = 0.0) {
		using namespace std;
		stringstream is(s, stringstream::in);
		read_newick<T>(is, tree, taxa_names, default_branchlength, true);
		return;
	}
	
	template<typename T> void read_newick(istream& is, phylo<T>& tree, vector<string>& taxa_names, double default_branchlength = 0.0) {
		using namespace std;
		read_newick<T>(is, tree, taxa_names, default_branchlength, true);
		return;
	}
	
	
	
}
#endif

