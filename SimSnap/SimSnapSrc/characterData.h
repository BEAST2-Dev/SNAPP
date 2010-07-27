/*
 *  characterData.h
 *  SingleSiteSorter
 *
 *  Created by David Bryant on 5/02/09.
 *  Copyright 2009 David Bryant
 *
 */
#ifndef CHARACTER_DATA_H
#define CHARACTER_DATA_H

using namespace std;

class Exception {
public: 
	Exception() : message("") {}
	Exception(string ex) : message(ex) {};
	string getMessage() {return message;}
private:
	string message;
};

/**
 class containing a single site Pattern... used for sorting.
 **/
template<typename BASE> class sitePattern {
private:
	vector<BASE> pattern;
	
public:
	int id;
	
	
	sitePattern() {
		pattern.clear();
		id = -1;
	}
	
	sitePattern(const sitePattern& x) {
		copy(x);
	} 
	
	sitePattern(const vector<BASE>& c, int id) {
		pattern.resize(c.size());
		std::copy(c.begin(),c.end(),pattern.begin());
		this->id = id;
	}
	
	~sitePattern() {
		pattern.clear();
	}
	
	bool operator<(const sitePattern& other) const {
		uint n = this->pattern.size();
		if (n!=other.pattern.size())
			throw Exception("Comparing sitePatterns of different sizes");
		uint i=0;
		while(i<n) {
			BASE f = this->pattern[i], s = other.pattern[i];
			if (f<s)
				return true;
			else if (f>s)
				return false;
			i++;
		}
		return false;
	}
	
	void copy(const sitePattern& other) {
		pattern.resize(other.pattern.size());
		std::copy(other.pattern.begin(),other.pattern.end(),pattern.begin());
		id = other.id;
	}
	
	sitePattern& operator=(const sitePattern& other)  {
		if (&other != this)
			copy(other);
		
		return *this;
	}
	
	void getPattern(vector<BASE>& c) {
		c.resize(pattern.size());
		std::copy(pattern.begin(),pattern.end(),c.begin());
	}
	
	bool operator==(const sitePattern& other) {
		uint n = this->pattern.size();
		if (n!=other.pattern.size())
			throw Exception("Comparing sitePatterns of different sizes");
		uint i=0;
		while(i<n) {
			if (this->pattern[i] != other.pattern[i])
				return false;
			i++;
		}
		return true;
	}
	
};


template<typename BASE> 
class characterData {
private:
	vector<vector<BASE> > patterns;
	vector<uint> site2pattern; //Which pattern is there at each site.
	vector<uint> freqs; 
public:
	
	characterData() {}
	
	/**
	 Takes the characters, sorts out unique characters, and sets up storage.
	 
	 @param const vector<vector<BASE> >& chars   Character data.
	 @param bool transposed		flag set to true if data is transposed (first index is site number) or not (first index is taxon number)
	 **/
	
	characterData(const vector<vector<BASE> >& chars, bool transposed=false) {
		setChars(chars, transposed);
	}
	
	~characterData() {
		for(uint i=0;i<patterns.size();i++)
			patterns[i].clear();
		patterns.clear();
		site2pattern.clear();
		freqs.clear();
	}
	
	
	void setChars(const vector<vector<BASE> >& chars, bool transposed) {
		
		uint ntax, nsites;
		if (transposed) {
			nsites = chars.size();
			ntax = chars[0].size();
		} else {
			ntax = chars.size();
			nsites = chars[0].size();
		}
		patterns.clear();
		freqs.clear();
		site2pattern.resize(nsites);
		
		vector< sitePattern<BASE> > theCharList;
		for(uint i=0;i<nsites;i++) {
			if (!transposed) {
				vector<uint> thisSite(ntax);
				for(uint j=0;j<ntax;j++)
					thisSite[j] = chars[j][i];
				theCharList.push_back(sitePattern<BASE>(thisSite,i));
			}
			else {
				//cerr<<"chars["<<i<<"] = ";
				//for(uint j=0;j<chars[i].size();j++)
				//	cerr<<" \t "<<chars[i][j];
				//cerr<<endl;
				sitePattern<BASE> theSite(chars[i],i);
				theCharList.push_back(theSite);
			}
		}
		
		
		sort(theCharList.begin(),theCharList.end());
		
		/*for(typename list<sitePattern<BASE> >::iterator i = theChars.begin(); i!=theChars.end();i++) {
		 vector<uint> pat;
		 (*i).getPattern(pat);
		 cout<<(*i).id<<"\t";
		 for(uint j=0;j<pat.size();j++)
		 cout<<pat[j]<<"\t";
		 cout<<endl;
		 }*/
		
		
		uint npatterns = 0;
		typename vector<sitePattern<BASE> >::iterator from = theCharList.end();
		for(typename vector<sitePattern<BASE> >::iterator i = theCharList.begin(); i!=theCharList.end();i++) {
			if ( from==theCharList.end() || !( (*i)==(*from))) {
				from = i;
				patterns.resize(npatterns+1);
				(*from).getPattern(patterns[npatterns]);
				freqs.resize(npatterns+1);
				freqs[npatterns]=1;
				site2pattern[(*i).id] = npatterns;
				npatterns++;
			} else {
				site2pattern[(*i).id]=npatterns-1;//cout<<(*i).id<<"->"<<npatterns-1<<endl;
				freqs[npatterns-1]++;
			}
		}
		theCharList.clear();
	}
	
	void clear() {
		for(uint i=0;i<patterns.size();i++)
			patterns[i].clear();
		patterns.clear();
		site2pattern.clear();
		freqs.clear();
	}
	
	uint getNSites() const {
		return site2pattern.size();
	}
	
	void getSite(uint id, vector<BASE>& site) const {
		uint pattern_id = site2pattern[id];
		site.resize(patterns[pattern_id].size());
		std::copy(patterns[pattern_id].begin(),patterns[pattern_id].end(),site.begin());
	}
	
	/** Get the base for a given taxon at a site 
	 **/
	BASE get(uint taxon, uint site) const {
		return patterns[site2pattern[site],taxon];
	}
	
	/**
	 Get the base for a give taxon in a given pattern
	 **/
	BASE getp(uint taxon, uint pattern) const {
		return patterns[pattern][taxon];
	}
	
	uint getNPatterns() const {
		return patterns.size();
	}
	
	void getPattern(uint id, vector<BASE>& pattern) const {
		pattern.resize(patterns[id].size());
		std::copy(patterns[id].begin(),patterns[id].end(),pattern.begin());
	}
	
	uint getPatternFreq(uint id) const {
		return freqs[id];
	}
	
	uint getSitePatternId(uint id) const {
		return site2pattern[id];
	}
	
	void printPatterns(ostream& os,string separator = "\t") const {
		for(uint i=0;i<patterns.size();i++) {
			os<<i<<"\t|\t";
			for(uint j=0;j<patterns[i].size();j++) {
				if (j>0)
					os<<separator;
				os<<patterns[i][j];
			}
			os<<"\n";
		}
		os<<"\n";
		//for(uint i = 0;i<site2pattern.size();i++)
		//	os<<i<<"\t"<<site2pattern[i]<<endl;
		//os<<"\n";
	}
	
	void printSites(ostream& os,string separator = "\t") const {
		for(uint i=0;i<site2pattern.size();i++) {
			os<<i<<"\t|\t";
			uint id = site2pattern[i];
			for(uint j=0;j<patterns[id].size();j++) {
				if (j>0)
					os<<separator;
				os<<patterns[id][j];
			}
			os<<"\n";
		}
	}
	
};

/**
 Read in the allele frequencies from a tab delimited file.
 
 The format of the file depends on whether the field transposed is set to false (the default) or true.
 
 if transposed = false, the format is
 
 <number-of-species> <number-of-characters>
 <species-name-1> <sample-size-species-1> f_11 f_12 f_13 ..... f_1m
 <species-name-2> <sample-size-species-2> f_21 f_22 f_23 ..... f_2m
 ...
 <species-name-n> <sample-size-species-n> f_n1 f_n2 f_n3 ..... f_nm
 
 
 if transposed = true, the format is
 
 <number-of-species> <number-of-characters>
 <species-name-1>  <species-name-2>  ...  <species-name-n>
 <sample-size-species-1>   <sample-size-species-2>   ...    <sample-size-species-n>
 f_11  f_21    ...    f_n1
 f_12  f_22    ...    f_n2
 
 ...
 
 f_1m  f_2m    ...    f_nm
 
 In both cases, <species-name-i> contains no whitespace, n is the number of species, m is the number of characters, and 
 f_ij is the frequency of the 1 allele for character j in species i.
 **/

template<typename BASE> void readCharactersTable(istream& is, characterData<BASE>& data, vector<string>& taxanames, vector<uint>& sampleSizes, bool transposed = false) {
	//FIrst two numbers are the number of species and the number of characters.
	uint ntax, nsites;
	is>>ntax;
	is>>nsites;
	
	if (ntax==0 || nsites == 0) {
		cerr<<"Error reading number of species or number of sites. Aborting."<<endl;
		exit(1);
	}
	vector<vector<uint> > counts;
	taxanames.clear();
	counts.resize(ntax);
	sampleSizes.resize(ntax);
	
	if (!transposed) {
		//Now the format is taxa name (without spaces or crap) followed by nsites integers.
		for(uint n=0;n<ntax;n++) {
			counts[n].resize(nsites);
			string name;
			is>>name;
			taxanames.push_back(name);
			is>>sampleSizes[n];
			for(uint i=0;i<nsites;i++) {
				is>>counts[n][i];
			}
		}
	}
	else {
		for(uint n=0;n<ntax;n++) {
			string name;
			is>>name;
			taxanames.push_back(name);
		}
		for(uint n=0;n<ntax;n++) {
			is>>sampleSizes[n];
			counts[n].resize(nsites);
		}
		
		for(uint i=0;i<nsites;i++) 
			for(uint n=0;n<ntax;n++) 
				is>>counts[n][i];
	}
	
	data.setChars(counts,false);
}


#endif

