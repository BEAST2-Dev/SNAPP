/**
 *  Collection of random generation utilities.
 * 
 *
 *  @author David Bryant
 *
 */

//TODO: Use a better random number generator.
//TODO: Rewrite this using a class, to incorporate better random number generator and get rid of inlined functions.


#ifndef RANDOM_H_INCLUDE
#define RANDOM_H_INCLUDE

#include"../global/stdIncludes.h"

namespace Phylib {
	
	inline void seed_random() {std::srand(getpid());}
	
	inline void seed_random(int r_seed) {std::srand(r_seed);}
	
	inline unsigned int random_num(unsigned int x) {
		if (x==0) 
			return 0;
		return std::rand()%x;
	}
	
	//gives a random double between 0 and 1
	inline double randu() {
		return ((double)(std::rand())/RAND_MAX);
	}
	
	
	/**
	 * Generates random double from exponential with given mean
	 *
	 * @param mean
	 * @return random double
	 */
	inline double random_exp(double mean) {
		return -mean * std::log(randu());
	}	
	
	/**
	 Generates random gaussian
	 **/
	inline double random_gaussian() {
		double x1, x2, w;
		
		do {
			x1 = 2.0 * randu() - 1.0;
			x2 = 2.0 * randu() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
		
		w = sqrt( (-2.0 * std::log( w ) ) / w );
		return x1*w;
		//y1 = x1 * w;
		//y2 = x2 * w;
	}
	
	
	/**
	 Generate a random number from the gamma distribution with parameter alpha (and beta =alpha)
	 
	 //TODO: Rewrite random number generators as class
	 **/
	inline double random_gamma(double alpha, double beta = 1.0) {
		if (alpha < 1.0) {
			double u = randu();
			double gamma = random_gamma(alpha + 1.0);
			return beta * gamma * std::pow(u, 1.0 / alpha);
		} else {
			double d, c, x, v, u;
			d = alpha - 1.0 / 3.0;
			c = 1.0 / std::sqrt(9.0 * d);
			for (; ;) {
				do {
					x = random_gaussian();
					v = 1.0 + c * x;
				} while (v <= 0.0);
				v = v * v * v;
				u = randu();
				if (u < 1.0 - 0.0331 * (x * x) * (x * x))
					return (beta * d * v);
				if (std::log(u) < 0.5 * x * x + d * (1.0 - v + std::log(v)))
					return (beta * d * v);
			}
		}
	}
	
	
	
}


#endif


