/*
 *  optimisation.h
 *  SingleSiteSorter
 *
 *  Created by David Bryant on 24/11/08.
 *  Copyright 2008 David Bryant
 *
 */


#ifndef OPTIMISATION_H
#define OPTIMISATION_H

#include"../global/stdIncludes.h"
#include "../utilities/phylibException.h"


namespace Phylib {
	
	using namespace std;
	
	template<typename DOUBLE> 
	class MultiDimFunction {
	public:
		virtual DOUBLE evaluate(const vector<DOUBLE>& x) = 0;
		DOUBLE evalAlongLine(const vector<DOUBLE>& x, const vector<DOUBLE>& s, DOUBLE alpha) {
			vector<DOUBLE> y(x.size());
			for(int i=0;i<x.size();i++)
				y[i] = x[i]+alpha*s[i];
			return evaluate(y);
		}
		uint getNumDimensions() { return ndim;}
		void setNumDimensions(uint dim) {ndim = dim;}
	private:
		uint ndim;
	};
	
	
	template<typename DOUBLE> DOUBLE simulatedAnnealing(MultiDimFunction<DOUBLE>& f, 	
														vector<DOUBLE>& x, 
														vector<DOUBLE>& lowerBounds,
														vector<DOUBLE>& upperBounds,	
														DOUBLE stepSize, 
														uint numSteps) {
		//TODO: 
		//Improved stopping conditions.
		DOUBLE qx = f.evaluate(x);
		uint n = x.size();
		vector<DOUBLE> xprime(n);
		DOUBLE Tn = 1.0;
		const uint nAcceptsPerTempChange = 20;
		const DOUBLE coolingRate = 0.9;
		
		
		
		//Form a list of moves for which there is room to move (allows us to set constants in some dimensions)
		vector<uint> moves;
		moves.clear();
		for(uint i=0;i<n;i++)
			if (lowerBounds[i]<upperBounds[i])
				moves.push_back(i);
		
		uint nAccepts = 0; //Number of accepts since last temperature change.
		for(int i=0;i<numSteps;i++) {
			
			//Generate next value.
			std::copy(x.begin(),x.end(),xprime.begin());
			//Choose a random dimension at uniform.
			uint m = moves[random_num(moves.size())];
			
			xprime[m] += (2.0*randu()-1.0)*stepSize;
			
			//Check if xprime is feasible
			if (xprime[m]> lowerBounds[m] && xprime[m] < upperBounds[m]) {
				DOUBLE qxprime = f.evaluate(xprime);
				DOUBLE alpha = std::exp(-(qxprime - qx)/Tn);
				
				cerr<<"Trying x'=[";
				for(int j=0;j<n;j++) {
					cerr<<xprime[j];
					if (j<n-1)
						cerr<<", ";
				}
				cerr<<"]\tqx = "<<qxprime<<"\tTn ="<<Tn<<endl;
				
				if (randu()<alpha) {
					//Accept xprime.
					qx = qxprime;
					x[m] = xprime[m];
					cerr<<"Accept"<<endl;
					nAccepts++;
				} else {
					cerr<<"Reject"<<endl;
				}
				if (nAccepts>=nAcceptsPerTempChange) {
					Tn *= coolingRate;
					nAccepts=0;
				}
			}
		}
		
		return qx;
		
		
	}
	
	/****
	 Bracket the minimum
	 
	 Extends [a,c] so that it contains a point b for which f(b)<f(a) and f(b) < f(c). 
	 This guarantees that the minimum is contained in the interval [a,c]
	 which are then used to start a line search.
	 
	 Here f(alpha) is actually F[x + alpha*s] where x s are vectors.
	 @param x Current point; start of line search.
	 @param s Direction of line search
	 @param a,c initial interval. Replaced by bracketting interval.
	 OUTPUT:
	 @param fa value of f(a) at termination.
	 @param fc value of f(c) at termination.
	 Returns false if no bracket can be found (e.g. function unbounded below in this direction)
	 
	 *****/
	
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
	
	
	
	template<typename DOUBLE> bool bracketLinesearch(MultiDimFunction<DOUBLE>& f,
													 const vector<DOUBLE>& x,
													 const vector<DOUBLE>& s,
													 DOUBLE& a, 
													 DOUBLE& c, 
													 DOUBLE& fa, 
													 DOUBLE& fc)
	{
		const DOUBLE GOLD =   (1.0 + std::sqrt(5))/2.0;;
		const DOUBLE TINY = 10e-20;
		const DOUBLE GLIMIT = 100.0;
		DOUBLE b = c;
		
		DOUBLE fb; 
		DOUBLE ulim,u,r,q,fu;
		DOUBLE tmp;
		
		
		fa=f.evalAlongLine(x, s, a); //fa = f(a)
		fb=f.evalAlongLine(x, s, b); //fb = f(b)
		
		if (fb > fa) {
			tmp = a; a=b; b=tmp;
			tmp = fa; fa = fb; fb = tmp;
		}
		
		c=(b)+GOLD*(b-a); //Initial guess for c (outside the interval)
		fc=f.evalAlongLine(x, s, c); 
		while (fb > fc) {
			
			
			r=(b-a)*(fb-fc);
			q=(b-c)*(fb-fa);
			u=(b)-((b-c)*q-(b-a)*r)/(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
			ulim=(b)+GLIMIT*(c-b);
			if ((b-u)*(u-c) > 0.0) {
				fu=f.evalAlongLine(x, s, u); //fu = f(u)
				if (fu < fc) {
					a=(b);
					b=u;
					fa=fb;
					fb=fu;
					return true;
				} else if (fu > fb) {
					c=u;
					fc=fu;
					return false;
				}
				u=(c)+GOLD*(c-b);
				fu=f.evalAlongLine(x, s, u); //fu = f(u)
			} else if ((c-u)*(u-ulim) > 0.0) {
				fu=f.evalAlongLine(x, s, u); //fu = f(u)
				if (fu < fc) {
					SHFT(b,c,u,c+GOLD*(c-b))
					SHFT(fb,fc,fu,f.evalAlongLine(x, s, u))
				}
			} else if ((u-ulim)*(ulim-c) >= 0.0) {
				u=ulim;
				fu=f.evalAlongLine(x, s, u); //fu = f(u)
			} else {
				u=(c)+GOLD*(c-b);
				fu=f.evalAlongLine(x, s, u); //fu = f(u)
			}
			SHFT(a,b,c,u)
			SHFT(fa,fb,fc,fu)
		}
		return false;
	}
	
	template<typename DOUBLE> DOUBLE brent(MultiDimFunction<DOUBLE>& f,const vector<DOUBLE>& X, const vector<DOUBLE>& S, DOUBLE a, DOUBLE c, DOUBLE fa, DOUBLE fc, DOUBLE tol, DOUBLE& xmin)
	{
		int iter;
		DOUBLE d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
		DOUBLE e=0.0;
		
		const DOUBLE ITMAX = 100;
		const DOUBLE GOLD = (3.0 - std::sqrt(5))/2.0;
		const DOUBLE ZEPS = 1.0e-10;
		
		
		a=((a < c) ? a : c);
		DOUBLE b=((a > c) ? a : c);
		
		x=w=v=b;
		fw=fv=fx=f.evalAlongLine(X, S, x);
		for (iter=1;iter<=ITMAX;iter++) {
			xm=0.5*(a+b);
			tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
			if (std::fabs(x-xm) <= (tol2-0.5*(b-a))) {
				xmin=x;
				cerr<<"xmin = "<<xmin<<endl;
				return fx;
			}
			if (fabs(e) > tol1) {
				r=(x-w)*(fx-fv);
				q=(x-v)*(fx-fw);
				p=(x-v)*q-(x-w)*r;
				q=2.0*(q-r);
				if (q > 0.0) p = -p;
				q=fabs(q);
				etemp=e;
				e=d;
				if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
					d=GOLD*(e=(x >= xm ? a-x : b-x));
				else {
					d=p/q;
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				}
			} else {
				d=GOLD*(e=(x >= xm ? a-x : b-x));
			}
			u=(std::fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
			fu=f.evalAlongLine(X, S, u); //fu = f(u)
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				SHFT(v,w,x,u)
				SHFT(fv,fw,fx,fu)
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					v=w;
					w=u;
					fv=fw;
					fw=fu;
				} else if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
			}
		}
		throw PhylibException("Too many iterations in BRENT");
		return(10e20);
	}
	
	
	template<typename DOUBLE> bool linesearch(MultiDimFunction<DOUBLE>& f,const vector<DOUBLE>& x, const vector<DOUBLE>& s, DOUBLE& alpha, DOUBLE& fmin) {
		
		const DOUBLE TOL = 2.0e-4;
		DOUBLE a = 0.0;
		DOUBLE c = 1.0; //Initial interval. 
		
		//			DOUBLE alphaMin, alphaMax;			
		
		DOUBLE fa,fc;
		
		//FInd a,b,c so that f(a)>
		if (!bracketLinesearch(f,x, s, a, c, fa, fc)) {
			alpha = 0;
			fmin = fa;
			return false;
		}
		
		fmin = brent<DOUBLE>(f,x, s, a, c, fa, fc, TOL, alpha);	
		return 	true;
	}	
	
	template<typename DOUBLE> DOUBLE simpleSearch(MultiDimFunction<DOUBLE>& f, vector<DOUBLE>& x, int numSteps) {
		int n=f.getNumDimensions();
		vector<DOUBLE> s(n);
		std::fill(s.begin(),s.end(),0.0);
		DOUBLE fLineMin = 10e20;
		
		for (int i=0;i<numSteps;i++) {
			for (int j=0;j<n;j++) {
				s[j]=1.0;
				DOUBLE alpha=0.5;
				DOUBLE fLineMin;
				bool success = linesearch(f, x, s,alpha, fLineMin);
				if (!success) {
					s[j]=-1.0;
					alpha = 0.5;
					success = linesearch(f, x, s, alpha, fLineMin);
				}
				x[j] = x[j] + alpha*s[j];
				
				cerr<<"x=[";
				for(int k=0;k<n;j++) {
					cerr<<x[k];
					if (k<n-1)
						cerr<<", ";
				}
				cerr<<"]\tf(x) = "<<fLineMin<<endl;
				
				
				
				s[j]=0.0;
			}
		}
		return fLineMin;
	}
	
	
	
}


#endif

