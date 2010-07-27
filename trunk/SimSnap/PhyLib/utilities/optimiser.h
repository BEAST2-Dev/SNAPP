/*
 *  Optimiser.h
 * 
 *
 *  Created by David Bryant on 10/11/08.
 *
 * Generic class to apply multidimensional optimisation routines
 */

#ifndef OPTIMISER_H
#define OPTIMISER_H

#include"../global/stdIncludes.h"
#include "../utilities/phylibException.h"

//TODO: Add bool option to say constrained or unconstrained.


namespace Phylib {
	
	using namespace std;
	
	
	template<typename DOUBLE> 
	class MultiDimFunction {
	public:
		virtual DOUBLE evaluate(const vector<DOUBLE>& x) = 0;
		uint getNumDimensions() { return ndim;}
	private:
		uint ndim;
	};
	
	template<typename DOUBLE> 
	class praxis {
		
		
		/**
		 * methods for minimization of a real-valued function of
		 * several variables without using derivatives
		 * - algorithm is based on Brent's modification of a conjugate direction
		 *   search method proposed by Powell (praxis)
		 *   Richard P. Brent.  1973.   Algorithms for finding zeros and extrema
		 */
		
	public:
		
		// Variables that control aspects of the inner workings of the
		// minimization algorithm. Setting them is optional, they
		// are all set to some reasonable default values given below.
		
		/**
		 *  controls the printed output from the routine
		 *  (0 -> no output, 1 -> print only starting and final values,
		 *   2 -> detailed map of the minimization process,
		 *   3 -> print also eigenvalues and vectors of the
		 *   search directions), the default value is 0
		 */
		const int prin = 0;
		
		/**
		 * step is a steplength parameter and should be set equal
		 * to the expected distance from the solution.
		 * exceptionally small or large values of step lead to
		 * slower convergence on the first few iterations
		 * the default value for step is 1.0
		 */
		const double step = 1.0;
		
		/**
		 * scbd is a scaling parameter. 1.0 is the default and
		 * indicates no scaling. if the scales for the different
		 * parameters are very different, scbd should be set to
		 * a value of about 10.0.
		 */
		const double scbd = 1.0;
		
		/**
		 * illc  should be set to true
		 * if the problem is known to
		 * be ill-conditioned. the default is false. this
		 * variable is automatically set, when the problem
		 * is found to to be ill-conditioned during iterations.
		 */
		const bool illc = false;
		
		const DOUBLE EPSILON = 2.220446049250313E-16;
		const DOUBLE SQRT_EPSILON = 1.4901161193847656E-8;
		const DOUBLE SQRT_SQRT_EPSILON = 1.220703125E-4;
		
		void optimize(MultiDimFunction& f, vector<DOUBLE>&  xvector, double tolfx, double tolx, const vector<DOUBLE>& lowerBounds, const vector<DOUBLE>& upperBounds /*, MinimiserMonitor monitor*/)
		{
			t = tolx;
			
			MultiDimFunction& fun = f;
			x = xvector;
			
			checkBounds(x, lowerBounds, upperBounds);
			h = step;
			
			dim = fun.getNumArguments();;
			
			d = new double[dim];
			y = new double[dim];
			z = new double[dim];
			q0 = new double[dim];
			q1 = new double[dim];
			v = new double[dim][dim];
			tflin = new double[dim];
			
			small = EPSILON*EPSILON;
			vsmall = small*small;
			large = 1.0/small;
			vlarge = 1.0/vsmall;
			ldfac = (illc ? 0.1 : 0.01);
			nl = kt = 0;
			numFun = 1;
			fx = fun.evaluate(x);
			
			stopCondition(fx, x, tolfx, tolx, true);
			
			qf1 = fx;
			t2 = small + std::abs(t);
			t = t2;
			dmin = small;
			
			if (h < 100.0*t) h = 100.0*t;
			ldt = h;
			for (i = 0; i < dim; i++)
			{
				for (j = 0; j < dim; j++)
				{
					v[i][j] = (i == j ? 1.0 : 0.0);
				}
			}
			d[0] = 0.0;
			qd0 = 0.0;
			for (i = 0; i < dim; i++)
				q1[i] = x[i];
			
			if (prin > 1)
			{
				cout<<"\n------------- enter function praxis -----------\n"<<endl;
				cout<<"... current parameter settings ..."<<endl;
				cout<<"... scaling ... " << scbd<<endl;
				cout<<"...   tolx  ... "<<t<<endl;
				cout<<"...  tolfx  ... "<<tolfx<<endl;
				cout<<"... maxstep ... "<<h<<endl;
				cout<<"...   illc  ... "<<illc<<endl;
				cout<<"... maxFun  ... "<<maxFun<<endl;
			}
			if (prin > 0)
				cout<<endl;
			
			while(true)
			{
				sf = d[0];
				s = d[0] = 0.0;
				
				/* minimize along first direction */
				
				min1 = d[0];
				min2 = s;
				min(fun,0, 2, fx, false, lowerBounds, upperBounds);
				d[0] = min1;
				s = min2;
				
				if (s <= 0.0)
					for (i = 0; i < dim; i++)
					{
						v[i][0] = -v[i][0];
					}
				if ((sf <= (0.9 * d[0])) || ((0.9 * sf) >= d[0]))
					for (i=1; i < dim; i++)
						d[i] = 0.0;
				
				bool gotoFret = false;
				for (k=1; k < dim; k++)
				{
					for (i=0; i< dim; i++)
					{
						y[i] = x[i];
					}
					sf = fx;
					illc = illc || (kt > 0);
					
					bool gotoNext;
					do
					{
						kl = k;
						df = 0.0;
						if (illc)
						{        /* random step to get off resolution valley */
							for (i=0; i < dim; i++)
							{
								z[i] = (0.1 * ldt + t2 * Math.pow(10.0,(double)kt)) * (MathUtils.nextDouble() - 0.5);
								s = z[i];
								for (j=0; j < dim; j++)
								{
									x[j] += s * v[j][i];
								}
							}
							
							checkBounds(x, lowerBounds, upperBounds);
							
							fx = fun.evaluate(x);
							numFun++;
						}
						
						/* minimize along non-conjugate directions */
						for (k2=k; k2 < dim; k2++)
						{
							sl = fx;
							s = 0.0;
							
							min1 = d[k2];
							min2 = s;
							min(fun,k2, 2, fx, false, lowerBounds, upperBounds);
							d[k2] = min1;
							s = min2;
							
							if (illc)
							{
								double szk = s + z[k2];
								s = d[k2] * szk*szk;
							}
							else
								s = sl - fx;
							if (df < s)
							{
								df = s;
								kl = k2;
							}
						}
						
						if (!illc && (df < std::abs(100.0 * EPSILON * fx)))
						{
							illc = true;
							gotoNext = true;
						}
						else
							gotoNext = false;
					} while (gotoNext);
					
					
					
					if ((k == 1) && (prin > 1))
						vecprint("\n... New Direction ...", d);
					/* minimize along conjugate directions */
					for (k2=0; k2<=k-1; k2++)
					{
						s = 0.0;
						
						min1 = d[k2];
						min2 = s;
						min(fun,k2, 2, fx, false, lowerBounds, upperBounds);
						d[k2] = min1;
						s = min2;
					}
					f1 = fx;
					fx = sf;
					lds = 0.0;
					for (i=0; i<dim; i++)
					{
						sl = x[i];
						x[i] = y[i];
						y[i] = sl - y[i];
						sl = y[i];
						lds = lds + sl*sl;
					}
					checkBounds(x, lowerBounds, upperBounds);
					
					lds = std::sqrt(lds);
					if (lds > small)
					{
						for (i=kl-1; i>=k; i--)
						{
							for (j=0; j < dim; j++)
								v[j][i+1] = v[j][i];
							d[i+1] = d[i];
						}
						d[k] = 0.0;
						for (i=0; i < dim; i++)
							v[i][k] = y[i] / lds;
						
						min1 = d[k];
						min2 = lds;
						min(fun,k, 4, f1, true, lowerBounds,upperBounds);
						d[k] = min1;
						lds = min2;
						
						if (lds <= 0.0)
						{
							lds = -lds;
							for (i=0; i< dim; i++)
								v[i][k] = -v[i][k];
						}
					}
					ldt = ldfac * ldt;
					if (ldt < lds)
						ldt = lds;
					if (prin > 1)
						print();
					/*if(monitor!=null) {
					 
					 monitor.newMinimum(fx,x);
					 }*/
					
					if(stopCondition(fx, x, tolfx, tolx, false))
					{
						kt++;
					}
					else
					{
						kt = 0;
					}
					if (kt > 1)
					{
						gotoFret = true;
						break;
					}
					
				}
				
				if (gotoFret) break;
				
				
				/*  try quadratic extrapolation in case    */
				/*  we are stuck in a curved valley        */
				quadr(fun, lowerBounds, upperBounds);
				dn = 0.0;
				for (i=0; i < dim; i++)
				{
					d[i] = 1.0 / std::sqrt(d[i]);
					if (dn < d[i])
						dn = d[i];
				}
				if (prin > 2)
					matprint("\n... New Matrix of Directions ...",v);
				for (j=0; j < dim; j++)
				{
					s = d[j] / dn;
					for (i=0; i < dim; i++)
						v[i][j] *= s;
				}
				if (scbd > 1.0)
				{       /* scale axis to reduce condition number */
					s = vlarge;
					for (i=0; i < dim; i++)
					{
						sl = 0.0;
						for (j=0; j < dim; j++)
							sl += v[i][j]*v[i][j];
						z[i] = std::sqrt(sl);
						if (z[i] < SQRT_SQRT_EPSILON)
							z[i] = SQRT_SQRT_EPSILON;
						if (s > z[i])
							s = z[i];
					}
					for (i=0; i < dim; i++)
					{
						sl = s / z[i];
						z[i] = 1.0 / sl;
						if (z[i] > scbd)
						{
							sl = 1.0 / scbd;
							z[i] = scbd;
						}
					}
				}
				for (i=1; i < dim; i++)
					for (j=0; j<=i-1; j++)
					{
						s = v[i][j];
						v[i][j] = v[j][i];
						v[j][i] = s;
					}
				minfit(dim, EPSILON, vsmall, v, d);
				if (scbd > 1.0)
				{
					for (i=0; i < dim; i++)
					{
						s = z[i];
						for (j=0; j < dim; j++)
							v[i][j] *= s;
					}
					for (i=0; i < dim; i++)
					{
						s = 0.0;
						for (j=0; j < dim; j++)
							s += v[j][i]*v[j][i];
						s = std::sqrt(s);
						d[i] *= s;
						s = 1.0 / s;
						for (j=0; j < dim; j++)
							v[j][i] *= s;
					}
				}
				for (i=0; i < dim; i++)
				{
					if ((dn * d[i]) > large)
						d[i] = vsmall;
					else if ((dn * d[i]) < small)
						d[i] = vlarge;
					else
						d[i] = Math.pow(dn * d[i],-2.0);
				}
				sort();               /* the new eigenvalues and eigenvectors */
				dmin = d[dim-1];
				if (dmin < small)
					dmin = small;
				illc = (SQRT_EPSILON * d[0]) > dmin;
				if ((prin > 2) && (scbd > 1.0))
					vecprint("\n... Scale Factors ...",z);
				if (prin > 2)
					vecprint("\n... Eigenvalues of A ...",d);
				if (prin > 2)
					matprint("\n... Eigenvectors of A ...",v);
				
				if ((maxFun > 0) && (nl > maxFun))
				{
					if (prin > 0)
						cout<<"\n... maximum number of function calls reached ..."<<endl;
					break;
				}
			}
			
			if (prin > 0)
			{
				vecprint("\n... Final solution is ...", x);
				cout<<"\n... Function value reduced to " << fx << " ..."<<endl;
				cout<<"... after " << numFun << " function calls."<<endl;
			}
			
			//return (fx);
		}
		
		
		//
		// Private stuff
		//
		
		// some global variables
		private int i, j, k, k2, nl, kl, kt;
		private double s, sl, dn, dmin,
		fx, f1, lds, ldt, sf, df,
		qf1, qd0, qd1, qa, qb, qc, small, vsmall, large,
		vlarge, ldfac, t2;
		
		// need to be initialised
		private vector<DOUBLE>  d;
		private vector<DOUBLE>  y;
		private vector<DOUBLE>  z;
		private vector<DOUBLE>  q0;
		private vector<DOUBLE>  q1;
		private vector< vector<DOUBLE> > v;
		
		private vector<DOUBLE>  tflin;
		
		private int dim;
		private vector<DOUBLE>  x;
		private MultiDimFunction& fun;
		
		// these will be set by praxis to the global control parameters
		private double h, t;
		
		// sort d and v in descending order
			
	template<typename DOUBLE> 
	class optimiser };
	
	
	{
		
	private:
		vector<DOUBLE> lowerBounds; /* Defines lower bound in each dimension*/
		vector<DOUBLE> upperBounds; /* Defines upper bound in each dimension*/
		
	public:
		~optimiser() {}
		optimiser() {
			lowerBounds.clear();
			upperBounds.clear();
		}
		
		virtual DOUBLE evaluate(const vector<DOUBLE>& x) = 0; //Evaluate the function at a point. Assumed positive (???).Must be overloaded.
		
		/**
		 Return current dimension. This will be 0 if the optimiser has not been initialised.
		 **/
		uint getDimension() {
			return lowerBounds.size();
		}
		
		/**
		 Set upper and lower bounds. Checks that dimensions are ok and that lower bounds do not exceed upper bounds.
		 @param lower vector of lower bounds
		 @param upper vector of upper bounds.
		 */
		void setBounds(const vector<DOUBLE>& lower, const vector<DOUBLE>& upper) {
			if (lower.size()!=upper.size())
				throw PhylibException("Lower bounds and upper bounds have different dimensions");
			
			uint n = lower.size();
			
			lowerBounds.resize(n);
			std::copy(lower.begin(),lower.end(),lowerBounds.begin());
			
			upperBounds.resize(n);
			std::copy(upper.begin(),upper.end(),upperBounds.begin());
			
			for(int i=0;i<n;i++)
				if (lowerBounds[i]>upperBounds[i])
					throw PhylibException("Lower bound exceeds upper bounds: infeasible region");
		}
		
		/**
		 Set all of the dimensions to the same lower and upper bound. Throws exception if lower > upper
		 */
		void setBounds(uint n, DOUBLE lower, DOUBLE upper) {
			
			if (lower>upper)
				throw PhylibException("Lower bound exceeds upper bounds: infeasible region");
			lowerBounds.resize(n);
			upperBounds.resize(n);
			std::fill(lowerBounds.begin(),lowerBounds.end(),lower);
			std::fill(upperBounds.begin(),upperBounds.end(),upper);
		}
		
		
		
		/**
		 Optimise using a simple simulated annealing algorithm. We use a simple step proposal 
		 distribution: x' = x + u where u is uniform on [-stepSize,stepSize]^n. 
		 The annealing schedule is T(n) = (coolingRate)^n  so we definitely need 0<coolingRate<1
		 
		 We follow the algorithm on page 8-15 of Nicholls, Fox, Tan, where evaluate(x) returns q(x).
		 
		 Assumes evaluate(x) > 0 for all x. 
		 
		 @param x State. Used as initial state (assumed feasible), and to return final state.
		 @param stepSize. Step size for proposal function
		 @param coolingRate. Parameter for annealing schedule.
		 @returns value at minimum (of q(x) := evaluate(x) ).
		 **/
		DOUBLE simulatedAnnealing(vector<DOUBLE>& x, double stepSize, double coolingRate, int numSteps) { //Minimise by simulated annealing
			
			//TODO: 
			//Improved stopping conditions.
			DOUBLE qx = evaluate(x);
			uint n = x.size();
			vector<DOUBLE> xprime(n);
			DOUBLE Tn = 1.0;
			const uint nAcceptsPerTempChange = 10;
			
			
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
					DOUBLE qxprime = evaluate(xprime);
					DOUBLE alpha = std::exp(-(qxprime - qx)/Tn);
					
					cerr<<"Trying x'=[";
					for(int j=0;j<n;j++) {
						cerr<<xprime[j];
						if (j<n-1)
							cerr<<", ";
					}
					cerr<<"]\tqx = "<<qxprime<<"\t"<<endl;
					
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
		
		
	private:		
		/*****
		 THE FOLLOWING TWO ROUTINES ARE ADAPTED FROM NUMERICAL RECIPES in C.
		 ***********/
		
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
		
		/**
		 Evaluate the function at x+alpha(s
		 **/
		DOUBLE evaluateLinesearch(const vector<DOUBLE>& x, const vector<DOUBLE>& s, DOUBLE alpha) {
			vector<DOUBLE> y(x.size());
			for(uint i=0;i<x.size();i++)
				y[i] = x[i] + alpha*s[i];
			return evaluate(y);
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
		
		bool bracketMinimum(const vector<DOUBLE>& x, const vector<DOUBLE>& s, DOUBLE& a, DOUBLE& c, DOUBLE& fa, DOUBLE& fc)
		{
			const DOUBLE GOLD =   (1.0 + std::sqrt(5))/2.0;;
			const DOUBLE TINY = 10e-20;
			const DOUBLE GLIMIT = 100.0;
			DOUBLE b = c;
			
			DOUBLE fb; 
			DOUBLE ulim,u,r,q,fu;
			DOUBLE tmp;
			
			fa=evaluateLinesearch(x, s, a); //fa = f(a)
			fb=evaluateLinesearch(x, s, b); //fb = f(b)
			
			if (fb > fa) {
				tmp = a; a=b; b=tmp;
				tmp = fa; fa = fb; fb = tmp;
			}
			
			c=(b)+GOLD*(b-a); //Initial guess for c (outside the interval)
			fc=evaluateLinesearch(x, s, c); 
			while (fb > fc) {
				
				
				r=(b-a)*(fb-fc);
				q=(b-c)*(fb-fa);
				u=(b)-((b-c)*q-(b-a)*r)/(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
				ulim=(b)+GLIMIT*(c-b);
				if ((b-u)*(u-c) > 0.0) {
					fu=evaluateLinesearch(x, s, u); //fu = f(u)
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
					fu=evaluateLinesearch(x, s, u); //fu = f(u)
				} else if ((c-u)*(u-ulim) > 0.0) {
					fu=evaluateLinesearch(x, s, u); //fu = f(u)
					if (fu < fc) {
						SHFT(b,c,u,c+GOLD*(c-b))
						SHFT(fb,fc,fu,evaluateLinesearch(x, s, u))
					}
				} else if ((u-ulim)*(ulim-c) >= 0.0) {
					u=ulim;
					fu=evaluateLinesearch(x, s, u); //fu = f(u)
				} else {
					u=(c)+GOLD*(c-b);
					fu=evaluateLinesearch(x, s, u); //fu = f(u)
				}
				SHFT(a,b,c,u)
				SHFT(fa,fb,fc,fu)
			}
			return false;
		}
		
		DOUBLE brent(const vector<DOUBLE>& X, const vector<DOUBLE>& S, DOUBLE a, DOUBLE c, DOUBLE fa, DOUBLE fc, DOUBLE tol, DOUBLE& xmin)
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
			fw=fv=fx=evaluateLinesearch(X, S, x);
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
				fu=evaluateLinesearch(X, S, u); //fu = f(u)
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
		
		
		bool linesearch(const vector<DOUBLE>& x, const vector<DOUBLE>& s, DOUBLE& alpha, DOUBLE& fmin) {
			
			const DOUBLE TOL = 2.0e-4;
			DOUBLE a = 0.0;
			DOUBLE c = 1.0; //Initial interval. 
			
			//			DOUBLE alphaMin, alphaMax;			
			
			DOUBLE fa,fc;
			
			//FInd a,b,c so that f(a)>
			if (!bracketMinimum(x, s, a, c, fa, fc)) {
				alpha = 0;
				fmin = fa;
				return false;
			}
			
			fmin = brent(x, s, a, c, fa, fc, TOL, alpha);	
			return 	true;
		};
		
	public:
		DOUBLE brentPowell(vector<DOUBLE>& x) {
			
			
			
		
		
	};
	
	
	
	
	
#endif

