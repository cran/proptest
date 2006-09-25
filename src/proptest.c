#include <R.h>
#include <Rmath.h>
#include "proptest.h"

/*
	SMOOTH TEST OF THE PROPORTIONAL HAZARDS ASSUMPTION
	called by smoothproptest
*/

void lwy_smoothproptest(double *u, int *nev, int *nvar, double *djinv, double *vinv, int *alt, int *nalt,
		  int *nsim, int *tested, double *stat_d, double *stat_bic, double *logn,
		  double *pval_d_sim, double *pval_bic_sim, double *wts, double *pval_bic_w)
{
/*	tested: in R (indexed from 1) points to the start of the tested cov. (index of tested itself)
	in C: points to the start of the tested smooth functions (index of tested smooth)
	alt(nalt,d)
*/
/*
	Output:
	pval_d_sim		simulated p-value for fixed dimension test
	pval_bic_sim	simul. p-val for data driven test
	pval_w_sim		p-value obtained from mixture of chisq with simul. weigths (doesn't work)
	wt		those^ weights
*/ 
	double g;
	int d=*nvar-*tested;
	double *ugend = (double *) R_alloc(*nvar,sizeof(double));
	double *u2gend = (double *) R_alloc(d,sizeof(double));
	double *scoreg = (double *) R_alloc(*nalt,sizeof(double));
	double *p_vinv;
	double *work = (double *) R_alloc(d,sizeof(double));
	double ***u2gendalt = (double ***) R_alloc(*nalt,sizeof(double **)); /* scorevector (array of arrays of pointers) */
	int s;
	int *dimalt = (int *) R_alloc(*nalt,sizeof(int));
	int i,j,b;
	GetRNGstate();
	for (i=0; i<*nalt; i++) {
		dimalt[i]=0;
		for (j=0; j<d; j++) {
			dimalt[i] += alt[i+j**nalt];
		}
		u2gendalt[i] = (double **) R_alloc(dimalt[i],sizeof(double *));
		b=0;
		for (j=0; j<d; j++) {
			if (alt[i+j**nalt]>0) {
				u2gendalt[i][b] = u2gend + j;
				++b;
			}
		}
	}
	for (i=0; i<d; i++) {
		wts[i]=0.;
	}
	*pval_d_sim=0.;
	*pval_bic_sim=0.;
	*pval_bic_w=0.;
	for (b=0; b<*nsim; b++) {
		g = norm_rand();
		for (j=0; j<*nvar; j++) {
			ugend[j]=u[j**nev+0]*g;
		}
		for (i=1; i<*nev; i++) {
			g = norm_rand();
			for (j=0; j<*nvar; j++) {
				ugend[j] += u[j**nev+i]*g;
			}
		}
		for (j=0; j<d; j++) {
			u2gend[j] = ugend[*tested+j];
			for (i=0; i<*tested; i++) {
				u2gend[j] -= djinv[i*d+j]*ugend[i];
			}
		}
		p_vinv = vinv;
		s=0;
		for (i=0; i<*nalt; i++) {
			scorestat(u2gendalt[i],p_vinv,dimalt+i,scoreg+i,work);
			scoreg[i] -= dimalt[i]**logn;
			if (scoreg[i]>scoreg[s]) s=i;
			p_vinv += dimalt[i]*dimalt[i];
		}
		*pval_bic_sim += (*stat_bic<=(scoreg[s]+dimalt[s]**logn));
		wts[dimalt[s]-1] += 1.;
		*pval_d_sim += (*stat_d<=(scoreg[0]+d**logn));
	}
	*pval_bic_sim /= *nsim;
	*pval_d_sim /= *nsim;
	for (i=0; i<d; i++) {
		wts[i] /= *nsim;
		*pval_bic_w += wts[i]*pchisq(*stat_bic,i+1,0,0);
	}
	
	PutRNGstate();
}

void scorestat(double **u, double *v, int *k, double *stat, double *work)
{
/*
	computes quadratic score statistic u^T*v*u; result in *stat
*/
	int i,j;
	for (i=0; i<*k; i++) {
		work[i]=0.;
		for (j=0; j<*k; j++) {
			work[i] += v[j**k+i]**(u[j]); /* v[i,j]*u[j] */
		}
	}
	*stat =0.;
	for (i=0; i<*k; i++) {
		*stat += *(u[i])*work[i];
	}
}

void h_approx_sim(double *x, double *vinv1, double *vinv2, double *cholv, int *d,
		  double *logn, int *nsim, double *h)
/*
	two term H-approximation of the distribution function of T_S (with BIC with all subsets)
*/
{
	double *norm_iid = (double *) R_alloc(*d,sizeof(double));
	double *u2 = (double *) R_alloc(*d,sizeof(double));
	int nalt2=*d*(*d-1)/2;
	double ***u2alt = (double ***) R_alloc(nalt2,sizeof(double **)); /* array of arrays of pointers; for each alt array of pointers to u2 */
	int i,j,k;
	double m1,m2; /* max of one-dimensional and two-dimensional score statistics, respectively */
	double t;
	int *p_two = (int *) R_alloc(1,sizeof(int));
	*p_two=2;

	GetRNGstate();
	
	i=0;
	for (j=0; j<*d-1; j++) {
		for (k=j+1; k<*d; k++) {
			u2alt[i] = (double **) R_alloc(2,sizeof(double *));
			u2alt[i][0] = u2+j;
			u2alt[i][1] = u2+k;
			++i;
		}
	}
	
	*h=0.;
	for (i=0; i<*nsim; i++) {
		m1=0.;
		m2=0.;
		for (j=0; j<*d; j++) {
			norm_iid[j]=norm_rand();
			u2[j]=0.;
		}
		/* u2 = t(cholv) %*% norm_iid */
		for (j=0; j<*d; j++) {
			for (k=0; k<*d; k++) {
				u2[j] += cholv[j**d+k]*norm_iid[k];    /* u2[j] = cholv[k,j]*norm_iid[k] */
			}
		}
		for (j=0; j<*d; j++) {
			t = u2[j]*vinv1[j]*u2[j];
			if (t > m1) m1=t;
		}
		for (j=0; j<nalt2; j++) {
			scorestat(u2alt[j],vinv2+j*4,p_two,&t,norm_iid); /* using norm_iid as work array */
			if (t > m2) m2=t;
		}
		*h += ((m1<=*x)&&(m2-m1<=*logn)) + ((m2<=*x)&&(m2-m1>*logn));
	}
	*h /= (*nsim);

	PutRNGstate();
}

/*
	TEST OF THE PH ASSUMPTION BASED ON THE SCORE PROCESS (KS, CM, AD type)
	called from scoreproptest
*/

void lwy_scoreproptest(double *u, int *nev, int *nvar, double *imatrow, double *dsigma, double *sigma,
		       int *firstlast, int *nsim, int *tested, double *statks, double *statcm,
		       double *statad, double *pks, double *pcm, double *pad, int *nsim_plot,
		       double *scoreprocess_sim)
{
	double g; /* = (double *) R_alloc(*nev,sizeof(double)); */
	double *ug = (double *) R_alloc(*nev**nvar,sizeof(double));
	double *ugtested = ug+(--*tested)**nev;
	int i,j,b;
	
	GetRNGstate();
	
	--firstlast[0];
	--firstlast[1];
	*pks=0.;
	*pcm=0.;
	*pad=0.;
	for (b=0; b<*nsim; b++) {
		g = norm_rand();
		for (j=0; j<*nvar; j++) {
			ug[j**nev+0]=u[j**nev+0]*g; /* *l[0] */
		}
		for (i=1; i<*nev; i++) {
			g = norm_rand();
			for (j=0; j<*nvar; j++) {
				ug[j**nev+i]=ug[j**nev+i-1]+u[j**nev+i]*g; /* *l[i]; */
			}
		}
		for (i=0; i<*nev; i++) {
			for (j=0; j<*nvar; j++) {
				ugtested[i] -= imatrow[j**nev+i]*ug[j**nev+*nev-1];
			}
		}
		
		if (b<*nsim_plot) {
			for (i=0; i<*nev; i++) {
				scoreprocess_sim[i+b**nev] = ugtested[i]; /* realisations for plotting */
			}
		}
		
		*pks += (*statks<=ksstat(ugtested,nev));
		*pcm += (*statcm<=cmstat(ugtested,firstlast,dsigma));
		*pad += (*statad<=adstat(ugtested,firstlast,dsigma,sigma));
	}
	*pks /= *nsim;
	*pcm /= *nsim;
	*pad /= *nsim;
	
	PutRNGstate();
}

double ksstat(double *x, int *n)
{
	double y=0.;
	int i;
	for (i=0; i<*n; i++) {
		if (fabs(x[i])>y) {
			y=fabs(x[i]);
		}
	}
	return y;
}

double cmstat(double *x, int *firstlast, double *dsigma)
{
	double y=0.;
	int i;
	for (i=firstlast[0]; i<=firstlast[1]; i++) {
		y += x[i]*x[i]*dsigma[i];
	}
	return y;
}

double adstat(double *x, int *firstlast, double *dsigma, double *sigma)
{
	double y=0.;
	int i;
	for (i=firstlast[0]; i<firstlast[1]; i++) {
		y += x[i]*x[i]*dsigma[i]/sigma[i]/(1.-sigma[i]);
	}
	return y;
}


