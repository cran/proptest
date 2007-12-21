double ksstat(double *x, int *n);
double cmstat(double *x, int *firstlast, double *dsigma);
double adstat(double *x, int *firstlast, double *dsigma, double *sigma);

void lwy_scoreproptest(double *u, int *nev, int *nvar, double *l, double *imatrow, double *dsigma,
		       double *sigma, int *firstlast, int *nsim, int *tested, double *statks,
		       double *statcm, double *statad, double *pks, double *pcm, double *pad,
		       int *nsim_plot, double *scoreprocess_sim);

void lwy_scoreproptest_global(double *du, int *nev, int *nvar, double *l, double *imatimatendinv,
			      double *sigmaend, int *nsim, double *stat, double *p, int *nsim_plot,
			      double *testprocess_sim);

void scorestat(double **u, double *v, int *k, double *stat, double *work);

void lwy_smoothproptest(double *u, int *nev, int *nvar, double *djinv, double *vinv, int *alt, int *nalt,
		  int *nsim, int *tested, double *stat_d, double *stat_bic, double *logn,
		  double *pval_d_sim, double *pval_bic_sim, double *wts, double *pval_bic_w);

void h_approx_sim(double *x, double *vinv1, double *vinv2, double *cholv, int *d,
		  double *logn, int *nsim, double *h, double *hmax);

