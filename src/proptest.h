double ksstat(double *x, int *n);
double cmstat(double *x, int *firstlast, double *dsigma);
double adstat(double *x, int *firstlast, double *dsigma, double *sigma);

void lwy_scoreproptest(double *u, int *nev, int *nvar, double *imatrow, double *dsigma, double *sigma,
		       int *firstlast, int *nsim, int *tested, double *statks, double *statcm,
		       double *statad, double *pks, double *pcm, double *pad, int *nsim_plot,
		       double *scoreprocess_sim);


void scorestat(double **u, double *v, int *k, double *stat, double *work);

void lwy_smoothproptest(double *u, int *nev, int *nvar, double *djinv, double *vinv, int *alt, int *nalt,
		  int *nsim, int *tested, double *stat_d, double *stat_bic, double *logn,
		  double *pval_d_sim, double *pval_bic_sim, double *wts, double *pval_bic_w);

