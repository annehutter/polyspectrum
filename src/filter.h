#ifndef FILTER_H
#define FILTER_H

/*------------------------------------------------------------------------*/
/* FILTER WITH MINIMAL KBINWIDTH */
/*------------------------------------------------------------------------*/
fftw_complex *generate_kfilter(int nbins, int local_n0, int local_n0_start, double k, double binwidth, double boxsize);
void construct_kfilter(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double k, double binwidth, double boxsize);

/*------------------------------------------------------------------------*/
/* FILTER WITH DETERMINING MINIMUM & MAXIMUM KBINWIDTH */
/*------------------------------------------------------------------------*/
fftw_complex *generate_kfilter_with_determining_kbinwidth(int nbins, int local_n0, int local_n0_start, double *kmin, double *kmax, double k, double binwidth, double boxsize);
void construct_kfilter_with_determining_kbinwidth(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double *kmin, double *kmax, double k, double binwidth, double boxsize);

/*------------------------------------------------------------------------*/
/* FILTER WITH IMPOSING MINIMUM & MAXIMUM K-VALUES */
/*------------------------------------------------------------------------*/
fftw_complex *generate_kfilter_with_kbinwidth(int nbins, int local_n0, int local_n0_start, double kmin, double kmax, double binwidth, double boxsize);
void construct_kfilter_with_kbinwidth(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double kmin, double kmax, double binwidth, double boxsize);

/*------------------------------------------------------------------------*/
/* FILTERS AS USED IN WATKINSON ET AL. 2017, 2018 */
/*------------------------------------------------------------------------*/
fftw_complex *generate_kfilter_Watkinson(int nbins, int local_n0, int local_n0_start, double k, double binwidth, double boxsize);
void construct_kfilter_Watkinson(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double k, double binwidth, double boxsize);
fftw_complex *generate_kfilter_Watkinson_kn(int nbins, int local_n0, int local_n0_start, double kmin, double kmax, double binwidth, double boxsize);
void construct_kfilter_Watkinson_kn(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double kmin, double kmax, double binwidth, double boxsize);

/*------------------------------------------------------------------------*/
/* AUXILIARY FUNCTIONS */
/*------------------------------------------------------------------------*/
void get_corners(int mx, int my, int mz, double binwidth, double *mmin, double *mmax);
void calc_k3_min_and_max(double *k, double **kmin, double **kmax);
void calc_kn_min_and_max(int n, double *k, double **kmin, double **kmax);

#endif
