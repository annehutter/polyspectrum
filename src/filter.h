#ifndef FILTER_H
#define FILTER_H

fftw_complex *generate_kfilter(int nbins, int local_n0, int local_n0_start, double k, double boxsize);
void construct_kfilter(int nbins, int local_n0, int local_n0_start, fftw_complex *array, double k, double boxsize);

#endif
