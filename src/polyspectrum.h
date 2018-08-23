#ifndef POLYSPECTRUM_H
#define POLYSPECTRUM_H

fftw_complex *generate_values_polygons(int nbins, int local_n0, fftw_complex *kfilter, fftw_complex *fft_array);
fftw_complex *generate_num_polygons(int nbins, fftw_complex *kfilter);

double polyspectrum(int nbins, int local_n0, int local_n0_start, fftw_complex *fft_array, int n, double *k, double boxsize);

#endif
