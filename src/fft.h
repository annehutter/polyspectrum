#ifndef FFT_H
#define FFT_H

void fft_real_to_kspace(int nbins, fftw_complex *input, fftw_complex *output);
void fft_k_to_realspace(int nbins, fftw_complex *input, fftw_complex *output);

#endif
