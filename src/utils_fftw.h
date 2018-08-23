#ifndef UTILS_FFTW
#define UTILS_FFTW

fftw_complex *allocate_3D_array_fftw_complex(int nbins);

void product_3D_fftw_arrays(int nbins, int local_n0, fftw_complex *array1, fftw_complex *array2, fftw_complex *output);
void product3_3D_fftw_arrays(int nbins, int local_n0, fftw_complex *array1, fftw_complex *array2, fftw_complex *array3, fftw_complex *output);

fftw_complex sum_3D_fftw_array(int nbins, int local_n0, fftw_complex *array);

#endif
