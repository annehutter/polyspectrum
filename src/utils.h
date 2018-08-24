#ifndef UTILS_H
#define UTILS_H

#ifndef MAXLENGTH
#define MAXLENGTH 1024
#endif

/* ------------------------------------------------------------*/
/* RANDOM NUMBERS                                              */
/* ------------------------------------------------------------*/
double randdouble();

/* ------------------------------------------------------------*/
/* STRING ADDITION                                             */
/* ------------------------------------------------------------*/
char* concat(const char *s1, const char *s2);
char* concat2(const char *s1, const char *s2, const char *s3);
char* concat3(const char *s1, const char *s2, const char *s3, const char *s4);

/* ------------------------------------------------------------*/
/* ARRAY INITIALIZATION                                        */
/* ------------------------------------------------------------*/
void initialize_array_int(int nbins, int *array, int value);
void initialize_array_long_int(int nbins, long int *array, int value);
void initialize_array_float(int nbins, float *array, float value);
void initialize_array_double(int nbins, double *array, double value);

/* ------------------------------------------------------------*/
/* ARRAY ALLOCATION                                            */
/* ------------------------------------------------------------*/
int* allocate_array_int(int length, char *arrayname);
long int* allocate_array_long_int(int length, char *arrayname);
float* allocate_array_float(int length, char *arrayname);
double* allocate_array_double(int length, char *arrayname);

/* ------------------------------------------------------------*/
/* FILE EXISTENCE                                              */
/* ------------------------------------------------------------*/
int file_exist(char *filename);

#endif
