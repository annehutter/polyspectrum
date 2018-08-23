#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>

/* ------------------------------------------------------------*/
/* RANDOM NUMBERS                                              */
/* ------------------------------------------------------------*/
double randdouble()
{
    return rand()/((double)(RAND_MAX)+1);
} 

/* ------------------------------------------------------------*/
/* STRING ADDITION                                             */
/* ------------------------------------------------------------*/
char* concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1) + strlen(s2) + 1);//+1 for the null-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    
    return result;
}

char* concat2(const char *s1, const char *s2, const char *s3)
{
    char *result = malloc(strlen(s1) + strlen(s2) + strlen(s3) + 1);//+1 for the null-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    strcat(result, s3);
    
    return result;
}

/* ------------------------------------------------------------*/
/* ARRAY INITIALIZATION                                        */
/* ------------------------------------------------------------*/
void initialize_array_int(int nbins, int *array, int value)
{
    for(int i=0; i<nbins; i++)
    {
        array[i] = value;
    }
}

void initialize_array_long_int(int nbins, long int *array, int value)
{
    for(int i=0; i<nbins; i++)
    {
        array[i] = value;
    }
}

void initialize_array_float(int nbins, float *array, float value)
{
    for(int i=0; i<nbins; i++)
    {
        array[i] = value;
    }
}

void initialize_array_double(int nbins, double *array, double value)
{
    for(int i=0; i<nbins; i++)
    {
        array[i] = value;
    }
}

/* ------------------------------------------------------------*/
/* ARRAY ALLOCATION                                            */
/* ------------------------------------------------------------*/

int* allocate_array_int(int length, char *arrayname)
{
    int *tmp = malloc(length * sizeof(int));
    if(tmp == NULL)
    {
        fprintf(stderr, "allocate_array_int: could not allocate array %s\n", arrayname);
        exit(EXIT_FAILURE);
    }
    initialize_array_int(length, tmp, 0.);
    return tmp;
}

long int* allocate_array_long_int(int length, char *arrayname)
{
    long int *tmp = malloc(length * sizeof(long int));
    if(tmp == NULL)
    {
        fprintf(stderr, "allocate_array_int: could not allocate array %s\n", arrayname);
        exit(EXIT_FAILURE);
    }
    initialize_array_long_int(length, tmp, 0.);
    return tmp;
}

float* allocate_array_float(int length, char *arrayname)
{
    float *tmp = malloc(length * sizeof(float));
    if(tmp == NULL)
    {
        fprintf(stderr, "allocate_array_float: could not allocate array %s\n", arrayname);
        exit(EXIT_FAILURE);
    }
    initialize_array_float(length, tmp, 0.);
    return tmp;
}

double* allocate_array_double(int length, char *arrayname)
{
    double *tmp = malloc(length * sizeof(double));
    if(tmp == NULL)
    {
        fprintf(stderr, "allocate_array_double: could not allocate array %s\n", arrayname);
        exit(EXIT_FAILURE);
    }
    initialize_array_double(length, tmp, 0.);
    return tmp;
}

/* ------------------------------------------------------------*/
/* FILE EXISTENCE                                              */
/* ------------------------------------------------------------*/

int file_exist(char *filename)
{
    FILE *file;
    if((file = fopen (filename, "rt"))){
        fclose(file);
        return 1;
    }else{
        return 0;
    }
}
