#ifndef KVECTORS_H
#define KVECTORS_H

/* structure for kvectors: This strcuture defines the order of the polyspectrum, the k-vectors */

typedef struct
{
    int n;
    double *kpolygon;
    
    int numValues;
    double *theta;
    double *k;
} kvectors_t;

/* functions */

kvectors_t *initKvectors();
kvectors_t *read_params_to_kvectors(confObj_t simParam);

void deallocate_kvectors(kvectors_t *theseKvectors);

double *generate_cosTheta_values(int numValues);
double *generate_k_values_bispectrum(int numValues, double k1, double k2);
double *generate_k_values_powerspectrum(int grid_size, double box_size);

#endif
