#ifndef RUN_H
#define RUN_H

void run(confObj_t simParam, int size, int thisRank);
void save_polyspectrum(confObj_t simParam, int num, double *theta, double *k, double *polyspectrum);
void save_polyspectrum_numpolygons(confObj_t simParam, int num, double *theta, double *k, double *polyspectrum, double *numpolygons);

#endif
