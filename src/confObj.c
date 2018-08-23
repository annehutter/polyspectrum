/*
 *  confObj.c
 *  uvff
 *
 *  Created by 
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

/*--- Includes ----------------------------------------------------------*/
#include "confObj.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <inttypes.h>
#include <stdbool.h>
#include "xmem.h"

/*--- Defines for the Ini structure -------------------------------------*/


/*--- Prototypes of local functions -------------------------------------*/


/*--- Implementations of exported functios ------------------------------*/
extern confObj_t
confObj_new(parse_ini_t ini)
{
    confObj_t config;
    assert(ini != NULL);
    
    config = xmalloc(sizeof(struct confObj_struct));
    
    //reading mandatory stuff

    //Input grid
    getFromIni(&(config->grid_size), parse_ini_get_int32,
               ini, "gridsize", "Grid");
    getFromIni(&(config->box_size), parse_ini_get_double,
               ini, "boxsize", "Grid");
    getFromIni(&(config->gas_inputs_in_dp), parse_ini_get_int32,
               ini, "gasInputsInDoublePrecision", "Grid");
    getFromIni(&(config->ion_inputs_in_dp), parse_ini_get_int32,
               ini, "ionInputsInDoublePrecision", "Grid");
    getFromIni(&(config->density_file), parse_ini_get_string,
               ini, "densityFile", "Grid");
    getFromIni(&(config->ion_file), parse_ini_get_string,
               ini, "ionFile", "Grid");  
    
    //Output
    getFromIni(&(config->output_dir), parse_ini_get_string,
               ini, "output_dir", "Output");
    getFromIni(&(config->output_basename), parse_ini_get_string,
               ini, "output_basename", "Output");
    
    //Cosmology
    getFromIni(&(config->h), parse_ini_get_double,
               ini, "hubble_h", "Cosmology");
    getFromIni(&(config->omega_b), parse_ini_get_double,
               ini, "omega_b", "Cosmology");
    getFromIni(&(config->omega_m), parse_ini_get_double,
               ini, "omega_m", "Cosmology");
    getFromIni(&(config->omega_l), parse_ini_get_double,
               ini, "omega_l", "Cosmology");
    getFromIni(&(config->sigma8), parse_ini_get_double,
               ini, "sigma8", "Cosmology");
    getFromIni(&(config->Y), parse_ini_get_double,
               ini, "Y", "Cosmology");

    //Polyspectrum
    getFromIni(&(config->n), parse_ini_get_int32,
               ini, "n", "Polyspectrum");
    getFromIni(&(config->k1), parse_ini_get_double,
               ini, "k1", "Polyspectrum");
    getFromIni(&(config->k2), parse_ini_get_double,
               ini, "k2", "Polyspectrum");
    getFromIni(&(config->sweep_through_theta), parse_ini_get_int32,
               ini, "sweepThroughTheta", "Polyspectrum");
    getFromIni(&(config->num_values), parse_ini_get_int32,
               ini, "numValues", "Polyspectrum");
    getFromIni(&(config->theta), parse_ini_get_double,
               ini, "theta", "Polyspectrum");
    
    return config;
}

extern void
confObj_del(confObj_t *config)
{
    assert(config != NULL);
    assert(*config != NULL);
    
    //Input
    xfree((*config)->density_file);
    xfree((*config)->ion_file);
    
    xfree((*config)->output_dir);
    xfree((*config)->output_basename);

    xfree(*config);
    *config = NULL;
}

extern confObj_t
readConfObj(char *fileName) {
    confObj_t newConfObj;
    parse_ini_t ini;
    
    ini = parse_ini_open(fileName);
    if (ini == NULL) {
        fprintf(stderr, "FATAL:  Could not open %s for reading.\n",
                fileName);
        exit(EXIT_FAILURE);
    }
    
    newConfObj = confObj_new(ini);
    
    parse_ini_close(&ini);
    
    return newConfObj;
}

/*--- Implementations of local functions --------------------------------*/
