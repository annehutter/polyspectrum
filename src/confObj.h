/*
 *  confObj.h
 *  uvff
 *
 *  Created by 
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CONFOBJ_H
#define CONFOBJ_H

/*--- Includes ----------------------------------------------------------*/
#include "parse_ini.h"
#include <stdint.h>
#include <stdbool.h>


/*--- ADT handle --------------------------------------------------------*/
typedef struct confObj_struct *confObj_t;


/*--- Implemention of main structure ------------------------------------*/
struct confObj_struct {
    //Input grid
    int            grid_size;
    double         box_size;
    int            gas_inputs_in_dp;
    int            ion_inputs_in_dp;
    char           *density_file;
    char           *ion_file;
    
    //Cosmology
    double         h;
    double         omega_b;
    double         omega_m;
    double         omega_l;
    double         sigma8;
    double         Y;
    
    //Polyspectrum 
    char           *which_field;
    int            n;
    int            equilateral;
    double         k1;
    double         k2;
    int            num_values;
    double         theta;
    double         kbinwidth;
    
    //Output
    char           *output_dir;
    char           *output_basename;
    int            write_numpolygons;
};


/*--- Prototypes of exported functions ----------------------------------*/
extern confObj_t
readConfObj(char *fileName);

extern confObj_t
confObj_new(parse_ini_t ini);

extern void
confObj_del(confObj_t *config);


#endif
