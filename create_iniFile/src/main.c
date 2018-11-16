#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXLENGTH 1024

int main (int argc, /*const*/ char * argv[]) 
{ 
  char filename[MAXLENGTH];
  FILE *file = NULL;
  
  int gridsize = 0;
  float boxsize = 0.;
  int gasInputsInDoublePrecision = 0;
  int ionInputsInDoublePrecision = 0;
  char densityFile[MAXLENGTH];
  char ionFile[MAXLENGTH];
  
  float hubble_h = 0.7;
  float omega_b = 0.0456;
  float omega_m = 0.27;
  float omega_l = 0.73;
  float sigma8 = 0.82;
  float Y = 0.24;
  
  char whichField[MAXLENGTH];
  int n = 0;
  int equilateral = 0;
  float k1 = 0.;
  float k2 = 0.;
  int numValues = 1;
  float theta = 0.;
  
  char output_dir[MAXLENGTH];
  char output_basename[MAXLENGTH];
  
  for(int i=0; i<MAXLENGTH; i++)
  {
    densityFile[i] = '\0';
    ionFile[i] = '\0';
    output_dir[i] = '\0';
    output_basename[i] = '\0';
  }
  
  if(argc != 22)
  {
    printf("create_inifile: (C) Use at own risk...\nUSAGE: create_inifile <FILE>\n");
  }
  else
  {
    strcpy(filename, argv[1]);
    gridsize = atoi(argv[2]);
    boxsize = atof(argv[3]);
    gasInputsInDoublePrecision = atoi(argv[4]);
    ionInputsInDoublePrecision = atoi(argv[5]);
    strcpy(densityFile, argv[6]);
    strcpy(ionFile, argv[7]);
    
    hubble_h = atof(argv[8]);
    omega_b = atof(argv[9]);
    omega_m = atof(argv[10]);
    omega_l = atof(argv[11]);
    sigma8 = atof(argv[12]);
    Y = atof(argv[13]);
    
    strcpy(whichField, argv[14]);
    n = atoi(argv[15]);
    equilateral = atoi(argv[16]);
    k1 = atof(argv[17]);
    k2 = atof(argv[18]);
    numValues = atoi(argv[19]);
    theta = atof(argv[20]);
    
    strcpy(output_dir, argv[21]);
    strcpy(output_basename, argv[22]);
  }
  
  file = fopen(filename, "wt");
  
  fprintf(file, "[Grid]\n");
  fprintf(file, "gridsize = %d\n", gridsize);
  fprintf(file, "boxsize = %f\t# in h^-1 Mpc\n", boxsize);
  fprintf(file, "gasInputsInDoublePrecision = %d\n", gasInputsInDoublePrecision);
  fprintf(file, "ionInputsInDoublePrecision = %d\n", ionInputsInDoublePrecision);
  fprintf(file, "densityFile = %s\n", densityFile);
  fprintf(file, "ionFile = %s\n", ionFile);
  
  fprintf(file, "\n[Cosmology]\n");
  fprintf(file, "hubble_h = %f\n", hubble_h);
  fprintf(file, "omega_b = %f\n", omega_b);
  fprintf(file, "omega_m = %f\n", omega_m);
  fprintf(file, "omega_l = %f\n", omega_l);
  fprintf(file, "sigma8 = %f\n", sigma8);
  fprintf(file, "Y = %f\n", Y);
  
  fprintf(file, "\n[Polyspectrum]\n");
  fprintf(file, "whichField = %s\n", whichField);
  fprintf(file, "# Possibilities: DENS, XHII, XHI_DENS\n");
  fprintf(file, "n = %d\n", n);
  fprintf(file, "equilateral = %d\n", equilateral);
  fprintf(file, "k1 = %f\t# in h Mpc^-1\n", k1);
  fprintf(file, "k2 = %f\t# in h Mpc^-1\n", k2);
  fprintf(file, "numValues = %d\n", numValues);
  fprintf(file, "theta = %f\n", theta);
  
  fprintf(file, "\n[Output]\n");
  fprintf(file, "output_dir = %s\n", output_dir);
  fprintf(file, "output_basename = %s\n", output_basename);
  
  fclose(file);
  
  return 0;
}

