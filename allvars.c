#include "allvars.h"


struct io_header_1 header1, header;

int WhichSpectrum;


int SphereMode;
int *Local_nx_table;

FILE *FdTmp, *FdTmpInput;

int Nmesh, Nsample;

long long IDStart;



char GlassFile[500];
char FileWithInputSpectrum[500];
char     FileWithDelta[500];
char     FileWithDeltaDot[500];
int TileFac;

double Box;
int Seed;

long long TotNumPart;

int NumPart;

int *Slab_to_task;

int NTaskWithN;

struct part_data *P;

int Nglass;

double InitTime, InputTime;
double Redshift;
double MassTable[6];


char OutputDir[100], FileBase[100];
int NumFilesWrittenInParallel;


int ThisTask, NTask;

int Local_nx, Local_x_start;

int IdStart;

rfftwnd_mpi_plan Inverse_plan, Inverse_plan2;
fftw_real *Disp, *Velq, *Workspace, *Workspace2;
fftw_complex *Cdata, *Cdata2;
fftw_complex *delta_complx, *deltadot_complx;
fftw_complex     *DeltaDotField, *VelPrefac;

double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
double InputSpectrum_UnitLength_in_cm;
double G, Hubble;
double RhoCrit;

double Omega, OmegaLambda, OmegaRadiation, OmegaDM_2ndSpecies, Sigma8;
double OmegaBaryon, HubbleParam, Z_eq;
double ShapeGamma;
double PrimordialIndex;
double Dplus;			/* growth factor */


int ReNormalizeInputSpectrum;
