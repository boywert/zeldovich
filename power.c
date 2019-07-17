#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "allvars.h"
#include "proto.h"

#define WORKSIZE 100000

static double R8;
static double r_tophat;

static double AA, BB, CC;
static double nu;
static double Norm;


static int NPowerTable;

static struct pow_table
{
  double logk, logD;
}
 *PowerTable;


double PowerSpec(double k)
{
  double power = 0;

  switch (WhichSpectrum)
    {
    case 1:
      power = PowerSpec_EH(k);
      break;

    case 2:
      power = PowerSpec_Tabulated(k);
      break;

    default:
      power = PowerSpec_Efstathiou(k);
      break;
    }

  power *= pow(k, PrimordialIndex - 1.0);

  return power;
}


void read_power_table(void)
{
  FILE *fd;
  char buf[500];
  double k, p;


  sprintf(buf, FileWithInputSpectrum);

  if(!(fd = fopen(buf, "r")))
    {
      printf("can't read input spectrum in file '%s' on task %d\n", buf, ThisTask);
      FatalError(17);
    }

  NPowerTable = 0;
  do
    {
      if(fscanf(fd, " %lg %lg ", &k, &p) == 2)

	NPowerTable++;
      else
	break;
    }
  while(1);

  fclose(fd);


  if(ThisTask == 0)
    {
      printf("found %d rows in input spectrum table\n", NPowerTable);
      fflush(stdout);
    }


  PowerTable = malloc(NPowerTable * sizeof(struct pow_table));

  sprintf(buf, FileWithInputSpectrum);

  if(!(fd = fopen(buf, "r")))
    {
      printf("can't read input spectrum in file '%s' on task %d\n", buf, ThisTask);
      FatalError(18);
    }

  NPowerTable = 0;
  do
    {
      double p;

      if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
	{
	  PowerTable[NPowerTable].logk = k;
	  PowerTable[NPowerTable].logD = p;
	  NPowerTable++;
	}
      else
	break;
    }
  while(1);

  fclose(fd);

  qsort(PowerTable, NPowerTable, sizeof(struct pow_table), compare_logk);
}

int compare_logk(const void *a, const void *b)
{
  if(((struct pow_table *) a)->logk < (((struct pow_table *) b)->logk))
    return -1;

  if(((struct pow_table *) a)->logk > (((struct pow_table *) b)->logk))
    return +1;

  return 0;
}

void initialize_powerspectrum(void)
{
  double res;

  InitTime = 1 / (1 + Redshift);
  OmegaRadiation = Omega/(2 + Z_eq);
  OmegaLambda = 1.0 - Omega - OmegaRadiation - 1e-30;


  R8 = 8 * (3.085678e24 / UnitLength_in_cm);	/* 8 Mpc/h */

  Norm = 1.0;
      /* tabulated file is already at the initial redshift */
  Dplus = 1.0;
  
}

double PowerSpec_Tabulated(double k)
{
  double logk, logD, P, kold, u, dlogk, Delta2;
  int binlow, binhigh, binmid;

  kold = k;

  k *= (InputSpectrum_UnitLength_in_cm / UnitLength_in_cm);	/* convert to h/Mpc */

  logk = log10(k);

  if(logk < PowerTable[0].logk || logk > PowerTable[NPowerTable - 1].logk)
    return 0;

  binlow = 0;
  binhigh = NPowerTable - 1;

  while(binhigh - binlow > 1)
    {
      binmid = (binhigh + binlow) / 2;
      if(logk < PowerTable[binmid].logk)
	binhigh = binmid;
      else
	binlow = binmid;
    }

  dlogk = PowerTable[binhigh].logk - PowerTable[binlow].logk;

  if(dlogk == 0)
    FatalError(777);

  u = (logk - PowerTable[binlow].logk) / dlogk;

  logD = (1 - u) * PowerTable[binlow].logD + u * PowerTable[binhigh].logD;

  Delta2 = pow(10.0, logD);

  P = Norm * Delta2 / (4 * M_PI * kold * kold * kold);

  return P;
}


double PowerSpec_Efstathiou(double k)
{
  return Norm * k / pow(1 + pow(AA * k + pow(BB * k, 1.5) + CC * CC * k * k, nu), 2 / nu);
}

double PowerSpec_EH(double k)	/* Eisenstein & Hu */
{
  return Norm * k * pow(tk_eh(k), 2);
}


double tk_eh(double k)		/* from Martin White */
{
  double q, theta, ommh2, a, s, gamma, L0, C0;
  double tmp;
  double omegam, ombh2, hubble;

  /* other input parameters */
  hubble = HubbleParam;

  omegam = Omega;
  ombh2 = OmegaBaryon * HubbleParam * HubbleParam;

  if(OmegaBaryon == 0)
    ombh2 = 0.04 * HubbleParam * HubbleParam;

  k *= (3.085678e24 / UnitLength_in_cm);	/* convert to h/Mpc */

  theta = 2.728 / 2.7;
  ommh2 = omegam * hubble * hubble;
  s = 44.5 * log(9.83 / ommh2) / sqrt(1. + 10. * exp(0.75 * log(ombh2))) * hubble;
  a = 1. - 0.328 * log(431. * ommh2) * ombh2 / ommh2
    + 0.380 * log(22.3 * ommh2) * (ombh2 / ommh2) * (ombh2 / ommh2);
  gamma = a + (1. - a) / (1. + exp(4 * log(0.43 * k * s)));
  gamma *= omegam * hubble;
  q = k * theta * theta / gamma;
  L0 = log(2. * exp(1.) + 1.8 * q);
  C0 = 14.2 + 731. / (1. + 62.5 * q);
  tmp = L0 / (L0 + C0 * q * q);
  return (tmp);
}



double TopHatSigma2(double R)
{
  double result, abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  F.function = &sigma2_int;

  r_tophat = R;

  gsl_integration_qag(&F, 0, 500.0 * 1 / R,
		      0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return result;

  /* note: 500/R is here chosen as (effectively) infinity integration boundary */
}


double sigma2_int(double k, void *param)
{
  double kr, kr3, kr2, w, x;

  kr = r_tophat * k;
  kr2 = kr * kr;
  kr3 = kr2 * kr;

  if(kr < 1e-8)
    return 0;

  w = 3 * (sin(kr) / kr3 - cos(kr) / kr2);
  x = 4 * PI * k * k * w * w * PowerSpec(k);

  return x;
}


double GrowthFactor(double astart, double aend)
{
  return growth(aend) / growth(astart);
}


double growth(double a)
{
  double hubble_a;

  hubble_a = sqrt(Omega / (a * a * a) + (1 - Omega - OmegaLambda) / (a * a) + OmegaLambda);

  double result, abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);

  F.function = &growth_int;

  gsl_integration_qag(&F, 0, a, 0, 1.0e-8, WORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);

  gsl_integration_workspace_free(workspace);

  return hubble_a * result;
}

double DEBA18_prefac(double k, double a) {
  double y;
  double hubble_a;
  // double kscale = 0.12; Sten Delos et al 2018 - not used anymore - Boyd
  double ddot, d;
  double a_eq = 1.0/(1+Z_eq);
  double res;
  double T_cmb = 2.7255; // Fixsen 2009
  double N_eff = 3.046; // Mangano et al (2002,2005)
  double T_ref = 2.7;
  double Theta = T_cmb/T_ref;
  double f_v = 1.0 - 1.0/(N_eff*(7./8.)*pow(4./11.,4./3.) + 1.0);
  double I_2 = 0.594*(1 - 0.631*f_v + 0.284*f_v*f_v); // eq B14 Hu & Sukiyama 1996
  double k_eq  = 9.67e-2 * Omega * HubbleParam * HubbleParam * sqrt(1.0 - f_v) / (Theta*Theta); // 1/Mpc
  double logk_term, aHaEQ;
  k *= (3.085678e24 / UnitLength_in_cm) * HubbleParam; //change to 1/Mpc
  hubble_a =
    Hubble * sqrt(Omega / pow(a, 3) + OmegaRadiation/ pow(a, 4) +
		  (1 - Omega - OmegaLambda ) / pow(a, 2) + OmegaLambda);
  y = a/a_eq;
  
  aHaEQ = (1.0 + sqrt(1.0 + 8.0*(k/k_eq)*(k/k_eq)))/(4.0 * (k/k_eq)*(k/k_eq)); //eq B13 Hu & Sukiyama 1996
  logk_term = log(4.0 * I_2 * exp(-3.0) / aHaEQ);
  d = (logk_term - log(( sqrt(1 + y) + 1 )/( sqrt(1 + y) - 1)))*(y + 2.0/3.0) + 2.0*sqrt(1 + y);  //eq D3 Hu & Sukiyama 1996
  ddot = logk_term - log((sqrt(1+y) + 1)/(sqrt(1+y) - 1)) + 2./sqrt(1.+y) + 2./(3.*y*sqrt(1 + y));
  
  res = y * hubble_a * ddot / d * sqrt(a);
  //printf("y = %g, k = %g h/Mpc, ddot = %g, d = %g, hubble_a = %g, res = %g\n", y, k, ddot, d, hubble_a, res);
  return res;
}

double DplusDEBA18(double k, double astart, double aend) {
  double y1, y2;
  double d1, d2;
  double kscale = 0.12;
  double a_eq = 1./(1+Z_eq);
  k *= (3.085678e24 / UnitLength_in_cm); //change to h/Mpc
  y1 = astart/a_eq;
  y2 = aend/a_eq;
  d1 = (log(k/kscale) - log((sqrt(1+y1)+1)/(sqrt(1+y1)-1)))*(y1+2.0/3.0)+2.0*sqrt(1+y1);
  d2 = (log(k/kscale) - log((sqrt(1+y2)+1)/(sqrt(1+y2)-1)))*(y2+2.0/3.0)+2.0*sqrt(1+y2);
  return d1/d2;
}

double growth_int(double a, void *param)
{
  return pow(a / (Omega + (1 - Omega - OmegaLambda) * a + OmegaLambda * a * a * a), 1.5);
}


double F_Omega(double a)
{
  double omega_a;

  omega_a = Omega / (Omega + a * (1 - Omega - OmegaLambda) + a * a * a * OmegaLambda);

  return pow(omega_a, 0.6);
}
