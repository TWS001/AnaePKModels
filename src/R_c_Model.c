#include <math.h>
#include <errno.h>
#include <stdbool.h>

double lbmJames(double wt,double ht,bool MALE);




void  // Propofol(incl. ke0), Schnider et al. Anesthesiology 1998
PropSchnider(
		   double *Age,// the patient's age
		   double *Weight,// the patient's weight
		   double *Height,// the patient's height (in cm)
		   bool *Gender, // Lean body mass
		   double *ArrParameter
           )
{
	double LBM = lbmJames(*Weight, *Height, *Gender);

	double V1 = 4.27;
	double V2 = 18.9 - (0.391*(*Age-53));
	double V3 = 238;

    double Cl1= 1.89+ ((*Height-177)*0.0264) + ((*Weight-77)*0.0456) - ((LBM-59)*0.0681);
	double Cl2 = 1.29 +((*Age-53)*-0.024);
	double Cl3 = 0.836;
	ArrParameter[0]= V1;
	ArrParameter[1]= V2;
	ArrParameter[2]= V3;
	ArrParameter[3]= Cl1;
	ArrParameter[4]= Cl2;
	ArrParameter[5]= Cl3;
	ArrParameter[6]= 0.456;
}

void  // Propofol(incl. ke0), Schnider et al. Anesthesiology 1998
PropSchnider2(
		   double *Age,// the patient's age
		   double *Weight,// the patient's weight
		   double *Height,// the patient's height (in cm)
		   double *lbm, // Lean body mass
		   double *ArrParameter
           )
{

	double V1 = 4.27;
	double V2 = 18.9 - (0.391*(*Age-53));
	double V3 = 238;

    double Cl1= 1.89+((*Height-177)*0.0264)+ ((*Weight-77)*0.0456)-((*lbm-59)*0.0681);
	double Cl2 = 1.29 +((*Age-53)*-0.024);
	double Cl3 = 0.836;
	ArrParameter[0]= V1;
	ArrParameter[1]= V2;
	ArrParameter[2]= V3;
	ArrParameter[3]= Cl1;
	ArrParameter[4]= Cl2;
	ArrParameter[5]= Cl3;
	ArrParameter[6]= 0.456;
}


void  // Remifentanil(incl. ke0), Minto et al. (Age, LBM) Anesthesiology 199?
RemiMinto(
		   double *Age,// the patient's age
		   double *LBM, // Lean body mass
		   double *ArrParameter
           )
{
	
	double V1 = 5.1-(0.0201*(*Age-40))+(0.072*(*LBM-55));
	double V2 = 9.82-(0.0811*(*Age-40))+(0.108*(*LBM-55));
	double V3 = 5.42;

    double Cl1= 2.6-(0.0162*(*Age-40))+(0.0191*(*LBM-55));
	double Cl2 = 2.05-(0.0301*(*Age-40));
	double Cl3 = 0.076-(0.00113*(*Age-40));
	double ke0 = 0.595-0.007*(*Age-40);

	ArrParameter[0]= V1;
	ArrParameter[1]= V2;
	ArrParameter[2]= V3;
	ArrParameter[3]= Cl1;
	ArrParameter[4]= Cl2;
	ArrParameter[5]= Cl3;
	ArrParameter[6]= ke0;

}


void  // Propofol, Eleveld BJA 2018
PropEleveldFinal(
		  double *weight, // Enter the weight in kg
		  double *height, // Enter the height in cm
		  double *pma, // Enter postmenstrual age
		  double *age, // Enter age as years
		  bool *male,  // Enter gender (1=MALE)
		  bool *venous,  // Enter if concentration is venous (1)
		  bool *opiates,  // Enter if opiates are given (1)
		  bool *patient, // For patient enter 1 for volunteer 0
		  double *ArrParameter
)
{
   /* some setup of variables so we can reuse NONMEM code */
	/* fake theta and eta matrix so we can keep code similar to NONMEM code for clarity */
	/* put dummy in index 0 so 1-offset fortran array addressing can be emulated */
    const double WGT = *weight;
	const double AGE = *age;
	const double HGT = *height;
	const double PMA = *pma;


	const unsigned M1F2 = *male ? 1 : 2;
	/* const unsigned P1V2 = *patient ? 1 : 2; */
	const int A1V2 = *venous ? 2 : 1;
	const int TECH = *opiates ? 2 : 1;
#define EXP(a)		(exp(a))

	const double PKTHETA[] = { 0, /* offset for 1-based FORTRAN arrays */
		1.837860e+00,
		3.238730e+00,
		5.608800e+00,
		5.819830e-01,
		5.596720e-01,
		1.030460e-01,
		1.913070e-01,
		3.744220e+00,
		2.203300e+00,
		-1.563300e-02,
		-2.857090e-03,
		3.513130e+00,
		-1.381660e-02,
		4.223570e+00,
		7.420430e-01,
		2.656420e-01,
		3.498850e-01,
		-3.849270e-01,
	};
	const double PDTHETA[] = { 0, /* offset for 1-based FORTRAN arrays */
		1.124750e+00,
		-1.921620e+00,
		9.298240e+01,
		3.877210e-01,
		2.082830e+00,
		5.173990e-02,
		-6.348370e-03,
		2.156050e-01,
		6.386410e-01,
	};
	const double _ETA[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, };
#define ETA(a)		(_ETA[a])

/* PK mdoel -------------------------------------------------------------- */
#define THETA(a)	(PKTHETA[a])

/* Al-sallami FFM */
    const double HT2=(HGT/100.)*(HGT/100.);
    const double MATM=0.88+((1-0.88)/(1+pow(AGE/13.4, -12.7)));
    const double MATF=1.11+((1-1.11)/(1+pow(AGE/7.1, -1.1)));
    const double MATR=0.88+((1-0.88)/(1+pow(35./13.4, -12.7)));
    const double FFMM=MATM*42.92*(HT2)*WGT/(30.93*(HT2)+WGT);
    const double FFMF=MATF*37.99*(HT2)*WGT/(35.98*(HT2)+WGT);
    const double FFMR=MATR*42.92*(1.7*1.7)*70./(30.93*(1.7*1.7)+70.);
    const double FFM=FFMM*(2-M1F2) + FFMF*(M1F2-1);
    const double NFFM=FFM/FFMR;
/* maturation */
    const double DV1=1;
    const double DV2=1;
    const double DV3=1;
/* sigmoidal maturation of CL */
    const double PMW=PMA*52.;
    const double PMR=(35.+40./52.)*52.;
    const double ME50=EXP(THETA(8));
    const double MGAM=EXP(THETA(9));
    const double MCL=pow(PMW, MGAM)/(pow(PMW, MGAM)+pow(ME50, MGAM));
    const double RCL=pow(PMR, MGAM)/(pow(PMR, MGAM)+pow(ME50, MGAM));
    const double DCL=MCL/RCL;
    const double DQ2=1;
/* sigmoidal maturation of Q3 based on 40 weeks gestation */
    const double PMEW=AGE*52.+40.;
    const double PMER=35.*52.+40.;
    const double QE50=EXP(THETA(14));
    const double MQ3=PMEW/(PMEW+QE50);
    const double RQ3=PMER/(PMER+QE50);
    const double DQ3=MQ3/RQ3;
/* aging */
    const double KV1=1;
    const double KV2=EXP(THETA(10)*(AGE-35.));
    const double KV3=EXP(THETA(13)*(AGE)*(TECH-1));
    const double KCL=EXP(THETA(11)*(AGE)*(TECH-1));
    const double KQ2=1;
    const double KQ3=1;
/* covariate structure */
/* V1 scales sigmoid with weight */
    const double VV50=EXP(THETA(12));
    const double CV1=WGT/(WGT+VV50);
    const double RV1=70./(70.+VV50);
    const double M1 =(CV1/RV1) * KV1 * DV1;
    const double VCV1=(A1V2-1)*(1-CV1);
    const double V1 =EXP(THETA(1)+ETA(1)) * M1 * (1+VCV1*EXP(THETA(17)));
    const double M2 =(WGT/70.) * KV2 * DV2;
    const double V2 =EXP(THETA(2)+ETA(2)) * M2;
    const double M3 =(NFFM) * KV3 * DV3;
    const double V3 =EXP(THETA(3)+ETA(3)) * M3;
    const double M4 =pow(WGT/70., 0.75) * KCL * DCL;
    const double CL =EXP((2-M1F2)*THETA(4)+(M1F2-1)*THETA(15)+ETA(4)) * M4;
    const double RV2=EXP(THETA(2));
    const double M5 =pow(V2/RV2, 0.75) * KQ2 * DQ2;
    const double KM5=1+EXP(THETA(16))*(1-MQ3);
    const double Q2 =EXP(THETA(5)+ETA(5)+(A1V2-1)*THETA(18)) * M5 * KM5;
    const double RV3=EXP(THETA(3));
    const double M6 =pow(V3/RV3, 0.75) * KQ3 * DQ3;
    const double Q3 =EXP(THETA(6)+ETA(6)) * M6;
/* error model */
/*  const double RESV=THETA(7); unused */

/* PD mdoel -------------------------------------------------------------- */
#undef THETA
#define THETA(a)	(PDTHETA[a])

/* effect compartment */
    /* const double E50 =EXP(THETA(1)+ETA(1)+THETA(7)*(AGE-35.)); */
    const double TKE0=(2-A1V2)*THETA(2) + (A1V2-1)*THETA(8);
    const double KE0 =EXP(TKE0+ETA(2))*pow(WGT/70, -0.25);
/*  const double EMAX=THETA(3); unused */
/*  const double GAM =EXP(THETA(4)); unused */
/*  const double GAM1=EXP(THETA(9)); unused */
/* residual error */
/*  const double RESD=EXP(THETA(5))*EXP(ETA(3)); unused */
/*  const double ALAG1=15./60. + EXP(THETA(6)*(AGE))/60.; unused */

#undef EXP
#undef THETA



	/* construct return type */
	ArrParameter[0]= V1;
	ArrParameter[1]= V2;
	ArrParameter[2]= V3;
	ArrParameter[3]= CL ;
	ArrParameter[4]= Q2 ;
	ArrParameter[5]= Q3;
	ArrParameter[6]= KE0;

}


void  // Propofol, Marsh et al. (Weight) BJA ???
PropMarsh(
		  double *Weight, // the patient's weigh
		  double *ArrParameter
           )
{

	double V1 = *Weight*0.228;
	double V2 = V1*(0.112/0.055);
	double V3 = V1*(0.0419/0.0033);

	double Cl1 = V1*0.119;
	double Cl2 = V2*0.055;
	double Cl3 = V3*0.0033;

	ArrParameter[0]= V1;
	ArrParameter[1]= V2;
	ArrParameter[2]= V3;
	ArrParameter[3]= Cl1;
	ArrParameter[4]= Cl2;
	ArrParameter[5]= Cl3;

}

void  // Propofol, Diprifusor (pers. com: Iain Glen  k12=0.114), incl. ke0
Diprifusor(
		  double *Age, // Age, a dummy parameter in Diprifusor!
		  double *Weight, // the patient's weight
		  double *ArrParameter
	)
{
	double V1 = *Weight*0.228;
		if(*Age<16)V1=0.0;
	double V2 = V1*(0.114/0.055);
	double V3 = V1*(0.0419/0.0033);

	double Cl1 = V1*0.119;
	double Cl2 = V2*0.055;
	double Cl3 = V3*0.0033;

	ArrParameter[0]= V1;
	ArrParameter[1]= V2;
	ArrParameter[2]= V3;
	ArrParameter[3]= Cl1;
	ArrParameter[4]= Cl2;
	ArrParameter[5]= Cl3;
	ArrParameter[6]= 0.26;
}

void  // Scott Fentanyl weight scaled with ke0 0.147
FentScott(
		  double *Weight, // the patient's weight
		  double *ArrParameter
           )
{

	ArrParameter[0]= *Weight*12.7/70;
	ArrParameter[1]= *Weight*49.34479167/70;
	ArrParameter[2]= *Weight*296.8831169/70;
	ArrParameter[3]= *Weight*0.7112/70;
	ArrParameter[4]= *Weight*4.7371/70;
	ArrParameter[5]= *Weight*2.286/70;
	ArrParameter[6]= 0.147;

}

double lbmJames(double wt,double ht,bool MALE){
  double a = 1.07;
  double b = 148;
  if(MALE){
    a = 1.1;
    b = 128;
  }
  return(a * wt - b * pow(wt/ht,2));
}
