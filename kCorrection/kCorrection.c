#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream.h>
#include <fstream.h>


double integralX(double a);                        //a is e1 or E1
double equalgreat(double upperBound);              //upperBound is e2 or E2
double f_less(double currentE);                    //currentE is current E integration at for numerical integration

double integral_dIdalpha(double a);             //a is e1,E1, or x, and b is e2,E2, or upperBound
double  f_dIdalpha_greater(double upperBound);  //self explanitory (or see above)
double f_dIdalpha_less(double currentE);        //self explanitory (or see above)

double integral_dIdbeta(double b);              //b is e2 or E2.                
double f_dIdbeta_greater(double currentE);      //f_dIdbeta_less = 0;

double integral_dIdEo(double a);
double f_dIdEo_greater(double upperBound);
double f_dIdEo_less(double currentE);

double sigma_k(double XObserv, double XBolom, double k, double sig_alpha_x, double sig_beta_x, double sig_Epeak_x);

double f_DL(double x);
double integrate_DL(double z); //lowerbound=0, upperbound=z, epsilon=1.e-8

   
//-------Global paramters DEFINED BY USER!-------------

char grb[10];  //GRB (yymmdd)

double z;  //redshift
double sig_z; 

//bolometric bandpass
double e1; //keV  
double e2; //keV

//observed bandpass
double E1; //keV 
double E2; //keV

double alpha;  //low spectral index
double sig_alpha_l;
double sig_alpha_r;
double beta;   //high spectral index
double sig_beta_l;
double sig_beta_r;

double Epeak;  //Peak Energy in keV (1 keV = 0.001 MeV)
double sig_Epeak_l;
double sig_Epeak_r;

//precision
double epsilon=1.e-2;

double Eo;
double Ec;

//luminosity distance, from my Integrate.c
double DL;
double sig_DL;

//fluence
double S;
double sig_S_l;
double sig_S_r;


int main()
{

  //READ IN THE GLOBAL PARAMETERS FROM FILE!

  cout<<"THE INPUT FILE MUST BE A TEXT FILE, and there should be no blanks for any values.\n";
  cout<<"(1) The columns should be in this order: \n";
  cout<<"      grb,z,dz,S,dS+,dS-,e1,e2,alpha,dalpha+,dalpha-,beta,dbeta+,dbeta-,Epeak,dEpeak+,dEpeak- \n";
  cout<<"(2) S is the fluence, in units of 10^-7 ergs cm^-2. \n";
  cout<<"(3) e1 and e2 define the bolometric bandpass, and E1 and E2 define the observed bandpass (assumed to be 0.1 and 10000 keV respectively); in units of keV.\n";
  cout<<"(4) The values for alpha and beta should be positive. \n";
  cout<<"      VERY IMPORTANT: if beta is a CPL or PL, then beta should equal 0.\n";  
  cout<<"      This is a flag to the program to use a CPL or PL.\n";
  cout<<"(5) Epeak should be in units of keV.\n";
  cout<<"*** Do you still wish to continue (y/n)? ***\n";
  char choice;
  cin>>choice;

  if(choice=='y' || choice=='Y')
  {

	ifstream inFile;
 	char inputFilename[50];
  	cout<<"Please input the filename (including the directory it's in) that you wish to read in:\n";
  	cin>>inputFilename;

  	int rows; //how many rows are in the input file
  	cout<<"How many rows are in this file?  \n";
 	cin>>rows;

	cout<<"The output file of k values is called \"k_vals\", and the columns are in this order: grb, k, dk+, dk- \n";
  	cout<<"The output file of the isotropic energy is called \"isotropic_Energy\", and the columns are in this order: grb, E_iso, dE+, dE- \n";
        
	char outdir[50];
        
        cout<<"What directory do you want these in?\n";
        cin>>outdir;

	char outdircpy[50];
	strcpy(outdircpy,outdir);

	char s1[10]="/k_vals";
		printf("s1 = %s\n",s1);
	char s2[20]="/isotropic_Energy";
		printf("s2 = %s\n",s2);	
	
	char kout[50];
	char isoout[50];
	strcpy(kout, strcat(outdir,s1));
		 printf("kout = %s\n",kout);
	strcpy(isoout,strcat(outdircpy,s2));
		 printf("isoout = %s\n",isoout);
  	FILE *output=fopen(kout,"w");
  	FILE *isotropic=fopen(isoout,"w");

  	inFile.open(inputFilename, ios::in);

  	if(!inFile)
  	{
  	  cout<<"Can't open file\n";
  	  exit(1);
  	}
  
  	int num;
  	int i=0;
  	printf("starting \n");

  	while(i<rows)
        {
  	  printf("i=%d\n",i);
    
  	  inFile>>grb;
		printf("grb = %s\n",grb);
  	  inFile>>z;
                printf("z = %3.3e\n",z);
  	  inFile>>sig_z;
                printf("sig_z = %3.3e\n",sig_z);
  	  inFile>>S;
                printf("original S = %3.3e\n",S);
  	  S=S*pow(10,-7);
                printf("S = %3.3e \n",S);
  	  inFile>>sig_S_l;
                printf("original sig_S_l =  %3.3e\n",sig_S_l);
  	      sig_S_l=sig_S_l*pow(10,-7);
                printf("sig_S_l\n =  %3.3e\n",sig_S_l);
  	  inFile>>sig_S_r;
                printf("original sig_S_r =  %3.3e\n",sig_S_r);
		sig_S_r=sig_S_r*pow(10,-7);
                printf("sig_S_r =  %3.3e\n",sig_S_r);
  	  inFile>>e1;
                printf("e1 =  %3.3e\n",e1);
  	  inFile>>e2;
                printf("e2 =  %3.3e\n",e2);
  	  //inFile>>E1;
  	  //inFile>>E2;

  	  inFile>>alpha;
                printf("original alpha =  %3.3e\n",alpha);
		alpha=(-1)*alpha;
                printf("alpha =  %3.3e\n",alpha);
	  inFile>>sig_alpha_l;
                printf("sig_alpha_l =  %3.3e\n",sig_alpha_l);
 	  inFile>>sig_alpha_r;
                printf("sig_alpha_r =  %3.3e\n",sig_alpha_r);
	  inFile>>beta;
                printf("original beta =  %3.3e\n",beta);
		if(beta!=0){beta=(-1)*beta;}
		else {beta=-50.0;}
                printf("beta =  %3.3e\n",beta);
   	  inFile>>sig_beta_l;
                printf("sig_beta_l =  %3.3e\n",sig_beta_l);
	  inFile>>sig_beta_r;
                printf("sig_beta_r =  %3.3e\n",sig_beta_r);
  	  inFile>>Epeak;
                printf("Epeak =  %3.3e\n",Epeak);
  	  inFile>>sig_Epeak_l;
                printf("sig_Epeak_l =  %3.3e\n",sig_Epeak_l);
  	  inFile>>sig_Epeak_r;
                printf("sig_Epeak_r =  %3.3e\n",sig_Epeak_r);

  	  E1=(0.1/(1+z)); //keV
                printf("E1 =  %3.3e\n",E1);
  	  E2=(10000/(1+z)); //keV
                printf("E2 =  %3.3e\n",E2);
  	  Eo=(Epeak)/(2+alpha);
                printf("Eo =  %3.3e\n",Eo);
  	  Ec=(alpha-beta)*Eo;
                printf("Ec =  %3.3e\n",Ec );

    	  if(beta>-10.0)
          {     
			printf("beta>-10.0\n");
    		//-----OBSERVER "FLUENCE"------
    		double XObserv=0.0;
   		double upperBound=e2;
  		XObserv=integralX(e1)+equalgreat(upperBound);
 
    		//------BOLOMETRIC "FLUENCE"-----
    		double XBolom=0.0;
    		upperBound=E2;
    		XBolom=integralX(E1)+equalgreat(upperBound);


		//------Flux / Normalization--------
		double flux = S / XObserv;
		printf("flux=%3.3e\n",flux);


    	
		//-------SOLVE FOR k!-------
   		double k=XBolom/XObserv;
    		printf("************GRB %s*********\n",grb);
    		printf("*          k=%3.3e        *\n", k);
    		//-----error!--------------
    		double sig_k_l=sigma_k(XObserv, XBolom, k, sig_alpha_l, sig_beta_l, sig_Epeak_l);
    		printf("*    sigma_k_l=%3.3e        *\n",sig_k_l); 
    		double sig_k_r=sigma_k(XObserv, XBolom, k, sig_alpha_r, sig_beta_r, sig_Epeak_r);
    		printf("*    sigma_k_r=%3.3e        *\n",sig_k_r);
    		fprintf(output,"%s  %3.3e   %3.3e   %3.3e   \n",grb,k,sig_k_l,sig_k_r);
		
		printf("******now finding E_iso:*************************\n");
    	
		//------Luminosity Distance and it's error, from my Integrate.c------------
		double c = 3.e5;      //speed of light, units: [km/s]
  		double H_o = 73;       //Hubble Constant, units are [km s-1 Mps-1]
		DL = c * (1+z) / H_o;
		DL=DL*integrate_DL(z);  //cm (converted in integrate_DL)
		sig_DL=pow((DL/ (1+z))*sig_z, 2)+pow( ( (c*(1+z) / H_o)*(sig_z*integrate_DL(z)) ) ,   2);
		sig_DL=sqrt(sig_DL);	

		//------E_iso--------------------	
		double t = 4 * 3.1415 * pow(DL,2);
		t=t/(1+z);
        	double E_iso = S*t*k;

		double partA = pow((E_iso/k),2)*pow(sig_k_l,2);
		double partB = pow((E_iso/S),2) * pow(sig_S_l,2);
		double partC = pow(((E_iso)/(1+z)),2)*pow(sig_z,2);
		double partD = pow((E_iso*(2/DL)),2)*pow(sig_DL,2);
		double sig_E_iso_l= partA + partB + partC + partD;
		sig_E_iso_l=sqrt(sig_E_iso_l);

 		partA = pow((E_iso/k),2)*pow(sig_k_r,2);
        	partB = pow((E_iso/S),2) * pow(sig_S_r,2); 
        	partC = pow(((E_iso)/(1+z)),2)*pow(sig_z,2); 
        	partD = pow((E_iso*(2/DL)),2)*pow(sig_DL,2);
        	double sig_E_iso_r= partA + partB + partC + partD;
		sig_E_iso_r=sqrt(sig_E_iso_r);

		printf("*  %s  %3.3e   %3.3e   %3.3e	%3.3e *\n", grb, E_iso, sig_E_iso_l, sig_E_iso_r,flux);
		printf("**************************************\n");
		fprintf(isotropic, "%s %3.3e %3.3e %3.3e\n", grb, E_iso, sig_E_iso_l, sig_E_iso_r,flux);
		printf("end if beta>-10.0\n");
	        
          }//end if      
    	  else
    	  {
    		//-----OBSERVER "FLUENCE"------
    		double XObserv=0.0;  
    		double temp=Ec;
    		Ec=e2;
    		XObserv=integralX(e1); 
    
    		//------BOLOMETRIC "FLUENCE"-----
    		double XBolom=0.0;
    		Ec=E2;
    		XBolom=integralX(E1);

 		//------Flux / Normalization--------
        	double flux = S / XObserv;


    		//-------SOLVE FOR k!-------
    		Ec=temp;
    		double k=XBolom/XObserv;
    		printf("************GRB %s*********\n",grb);
    		printf("*          k=%3.3e        *\n", k);
    		//-----error!--------------
    		double sig_k_l=sigma_k(XObserv, XBolom, k, sig_alpha_l, sig_beta_l, sig_Epeak_l);
    		printf("*    sigma_k_l=%3.3e        *\n",sig_k_l);
    		double sig_k_r=sigma_k(XObserv, XBolom, k, sig_alpha_r, sig_beta_r, sig_Epeak_r);
    		printf("*    sigma_k_r=%3.3e        *\n",sig_k_r);
    		printf("*******************************\n");
    		fprintf(output,"%s	%3.3e	%3.3e	%3.3e	\n",grb,k,sig_k_l,sig_k_r);
	

                printf("******now finding E_iso:*************************\n");
                
                //------Luminosity Distance and it's error, from my Integrate.c------------
                double c = 3.e5;      //speed of light, units: [km/s]
                double H_o = 73;       //Hubble Constant, units are [km s-1 Mps-1]
                DL = c * (1+z) / H_o;
                DL=DL*integrate_DL(z);  //cm (converted in integrate_DL)
                sig_DL=pow((DL/ (1+z))*sig_z, 2)+pow( ( (c*(1+z) / H_o)*(sig_z*integrate_DL(z)) ) ,   2);
                sig_DL=sqrt(sig_DL);
                
                //------E_iso--------------------
                double t = 4 * 3.1415 * pow(DL,2);
                t=t/(1+z);
                double E_iso = S*t*k;
                
                double partA = pow((E_iso/k),2)*pow(sig_k_l,2);
                double partB = pow((E_iso/S),2) * pow(sig_S_l,2);
                double partC = pow(((E_iso)/(1+z)),2)*pow(sig_z,2);
                double partD = pow((E_iso*(2/DL)),2)*pow(sig_DL,2);
                double sig_E_iso_l= partA + partB + partC + partD;
                sig_E_iso_l=sqrt(sig_E_iso_l);
                
                partA = pow((E_iso/k),2)*pow(sig_k_r,2);
                partB = pow((E_iso/S),2) * pow(sig_S_r,2);
                partC = pow(((E_iso)/(1+z)),2)*pow(sig_z,2);
                partD = pow((E_iso*(2/DL)),2)*pow(sig_DL,2);
                double sig_E_iso_r= partA + partB + partC + partD;
                sig_E_iso_r=sqrt(sig_E_iso_r);
        
                printf("*  %s  %3.3e   %3.3e   %3.3e    %3.3e *\n", grb, E_iso, sig_E_iso_l, sig_E_iso_r,flux);
                printf("**************************************\n");
                fprintf(isotropic, "%s %3.3e %3.3e %3.3e\n", grb, E_iso, sig_E_iso_l, sig_E_iso_r,flux);
          }//end else
          i++;     
  	}//end while
  	printf("done\n");
  	fclose(output);
  	fclose(isotropic);
	}//end if yes
	else {exit(1);}

}//end main


/*
 * function to integrate for Luminosity Distance.
 */
double f_DL(double x)
{
    
  //WMAP3 Parameters:
  double Omega_m = 0.24;        //matter (WMAP3: 0.24)
  double Omega_lambda = 0.76;   //dark matter (WMAP3: 0.76)
  // printf("x=%15.3e\n",x);  
  //printf("integral: %15.3e\n",  1/( sqrt( Omega_m * pow((1+x),3) + Omega_lambda ) ));
  double val;
  val = 1+x;
  val = Omega_m*(val*val*val);
  val = val + Omega_lambda;
  val = sqrt(val);
  return 1/(val);
}

double integrate_DL(double z)
{
    double a=0;
    double b=z;
    double epsilon=1.e-8;
    double result;
    int i;
    int n;
    double h;
    double s;
    double s1;
    double s2;
    double s3;
    double x;
   
    s2 = 1;
    h = b-a;
    s = f_DL(a)+f_DL(b);
    do
    {
        s3 = s2;
        h = h/2;
        s1 = 0;
        x = a+h;
        do
        {
            s1 = s1+2*f_DL(x);
            x = x+2*h;
        }
        while(x<b);
        s = s+s1;
        s2 = (s+s1)*h/3;
        x = fabs(s3-s2)/15;
    }
   while(x>epsilon);
    result = s2; //Mpc
    result = result * 3.086e24;  //cm
    return result;
}

double integralX(double a)  
{
    double result;
    int i;
    int n;
    double h;
    double s;
    double s1;  
    double s2;
    double s3;
    double x;
    s2 = 1;
    h = Ec-a; 
    s = f_less(a) + f_less(Ec);

    do
    {
        s3 = s2;
        h = h/2;
        s1 = 0;
        x = a+h;
        do
        {
	    s1 = s1+2*f_less(x);
            x = x+2*h;
        }
        while(x<Ec);
        s = s+s1;
        s2 = (s+s1)*h/3;
        x = fabs(s3-s2)/15;
    }
    while(x>epsilon);
    result = s2;
     //printf("RESULT OF INTEGRAL: %3.3e\n",result);
    return result;
}

double f_less(double currentE)
{
  double val1=currentE/100; 
  val1=pow(val1,alpha); 
  double val2=currentE*(2+alpha)*(-1);
  val2=val2/Epeak;
  val2=exp(val2);
  double val=currentE*val1*val2;
	//printf("E<Ec:%3.3e\n",val);
  return val;
}

double equalgreat(double upperBound)
{
  double val=(alpha-beta)*Epeak;
  val=val/(100*(2+alpha));
  val=pow(val,(alpha-beta));
  val=val*exp(beta-alpha);
   //printf("ec:%3.3e\n",Ec);
   //printf("upperBound:%3.3e\n",upperBound);
  double boundB = pow((upperBound),(beta+2))/(beta+2);
  double boundA = pow((Ec),(beta+2))/(beta+2);

  val=val*(boundB-boundA); 
  val=val/pow(100,beta); 
   //printf("result of equalgreat:%3.3e\n",val);
  return val;
}
 
double f_dIdalpha_less(double currentE)
{   
   double val=currentE*pow((currentE/100),alpha)*exp((-currentE*(2+alpha))/Epeak)*log(currentE/100);
     //printf("val of f_dIdalpha_less:%3.3e\n",val);
   return val;
}

double  f_dIdalpha_greater(double upperBound)
{
   double val=equalgreat(upperBound)*log(((alpha-beta)*Eo)/100);
       //printf("result of f_dIdalpha_greater:%3.3e\n",val);
   return val;
}

double integral_dIdalpha(double a)
{
    double result;
    int i;
    int n;
    double h;
    double s;
    double s1;
    double s2;
    double s3;
    double x;
    s2 = 1;
    h = Ec-a;
    s = f_dIdalpha_less(a) + f_dIdalpha_less(Ec);
    do
    {
        s3 = s2;
        h = h/2;
        s1 = 0;
        x = a+h;
        do
        {
            s1 = s1+2*f_dIdalpha_less(x);
            x = x+2*h;
        }
        while(x<Ec);
        s = s+s1;
        s2 = (s+s1)*h/3;
        x = fabs(s3-s2)/15;
    }
    while(x>epsilon);
    result = s2;
      //printf("result of integral_dIdalpha: %3.3e\n",result);
    return result;
}

double f_dIdbeta_greater(double upperBound) //the  f_dIdbeta_less is 0
{
   double c=pow((Ec/100),(alpha-beta))*(1/pow(100,beta))*exp(beta-alpha);
      //printf("should be 88732 --> c:%3.3e\n",c);
   
   double d=log(Ec);
   double boundBd=pow(upperBound,(beta+2))/(beta+2);
   double boundAd=pow(Ec,(beta+2))/(beta+2);
   d=d*(boundBd-boundAd);
      //printf("should be .32008 --> d:%3.3e\n",d);

   double boundBf=log(upperBound)/(beta+2)*pow(upperBound,(beta+2));   
   double boundAf=log(Ec)/(beta+2)*pow(Ec,(beta+2));
   double f=boundBf-boundAf;
     //printf("should be 234874 --> f:%3.3e\n",f);

   double g=1/(beta+2);
   double boundBg=pow(upperBound,(beta+2))/(beta+2);
   double boundAg=pow(Ec,(beta+2))/(beta+2);
   g=g*(boundBg-boundAg);
     //printf("should be -.116778 --> g:%3.3e\n",g);

   double val=c*((f-g)-d);
    //printf("should be 2807.84 --> val of f_dIdbeta_greater: %3.3e\n",val);
   return val;
}

double f_dIdEo_less(double currentE) 
{
   double val=pow((currentE/100),alpha)*exp(-currentE/Eo)*pow((currentE/Eo),2);
      //printf("val of f_dIdEo_less:%3.3e\n",val);
   return val;
}

double f_dIdEo_greater(double upperBound)
{
   double val=pow((Ec/100),(alpha-beta))*exp(beta-alpha)*(1/pow(100,beta))*((alpha-beta)/Eo);
   double boundB=pow(upperBound,(beta+2))/(beta+2);
   double boundA=pow(Ec,(beta+2))/(beta+2);
   val=val*(boundB-boundA);
     //printf("result of f_dIdEo_greater:%3.3e\n",val);
   return val;
}

double integral_dIdEo(double a)
{
    double result;
    int i;
    int n;   
    double h;
    double s; 
    double s1;
    double s2;
    double s3;
    double x; 
    s2 = 1;  
    h = Ec-a;
    s = f_dIdEo_less(a) + f_dIdEo_less(Ec);
    do
    {
        s3 = s2;
        h = h/2;
        s1 = 0; 
        x = a+h;
        do
        {
            s1 = s1+2*f_dIdEo_less(x);
            x = x+2*h;
        }
        while(x<Ec);
        s = s+s1;
        s2 = (s+s1)*h/3;   
        x = fabs(s3-s2)/15;
    }
    while(x>epsilon);
    result = s2;
       //printf("result of integral_dIEo: %3.3e\n",result);  
    return result;
}


double sigma_k(double XObserv, double XBolom, double k, double sig_alpha_x, double sig_beta_x, double sig_Epeak_x)
{
    double result;
    if(beta>-5.0)
    {
    double partA=( (k/XBolom)*(integral_dIdalpha(E1) + f_dIdalpha_greater(E2))) ;
    partA=partA-((k/XObserv)*(integral_dIdalpha(e1)+f_dIdalpha_greater(e2)));
    partA=pow(partA,2)*pow(sig_alpha_x,2);

    double partB=((k/XBolom)*(f_dIdbeta_greater(E2)));
    partB=partB-((k/XObserv)*(f_dIdbeta_greater(e2)));
    partB=pow(partB,2)*pow(sig_beta_x,2);

    double partC=((k/XBolom)*(integral_dIdEo(E1)+f_dIdEo_greater(E2)));
    partC=partC-((k/XObserv)*(integral_dIdEo(e1)+f_dIdEo_greater(e2)));
    double sig_Eo=pow((1/(2+alpha)),2)*pow(sig_Epeak_x,2) + pow((-1/(2+alpha)),2)*pow(sig_alpha_x,2);
    sig_Eo = sqrt(sig_Eo);
    partC=pow(partC,2)*pow(sig_Eo,2);


    result=partA+partB+partC;
    result=sqrt(result);
    }

    else
    {
    double temp=Ec;

    Ec=E2;
    double partA=((k/XBolom)*(integral_dIdalpha(E1))) ;
    double partC=((k/XBolom)*(integral_dIdEo(E1)));

    Ec=e2;    
    partA=partA-((k/XObserv)*(integral_dIdalpha(e1)));
    partC=partC-((k/XObserv)*(integral_dIdEo(e1)));

    Ec=temp;
    partA=pow(partA,2)*pow(sig_alpha_x,2);
    double partB=0;
    double sig_Eo=pow((1/(2+alpha)),2)*pow(sig_Epeak_x,2) + pow((-1/(2+alpha)),2)*pow(sig_alpha_x,2);
    sig_Eo = sqrt(sig_Eo);
    partC=pow(partC,2)*pow(sig_Eo,2);
    
    
    result=partA+partB+partC;
    result=sqrt(result);    
    }
    return result;

}

