#include "mylib.h"

int main (int narg, char * argc[]) {

double L=1, diam=1.0, shell=0.5;
long int Nsphere,NsphereCond,NsphereIns, NsphereCondFin;

/*
cout << "Get Nsphere:";
cin>>Nsphere;


cout << "Get L:";
cin>>L;
*/
string suffix (argc[1]);

string strfile ("percFile.dat");
ofstream percFile;
percFile.open((strfile+suffix).c_str(), fstream::out);

string strfile1 ("condFile.dat");
ofstream condFile;
condFile.open((strfile1+suffix).c_str(), fstream::out);


//int Nsquare =(int)(0.98 *L/(diam));
//double Lsquare = L/Nsquare;
//clock_t start_all =start_time();
//penetrable_spheres(shell,diam,Nsphere, L);
double prob;
double eta2d=(diam+2*shell)*(diam+2*shell)*PI/4;
double eta3d=(diam+2*shell)*(diam+2*shell)*(diam+2*shell)*PI/6;
double phi3d=(diam)*(diam)*(diam)*PI/6;


/*
for (int i=Nsphere; i<Nsphere+120; i=i+4){
      prob=0;
      prob=impenetrable_spheres(shell,diam,i, L);
      percFile<<eta3d*i/(L*L*L)<<" "<<prob<<endl;
      //pi*0.55^2=0.95
      //pi/6=0.5236
}
*/
//getTime(start_all);
double num, shellin, shellfin;
double eta_in, eta_fin, phi2;
int numshell=NUM_SHELL, flag_s;

ofstream zrFile;
//ofstream zrFileCube;
ofstream zrFile2;
ofstream zrFile3;

/*string file1 ("grHomo_cc.dat");
zrFile.open((file1).c_str(), fstream::app);
zrFile<<"Realiz Nsph length Lato Nbin deltar "; 
for (int i=0;i<NBIN_GR;i++) {
     
     zrFile << i*DR_GR<<" ";
  } zrFile<<endl;


string file2 ("grHomo_angle.dat");
zrFile2.open((file2).c_str(), fstream::app);
zrFile2<<"Realiz Nsph length Lato Nbin deltar "; 
for (int i=0;i<NBIN_GR_ANGLE;i++) {
     
     zrFile2 << i*DR_GR_ANGLE<<" ";
  } zrFile2<<endl;*/

string file3 ("connectFile.dat");
zrFile3.open((file3+suffix).c_str(), fstream::app);
zrFile3<<"Realiz Nsph length Lato Nbin deltar "; 
for (int i=0;i<NBIN_MINDIST;i++) {
     
     zrFile3 << i*DR_MINDIST<<" ";
  } zrFile3<<endl;

/*
zrFileCube<<"Realiz Nsph Lato Nbin deltar "; 
for (int i=0;i<NBIN_GR;i++) {
     
     zrFileCube << i*DR_GR<<" ";
  } zrFileCube<<endl;
*/

/*   HpFile.close();
   ConFile.close();
   grFile.close();
    HpFileP.close();
   ConFileP.close();
   grFileP.close();*/

//zrFile.close();
//zrFile2.close();
zrFile3.close();




double diamCond=1.0;
double diamIns=4.0;

double volCond=PI/6*(diamCond+2*shell)*(diamCond+2*shell)*(diamCond+2*shell);
double volIns=PI/6*(diamIns)*(diamIns)*(diamIns);
double phiCond=PI/6*(diamCond)*(diamCond)*(diamCond);


NsphereCond=(int)(eta_in*L*L*L/phiCond);
NsphereIns=(int)(-log(1-phi2)*L*L*L/volIns);
NsphereCondFin=(int)(eta_fin*L*L*L/phiCond);

/*
0.2 4 0.1784
0.2 6 0.1555
0.2 8 0.138424499815182
0.2 10 0.1236
0.2 12 0.112669589334313
0.2 14 0.103101850654296
0.2 17 0.0911
0.2 20 0.0822636639194041
0.2 25 0.0702
0.2 30 0.0617933517384643
0.2 40 0.0495378902872627
*/


/*int nphi=11;
double length[]={10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
double phiVect[]={0.001, 0.003, 0.005, 0.01, 0.02, 0.035, 0.05, 0.09, 0.1236, 0.14, 0.16};
double NsphVect[]={500, 1000, 3000};
double LatoVect[]={(10+1)*5, (10+1)*5, (10+1)*5, (10+1)*5, (10+1)*5, (10+1)*5, (10+1)*5, (10+1)*5, (10+1)*5, (10+1)*5, (10+1)*5};
*/

/*int nphi=11;
double length[]={20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};
double phiVect[]={0.001, 0.003, 0.005, 0.008, 0.01, 0.02, 0.035, 0.05, 0.082, 0.09, 0.11};
double NsphVect[]={500, 1000, 3000};
double LatoVect[]={(20+1)*5, (20+1)*5, (20+1)*5, (20+1)*5, (20+1)*5, (20+1)*5, (20+1)*5, (20+1)*5, (20+1)*5, (20+1)*5, (20+1)*5};
*/
int nphi=11;
double length[]={5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
double phiVect[]={0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.125, 0.14, 0.155, 0.17, 0.19};
double NsphVect[]={500, 1000, 3000};
double LatoVect[]={(5+1)*5, (5+1)*5, (5+1)*5, (5+1)*5, (5+1)*5, (5+1)*5, (5+1)*5, (5+1)*5, (5+1)*5, (5+1)*5, (5+1)*5};








/*int nphi=5;
double length[]={10, 10, 10, 10, 10};
double phiVect[]={0.005, 0.01, 0.02, 0.035, 0.05,};
double NsphVect[]={500, 1000, 3000};
double LatoVect[]={(10+1)*5,(10+1)*5, (10+1)*5, (10+1)*5, (10+1)*5};
*/


double *percVect;
percVect=allocate_double_vector(REALIZ);

double *condVect;
condVect=allocate_double_vector(REALIZ);


//0.35 max conc for impenetrable spheres phi2=0
//0.2  max conc for segregation phi2=0.35
//~0.08 max conc for segregaton phi2=0.7
//getchar();


for (int i=0; i<nphi; i=i+1){

       /*if (phiVect[i]<5e-03) L=60;
       else if (phiVect[i]>=5e-2 && phiVect[i]<=1e-01) L=35;
       else if (phiVect[i]>1e-1 && phiVect[i]<2e-1) L=25;
       else L=15;  
       */
      double length1=length[i];
      double vspherocyl=PI*(length1*diamCond*diamCond/4 + diamCond*diamCond*diamCond/6);

      for (int l=0; l<REALIZ; l++)
            percVect[l]=0;

      for (int l=0; l<REALIZ; l++)
            condVect[l]=0;

      //int x=(int)(phiVect[i]*(6/PI)*(L*L*L));
     cerr<<endl<<"---------------------------------------------"<<endl; 
      //int x=NsphVect[i];
      L=LatoVect[i];
      //cerr<<"Nsphere "<<x<<endl;
      //double volume= vspherocyl*x/phiVect[i]; 
       //cerr<<"Volume " <<volume<<endl;
       //double num=1.0/3.0;
       //L=pow(volume, num);
       //cerr<<"Lato "<<L<<endl;
       //int Lato=(int)L;
       int x=(int) (phiVect[i]*(L*L*L)/vspherocyl);
       
       if (x<1000)
       {
            x=1000;
            double volume= vspherocyl*x/phiVect[i]; 
      
            double num=1.0/3.0;
            L=pow(volume, num);
       //cerr<<"Lato "<<L<<endl;
            int Lato=(int)L;      
       }
       
       
       
      cerr<<"Lato: "<<L<<" Nsphere: "<<x<<" Conc: "<<vspherocyl*x/(L*L*L)<<endl;
      double Tstar= 1/2.169;
      double lambda=0.05;
      
     //getchar();
         //impenetrable_spheres (shellin, shellfin,numshell, percVect , diamCond, x, L, condVect, Tstar, lambda);
      impenetrable_spherocylinder(numshell, percVect , diamCond, length1, x, L, condVect, Tstar, lambda, suffix);
      
      percFile<<x<<" "<< L<<" "<<phiVect[i]<<" "<<length1<<" ";
      condFile<<x<<" "<< L<<" "<<phiVect[i]<<" "<<length1<<" ";
       cerr<< L<<" "<<phiVect[i]<<" ";
       
       
        for (int l=0; l<REALIZ; l++){
              cerr<< 2*percVect[l]<<" ";
              percFile<<2*percVect[l]<<" ";
              condFile<<condVect[l]<<" ";

       }
       percFile<<endl;
       cerr <<endl;
       condFile<< endl;
       
}


      

DestroyDoubleVector(percVect);
DestroyDoubleVector(condVect);
condFile.close();
percFile.close();

//cout  << "Exit..."<<endl;
return 0;

}

