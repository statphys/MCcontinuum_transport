#include "mylib.h"

void penetrable_spheres (double shell,double diam, long int Nsphere, double Lato){

   int Nsquare =(int)(0.98 *Lato/(diam+shell));
   double Lsquare = Lato/Nsquare;
   long double randc=0.0;
   cerr << "Nsquare = "<< Nsquare << " Lsquare = " << Lsquare <<endl;
   cerr << "Nsquare*Lsquare = " <<  Nsquare*Lsquare << endl;
   sgenrand(time(NULL));

   double  **pos;
   pos = allocate_double_matrix (Nsphere+2, 6);
     for (int l=0; l<Nsphere+2; l++)
         for (int k=0; k<6; k++)
            pos[l][k]=0.0;

   ofstream myfile;
   myfile.open ("penetrable.dat");

   clock_t start_3 =start_time();


   for ( long int  i=1 ; i <= Nsphere; i++ ){

         pos[i][0]=(genrand())*Lato;
         pos[i][1]=(genrand())*Lato;
         pos[i][2]=(genrand())*Lato;

         myfile << pos[i][0] << " ";
         myfile << pos[i][1] << " ";
         myfile << pos[i][2] << " ";

         int Nx= (int)(pos[i][0]/Lsquare) +1;
         int Ny= (int)(pos[i][1]/Lsquare) +1;
         int Nz= (int)(pos[i][2]/Lsquare) +1;

         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;

         pos[i][3]=(double)Nx;
         pos[i][4]=(double)Ny;
         pos[i][5]=(double)Nx;

         myfile<< endl;

   }

   getTime(start_3);
   myfile.close();
   DestroyDoubleMatrix(pos, Nsphere+2);

}

void impenetrable_spheres (double shellin, double shellfin, int numshell, double *percVect ,double diam, long int Nsphere, double Lato, double *condVect, double Tstar, double lambda){

   int Nsquare =(int)(0.98*Lato/(diam));
   int Nbin_gr=NBIN_GR, Nbin_Hp=NBIN_HP, realiz=REALIZ;
   double deltar_gr=DR_GR, deltar_Hp=DR_HP;
   int count_perc=0;
   double prob_perc;
   int *clusterRec;
   int * clusterHK;
   double *hist_gr;
   double *countVect;
   hist_gr=allocate_double_vector(Nbin_gr);
   for (int l=0; l<Nbin_gr; l++)
            hist_gr[l]=0;
   
   double **pos;
   pos = allocate_double_matrix (Nsphere+2, 6);
   //memset(pos , 0 , Nsphere*6 * sizeof(double));
   for (int l=0; l<Nsphere+2; l++)
         for (int k=0; k<6; k++)
            pos[l][k]=0.0;

   int *NpartInSquare;
   NpartInSquare=allocate_int_vector(Nsquare*Nsquare*Nsquare);
   //memset(NpartInSquare , 0 , Nsquare*Nsquare*Nsquare * sizeof(* NpartInSquare));
   for (int l=0; l<Nsquare*Nsquare*Nsquare; l++)
            NpartInSquare[l]=0;

   int **PartInSquarePos;
   PartInSquarePos=allocate_int_matrix(Nsquare*Nsquare*Nsquare,MAX_SPHERES);
   //memset(PartInSquarePos , 0 , Nsquare*Nsquare*Nsquare*MAX_SPHERES * sizeof(** PartInSquarePos));
   for (int l=0; l<Nsquare*Nsquare*Nsquare; l++)
         for (int k=0; k<MAX_SPHERES; k++)
            PartInSquarePos[l][k]=0;


   int** connect1;
   double shell=1;
   double conc=(PI/6)*Nsphere/(Lato*Lato*Lato);
   //***** Place Particles *****//
   //clock_t start_3 =start_time();
   
   if (conc<=0.32){
     // cerr<<"Conc= "<<conc<<endl;
      //getchar();
      place_impenetrable_sphere(shell, diam, Nsphere, Lato, pos, NpartInSquare, PartInSquarePos);
      MC_impenetrable_pure(shell, diam, Nsphere, Lato,  pos, NpartInSquare, PartInSquarePos,connect1);
      cerr<<"End Pure MC"<<endl;
   }
   else{
      place_impenetrable_lattice(shell, diam, Nsphere, Lato, pos, NpartInSquare, PartInSquarePos);
      MC_impenetrable_pure(shell, diam, Nsphere, Lato,  pos, NpartInSquare, PartInSquarePos,connect1);
       cerr<<"End Pure MC"<<endl;
   }
   
   //getTime(start_3);
   cerr<<"Placed "<<Nsphere <<" cond spheres"<<endl;
   //print_filePos_impenetrable_beforeMC(Nsphere,  pos);
   
   
   //***** Equilibrate system *****//


   double passo=0;
   //clock_t start_4 =start_time();
   MC_impenetrable_sphere(shell, diam, Nsphere, Lato,  pos, NpartInSquare, PartInSquarePos,connect1, Tstar,  lambda, NSWEEP_SW_IN, passo);
   cerr<<"End SW MC after "<<NSWEEP_SW_IN <<" sweep "<<endl;
   

 for (int step=1; step<=realiz; step++){
   //cerr<<endl<<"step realization "<< step<<" ...................."<< endl;

   MC_impenetrable_sphere(shell, diam, Nsphere, Lato,  pos, NpartInSquare, PartInSquarePos,connect1, Tstar,  lambda, NSWEEP_SW, passo); 
   cerr<<"End SW MC after "<<NSWEEP_SW <<" sweep "<<endl;



   int** connect;
   connect=allocate_int_matrix(Nsphere,2);
   for (int l=0; l<Nsphere; l++)
         for (int k=0; k<2; k++)
            connect[l][k]=0;

   int *clusterRec;
   clusterRec=allocate_int_vector(Nsphere+2);
   for (int l=0; l<Nsphere+2; l++)
            clusterRec[l]=0;

   double conductance=0;

   ofstream pippo;
   //pippo.open("pipposotto.txt");
   
   

   cerr<<"FinishMC at step "<<step<< endl;
   //PrintAdjList(NsphereCond,AdjList);
   update_hist_gr(pos, Nsphere, hist_gr, Lato, deltar_gr, Nbin_gr);

  
   double s=0;
   double deltac_D=0;
   if (conc<0.005){
     s=3.5;
   }   
   else if (conc>=0.005 && conc<0.01){
       s=2.2;
   }
   else  if (conc>=0.01 && conc<0.1){
      s=1.4;
      }
   else {

      deltac_D=1.65*(1-conc)*(1-conc)*(1-conc)/(12*conc*(2-conc));
      //cerr<<"deltac == "<<deltac_D<<endl;
      deltac_D=deltac_D*2*diam/XI;
      //cerr<<"deltac == "<<deltac_D<<endl;
      deltac_D=log10(exp(-deltac_D))-6;
      //cerr<<"deltac == "<<deltac_D<<endl;
      deltac_D=log(pow(10,deltac_D));
      //cerr<<"deltac == "<<deltac_D<<endl;
      deltac_D=deltac_D*XI/(2*diam);
      s=-deltac_D/2;
        
      
   }
   //cerr<<"Conc  == "<<conc<<endl; 
   
  int span_cluster_index=0;
  
   int ** AdjListFirst=createAdj(s, diam, Nsphere, Lato, pos, connect);
   
   ClusterRec(clusterRec, Nsphere, AdjListFirst, pippo);
   
   span_cluster_index=percolation(Nsphere, connect, clusterRec, percVect, 1, realiz);
   DestroyIntMatrix(AdjListFirst, Nsphere+2);
   shell=s;
   if (span_cluster_index==0){
       double shellTemp=s;
       while (span_cluster_index==0){

          for (int l=0; l<Nsphere+2; l++)
             clusterRec[l]=0;
          
          shellTemp=shellTemp+s;
          
          int ** AdjListSecond=createAdj(shellTemp, diam, Nsphere, Lato, pos, connect);
          
          ClusterRec(clusterRec, Nsphere, AdjListSecond, pippo);
          span_cluster_index=percolation(Nsphere, connect, clusterRec, percVect, 1, realiz);
          DestroyIntMatrix(AdjListSecond, Nsphere+2);
        cerr<<"shell incremented = "<<shellTemp<<endl;

       }

      shell=shellTemp;

   }
   if (span_cluster_index!=0){
  double shellnew;
  double shell_fin=shell, shell_in=0;
   while( (shell_fin-shell_in)/shell_in>=EPSLON_SHELL )
   {
       
       for (int l=0; l<Nsphere+2; l++)
            clusterRec[l]=0;

       shellnew=(shell_fin+shell_in)/2;
      // cerr<<"shell actual"<<shellnew<<endl;
       int ** AdjList=createAdj(shellnew, diam, Nsphere, Lato, pos, connect);
       ClusterRec(clusterRec, Nsphere, AdjList, pippo);
       span_cluster_index=percolation (Nsphere, connect, clusterRec, countVect, 1, realiz);
       DestroyIntMatrix(AdjList, Nsphere+2);
       
       if (span_cluster_index== 0){
          shell_in=shellnew;
          //cerr<<"shell in = "<<shell_in<<endl;
        
       }

       else {
         shell_fin=shellnew;
         //cerr<<"shell finale = "<<shell_fin<<endl;
         
       }
   }

   double shelltrue = shellnew+EPSLON_SHELL*shell_in;
   cerr<<"D_c/2 = "<<shelltrue<<endl;
    double min_cond=exp(-2*2*shellnew/XI);
    //cerr<<" min  cond = "<<min_cond<<endl;
    //cerr<<" log min  cond = "<<log10(min_cond)<<endl;
     int ordine;
    if (log10(min_cond)<=-50)
          ordine=15;
    else if (log10(min_cond)>-50 && log10(min_cond)<=-10)
           ordine=10;
    else if (log10(min_cond)>-10) ordine=5;



      //cerr<<(log10(min_cond)-ordine)<<endl;
    shellnew=-(log10(min_cond)-ordine);
    shellnew=log(pow(10, shellnew)); 
    shellnew=shellnew*XI/4; 
    
    
    
     min_cond=exp(-2*2*shellnew/XI);

     cerr<<"Shell Augmented "<<shellnew<<endl;
     cerr<<" min  cond = "<<min_cond<<endl;
  

    double toll;

    double toll1=log10(exp(-2*shellnew*2*diam/XI));
    toll1=toll1/2;
    //cerr<<toll1<<endl;
    toll=pow(10, toll1) ;
    
    //cerr<<"toll "<<toll<<endl;
    

    

   for (int l=0; l<Nsphere+2; l++)
         clusterRec[l]=0;

   int ** AdjList=createAdj(shellnew, diam, Nsphere, Lato, pos, connect);
   ClusterRec(clusterRec, Nsphere, AdjList, pippo);
   span_cluster_index=percolation (Nsphere, connect, clusterRec, countVect, 1, realiz);
    
   vector <double>  cond_tip;


     
        //cerr<<"*********************** CONDUCTION **********************"<<endl;
        //clock_t start_5=start_time();
        // calcConductivity(diam, Nsphere, Lato, pos, clusterRec, AdjList, connect, span_cluster_index, &conductance, toll, min_cond);
        //getTime(start_5);
        
        condVect[step-1]=conductance;
        percVect[step-1]=shelltrue;
        //cerr<<"*********************** END  CONDUCTION **********************"<<endl;
         DestroyIntMatrix(AdjList, Nsphere+2);
     
   }
   /*else{
         percVect[step-1]=0;
         condVect[step-1]=0;
     
     }*/
   

   //pippo.close();
  
  // DestroyIntMatrix(AdjList, Nsphere+2);
   DestroyIntMatrix(connect, Nsphere);
   DestroyIntVector(clusterRec);
   
   //DestroyIntVector(clusterHK);
  }

  normal_write_gr(Nsphere, hist_gr, Lato, deltar_gr, Nbin_gr, realiz);
  //normal_write_Hp(Nsphere, hist_Hp, deltar_Hp, Nbin_Hp, realiz);
  
   DestroyDoubleMatrix(pos, Nsphere+2);
   DestroyIntMatrix(PartInSquarePos, Nsquare*Nsquare*Nsquare);
   DestroyIntVector(NpartInSquare);
   DestroyDoubleVector(hist_gr);
  //DestroyDoubleVector(hist_Hp);

  
}

void  segregation (double shellin, double shellfin, int numshell, double *percVect, double diamCond, double diamIns, long int NsphereCond,long int NsphereIns, double Lato){

   //getchar();
   int NsquareCond =(int)(0.98 *Lato/(diamCond));
   int NsquareIns=(int)(0.98 *Lato/(diamCond+diamIns));
   int realiz=REALIZ;

   ofstream pippo;
   //pippo.open("pipposotto.txt");

   double count_perc=0;
   double prob_perc;
   //int * clusterHK;
   int *clusterRec;

   int NsphereCond_new;


  for (int step=1; step<=realiz; step++){
    //cerr<<endl<<"step realization "<< step<<" ...................."<< endl;

   double  **posCond;
   posCond = allocate_double_matrix (NsphereCond+2, 6);
   //memset(pos , 0 , Nsphere*6 * sizeof(double));
   for (int l=0; l<NsphereCond+2; l++)
         for (int k=0; k<6; k++)
            posCond[l][k]=0.0;

   double  **posIns;
   posIns = allocate_double_matrix (NsphereIns+2, 6);
   //memset(pos , 0 , Nsphere*6 * sizeof(double));
   for (int l=0; l<NsphereIns+2; l++)
         for (int k=0; k<6; k++)
            posIns[l][k]=0.0;


   int *NpartInSquareCond;
   NpartInSquareCond=allocate_int_vector(NsquareCond*NsquareCond*NsquareCond);
   for (int l=0; l<NsquareCond*NsquareCond*NsquareCond; l++)
            NpartInSquareCond[l]=0;

   int *NpartInSquareIns;
   NpartInSquareIns=allocate_int_vector(NsquareIns*NsquareIns*NsquareIns);
   for (int l=0; l<NsquareIns*NsquareIns*NsquareIns; l++)
            NpartInSquareIns[l]=0;

   int **PartInSquarePosCond;
   PartInSquarePosCond=allocate_int_matrix(NsquareCond*NsquareCond*NsquareCond,MAX_SPHERES);
   //memset(PartInSquarePos , 0 , Nsquare*Nsquare*Nsquare*MAX_SPHERES * sizeof(** PartInSquarePos));
   for (int l=0; l<NsquareCond*NsquareCond*NsquareCond; l++)
         for (int k=0; k<MAX_SPHERES; k++)
            PartInSquarePosCond[l][k]=0;

   int **PartInSquarePosIns;
   PartInSquarePosIns=allocate_int_matrix(NsquareIns*NsquareIns*NsquareIns,MAX_SPHERES);
   //memset(PartInSquarePos , 0 , Nsquare*Nsquare*Nsquare*MAX_SPHERES * sizeof(** PartInSquarePos));
   for (int l=0; l<NsquareIns*NsquareIns*NsquareIns; l++)
         for (int k=0; k<MAX_SPHERES; k++)
            PartInSquarePosIns[l][k]=0;

   int ** AdjList;
   AdjList=allocate_int_matrix(NsphereCond+2, MAX_SPHERES);
   for (int l=0; l<NsphereCond+2; l++)
         for (int k=0; k<MAX_SPHERES; k++)
            AdjList[l][k]=0;

   int** connect;
   connect=allocate_int_matrix(NsphereCond,2);
   for (int l=0; l<NsphereCond; l++)
         for (int k=0; k<2; k++)
            connect[l][k]=0;

    int *clusterRec;
   clusterRec=allocate_int_vector(NsphereCond+2);
   for (int l=0; l<NsphereCond+2; l++)
            clusterRec[l]=0;

   int *clusterHK;
   clusterHK=allocate_int_vector(NsphereCond+2);
   for (int l=0; l<NsphereCond+2; l++)
            clusterHK[l]=0;



   clock_t start_3 =start_time();

   double shell=1;
  place_segregation(shell, diamCond, diamIns, NsphereCond, NsphereIns, Lato, AdjList, posCond, posIns, NpartInSquareCond, NpartInSquareIns, PartInSquarePosCond, PartInSquarePosIns);

  //  place_segregation_lattice(shell, diamCond, diamIns, NsphereCond, NsphereIns, Lato, AdjList, posCond, posIns, NpartInSquareCond, NpartInSquareIns, PartInSquarePosCond, PartInSquarePosIns);

  getTime(start_3);
  cerr<<"Placed "<<NsphereCond <<" cond spheres and "<<NsphereIns<<" Ins spheres"<<endl;

   clock_t start_4 =start_time();

   MC_segregation(shell, diamCond, diamIns, NsphereCond, NsphereIns, Lato, AdjList, posCond, posIns, NpartInSquareCond, NpartInSquareIns, PartInSquarePosCond, PartInSquarePosIns, connect);


   getTime(start_4);
   DestroyIntMatrix(PartInSquarePosCond, NsquareCond*NsquareCond);
   DestroyIntVector(NpartInSquareCond);

   DestroyDoubleMatrix(posIns, NsphereIns+2);
   DestroyIntMatrix(PartInSquarePosIns, NsquareIns*NsquareIns);
   DestroyIntVector(NpartInSquareIns);

   cerr<<"FinishMC"<<endl;
   //PrintAdjList(NsphereCond,AdjList);


   double x;
   double point=((shellfin-shellin)/numshell);

   for(int s=0; s<=numshell; s++){
      x=shellin+s*point;

      for (int l=0; l<NsphereCond+2; l++)
            clusterRec[l]=0;

      /*for (int l=0; l<NsphereCond+2; l++)
            clusterHK[l]=0; */
      cerr<<"Actual shell = "<<x<<endl;

      //createAdj(x, diamCond, NsphereCond, Lato, posCond, connect, AdjList);
      //PrintAdjList(NsphereCond,AdjList);
      //getchar();

      //PrintAdjList(NsphereCond,AdjList);
      //getchar();
     //***** Find Clusters *****//
     clock_t start_hk =start_time();
     start_hk =start_time();
     //cerr<<"###### Run HK ######"<<endl;
     ClusterHK(clusterHK, NsphereCond, AdjList, pippo);
     getTime(start_hk);
     //PrintAdjList(NsphereCond,AdjList);
      //getchar();
     //percolation (NsphereCond_new, connect, clusterHK, percVect, s, realiz);
     //createAdj2(x, diamCond, NsphereCond, Lato, posCond, connect, AdjList);

     clock_t start_c =start_time();
     start_c =start_time();
     cerr<<"###### Run Rec ######"<<endl;
     ClusterRec(clusterRec, NsphereCond, AdjList, pippo);
     getTime(start_c);
      //PrintAdjList(NsphereCond,AdjList);
      //getchar();

     //update_hist_gr(pos, Nsphere, hist_gr, Lato, deltar_gr, Nbin_gr);
     //update_hist_Hp(pos, AdjList, Nsphere, Lato, hist_Hp, deltar_Hp, Nbin_Hp);
     percolation (NsphereCond, connect, clusterRec, percVect, s, realiz);
     //cerr<<"step "<<step<<endl;
     /*for (int l=0; l<numshell; l++)
            cerr<< percVect[l]<<" ";


    cerr<<endl;*/

   }


   // get free mamory
   DestroyDoubleMatrix(posCond, NsphereCond+2);

   DestroyIntVector(clusterRec);
   DestroyIntVector(clusterHK);

   DestroyIntMatrix(AdjList, NsphereCond+2);
   DestroyIntMatrix(connect, NsphereCond);
   //pippo.close();

  }


  //cerr<<"count perc= "<<count_perc<<endl;
  //cerr<<"---------------------------------------------Prob =="<<prob_perc<<endl;

   for(int i=0; i<=numshell; i++)
     percVect[i]=percVect[i]/realiz;


}




void impenetrable_spherocylinder(int numshell, double *percVect ,double diam, double length, long int Nsphere, double Lato, double *condVect, double Tstar, double lambda, string suffix){

   int Nsquare =(int)(0.99*Lato/(length+diam));
   if (Nsquare<2){
      cerr<<"Too few squares "<<endl;
      exit(1);
      
   }
   
   int Nbin_gr=NBIN_GR, Nbin_Hp=NBIN_HP, Nbin_gr_angle=NBIN_GR_ANGLE, Nbin_mindist=NBIN_MINDIST, realiz=REALIZ;
   double deltar_gr=DR_GR, deltar_Hp=DR_HP, deltar_gr_angle=DR_GR_ANGLE, deltar_mindist=DR_MINDIST;
   int count_perc=0;
   double prob_perc;
   int *clusterRec;
   int * clusterHK;
   double *hist_gr, *hist_gr_angle, *hist_mindist;
   double *countVect;
   
   /*hist_gr=allocate_double_vector(Nbin_gr);
   for (int l=0; l<Nbin_gr; l++)
            hist_gr[l]=0;
   
   hist_gr_angle=allocate_double_vector(Nbin_gr_angle);
   for (int l=0; l<Nbin_gr_angle; l++)
            hist_gr_angle[l]=0;*/
   
   hist_mindist=allocate_double_vector(Nbin_mindist);
   for (int l=0; l<Nbin_mindist; l++)
            hist_mindist[l]=0;
            
            
   double **pos;
   pos = allocate_double_matrix (Nsphere+2, 6);
   //memset(pos , 0 , Nsphere*6 * sizeof(double));
   for (int l=0; l<Nsphere+2; l++)
         for (int k=0; k<6; k++)
            pos[l][k]=0.0;

   double **direct;
   direct = allocate_double_matrix (Nsphere+2, 3);
   //memset(pos , 0 , Nsphere*6 * sizeof(double));
   for (int l=0; l<Nsphere+2; l++)
         for (int k=0; k<3; k++)
            direct[l][k]=0.0;

   int *NpartInSquare;
   NpartInSquare=allocate_int_vector(Nsquare*Nsquare*Nsquare);
   //memset(NpartInSquare , 0 , Nsquare*Nsquare*Nsquare * sizeof(* NpartInSquare));
   for (int l=0; l<Nsquare*Nsquare*Nsquare; l++)
            NpartInSquare[l]=0;

   int **PartInSquarePos;
   PartInSquarePos=allocate_int_matrix(Nsquare*Nsquare*Nsquare,MAX_SPHERES);
   //memset(PartInSquarePos , 0 , Nsquare*Nsquare*Nsquare*MAX_SPHERES * sizeof(** PartInSquarePos));
   for (int l=0; l<Nsquare*Nsquare*Nsquare; l++)
         for (int k=0; k<MAX_SPHERES; k++)
            PartInSquarePos[l][k]=0;


   int** connect1;
   double shell=1;
   double vspherocyl=PI*(length*diam*diam/4 + diam*diam*diam/6);
   double conc=vspherocyl*Nsphere/(Lato*Lato*Lato);
   
   //***** Place Particles *****//
   //clock_t start_3 =start_time();
   cerr<<"CIAOOOOO conc "<<conc<<endl;
   
      place_impenetrable_spherocylinder(shell, diam,length,  Nsphere, Lato, pos, direct,  NpartInSquare,  PartInSquarePos);
    //print_filePos_impenetrable_beforeMC(Nsphere,  pos);
    //print_fileDir_impenetrable_beforeMC(Nsphere,  direct);
   
   //getchar();
   double passo=0, ds=0;
     MC_spherocyl_pure(diam, Nsphere, Lato, length, pos, direct, NpartInSquare, PartInSquarePos, NSWEEP_SCY_IN, passo, ds);
     cerr<<"End SW MC after "<<NSWEEP_SCY_IN <<" sweep "<<endl;
   //getTime(start_3);
   cerr<<"Placed "<<Nsphere <<" cond spheres"<<endl;
   //print_filePos_impenetrable_afterMC(Nsphere,  pos);
   //print_fileDir_impenetrable_afterMC(Nsphere,  direct);
   
   //***** Equilibrate system *****//


 
 for (int step=1; step<=realiz; step++){
   //cerr<<endl<<"step realization "<< step<<" ...................."<< endl;
   MC_spherocyl_pure(diam, Nsphere, Lato, length, pos, direct, NpartInSquare, PartInSquarePos, NSWEEP_SCY, passo, ds);
  // MC_impenetrable_sphere(shell, diam, Nsphere, Lato,  pos, NpartInSquare, PartInSquarePos,connect1, Tstar,  lambda, NSWEEP_SW, passo); 
  

   int** connect;
   connect=allocate_int_matrix(Nsphere,2);
   for (int l=0; l<Nsphere; l++)
         for (int k=0; k<2; k++)
            connect[l][k]=0;

   int *clusterRec;
   clusterRec=allocate_int_vector(Nsphere+2);
   for (int l=0; l<Nsphere+2; l++)
            clusterRec[l]=0;

   double conductance=0;

   ofstream pippo;
   //pippo.open("pipposotto.txt");
   
   

   cerr<<"FinishMC at step "<<step<< endl;
   //PrintAdjList(NsphereCond,AdjList);
   //update_hist_gr(pos, Nsphere, hist_gr, Lato, deltar_gr, Nbin_gr);
   //update_hist_gr_angle(direct, pos, Nsphere, hist_gr_angle, Lato, deltar_gr_angle, Nbin_gr_angle);
   update_hist_gr_mindist(direct, pos, Nsphere, length, hist_mindist, Lato, deltar_mindist, Nbin_mindist);
   
   double s=diam/2;
   double deltac_D=0;
   
   
  int span_cluster_index=0;
  
   int ** AdjListFirst=createAdjcylinder(s, diam, length, Nsphere, Lato, pos, direct, connect);
   
   ClusterRec(clusterRec, Nsphere, AdjListFirst, pippo);
   
   span_cluster_index=percolation(Nsphere, connect, clusterRec, percVect, 1, realiz);
   DestroyIntMatrix(AdjListFirst, Nsphere+2);
   shell=s;
   if (span_cluster_index==0){
       double shellTemp=s;
       while (span_cluster_index==0){

          for (int l=0; l<Nsphere+2; l++)
             clusterRec[l]=0;
          
          shellTemp=shellTemp+s;
          
          int ** AdjListSecond=createAdjcylinder(shellTemp, diam,length, Nsphere, Lato, pos,direct, connect);
          
          ClusterRec(clusterRec, Nsphere, AdjListSecond, pippo);
          span_cluster_index=percolation(Nsphere, connect, clusterRec, percVect, 1, realiz);
          DestroyIntMatrix(AdjListSecond, Nsphere+2);
        cerr<<"shell incremented = "<<shellTemp<<endl;

       }

      shell=shellTemp;

   }
   if (span_cluster_index!=0){
  double shellnew;
  double shell_fin=shell, shell_in=0;
   while( (shell_fin-shell_in)/shell_in>=EPSLON_SHELL )
   {
       
       for (int l=0; l<Nsphere+2; l++)
            clusterRec[l]=0;

       shellnew=(shell_fin+shell_in)/2;
      // cerr<<"shell actual"<<shellnew<<endl;
       int ** AdjList=createAdjcylinder(shellnew, diam,length, Nsphere, Lato, pos, direct, connect);
       ClusterRec(clusterRec, Nsphere, AdjList, pippo);
       span_cluster_index=percolation (Nsphere, connect, clusterRec, countVect, 1, realiz);
       DestroyIntMatrix(AdjList, Nsphere+2);
       
       if (span_cluster_index== 0){
          shell_in=shellnew;
          //cerr<<"shell in = "<<shell_in<<endl;
        
       }

       else {
         shell_fin=shellnew;
         //cerr<<"shell finale = "<<shell_fin<<endl;
         
       }
   }

   double shelltrue = shellnew+EPSLON_SHELL*shell_in;
   cerr<<"D_c = "<<2*shelltrue<<endl;
   //getchar();
   
    double min_cond=exp(-2*2*shelltrue/XI);
    double condCPA=min_cond;
    cerr<<" min  cond = "<<min_cond<<endl;
    //cerr<<" log min  cond = "<<log10(min_cond)<<endl;
     int ordine;
    if (log10(min_cond)<=-50)
          ordine=20;
    else if (log10(min_cond)>-50 && log10(min_cond)<=-10)
           ordine=20;
    else if (log10(min_cond)>-10) ordine=5;



      //cerr<<(log10(min_cond)-ordine)<<endl;
    shellnew=-(log10(min_cond)-ordine);
    shellnew=log(pow(10, shellnew)); 
    shellnew=shellnew*XI/4; 
    
    
    
     min_cond=exp(-2*2*shellnew/XI);

     cerr<<"Shell Augmented "<<2*shellnew<<endl;
     cerr<<" min  cond = "<<min_cond<<endl;
  

    double toll;

    double toll1=log10(exp(-2*shellnew*2*diam/XI));
    toll1=toll1/2;
    //cerr<<toll1<<endl;
    toll=pow(10, toll1) ;
    
    //cerr<<"toll "<<toll<<endl;
    

    

   for (int l=0; l<Nsphere+2; l++)
         clusterRec[l]=0;

   int ** AdjList=createAdjcylinder(shellnew, diam,length, Nsphere, Lato, pos, direct, connect);
   ClusterRec(clusterRec, Nsphere, AdjList, pippo);
   span_cluster_index=percolation (Nsphere, connect, clusterRec, countVect, 1, realiz);

    //if (span_cluster_index!=0) cerr<<"Ho percolato"<<endl;
  //PrintAdjList(Nsphere,AdjList);
   //getchar();
     
        //cerr<<"*********************** CONDUCTION **********************"<<endl;
        //clock_t start_5=start_time();
     
    calcConductivityCyl(diam,  length, Nsphere, Lato, pos, direct, clusterRec, AdjList, connect,  span_cluster_index, &conductance, toll, condCPA);
        //getTime(start_5);
        
    condVect[step-1]=conductance;
        
        //cerr<<"*********************** END  CONDUCTION **********************"<<endl;
    DestroyIntMatrix(AdjList, Nsphere+2);
    percVect[step-1]=shelltrue;
   }
   //else{
        // percVect[step-1]=0;
      //   condVect[step-1]=0;
     
     //}
   

   //pippo.close();
  
  // DestroyIntMatrix(AdjList, Nsphere+2);
   DestroyIntMatrix(connect, Nsphere);
   DestroyIntVector(clusterRec);
   
   //DestroyIntVector(clusterHK);
  }
  normal_write_gr_mindist(Nsphere, length, hist_mindist, Lato, deltar_mindist, Nbin_mindist, realiz, suffix);
  //normal_write_cc(Nsphere, length, hist_gr, Lato, deltar_gr, Nbin_gr, realiz);
  //normal_write_gr_angle(Nsphere, length, hist_gr_angle, Lato, deltar_gr_angle, Nbin_gr_angle, realiz);
   
 cerr<<"CIOAAAAAA"<<endl;
 
   DestroyDoubleMatrix(pos, Nsphere+2);  
   DestroyDoubleMatrix(direct, Nsphere+2);
   DestroyIntMatrix(PartInSquarePos, Nsquare*Nsquare*Nsquare);
   DestroyIntVector(NpartInSquare);
   //DestroyDoubleVector(hist_gr);
   //DestroyDoubleVector(hist_gr_angle);
   DestroyDoubleVector(hist_mindist);
  
}





clock_t start_time(){

   clock_t start ;
   start = clock() ;
   return start;
}

void getTime ( clock_t start )
{
    int exec_time_us;
    clock_t stop;

    stop=clock();
    exec_time_us = (stop - start) * 1000/CLOCKS_PER_SEC ;

    //printf( "Execution time = %d milliseconds\n", exec_time_us ) ;
    if (exec_time_us <1000)
       cerr<<"Execution time =" << (double)exec_time_us<< "millisecons"<< endl;
       //printf( "Execution time = %f millisecons\n", (double)exec_time_us );
    else if (exec_time_us >1000 && exec_time_us < 60000)
       cerr<<"Execution time =" << (double)exec_time_us/1000<< "seconds"<< endl;
       //printf( "Execution time = %f seconds\n", (double)exec_time_us/1000 ) ;
    else
       cerr<<"Execution time =" << (double)exec_time_us/1000/60<< "minutes"<< endl;
       //printf( "Execution time = %f minutes\n", (double)exec_time_us/1000/60 ) ;
}

double ** allocate_double_matrix(int row, int col)
{


         int i,j,n;
         double **mat;
       mat =  (double **)calloc(row, sizeof(double *));
       mat[0] = (double*) calloc(row*col, sizeof(double));



         for( i=1; i<row; i++)
           mat[i] = &mat[0][i*col];


        if (mat==NULL)
            cerr<< "Unable to allocate double matrix memory"<<endl;

        return mat;
}

int ** allocate_int_matrix(int row, int col)
{

         int i,j,n;
         int **mat;
       mat =  (int**)calloc(row, sizeof(double *));
       mat[0] = (int*) calloc(row*col, sizeof(double));



         for( i=1; i<row; i++)
           mat[i] = &mat[0][i*col];


        if (mat==NULL)
            cerr<< "Unable to allocate double matrix memory"<<endl;

        return mat;
}

double * allocate_double_vector(int n){

        double *v;

        v = (double*)calloc(n, sizeof(double));
        if (v==NULL)
            cerr<< "Unable to allocate double vector memory"<<endl;

        return v;


}

int * allocate_int_vector(int n){

        int *v;

        v = (int*)calloc(n, sizeof(int));
        if (v==NULL)
            cerr<< "Unable to allocate int vector memory"<<endl;

        return v;


}

void DestroyDoubleMatrix(double ** mat, int nrow)
{
	int i;

	if (mat != NULL) {
	   for( i=1; i<nrow; i++)
               mat[i]=NULL;

           free(mat[0]);

           free(mat);

	}


      //cerr<< "Deallocating double matrix..." << endl;


}



void DestroyIntMatrix(int ** mat, int nrow)
{
	int i;

	if (mat != NULL) {
	   for( i=1; i<nrow; i++)
               mat[i]=NULL;

           free(mat[0]);
           free(mat);

	}


      //cerr<< "Deallocating int matrix..." << endl;


}

void DestroyDoubleVector(double * vet)
{
	if (vet != NULL) {
		free(vet);
	}
        //cerr<< "Deallocating double vector..." << endl;
}

void DestroyIntVector(int * vet)
{
	if (vet != NULL) {
		free(vet);
	}
        //cerr<< "Deallocating int vector..." << endl;
}



/* A C-program for MT19937: Real number version                */
/* genrand() generates one pseudorandom real number (double)   */
/* which is uniformly distributed on [0,1]-interval, for each  */
/* call. sgenrand(seed) set initial values to the working area */
/* of 624 words. Before genrand(), sgenrand(seed) must be      */
/* called once. (seed is any 32-bit integer except for 0).     */
/* Integer generator is obtained by modifying two lines.       */
/*   Coded by Takuji Nishimura, considering the suggestions by */
/* Topher Cooper and Marc Rieffel in July-Aug. 1997. Comments  */
/* should be addressed to: matumoto@math.keio.ac.jp            */

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
/* Tempering parameters */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y) (y >>  11)
#define TEMPERING_SHIFT_S(y) (y <<  7)
#define TEMPERING_SHIFT_T(y) (y <<  15)
#define TEMPERING_SHIFT_L(y) (y >>  18)
static unsigned long mt[N]; /* the array for the state vector */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */
/* initializing the array with a NONZERO seed */
void sgenrand(unsigned long seed)
{
    /* setting initial seeds to mt[N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}
double /* generating reals */
/* unsigned long */ /* for integer generation */
genrand()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A for x=0,1 */
    if (mti >= N) { /* generate N words at one time */
        int kk;
        if (mti == N+1) /* if sgenrand() has not been called, */
            sgenrand(4357); /* a default initial seed is used */
        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];
        mti = 0;
    }
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    return ( (double)y / (unsigned long)0xffffffff ); /* reals */
    /* return y; */ /* for integer generation */
}

void findClusterRec(int * cluster,  int Nsphere, int** AdjList, int col){

   int count=1;
   for (int i=1; i<=Nsphere; i++){
      //cerr<<"cluster["<<i<<"] = "<<cluster[i]<<endl;
      if (cluster[i]==0)
         cluster[i]=count;
         //cerr<<"Adesso cluster["<<i<<"] = "<<cluster[i]<<endl;
         growClusterRec(cluster, Nsphere, AdjList, col, i,count);
      count++;
   }
   //return cluster;
}

void growClusterRec(int * cluster,  int Nsphere, int** AdjList, int col, int x, int count){
    int a;
    int j=0;
    //cerr<<"entro in grow"<<endl;
    while(AdjList[x][j]!=0){
        //cerr<<"cluster["<<AdjList[x][j]<<"] = "<<cluster[AdjList[x][j]]<<endl;
        if (cluster[AdjList[x][j]]==0){
           //cerr<<"sono dentro l' if"<<endl;
           //cerr<<"AdjList["<<x<<"]["<<j<<"] = "<<AdjList[x][j]<<endl;
           cluster[AdjList[x][j]]=count;
           growClusterRec(cluster, Nsphere, AdjList, col, AdjList[x][j], count);
           j++;
          //getchar();
        }
       else j++;
        //getchar();
    }

}

void adjust_group(int * cluster, int Nsphere){


           int count=1;
           int numGroup= Ngroup_search(cluster, Nsphere) ;

           for(int i=1; i<=numGroup; i++){

                  if (bool_search(cluster, Nsphere, i)==FALSE){
                      //int max=*max_element(cluster+1, cluster+Nsphere+1);
                      int max = find_max_el(cluster, Nsphere);

                      for(int j=1; j<=Nsphere; j++)
                            if (cluster[j]==max)
                                cluster[j]=i;

                      count++;
                      }

            }


}

int find_max_el(int *cluster, int Nsphere){

   int max=0;
   for (int i=1; i<=Nsphere; i++)
      if (cluster[i]>max)
         max=cluster[i];

    return max;

}

int Ngroup_search(int *vet, int n){

        int count=1, max=vet[1];

        for (int i=2; i<=n; i++)
            if (vet[i]>vet[i-1] && vet[i]>max)
            {
                max=vet[i];
                count++;

            }
         return count;
}

int bool_search(int *vet, int n, int num){
     int i;

      for(i=1; i<=n; i++)
           if (vet[i]==num)
                 return TRUE;
      return FALSE;

}

void ClusterRec(int *clusterRec, int Nsphere, int **AdjList, std::ofstream& pippo){

   //ofstream clusterFile;

   //clusterFile.open("cluster3dRec.dat");
   //int *cluster= clusterRec;

   findClusterRec(clusterRec, Nsphere, AdjList, MAX_SPHERES);
   //for (int l=1; l<=Nsphere; l++)
            //pippo<<cluster[l]<<" ";
   //pippo<<endl;
   //int max=*max_element(clusterRec+1, clusterRec+Nsphere+1);

   adjust_group(clusterRec, Nsphere);

   /*for (int l=1; l<=Nsphere; l++){
            //pippo<<cluster[l]<<" ";
            clusterFile<<clusterRec[l]<<" ";
            }
   */
   //pippo<<endl;

   //DestroyIntVector(cluster);
   //clusterFile.close();


}


void PrintAdjList(int Nsphere, int **AdjList){

      ofstream AdjFile;

           AdjFile.open("Adj.txt");
           for(int ic=1; ic<=Nsphere; ic++){
               AdjFile<<ic<<" ";
               //pippo<< ic <<" -- ";
                for(int l=0; l<MAX_SPHERES; l++){
                  if (AdjList[ic][l]==0)
                         break;
                  //pippo<< AdjList[ic][l]<< " ";
                 AdjFile << AdjList[ic][l]<<" ";

              }
              AdjFile<<endl;
              //pippo<< endl;
             }
          AdjFile.close();

}

void ClusterHK(int *clusterHK, int Nsphere, int **AdjList, std::ofstream& pippo){
    ofstream clusterFile;

   //clusterFile.open("cluster3dHK.dat");
   int *Neigh, *NodeL=clusterHK;

   Neigh=allocate_int_vector(Nsphere+2);

   for (int l=0; l<Nsphere+2; l++){
            Neigh[l]=0;
           // NodeL[l]=0;
   }

   int count=0;
   for(int i=1; i<=Nsphere; i++){
      int m;

      for(m=0; m<MAX_SPHERES;m++){
         //pippo<<AdjList[i][m]<<" ";
         if (AdjList[i][m]==0)
              // pippo<<endl;
              break;
        }

      Neigh[i]=m;
      if (count_zeros(NodeL, i, Neigh[i], AdjList)==0){
         count++;
         //pippo<<Neigh[i]<<endl;
         NodeL[i]=count;
         //pippo<<"Inizio NodeL["<<i<<"] = "<< NodeL[i]<< endl;

        for(int n=0;n<Neigh[i];n++){
            NodeL[AdjList[i][n]]=count;
           // pippo<<"NodeL["<<AdjList[i][n]<<"] = "<< NodeL[i]<< " ";

          }
        //pippo<<endl;
        //getchar();
      }
      else{
         //pippo<<"Inizio gia etichettato NodeL["<<i<<"] = "<< NodeL[i]<< endl;
         count++;
         int minimum=count;
         for (int j=0; j<Neigh[i]; j++)
            if (NodeL[AdjList[i][j]]>0){
                  minimum = min(minimum, NodeL[AdjList[i][j]]);
                  //NodeL[AdjList[i][j]]=minimum;
              }

         //pippo<<"il minimo Ã© "<<minimum<<endl;
         NodeL[i]=minimum;
         for (int j=0; j<Neigh[i]; j++){
            int temp = NodeL[AdjList[i][j]];
            for(int lc=1; lc<=Nsphere; lc++){

                 if (NodeL[lc]==temp && NodeL[lc]!=0){
                     //pippo<<"node["<<lc<<"] = "<<NodeL[AdjList[i][j]]<<" ";
                     NodeL[lc]=minimum;
                     //pippo<<"adesso node["<<lc<<"] = "<<NodeL[lc]<<endl;
                 }
              //pippo<<"node["<<lc<<"] != "<<NodeL[AdjList[i][j]]<<endl;
             }
             NodeL[AdjList[i][j]]=minimum;

          }

        }

    // print_int_vect(NodeL, Nsphere, pippo );
     //getchar();

   }

   adjust_group(NodeL, Nsphere) ;

   for (int l=1; l<=Nsphere; l++){
            //pippo<<NodeL[l]<<" ";
            clusterFile<<NodeL[l]<<" ";
            }

   //DestroyIntVector(cluster);
   DestroyIntVector(Neigh);
   //DestroyIntVector(NodeL);
   //clusterFile.close();

   //return NodeL;
}


int count_zeros(int *vet, int x, int m, int **AdjList){

     int count=0;

     for (int j=0; j<m; j++)
        if (vet[AdjList[x][j]]!=0)
           count++;

     return count;


}

void print_int_vect(int * vect, int n,  std::ofstream& pippo ){

      for (int i=1; i<=n; i++)
           cerr<< vect[i]<<" ";
      cerr<<endl;


}

void print_filePos_impenetrable_beforeMC(long int Nsphere,  double **pos){

    ofstream myfile;
    myfile.open ("impenetrable_noMC.dat");

    for(int i=1; i<=Nsphere; i++){

         myfile<<pos[i][0]<<" ";
         myfile<<pos[i][1]<<" ";
         myfile<<pos[i][2]<<endl;

     }


    myfile.close();
}


void print_fileDir_impenetrable_beforeMC(long int Nsphere,  double **pos){

    ofstream myfile;
    myfile.open ("direct_noMC.dat");

    for(int i=1; i<=Nsphere; i++){

         myfile<<pos[i][0]<<" ";
         myfile<<pos[i][1]<<" ";
         myfile<<pos[i][2]<<endl;

     }


    myfile.close();
}


void print_filePos_impenetrable_afterMC(long int Nsphere,  double **pos){

    ofstream myfile;
    myfile.open ("impenetrable_dopoMC.dat");

    for(int i=1; i<=Nsphere; i++){

         myfile<<pos[i][0]<<" ";
         myfile<<pos[i][1]<<" ";
         myfile<<pos[i][2]<<endl;

     }


    myfile.close();
}

void print_fileDir_impenetrable_afterMC(long int Nsphere,  double **pos){

    ofstream myfile;
    myfile.open ("direct_dopoMC.dat");

    for(int i=1; i<=Nsphere; i++){

         myfile<<pos[i][0]<<" ";
         myfile<<pos[i][1]<<" ";
         myfile<<pos[i][2]<<endl;
       }


    myfile.close();
}
  

void place_impenetrable_sphere(double shell,double diam, long int Nsphere, double Lato, double **pos, int * NpartInSquare, int ** PartInSquarePos){


   int Nsquare =(int)(0.98*Lato/(diam));
   double Lsquare = Lato/Nsquare;
   long double randc=0.0;
   //cerr << "Nsquare = "<< Nsquare << " Lsquare = " << Lsquare <<endl;
   //cerr << "Nsquare*Lsquare = " <<  Nsquare*Lsquare << endl;

   sgenrand(time(NULL));

   ofstream myfile;//AdjFile;
   //myfile.open ("impenetrable3d.dat");


   int flag=0;
   for ( long int  i=1 ; i <= Nsphere; i++ ){
         long int numSquare=0;

         //getchar();
         //randc = (genrand())*Lato;
         //cerr<<endl<<"I'm here...i = "<<i<<endl;
         pos[i][0]=(genrand())*Lato;
         pos[i][1]=(genrand())*Lato;
         pos[i][2]=(genrand())*Lato;



         int Nx= (int)(pos[i][0]/Lsquare) +1;
         int Ny= (int)(pos[i][1]/Lsquare) +1;
         int Nz= (int)(pos[i][2]/Lsquare) +1;

         //debug

        //numSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);
        //int numSquareTemp=numSquare;

        //cerr<< "--- New position = ---"<< numSquare<< endl;
        //cerr<< "(x,y,z) = "<< "("<< pos[i][0]<<", "<<pos[i][1]<<", "<<pos[i][2]<<")";
         //se pos[i][1]==L ==> Nx=Nsquare i.e. choose the last little cube
         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;

         int count =0;
         for (int n1 =Nx-1; n1<=Nx+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
              for (int n2 =Ny-1; n2<=Ny+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =Nz-1; n3<=Nz+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }
                      count++;
                      //find position in numerized squares from 1 to Nsquare^3
                      numSquare=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)
                      if ( NpartInSquare[numSquare-1]!=0){

                        //cerr<<"numSquare = "<< numSquare<< endl;

                        //cerr<<"There are "<< NpartInSquare[numSquare-1] <<" neighbors in square = "<< numSquare <<endl;
                        for (long int j=1; j<=NpartInSquare[numSquare-1]; j++){

                           //find the index of the neighboring particle
                           int k=PartInSquarePos[numSquare-1][j];
                           //cerr << "Pos vicino k = "<< k<< endl;

                           double dist=0.0;
                           //find distance of part i from each particle in the around squares
                           dist=pow((pos[i][0]-pos[k][0]+Lato*k1),2)+pow((pos[i][1]-pos[k][1]+Lato*k2),2)
                                +pow((pos[i][2]-pos[k][2]+Lato*k3),2);

                           dist=sqrt(dist);
                           //cerr<<"---dist = " <<dist<<" diam = "<<diam<<endl;
                           if (dist<diam){ //if sphere i overlaps sphere k, go out
                              flag=1;
                              //cerr<<"Overlap..!!"<<endl;
                              //getchar();
                              break;
                              }
                          }
                      }
                     if (flag==1) break;

                   }//end Nz for
                 if (flag==1) break;
                }//end Ny for
               if (flag==1) break;
             }//end Nx for
         if (flag==1) {i--; flag=0; continue;};

         //cerr << "Count = "<< count;

         //cerr<<" Placing sphere...i = "<<i<<endl;
         //myfile << pos[i][0] << " ";
         //myfile << pos[i][1] << " ";
         //myfile << pos[i][2] << " ";

         pos[i][3]=(double)Nx;
         pos[i][4]=(double)Ny;
         pos[i][5]=(double)Nz;

         numSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);
         //cerr<<"Before there were N part = "<< NpartInSquare[numSquare-1] << " in square = "<< numSquare <<endl;
         NpartInSquare[numSquare-1]++;
         //cerr<<"Now there are N part = "<< NpartInSquare[numSquare-1] << " in square = "<< numSquare <<endl;
         if (NpartInSquare[numSquare-1]>MAX_SPHERES){
            cerr<<"Too many spheres in square # "<< numSquare<<endl;
            exit(1);
            }

         PartInSquarePos[numSquare-1][NpartInSquare[numSquare-1]]=i;

         //cerr<< "PartPos["<<numSquare-1<<"]["<< NpartInSquare[numSquare-1] <<"] = "<<
             //           PartInSquarePos[numSquare-1][NpartInSquare[numSquare-1]]<<endl;

        // myfile<< endl;


   } //end for i-Nsphere

   //myfile.close();
}

void place_impenetrable_lattice(double shell,double diam,  long int Nsphere,  double Lato, double **pos,  int * NpartInSquare, int ** PartInSquarePos){

    

   int c;
   int Nsquare =(int)(0.98 *Lato/(diam));
 
   double conc=Nsphere*(PI/6)*(diam)*(diam)*(diam)/(Lato*Lato*Lato);
   
   

   double Lsquare = Lato/Nsquare;
 
   srand(time(NULL));

   if (conc>=PI/6 ){
      cerr<<"conc eff ="<<conc<<endl;
     
      cerr<<"concentration is too high for a square Lattice"<<endl;
      exit(1);
   }
   
   //getchar();

   //Firstly Insulating particles are placed, saving position and cube position in PosIns
   int count=0;
   ofstream pippo1;
   //pippo1.open("segregation3d_Ins.dat");
   

    //Conducting particles are added, by saving position and cube position in PosCond
     //For each added parlicle check if it intersects the previous placed  conducting
     //particles
   
   int flag=0;
   for ( long int  i=1 ; i <= Nsphere; i++ ){
      try{
         
         long int numSquare=0;

         //cerr<<endl<<"...Placing Cond particle i = "<<i<<endl;


         pos[i][0]=(double)(rand()%((int)Lato))+diam/2;
         pos[i][1]=(double)(rand()%((int)Lato))+diam/2;
         pos[i][2]=(double)(rand()%((int)Lato))+diam/2;

         int Nx= (int)(pos[i][0]/Lsquare) +1;
         int Ny= (int)(pos[i][1]/Lsquare) +1;
         int Nz= (int)(pos[i][2]/Lsquare) +1;
    
        
         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;

         numSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);
         
         for (int n1 =Nx-1; n1<=Nx+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =Ny-1; n2<=Ny+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                  for (int n3 =Nz-1; n3<=Nz+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquare=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)
                      if ( NpartInSquare[numSquare-1]!=0){


                        for (long int j=1; j<=NpartInSquare[numSquare-1]; j++){
                          //Preparare aperitivo
                           //find the index of the neighboring particle
                           int k=PartInSquarePos[numSquare-1][j];

                           double dist=0.0;
                           //find distance of part i from each particle in the around squares
                        dist=pow((pos[i][0]-pos[k][0]+Lato*k1),2)+pow((pos[i][1]-pos[k][1]+Lato*k2),2)+
                                pow((pos[i][2]-pos[k][2]+Lato*k3),2);

                           dist=sqrt(dist);
                           //cerr<<"---dist = " <<dist<<" diam = "<<diam<<endl;
                           if (dist<diam){ //if sphere i overlaps sphere k, go out
                              flag=1;
                              /*cerr<<"---dist = " <<dist<<" diam = "<<diam<<endl;
                              cerr<<"Overlap..!!"<<endl;
                              getchar();*/
                              break;
                              }
                          }
                      }
                     if (flag==1) break;

                   }//end Nz for
                 if (flag==1) break;
                }//end Ny for
               if (flag==1) break;
             }//end Nx for
         if (flag==1) {i--; flag=0;   continue;};


         pos[i][3]=(double)Nx;
         pos[i][4]=(double)Ny;
         pos[i][5]=(double)Nz;
        
         numSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);

         NpartInSquare[numSquare-1]++;
         //cerr<<"Now there are N part = "<< NpartInSquare[numSquare-1] << " in square = "<< numSquare <<endl;

         if (NpartInSquare[numSquare-1]>MAX_SPHERES){
            cerr<<"Too many spheres in square # "<< numSquare<<endl;
            exit(1);
            }
        
         PartInSquarePos[numSquare-1][NpartInSquare[numSquare-1]]=i;
       // cerr<<"Placed "<<i<<endl;
         //myfile<< endl;
   }
     catch (exception& e)
  {
    cerr << "exception caught: " << e.what() << endl;
  }



   } //end for i-Nsphere
   


}

void place_impenetrable_lattice2(double shell,double diam, long int Nsphere, double Lato, double **pos, int * NpartInSquare, int ** PartInSquarePos){

    
   int Nsquare =(int)(0.98*Lato/(diam));
   double Lsquare = Lato/Nsquare;
   long double randc=0.0;
   cerr << "Nsquare = "<< Nsquare << " Lsquare = " << Lsquare <<endl;
   cerr << "Nsquare*Lsquare = " <<  Nsquare*Lsquare << endl;
   
   double conc=Nsphere*(PI/6)/(Lato*Lato*Lato);
   double num=1.0/3.0;
   double Lgrid=pow((PI/6)/conc,num);

   int NsquareGrid =(int)(Lato/(Lgrid));
   cerr<< "NsquareGrid = "<<NsquareGrid<<endl;
   cerr<< "Lgrid = "<<Lgrid<<endl;
   cerr << "NsquareGrid*Lgrid = " << Lgrid*NsquareGrid << endl;
   NsquareGrid++;                                                                              
   cerr << "NsquareGrid*Lgrid = " << Lgrid*NsquareGrid << endl;
   cerr<< "Nsphere = "<<Nsphere<<" part add = "<<NsquareGrid*NsquareGrid*NsquareGrid<<endl;
   //cerr<<Lgrid*NsquareGrid<<endl;   
  //getchar();
   cerr<<(Lato/NsquareGrid)<<endl;


   Lgrid=(Lato/NsquareGrid);
   cerr << "NsquareGrid*Lgrid = " <<  Lgrid*NsquareGrid << endl;
   cerr<<conc<<endl;
   cerr<<"NsquareGrid = "<< NsquareGrid<<endl;
   cerr<<"New Lgrrid"<<Lgrid<<endl;
   //cerr<<Nsphere<<endl;
   //cerr<<conc*Lato*Lato*Lato;
   
  /* NsquareGrid=(int)pow(conc*(Lato*Lato*Lato)*6/PI, num);
   
   cerr<< "----NsquareGrid = "<<NsquareGrid<<endl;
   
   Lgrid=(Lato/NsquareGrid);
   cerr<< "----Lgrid = "<<Lgrid<<endl;
   */
   if (conc>=PI/6 || Lgrid <= 1){
      cerr<<"Effective concentration is too high for a square Lattice"<<endl;
      
      exit(1);
   }
     
  
   int count=0;
   ofstream pippo1;
   srand(time(NULL));

     //Secondly Conducting particles are added, saving position and cube position in PosCond
     //For each added parlicle check if it intersects the previous placed insulating particles and the other conducting
     //particles
  
   int count_pos=1;
   for ( int ii=0 ; ii <= NsquareGrid-1; ii++ ){
      for ( int jj=0 ; jj <= NsquareGrid-1; jj++ ){
         for ( int kk=0 ; kk <= NsquareGrid-1; kk++ ){

              if (count_pos>Nsphere){
                 //cerr<<"SUPERATO Nsphere cpos = "<<count_pos<<endl;

                 
                 //choose a previously indexed sphere considering also the current particle
                 int partic=rand()%(Nsphere+1)+1;//choose random num between 1 and Nsphere+1

                 //cerr<<"scelgo la partic "<<partic<<endl;

                 if (partic == Nsphere+1)
                      break;
                 else {
                    /*delete old information*/
                    int oldSquare = 
                    +((int)pos[partic][3]+Nsquare*((int)pos[partic][4]-1)+Nsquare*Nsquare*((int)pos[partic][5]-1));

                    int k=0; 
                    for (long int j=1; j<=NpartInSquare[oldSquare-1]-1; j++){
                       k++;
                       if (PartInSquarePos[oldSquare-1][j]==partic) k++;
                       PartInSquarePos[oldSquare-1][j]=PartInSquarePos[oldSquare-1][k];
                    }

                    PartInSquarePos[oldSquare-1][NpartInSquare[oldSquare-1]]=0;
                    NpartInSquare[oldSquare-1]--;

                    /*update with new information*/
                  
/////////////////////////////////////////////////////////////////////////////
                    double x = ii*Lgrid;
                    double y = jj*Lgrid;
                    double z = kk*Lgrid;
  
                    pos[partic][0]=x;
                    pos[partic][1]=y;
                    pos[partic][2]=z;  
              
                    int Nx= (int)(pos[partic][0]/Lsquare) +1;
                    int Ny= (int)(pos[partic][1]/Lsquare) +1;
                    int Nz= (int)(pos[partic][2]/Lsquare) +1;

                    if (Nx==Nsquare+1) Nx=Nsquare;
                    if (Ny==Nsquare+1) Ny=Nsquare;
                    if (Nz==Nsquare+1) Nz=Nsquare;

                    pos[partic][3]=(double)Nx;
                    pos[partic][4]=(double)Ny;
                    pos[partic][5]=(double)Nz; 
                    
                    int newSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);   
                    NpartInSquare[newSquare-1]++;

                    if (NpartInSquare[newSquare-1]>MAX_SPHERES){
                       cerr<<"Too many spheres in square # "<< newSquare<<endl;
                       exit(1);     
                    }
      
                    PartInSquarePos[newSquare-1][NpartInSquare[newSquare-1]]=partic;


                 }


 
              }//end if
              else{     
                 double x = ii*Lgrid;
                 double y = jj*Lgrid;
                 double z = kk*Lgrid;
                 //cerr<<"(x,,y,z) = "<<x<<", "<<y<<", "<<z<<endl;
                 //getchar();
              
                 pos[count_pos][0]=x;
                 pos[count_pos][1]=y;
                 pos[count_pos][2]=z;  
              
                 int Nx= (int)(pos[count_pos][0]/Lsquare) +1;
                 int Ny= (int)(pos[count_pos][1]/Lsquare) +1;
                 int Nz= (int)(pos[count_pos][2]/Lsquare) +1;

                 if (Nx==Nsquare+1) Nx=Nsquare;
                 if (Ny==Nsquare+1) Ny=Nsquare;
                 if (Nz==Nsquare+1) Nz=Nsquare;
        

                pos[count_pos][3]=(double)Nx;
                pos[count_pos][4]=(double)Ny;
                pos[count_pos][5]=(double)Nz;
              
                int numSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);
         
                NpartInSquare[numSquare-1]++;
         
                if (NpartInSquare[numSquare-1]>MAX_SPHERES){
                    cerr<<"Too many spheres in square # "<< numSquare<<endl;
                    exit(1);     
                }
               // cerr<<"---Placed "<<i<<endl;
               PartInSquarePos[numSquare-1][NpartInSquare[numSquare-1]]=count_pos;
            }//end else
   
         count_pos++;
                                
   
         }//end ii  
         
      }//end jj
      
   }//end kk
    
  cerr<<"pos_cond "<< count_pos << endl;
   
}


void MC_impenetrable_pure(double shell,double diam, long int Nsphere, double Lato,  double **pos, int * NpartInSquare, int ** PartInSquarePos, int **connect){

   int Nsquare =(int)(0.98*Lato/(diam));
   double Lsquare = Lato/Nsquare;

   cerr<<"****************I'm beginning Monte Carlo*****************"<<endl;
   //getchar();
   sgenrand(time(NULL));

   //ofstream distFile;
   //distFile.open("distFile.dat");

   double move[3];

   double gf=(1-PHI_F/2)/((1-PHI_F)*(1-PHI_F)*(1-PHI_F));
   double phi1=(Nsphere/(Lato*Lato*Lato))*(PI/6)*(diam*diam*diam);
   double passo;

   if (phi1<PHI_F) //freezing paking fraction phi_f=0.49 for 3d
        passo= 1+(1-phi1)*(1-phi1)*(1-phi1)/(24*phi1*(1-phi1/2));
   else
        passo= 1+ (PHI_P - phi1)/(24*gf*phi1*(PHI_P-PHI_F));


   cerr<<"Passo = "<<passo-diam<<endl;
   //int step_count=0;
   for (int step=1; step<=NSWEEP; step++){
     //double mean_dist=adaptMC(shell, diam, Nsphere, Lato, pos, NpartInSquare, PartInSquarePos, distFile);

      double f=2*(passo-diam);
      int flag=0;
      //double sumNeigh=0;
      //int countDist=0;

      for ( long int  i=1 ; i <= Nsphere; i++ ){

         long int numSquare=0;

         move[0]=pos[i][0]+ f*(genrand()-0.5)*2;
         move[1]=pos[i][1]+ f*(genrand()-0.5)*2;
         move[2]=pos[i][2]+ f*(genrand()-0.5)*2;

         for(int j=0; j<3; j++){

            if(move[j]>Lato)
               move[j]=move[j]-Lato;
            if(move[j]<0.0)
               move[j]=move[j]+Lato;
         }

         int Nx= (int)(move[0]/Lsquare) +1;
         int Ny= (int)(move[1]/Lsquare) +1;
         int Nz= (int)(move[2]/Lsquare) +1;

         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;

         int count =0;


         for (int n1 =Nx-1; n1<=Nx+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
              for (int n2 =Ny-1; n2<=Ny+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =Nz-1; n3<=Nz+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }
                      count++;
                      //find position in numerized squares from 1 to Nsquare^3
                      numSquare=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);

                  for (long int j=1; j<=NpartInSquare[numSquare-1]; j++){
                     double  delta=0.0;

                     int k=PartInSquarePos[numSquare-1][j];
                     if (k==i){continue;}

                     delta=pow((move[0]-pos[k][0]+Lato*k1),2)+pow((move[1]-pos[k][1]+Lato*k2),2)
                           +pow((move[2]-pos[k][2]+Lato*k3),2);

                     delta=sqrt(delta);

                     if (delta<diam){
                        flag=1;

                        //cerr<<"Overlap..part i = "<<i<<" !! "<<delta<<" < "<< diam<<endl;
                        break;
                     }
                     

                  }/*end for j particles*/

               if (flag==1) break;

               }//end Nz for
            if (flag==1) break;
           }//end Ny for
         if (flag==1) break;
        }//end Nx for
         if (flag==1) {flag=0; continue;};

         
         /*delete old information*/
         int oldSquare = ((int)pos[i][3]+Nsquare*((int)pos[i][4]-1)+Nsquare*Nsquare*((int)pos[i][5]-1));
         int k=0;
         for (long int j=1; j<=NpartInSquare[oldSquare-1]-1; j++){
               k++;
               if (PartInSquarePos[oldSquare-1][j]==i) k++;
               PartInSquarePos[oldSquare-1][j]=PartInSquarePos[oldSquare-1][k];
         }

         PartInSquarePos[oldSquare-1][NpartInSquare[oldSquare-1]]=0;
         NpartInSquare[oldSquare-1]--;

         /*update with new information*/
         pos[i][0]=move[0];
         pos[i][1]=move[1];
         pos[i][2]=move[2];

         pos[i][3]=(double)Nx;
         pos[i][4]=(double)Ny;
         pos[i][5]=(double)Nz;

         int newSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);
         NpartInSquare[newSquare-1]++;

         if (NpartInSquare[newSquare-1]>MAX_SPHERES){
            cerr<<"Too many spheres in square # "<< newSquare<<endl;
            exit(1);
         }

         PartInSquarePos[newSquare-1][NpartInSquare[newSquare-1]]=i;

      } /*end for i-Nsphere*/

   }/*end for sweep*/

   //distFile.close();


}//end function



void MC_impenetrable_sphere(double shell,double diam, long int Nsphere, double Lato,  double **pos, int * NpartInSquare, int ** PartInSquarePos, int **connect, double Tstar, double lambda, int sweep, double &pas){

   int Nsquare =(int)(0.98*Lato/(diam));
   double Lsquare = Lato/Nsquare;
   double cc=0, a, b, c;
   cerr<<"****************I'm beginning Monte Carlo*****************"<<endl;
   //getchar();
   sgenrand(time(NULL));
      

   //ofstream distFile;
   //distFile.open("distFile.dat");

   double move[3];

   double gf=(1-PHI_F/2)/((1-PHI_F)*(1-PHI_F)*(1-PHI_F));
   double phi1=(Nsphere/(Lato*Lato*Lato))*(PI/6)*(diam*diam*diam);
   
   double passo;
   cerr<<"pass "<<pas<<endl; 
  
   if (phi1<PHI_F) //freezing paking fraction phi_f=0.49 for 3d
        passo= 1+(1-phi1)*(1-phi1)*(1-phi1)/(24*phi1*(1-phi1/2));
   else
        passo= 1+ (PHI_P - phi1)/(24*gf*phi1*(PHI_P-PHI_F));

    double  f=(passo-diam)/3;
   //cerr<<"Passo = "<<passo-diam<<endl;
   
   double count_accept=0;
   
   for (int step=1; step<=sweep; step++){
     //double mean_dist=adaptMC(shell, diam, Nsphere, Lato, pos, NpartInSquare, PartInSquarePos, distFile);
  
   
    if (pas==0){
      if(step%IRATIO == 0 && f<Lato/2){
           //cerr<<"STEP = "<<step<< ", IRATIO "<<IRATIO<<endl;
          // ADJUST MAXIMUM DISPLACEMENT 

              double RATIO = (double)count_accept / (Nsphere * IRATIO);
            //cerr<<"count = "<<count_accept <<"/"<< (Nsphere * IRATIO)<<" Ratio = "<<RATIO<<endl; 
            //cerr<<" Old passo = "<<f<<endl;
              if( RATIO > 0.5 )
                   f=f + f*20/100;
              else
                 f=f - f*20/100;
                
             // cerr<<"New passo = "<< f<<endl;
              //getchar();  
              count_accept = 0.0;
              
          }
      
      }
      else {
      f=2*pas;
      
      }
      int flag=0;
      //double sumNeigh=0;
      //int countDist=0;

      for ( long int  i=1 ; i <= Nsphere; i++ ){

         long int numSquare=0;
         
         double sumEnergyNew=0;
         double sumEnergyOld=0;
         double diff;
         
         //cerr<<"-----------------------------Sphere i = "<<i<<endl;
    
         move[0]=pos[i][0]+ f*(genrand()-0.5)*2;
         move[1]=pos[i][1]+ f*(genrand()-0.5)*2;
         move[2]=pos[i][2]+ f*(genrand()-0.5)*2;

         for(int j=0; j<3; j++){

            if(move[j]>Lato)
               move[j]=move[j]-Lato;
            if(move[j]<0.0)
               move[j]=move[j]+Lato;
         }

         int Nx= (int)(move[0]/Lsquare) +1;
         int Ny= (int)(move[1]/Lsquare) +1;
         int Nz= (int)(move[2]/Lsquare) +1;

         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;

         int count =0;


         for (int n1 =Nx-1; n1<=Nx+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
              for (int n2 =Ny-1; n2<=Ny+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =Nz-1; n3<=Nz+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }
                      count++;
                      //find position in numerized squares from 1 to Nsquare^3
                      numSquare=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);

                  for (long int j=1; j<=NpartInSquare[numSquare-1]; j++){
                     double  delta=0.0;
                                         
                     int k=PartInSquarePos[numSquare-1][j];
                     if (k==i){continue;}
                     
                     
                     delta=pow((move[0]-pos[k][0]+Lato*k1),2)+pow((move[1]-pos[k][1]+Lato*k2),2)
                           +pow((move[2]-pos[k][2]+Lato*k3),2);

                     delta=sqrt(delta);

                     if (delta<diam){
                        flag=1;

                        //cerr<<"Overlap..part i = "<<i<<" !! "<<delta<<" < "<< diam<<endl;
                        break;
                     }
                     else if (delta>=diam && delta<= diam + lambda){
                       
                                             
                        sumEnergyNew=sumEnergyNew-1/Tstar;
                     }

                  }/*end for j particles*/
                  if (flag==1) break;
                  
                

               }//end Nz for
            if (flag==1) break;
           }//end Ny for
         if (flag==1) break;
        }//end Nx for
         if (flag==1) {flag=0; continue;};

         //cerr << "-------------- No overlap " <<endl;
         // Calculate the energy relative to the old position
         int Nxold= (int)(pos[i][0]/Lsquare) +1;
         int Nyold= (int)(pos[i][1]/Lsquare) +1;
         int Nzold= (int)(pos[i][2]/Lsquare) +1;
    
        
         if (Nxold==Nsquare+1) Nxold=Nsquare;
         if (Nyold==Nsquare+1) Nyold=Nsquare;
         if (Nzold==Nsquare+1) Nzold=Nsquare;

         numSquare=Nxold+Nsquare*(Nyold-1)+Nsquare*Nsquare*(Nzold-1);
         
         
         for (int n1 =Nxold-1; n1<=Nxold+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =Nyold-1; n2<=Nyold+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                  for (int n3 =Nzold-1; n3<=Nzold+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquare=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)
                      if ( NpartInSquare[numSquare-1]!=0){


                        for (long int j=1; j<=NpartInSquare[numSquare-1]; j++){
                        
                           int k=PartInSquarePos[numSquare-1][j];
                           if (k==i){continue;}
                           double dist0=0.0;
                           //find distance of part i from each particle in the around squares
                           dist0=pow((pos[i][0]-pos[k][0]+Lato*k1),2)+pow((pos[i][1]-pos[k][1]+Lato*k2),2)+
                                pow((pos[i][2]-pos[k][2]+Lato*k3),2);

                           dist0=sqrt(dist0);
                          // cerr<<"---dist0 = " <<dist0<<" diam = "<<diam<<endl;
                           if (dist0<diam){ //if sphere i overlaps sphere k, go out
                              
                              cerr<<"Overlap..!! BIG PROBLEM"<<endl;
                              exit(1);
                             
                              break;
                              }
                              
                           if (dist0>=diam && dist0<= diam + lambda){
                       
                                sumEnergyOld=sumEnergyOld-1/Tstar;
                            }  
                              
                          }
                      }
                    

                   }//end Nz for
               
                }//end Ny for
               
             }//end Nx for
        
         
         
           
         diff=sumEnergyNew-sumEnergyOld;
         //cerr<<"Diff = "<< sumEnergyNew <<" - "<<sumEnergyOld<<" = "<<diff<<endl;
         cc=genrand();
         
         if (diff<=0 || cc < exp(-diff)){
             
           /*  if (diff<=0) cerr <<"Diff < = 0"<<endl;
             else {cerr<<"Prob"<<endl;
             cerr<<"Enew -Eold --> "<<sumEnergyNew<<" - "<<sumEnergyOld<<" = "<<diff<<endl;
             cerr<<cc<<" < "<<exp(-diff)<<endl; }
             
                     
            cerr<<"Accept and update movement of the sphere...i = "<<i<<endl;
            
            getchar();*/
            count_accept=count_accept+1;
            
           // getchar();
            int oldSquare = ((int)pos[i][3]+Nsquare*((int)pos[i][4]-1)+Nsquare*Nsquare*((int)pos[i][5]-1));
            int k=0;
            for (long int j=1; j<=NpartInSquare[oldSquare-1]-1; j++){
                k++;
                if (PartInSquarePos[oldSquare-1][j]==i) k++;
                PartInSquarePos[oldSquare-1][j]=PartInSquarePos[oldSquare-1][k];
            }

            PartInSquarePos[oldSquare-1][NpartInSquare[oldSquare-1]]=0;
            NpartInSquare[oldSquare-1]--;

            /*update with new information*/
            pos[i][0]=move[0];
            pos[i][1]=move[1];
            pos[i][2]=move[2];

            pos[i][3]=(double)Nx;
            pos[i][4]=(double)Ny;
            pos[i][5]=(double)Nz;

            int newSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);
            NpartInSquare[newSquare-1]++;

            if (NpartInSquare[newSquare-1]>MAX_SPHERES){
               cerr<<"Too many spheres in square # "<< newSquare<<endl;
               exit(1);
            }

            PartInSquarePos[newSquare-1][NpartInSquare[newSquare-1]]=i;
        
        
        }//end if sumEnergy  
        //else{ cerr<<"Rejected"<<endl; getchar();} 
        
      } /*end for i-Nsphere*/

   }/*end for sweep*/

   if (pas==0) pas=f;
   //distFile.close();


}//end function


void place_segregation(double shell,double diamCond, double diamIns, long int NsphereCond, long int NsphereIns, double Lato, int ** AdjList, double **posCond, double **posIns, int * NpartInSquareCond, int *NpartInSquareIns, int ** PartInSquarePosCond, int ** PartInSquarePosIns){

   int c;
   int NsquareCond =(int)(0.98 *Lato/(diamCond));
   int NsquareIns=(int)(0.98 *Lato/(diamCond+diamIns));

   //cerr <<"N square Ins = "<<NsquareIns<<endl;
   //cerr << "N square Cond = "<< NsquareCond<<endl;

   double LsquareCond = Lato/NsquareCond;
   double LsquareIns = Lato/NsquareIns;
   sgenrand(time(NULL));
   //cerr <<"L square Ins = "<<LsquareIns<<endl;
   //cerr << "L square Cond = "<< LsquareCond<<endl;
//Firstly Insulating particles are placed, saving position and cube position in PosIns
   int count=0;
   for ( long int  i=1 ; i <= NsphereIns; i++ ){


         posIns[i][0]=(genrand())*Lato;
         posIns[i][1]=(genrand())*Lato;
         posIns[i][2]=(genrand())*Lato;
         /*
         myfile << posIns[i][0] << " ";
         myfile << posIns[i][1] << " ";
         myfile << posIns[i][2] << " ";
         */
         int NxIns= (int)(posIns[i][0]/LsquareIns) +1;
         int NyIns= (int)(posIns[i][1]/LsquareIns) +1;
         int NzIns= (int)(posIns[i][2]/LsquareIns) +1;

         if (NxIns==NsquareIns+1) NxIns=NsquareIns;
         if (NyIns==NsquareIns+1) NyIns=NsquareIns;
         if (NzIns==NsquareIns+1) NzIns=NsquareIns;

         posIns[i][3]=(double)NxIns;
         posIns[i][4]=(double)NyIns;
         posIns[i][5]=(double)NzIns;
         //myfile<<endl;

         long int numSquareIns=NxIns+NsquareIns*(NyIns-1)+NsquareIns*NsquareIns*(NzIns-1);;
         NpartInSquareIns[numSquareIns-1]++;
         PartInSquarePosIns[numSquareIns-1][NpartInSquareIns[numSquareIns-1]]=i;
       count++;
   }
   //myfile.close();
   //cerr <<"count = "<< count;
   //cerr<<endl<<"---Placed "<< NsphereIns<<" Insulating Particles---"<<endl;
   //int a;
   //cin>>a;

   //myfile.open ("segregation3d_Cond.dat");

     //Secondly Conducting particles are added, saving position and cube position in PosCond
     //For each added parlicle check if it intersects the previous placed insulating particles and the other conducting
     //particles
   int flagIns=0;
   int flagCond=0;
   for ( long int  i=1 ; i <= NsphereCond; i++ ){
      try{
         long int numSquareIns=0;
         long int numSquareCond=0;

         //cerr<<endl<<"...Placing Cond particle i = "<<i<<endl;


         posCond[i][0]=(genrand())*Lato;
         posCond[i][1]=(genrand())*Lato;
         posCond[i][2]=(genrand())*Lato;

         int NxCond= (int)(posCond[i][0]/LsquareCond) +1;
         int NyCond= (int)(posCond[i][1]/LsquareCond) +1;
         int NzCond= (int)(posCond[i][2]/LsquareCond) +1;

         int NxIns= (int)(posCond[i][0]/LsquareIns) +1;
         int NyIns= (int)(posCond[i][1]/LsquareIns) +1;
         int NzIns= (int)(posCond[i][2]/LsquareIns) +1;

        numSquareCond=NxCond+NsquareCond*(NyCond-1)+NsquareCond*NsquareCond*(NzCond-1);
        numSquareIns=NxIns+NsquareIns*(NyIns-1)+NsquareIns*NsquareIns*(NzIns-1);

        //cerr<<" I am in square ="<< numSquareIns<< endl;
         //se pos[i][1]==L ==> Nx=Nsquare i.e. choose the last little cube
         if (NxCond==NsquareCond+1) NxCond=NsquareCond;
         if (NyCond==NsquareCond+1) NyCond=NsquareCond;
         if (NzCond==NsquareCond+1) NzCond=NsquareCond;

         if (NxIns==NsquareIns+1) NxIns=NsquareIns;
         if (NyIns==NsquareIns+1) NyIns=NsquareIns;
         if (NzIns==NsquareIns+1) NzIns=NsquareIns;


         for (int n1 =NxIns-1; n1<=NxIns+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=NsquareIns;
                k1=1;
               }
             if (n1>NsquareIns){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =NyIns-1; n2<=NyIns+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=NsquareIns;
                     k2=1;
                   }
                  if (n2>NsquareIns){
                     m2=1;
                     k2=-1;
                   }
                  for (int n3 =NzIns-1; n3<=NzIns+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=NsquareIns;
                        k3=1;
                       }
                      if (n3>NsquareIns){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquareIns=m1+NsquareIns*(m2-1)+NsquareIns*NsquareIns*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)
                      if ( NpartInSquareIns[numSquareIns-1]!=0){

                        //cerr<<"numSquare = "<< numSquare<< endl;

                        //cerr<<"There are "<< NpartInSquareIns[numSquareIns-1] <<" neighbors in square ="<<
                        //numSquareIns<< endl;
                        for (long int j=1; j<=NpartInSquareIns[numSquareIns-1]; j++){

                           //find the index of the neighboring particle
                           int k=PartInSquarePosIns[numSquareIns-1][j];

                           double distIns=0.0;
                           //find distance of part i from each particle in the around squares
                           distIns=pow((posCond[i][0]-posIns[k][0]+Lato*k1),2)+pow((posCond[i][1]-posIns[k][1]+Lato*k2),2)+
                                pow((posCond[i][2]-posIns[k][2]+Lato*k3),2);
                           distIns=sqrt(distIns);

                           if (distIns<diamIns/2+diamCond/2){ //if sphere i overlaps sphere k, go out
                              flagIns=1;
                              //cerr<<"Overlap..!!"<<endl;
                              //getchar();
                              break;
                              }
                          }
                      }
                     if (flagIns==1) break;

                   }//end Nz for
                 if (flagIns==1) break;
                }//end Ny for
               if (flagIns==1) break;
             }//end Nx for
         if (flagIns==1) {i--; flagIns=0; continue;};

         //if (Nz==Nsquare+1) Nz==Nsquare;
        // cin>>a;
         for (int n1 =NxCond-1; n1<=NxCond+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=NsquareCond;
                k1=1;
               }
             if (n1>NsquareCond){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =NyCond-1; n2<=NyCond+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=NsquareCond;
                     k2=1;
                   }
                  if (n2>NsquareCond){
                     m2=1;
                     k2=-1;
                   }
                  for (int n3 =NzCond-1; n3<=NzCond+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=NsquareCond;
                        k3=1;
                       }
                      if (n3>NsquareCond){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquareCond=m1+NsquareCond*(m2-1)+NsquareCond*NsquareCond*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)
                      if ( NpartInSquareCond[numSquareCond-1]!=0){


                        for (long int j=1; j<=NpartInSquareCond[numSquareCond-1]; j++){
                          //Preparare aperitivo
                           //find the index of the neighboring particle
                           int k=PartInSquarePosCond[numSquareCond-1][j];

                           double distCond=0.0;
                           //find distance of part i from each particle in the around squares
                        distCond=pow((posCond[i][0]-posCond[k][0]+Lato*k1),2)+pow((posCond[i][1]-posCond[k][1]+Lato*k2),2)+
                                pow((posCond[i][2]-posCond[k][2]+Lato*k3),2);
                           //inviare il mail!!!!!!!!!!!!!!
                           distCond=sqrt(distCond);
                           //cerr<<"---dist = " <<dist<<" diam = "<<diam<<endl;
                           if (distCond<diamCond){ //if sphere i overlaps sphere k, go out
                              flagCond=1;
                              //cerr<<"Overlap..!!"<<endl;
                              //getchar();
                              break;
                              }
                          }
                      }
                     if (flagCond==1) break;

                   }//end Nz for
                 if (flagCond==1) break;
                }//end Ny for
               if (flagCond==1) break;
             }//end Nx for
         if (flagCond==1) {i--; flagCond=0; continue;};

         //cerr << "Count = "<< count;
         /*
         cerr<<" Placing sphere...i = "<<i<<endl;
         myfile << posCond[i][0] << " ";
         myfile << posCond[i][1] << " ";
         myfile << posCond[i][2] << " ";
         //myfile << pos[i][2] << " ";
         */
         posCond[i][3]=(double)NxCond;
         posCond[i][4]=(double)NyCond;
         posCond[i][5]=(double)NzCond;
         //pos[i][5]=(double)Nz;

         //numSquareIns=NxIns+NsquareIns*(NyIns-1)+NsquareIns*NsquareIns*(NzIns-1);
         numSquareCond=NxCond+NsquareCond*(NyCond-1)+NsquareCond*NsquareCond*(NzCond-1);

         NpartInSquareCond[numSquareCond-1]++;
         //cerr<<"Now there are N part = "<< NpartInSquare[numSquare-1] << " in square = "<< numSquare <<endl;

         if (NpartInSquareCond[numSquareCond-1]>MAX_SPHERES){
            cerr<<"Too many spheres in square # "<< numSquareCond<<endl;
            exit(1);
            }
        // cerr<<"---Placed "<<i<<endl;
         PartInSquarePosCond[numSquareCond-1][NpartInSquareCond[numSquareCond-1]]=i;
         //cerr<<"Placed "<<i<<endl;

         //myfile<< endl;
   }
     catch (exception& e)
  {
    cerr << "exception caught: " << e.what() << endl;
  }



   } //end for i-Nsphere
   //cerr<<endl<<"---Placed "<<NsphereCond<<" Conducting Particles---"<<endl;

}

int place_segregation_lattice(double shell,double diamCond, double diamIns, long int NsphereCond, long int NsphereIns, double Lato, int ** AdjList, double **posCond, double **posIns, int * NpartInSquareCond, int *NpartInSquareIns, int ** PartInSquarePosCond, int ** PartInSquarePosIns){

   int c;
   int NsquareCond =(int)(0.98 *Lato/(diamCond));
   int NsquareIns=(int)(0.98 *Lato/(diamCond+diamIns));

   double volIns=PI/6*(diamIns)*(diamIns)*(diamIns);
   double phi2=1-exp(-volIns*NsphereIns/(Lato*Lato*Lato));
   double espon=(1+diamCond/diamIns)*(1+diamCond/diamIns)*(1+diamCond/diamIns);
   double phi_avail=pow(1-phi2,espon);

   double phi1=NsphereCond*(PI/6)*(diamCond)*(diamCond)*(diamCond)/(Lato*Lato*Lato);
   double conc_eff=phi1/phi_avail;

   double num=1.0/3.0;

   double Lgrid=pow(PI/6/conc_eff,num);

   int NsquareGrid =(int)(Lato/(Lgrid));
   cerr<< "NsquareGrid = "<<NsquareGrid<<endl;
   cerr<< "Lgrid = "<<Lgrid<<endl;
   cerr<<"conc_eff = "<<conc_eff<<" conc iniz  ="<<phi1<<endl;
   cerr << "NsquareGrid*Lgrid = " <<  Lgrid*NsquareGrid << endl;
   //cerr<<Lgrid*NsquareGrid<<endl;
   //getchar();
   cerr<<(Lato/NsquareGrid)<<endl;
   Lgrid=(Lato/NsquareGrid);
   cerr << "NsquareGrid*Lgrid = " <<  Lgrid*NsquareGrid << endl;



   double LsquareCond = Lato/NsquareCond;
   double LsquareIns = Lato/NsquareIns;
   sgenrand(time(NULL));
   //cerr <<"L square Ins = "<<LsquareIns<<endl;
   //cerr << "L square Cond = "<< LsquareCond<<endl;

   if (conc_eff>=PI/6 || Lgrid <= 1){
      cerr<<"Effective concentration is too high for a square Lattice"<<endl;
      return 0;
   }

   getchar();

   //Firstly Insulating particles are placed, saving position and cube position in PosIns
   int count=0;
   ofstream pippo1;
   //pippo1.open("segregation3d_Ins.dat");
   for ( long int  i=1 ; i <= NsphereIns; i++ ){


         posIns[i][0]=(genrand())*Lato;
         posIns[i][1]=(genrand())*Lato;
         posIns[i][2]=(genrand())*Lato;
         /*
         myfile << posIns[i][0] << " ";
         myfile << posIns[i][1] << " ";
         myfile << posIns[i][2] << " ";
         */
         int NxIns= (int)(posIns[i][0]/LsquareIns) +1;
         int NyIns= (int)(posIns[i][1]/LsquareIns) +1;
         int NzIns= (int)(posIns[i][2]/LsquareIns) +1;

         if (NxIns==NsquareIns+1) NxIns=NsquareIns;
         if (NyIns==NsquareIns+1) NyIns=NsquareIns;
         if (NzIns==NsquareIns+1) NzIns=NsquareIns;

         posIns[i][3]=(double)NxIns;
         posIns[i][4]=(double)NyIns;
         posIns[i][5]=(double)NzIns;
       //  pippo1<<posIns[i][0]<<" ";
        // pippo1<<posIns[i][1]<<" ";
         //pippo1<<posIns[i][2]<<" ";
         //pippo1<<endl;

         //myfile<<endl;

         long int numSquareIns=NxIns+NsquareIns*(NyIns-1)+NsquareIns*NsquareIns*(NzIns-1);;
         NpartInSquareIns[numSquareIns-1]++;
         PartInSquarePosIns[numSquareIns-1][NpartInSquareIns[numSquareIns-1]]=i;
       count++;
   }
   //myfile.close();
   //cerr <<"count = "<< count;
   cerr<<endl<<"---Placed "<< NsphereIns<<" Insulating Particles---"<<endl;


     //Secondly Conducting particles are added, saving position and cube position in PosCond
     //For each added parlicle check if it intersects the previous placed insulating particles and the other conducting
     //particles
   int flagIns=0;
   int flag_pos=0;
   int count_pos=1;
   for ( int ii=0 ; ii <= NsquareGrid-1; ii++ ){
      for ( int jj=0 ; jj <= NsquareGrid-1; jj++ ){
         for ( int kk=0 ; kk <= NsquareGrid-1; kk++ ){

              int flagCond=0;
              double x = ii*Lgrid;
              double y = jj*Lgrid;
              double z = kk*Lgrid;
              //cerr<<"(x,,y,z) = "<<x<<", "<<y<<", "<<z<<endl;
              //getchar();
              int NxIns= (int)(x/LsquareIns) +1;
              int NyIns= (int)(y/LsquareIns) +1;
              int NzIns= (int)(z/LsquareIns) +1;

              long int numSquareIns=0;

              for (int n1 =NxIns-1; n1<=NxIns+1; n1++){
               int m1=n1;
               int k1=0;
               if (n1<1){ //i.e. Nx==1
                  m1=NsquareIns;
                  k1=1;
                 }
               if (n1>NsquareIns){ //i.e. Nx==Nsquare
                 m1=1;
                 k1=-1;
               }
               for (int n2 =NyIns-1; n2<=NyIns+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=NsquareIns;
                     k2=1;
                   }
                  if (n2>NsquareIns){
                     m2=1;
                     k2=-1;
                   }
                  for (int n3 =NzIns-1; n3<=NzIns+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=NsquareIns;
                        k3=1;
                       }
                      if (n3>NsquareIns){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquareIns=m1+NsquareIns*(m2-1)+NsquareIns*NsquareIns*(m3-1);
                      if ( NpartInSquareIns[numSquareIns-1]==0){
                               posCond[count_pos][0]=x;
                               posCond[count_pos][1]=y;
                               posCond[count_pos][2]=z;
                               //cerr<<"posCond["<<count_pos<<"][0]= "<<x<<endl;
                               flagCond=1;
                       }


                      else if ( NpartInSquareIns[numSquareIns-1]!=0){


                        for (long int j=1; j<=NpartInSquareIns[numSquareIns-1]; j++){

                           //find the index of the neighboring particle
                           int k=PartInSquarePosIns[numSquareIns-1][j];

                           double distIns=0.0;
                           //find distance of part i from each particle in the around squares
                           distIns=pow((x-posIns[k][0]+Lato*k1),2)+pow((y-posIns[k][1]+Lato*k2),2)+
                                pow((z-posIns[k][2]+Lato*k3),2);
                           distIns=sqrt(distIns);

                           if (distIns>diamIns/2+diamCond/2){ //if sphere i overlaps sphere k, go out
                               posCond[count_pos][0]=x;
                               posCond[count_pos][1]=y;
                               posCond[count_pos][2]=z;
                               //cerr<<"posCond["<<count_pos<<"][0]= "<<x<<endl;
                               flagCond=1;

                            }
                           else {
                                //cerr<<"Overlap..."<<endl;
                                 flagIns=1;

                                 break;
                                }
                          }//for NpartInSquareIns
                      }//if NpartInSquareIns
                    if (flagIns==1) break;

                   }//end Nz for
                 if (flagIns==1) break;
                }//end Ny for
               if (flagIns==1) break;
             }//end Nx for
            if (flagIns==1) { flagIns=0; continue;};


             if (count_pos==NsphereCond) {
                   cerr<<"-----"<<NsphereCond;
                   flag_pos=1;
                   count_pos++;
                   break;
             }
             if (flagCond==1) count_pos++;


         }//end ii
         if (flag_pos==1) {break;}
      }//end jj
      if (flag_pos==1) {break;}
   }//end kk

  cerr<<"pos_cond "<< count_pos << endl;
   for ( long int  i=1 ; i <= count_pos-1; i++ ){
         long int numSquare=0;

         int Nx= (int)(posCond[i][0]/LsquareCond) +1;
         int Ny= (int)(posCond[i][1]/LsquareCond) +1;
         int Nz= (int)(posCond[i][2]/LsquareCond) +1;

         if (Nx==NsquareCond+1) Nx=NsquareCond;
         if (Ny==NsquareCond+1) Ny=NsquareCond;
         if (Nz==NsquareCond+1) Nz=NsquareCond;


         posCond[i][3]=(double)Nx;
         posCond[i][4]=(double)Ny;
         posCond[i][5]=(double)Nz;


         numSquare=Nx+NsquareCond*(Ny-1)+NsquareCond*NsquareCond*(Nz-1);

         NpartInSquareCond[numSquare-1]++;

         if (NpartInSquareCond[numSquare-1]>MAX_SPHERES){
            cerr<<"Too many spheres in square # "<< numSquare<<endl;
            exit(1);
            }
        // cerr<<"---Placed "<<i<<endl;
         PartInSquarePosCond[numSquare-1][NpartInSquareCond[numSquare-1]]=i;
         //Update connect


   } //end for i-Nsphere



   if (count_pos-1<NsphereCond){

   int flagIns=0;
   int flagCond=0;
   for ( long int  i=count_pos ; i <= NsphereCond; i++ ){
      try{
         long int numSquareIns=0;
         long int numSquareCond=0;


            posCond[i][0]=(genrand())*Lato;
            posCond[i][1]=(genrand())*Lato;
            posCond[i][2]=(genrand())*Lato;

         int NxCond= (int)(posCond[i][0]/LsquareCond) +1;
         int NyCond= (int)(posCond[i][1]/LsquareCond) +1;
         int NzCond= (int)(posCond[i][2]/LsquareCond) +1;

         int NxIns= (int)(posCond[i][0]/LsquareIns) +1;
         int NyIns= (int)(posCond[i][1]/LsquareIns) +1;
         int NzIns= (int)(posCond[i][2]/LsquareIns) +1;

        numSquareCond=NxCond+NsquareCond*(NyCond-1)+NsquareCond*NsquareCond*(NzCond-1);
        numSquareIns=NxIns+NsquareIns*(NyIns-1)+NsquareIns*NsquareIns*(NzIns-1);

         //se pos[i][1]==L ==> Nx=Nsquare i.e. choose the last little cube
         if (NxCond==NsquareCond+1) NxCond=NsquareCond;
         if (NyCond==NsquareCond+1) NyCond=NsquareCond;
         if (NzCond==NsquareCond+1) NzCond=NsquareCond;

         if (NxIns==NsquareIns+1) NxIns=NsquareIns;
         if (NyIns==NsquareIns+1) NyIns=NsquareIns;
         if (NzIns==NsquareIns+1) NzIns=NsquareIns;


         for (int n1 =NxIns-1; n1<=NxIns+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=NsquareIns;
                k1=1;
               }
             if (n1>NsquareIns){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =NyIns-1; n2<=NyIns+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=NsquareIns;
                     k2=1;
                   }
                  if (n2>NsquareIns){
                     m2=1;
                     k2=-1;
                   }
                  for (int n3 =NzIns-1; n3<=NzIns+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=NsquareIns;
                        k3=1;
                       }
                      if (n3>NsquareIns){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquareIns=m1+NsquareIns*(m2-1)+NsquareIns*NsquareIns*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)
                      if ( NpartInSquareIns[numSquareIns-1]!=0){


                        //numSquareIns<< endl;
                        for (long int j=1; j<=NpartInSquareIns[numSquareIns-1]; j++){

                           //find the index of the neighboring particle
                           int k=PartInSquarePosIns[numSquareIns-1][j];

                           double distIns=0.0;
                           //find distance of part i from each particle in the around squares
                           distIns=pow((posCond[i][0]-posIns[k][0]+Lato*k1),2)+pow((posCond[i][1]-posIns[k][1]+Lato*k2),2)+
                                pow((posCond[i][2]-posIns[k][2]+Lato*k3),2);
                           distIns=sqrt(distIns);

                           if (distIns<diamIns/2+diamCond/2){ //if sphere i overlaps sphere k, go out
                              flagIns=1;
                              //cerr<<"Overlap..!!"<<endl;
                              //getchar();
                              break;
                              }
                          }
                      }
                     if (flagIns==1) break;

                   }//end Nz for
                 if (flagIns==1) break;
                }//end Ny for
               if (flagIns==1) break;
             }//end Nx for
         if (flagIns==1) {i--; flagIns=0; continue;};

         //if (Nz==Nsquare+1) Nz==Nsquare;
        // cin>>a;
         for (int n1 =NxCond-1; n1<=NxCond+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=NsquareCond;
                k1=1;
               }
             if (n1>NsquareCond){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =NyCond-1; n2<=NyCond+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=NsquareCond;
                     k2=1;
                   }
                  if (n2>NsquareCond){
                     m2=1;
                     k2=-1;
                   }
                  for (int n3 =NzCond-1; n3<=NzCond+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=NsquareCond;
                        k3=1;
                       }
                      if (n3>NsquareCond){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquareCond=m1+NsquareCond*(m2-1)+NsquareCond*NsquareCond*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)
                      if ( NpartInSquareCond[numSquareCond-1]!=0){


                        for (long int j=1; j<=NpartInSquareCond[numSquareCond-1]; j++){

                           //find the index of the neighboring particle
                           int k=PartInSquarePosCond[numSquareCond-1][j];

                           double distCond=0.0;
                           //find distance of part i from each particle in the around squares
                        distCond=pow((posCond[i][0]-posCond[k][0]+Lato*k1),2)+pow((posCond[i][1]-posCond[k][1]+Lato*k2),2)+
                                pow((posCond[i][2]-posCond[k][2]+Lato*k3),2);

                           distCond=sqrt(distCond);
                           //cerr<<"---dist = " <<dist<<" diam = "<<diam<<endl;
                           if (distCond<diamCond){ //if sphere i overlaps sphere k, go out
                              flagCond=1;
                              //cerr<<"Overlap..!!"<<endl;
                              //getchar();
                              break;
                              }
                          }
                      }
                     if (flagCond==1) break;

                   }//end Nz for
                 if (flagCond==1) break;
                }//end Ny for
               if (flagCond==1) break;
             }//end Nx for
         if (flagCond==1) {i--; flagCond=0; continue;};
         //cerr<<"---Placed "<<i<<endl;
         //cerr << "Count = "<< count;

         posCond[i][3]=(double)NxCond;
         posCond[i][4]=(double)NyCond;
         posCond[i][5]=(double)NzCond;


         //numSquareIns=NxIns+NsquareIns*(NyIns-1)+NsquareIns*NsquareIns*(NzIns-1);
         numSquareCond=NxCond+NsquareCond*(NyCond-1)+NsquareCond*NsquareCond*(NzCond-1);

         NpartInSquareCond[numSquareCond-1]++;
         //cerr<<"Now there are N part = "<< NpartInSquare[numSquare-1] << " in square = "<< numSquare <<endl;

         if (NpartInSquareCond[numSquareCond-1]>MAX_SPHERES){
            cerr<<"Too many spheres in square # "<< numSquareCond<<endl;
            exit(1);
            }
         //cerr<<"---Placed "<<i<<endl;
         PartInSquarePosCond[numSquareCond-1][NpartInSquareCond[numSquareCond-1]]=i;


   }
     catch (exception& e)
  {
    cerr << "exception caught: " << e.what() << endl;
  }



   } //end for i-Nsphere


   }



      if (posCond[NsphereCond][0]==0 && posCond[NsphereCond][1]==0 && posCond[NsphereCond][2]==0)
        {
            cerr<<"There are not Nspherecond Spheres #"<<endl;
            exit(1);
            }


     // pippo.close();

      return count_pos-1;


}


void MC_segregation(double shell,double diamCond, double diamIns, long int NsphereCond, long int NsphereIns, double Lato, int ** AdjList, double **posCond, double **posIns, int * NpartInSquareCond, int *NpartInSquareIns, int ** PartInSquarePosCond, int ** PartInSquarePosIns, int **connect){

   int NsquareCond =(int)(0.98 *Lato/(diamCond));
   int NsquareIns=(int)(0.98 *Lato/(diamCond+diamIns));

   double LsquareCond = Lato/NsquareCond;
   double LsquareIns = Lato/NsquareIns;
   ofstream percFile;
   //cerr<<"****************I'm beginning Monte Carlo*****************"<<endl;
   //getchar();
   sgenrand(time(NULL));

   double move[3];
   double volIns=PI/6*(diamIns)*(diamIns)*(diamIns);
   double phi2=1-exp(-volIns*NsphereIns/(Lato*Lato*Lato));
   double espon=(1+diamCond/diamIns)*(1+diamCond/diamIns)*(1+diamCond/diamIns);

   double phi_avail=pow(1-phi2,espon);
   double V_avail=phi_avail*(Lato*Lato*Lato);
   double phi_eff=(NsphereCond/(V_avail))*(PI/6)*(diamCond*diamCond*diamCond);

   double gf=(1-PHI_F/2)/((1-PHI_F)*(1-PHI_F)*(1-PHI_F));
   double passo;
   if (phi_eff<PHI_F)
        passo= 1+(1-phi_eff)*(1-phi_eff)*(1-phi_eff)/(24*phi_eff*(1-phi_eff/2));
   else
        passo= 1+ (PHI_P - phi_eff)/(24*gf*phi_eff*(PHI_P-PHI_F));


   //double mean_dist=adaptMC(shell, diamCond, NsphereCond, Lato, posCond, NpartInSquareCond, PartInSquarePosCond,  //percFile);
   int step_count=0;
   //cerr<<"MeanDist = "<<mean_dist<<endl;
   for (int step=1; step<=NSWEEP; step++){
      //cerr<<step<<"xxxxxxxxxxxxxx";
      double f=(passo-diamCond);
      int flag=0;
      //double sumNeigh=0;
      //int countDist=0;
      //cerr <<step<<endl;
      int flagIns=0;
      int flagCond=0;
      for ( long int  i=1 ; i <= NsphereCond; i++ ){
         try{
            long int numSquareIns=0;
            long int numSquareCond=0;

            move[0]=posCond[i][0]+ f*(genrand()-0.5)*2/2;
            move[1]=posCond[i][1]+ f*(genrand()-0.5)*2/2;
            move[2]=posCond[i][2]+ f*(genrand()-0.5)*2/2;

            for(int j=0; j<3; j++){

               if(move[j]>Lato)
                  move[j]=move[j]-Lato;
               if(move[j]<0.0)
                  move[j]=move[j]+Lato;
            }

            int NxCond= (int)(posCond[i][0]/LsquareCond) +1;
            int NyCond= (int)(posCond[i][1]/LsquareCond) +1;
            int NzCond= (int)(posCond[i][2]/LsquareCond) +1;

            int NxIns= (int)(posCond[i][0]/LsquareIns) +1;
            int NyIns= (int)(posCond[i][1]/LsquareIns) +1;
            int NzIns= (int)(posCond[i][2]/LsquareIns) +1;

            if (NxCond==NsquareCond+1) NxCond=NsquareCond;
            if (NyCond==NsquareCond+1) NyCond=NsquareCond;
            if (NzCond==NsquareCond+1) NzCond=NsquareCond;

            if (NxIns==NsquareIns+1) NxIns=NsquareIns;
            if (NyIns==NsquareIns+1) NyIns=NsquareIns;
            if (NzIns==NsquareIns+1) NzIns=NsquareIns;


            int count =0;


            for (int n1 =NxIns-1; n1<=NxIns+1; n1++){
                int m1=n1;
                int k1=0;
                if (n1<1){ //i.e. Nx==1
                   m1=NsquareIns;
                   k1=1;
                  }
                if (n1>NsquareIns){ //i.e. Nx==Nsquare
                   m1=1;
                   k1=-1;
                  }
                for (int n2 =NyIns-1; n2<=NyIns+1; n2++){
                     int m2=n2;
                     int k2=0;
                     if (n2<1){
                        m2=NsquareIns;
                        k2=1;
                      }
                     if (n2>NsquareIns){
                        m2=1;
                        k2=-1;
                      }
                     for (int n3 =NzIns-1; n3<=NzIns+1; n3++){
                       int m3=n3;
                       int k3=0;
                       if (n3<1){
                          m3=NsquareIns;
                          k3=1;
                       }
                       if (n3>NsquareIns){
                          m3=1;
                          k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                       numSquareIns=m1+NsquareIns*(m2-1)+NsquareIns*NsquareIns*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)

                          for (long int j=1; j<=NpartInSquareIns[numSquareIns-1]; j++){

                           //find the index of the neighboring particle
                             int k=PartInSquarePosIns[numSquareIns-1][j];
                             if (k==i){continue;}

                             double distIns=0.0;
                           //find distance of part i from each particle in the around squares
                             distIns=pow((move[0]-posIns[k][0]+Lato*k1),2)+pow((move[1]-posIns[k][1]+Lato*k2),2)+
                                pow((move[2]-posIns[k][2]+Lato*k3),2);
                             distIns=sqrt(distIns);

                             if (distIns<diamIns/2+diamCond/2){ //if sphere i overlaps sphere k, go out
                                flagIns=1;
                                //cerr<<"Overlap..Ins!!"<<endl;
                              //getchar();
                                break;
                              }
                          }

                      if (flagIns==1) break;

                   }//end Nz for
                   if (flagIns==1) break;
                 }//end Ny for
                 if (flagIns==1) break;
               }//end Nx for
               if (flagIns==1) { flagIns=0; continue;};

         //if (Nz==Nsquare+1) Nz==Nsquare;
        // cin>>a;
            for (int n1 =NxCond-1; n1<=NxCond+1; n1++){
               int m1=n1;
               int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=NsquareCond;
                k1=1;
               }
             if (n1>NsquareCond){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =NyCond-1; n2<=NyCond+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=NsquareCond;
                     k2=1;
                   }
                  if (n2>NsquareCond){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =NzCond-1; n3<=NzCond+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=NsquareCond;
                        k3=1;
                       }
                      if (n3>NsquareCond){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquareCond=m1+NsquareCond*(m2-1)+NsquareCond*NsquareCond*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)

                        for (long int j=1; j<=NpartInSquareCond[numSquareCond-1]; j++){

                           //find the index of the neighboring particle
                           int k=PartInSquarePosCond[numSquareCond-1][j];
                           if (k==i){continue;}
                           double distCond=0.0;
                           //find distance of part i from each particle in the around squares
                          distCond=pow((move[0]-posCond[k][0]+Lato*k1),2)+pow((move[1]-posCond[k][1]+Lato*k2),2)+
                                pow((move[2]-posCond[k][2]+Lato*k3),2);

                           distCond=sqrt(distCond);

                           //cerr<<"---dist = " <<dist<<" diam = "<<diam<<endl;
                           if (distCond<diamCond){ //if sphere i overlaps sphere k, go out
                              flagCond=1;
                              //cerr<<"Overlap..!!"<<endl;
                              //getchar();
                              break;
                              }
                            }

                        if (flagCond==1) break;

                      }//end Nz for
                     if (flagCond==1) break;
                   }//end Ny for
                   if (flagCond==1) break;
                }//end Nx for
                if (flagCond==1) { flagCond=0; continue;};

         /*delete old information*/
         int oldSquareCond = ((int)posCond[i][3]+NsquareCond*((int)posCond[i][4]-1)+NsquareCond*NsquareCond*((int)posCond[i][5]-1));
         int k=0;
         for (long int j=1; j<=NpartInSquareCond[oldSquareCond-1]-1; j++){
               k++;
               if (PartInSquarePosCond[oldSquareCond-1][j]==i) k++;
               PartInSquarePosCond[oldSquareCond-1][j]=PartInSquarePosCond[oldSquareCond-1][k];
         }

         PartInSquarePosCond[oldSquareCond-1][NpartInSquareCond[oldSquareCond-1]]=0;
         NpartInSquareCond[oldSquareCond-1]--;

         /*update with new information*/

            posCond[i][0]=move[0];
            posCond[i][1]=move[1];
            posCond[i][2]=move[2];

            posCond[i][3]=(double)NxCond;
            posCond[i][4]=(double)NyCond;
            posCond[i][5]=(double)NzCond;
         //pos[i][5]=(double)Nz;
            int newSquareCond=NxCond+NsquareCond*(NyCond-1)+NsquareCond*NsquareCond*(NzCond-1);
         NpartInSquareCond[newSquareCond-1]++;

         if (NpartInSquareCond[newSquareCond-1]>MAX_SPHERES){
            cerr<<"Too many spheres in square # "<< newSquareCond<<endl;
            exit(1);
         }

         PartInSquarePosCond[newSquareCond-1][NpartInSquareCond[newSquareCond-1]]=i;

      }//end try
      catch (exception& e)
      {
        cerr << "exception caught: " << e.what() << endl;
      }

     } //end for i-Nsphere

   }//end step


}


void place_impenetrable_spherocylinder(double shell,double diam, double length, long int Nsphere, double Lato, double **pos, double **direct, int * NpartInSquare, int ** PartInSquarePos){


   int Nsquare =(int)(0.99*Lato/(length+diam));
   double Lsquare = Lato/Nsquare;
   long double randc=0.0;
   cerr << "Nsquare = "<< Nsquare << " Lsquare = " << Lsquare <<endl;
   cerr << "Nsquare*Lsquare = " <<  Nsquare*Lsquare << endl;
   cerr<<"Length "<<length<<endl;
   sgenrand(time(NULL));

   ofstream myfile;//AdjFile;
   //myfile.open ("impenetrable3d.dat");


   int flag=0;
   for ( long int  i=1 ; i <= Nsphere; i++ ){
         long int numSquare=0;

         //getchar();
         //randc = (genrand())*Lato;
         //cerr<<endl<<"I'm here...i = "<<i<<endl;
         pos[i][0]=(genrand())*Lato;
         pos[i][1]=(genrand())*Lato;
         pos[i][2]=(genrand())*Lato;
         
         direct[i][0]=(genrand())*2-1;
         direct[i][1]=(genrand())*2-1;
         direct[i][2]=(genrand())*2-1;
         
         double dir[]={direct[i][0], direct[i][1], direct[i][2]};
         direct[i][0]= direct[i][0]/VecLength(dir);    
         direct[i][1]=direct[i][1]/VecLength(dir); 
         direct[i][2]=direct[i][2]/VecLength(dir); 
         
         
         int Nx= (int)(pos[i][0]/Lsquare) +1;
         int Ny= (int)(pos[i][1]/Lsquare) +1;
         int Nz= (int)(pos[i][2]/Lsquare) +1;

         //debug

        numSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);
        //int numSquareTemp=numSquare;

        //cerr<< "--- New position = ---"<< numSquare<< endl;
        //cerr<< "(x,y,z) = "<< "("<< pos[i][0]<<", "<<pos[i][1]<<", "<<pos[i][2]<<")";
        
        
         //se pos[i][1]==L ==> Nx=Nsquare i.e. choose the last little cube
         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;

         int count =0;
         for (int n1 =Nx-1; n1<=Nx+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
              for (int n2 =Ny-1; n2<=Ny+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =Nz-1; n3<=Nz+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }
                      count++;
                      //find position in numerized squares from 1 to Nsquare^3
                      numSquare=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)
                      if ( NpartInSquare[numSquare-1]!=0){

                        //cerr<<"numSquare = "<< numSquare<< endl;

                        //cerr<<"There are "<< NpartInSquare[numSquare-1] <<" neighbors in square = "<< numSquare <<endl;
                        for (long int j=1; j<=NpartInSquare[numSquare-1]; j++){

                           //find the index of the neighboring particle
                           int k=PartInSquarePos[numSquare-1][j];
                           //cerr << "Pos vicino k = "<< k<< endl;

                           double centredist=0.0;
                           //find distance of part i from each particle in the around squares
                           double rx=(pos[i][0]-pos[k][0]+Lato*k1);
                           double ry=(pos[i][1]-pos[k][1]+Lato*k2);
                           double rz=(pos[i][2]-pos[k][2]+Lato*k3);
                           
                           centredist=rx*rx+ry*ry+rz*rz;
                           if (centredist<=(length+diam)*(length+diam)){
                           
                              //double vr[]={fabs(rx), fabs(ry), fabs(rz)};
                              double vr[]={rx, ry, rz};
                              double w1[]={direct[i][0], direct[i][1], direct[i][2]};
                              double w2[]={direct[k][0], direct[k][1], direct[k][2]};
                           
                      
                              double dist = dist2_rods(centredist, vr, w1, w2, length/2, length/2);
                           
                           //dist=sqrt(dist);
                           //cerr<<"---dist = " <<dist<<" diam = "<<diam<<endl;
                           
                           
                              if (dist<diam*diam){ //if sphere i overlaps sphere k, go out
                                 flag=1;
                                 //cerr<<"Overlap..!! dist=  " << dist<<endl;
                                 //getchar();
                                 break;
                                 }
                            }//end if centredist
                         }//end forNpartInSquare
                      }
                     if (flag==1) break;

                   }//end Nz for
                 if (flag==1) break;
                }//end Ny for
               if (flag==1) break;
             }//end Nx for
         if (flag==1) {i--; flag=0; continue;};

         //cerr << "Count = "<< count;

         //cerr<<" Placing cylinder...i = "<<i<<endl;
         //myfile << pos[i][0] << " ";
         //myfile << pos[i][1] << " ";
         //myfile << pos[i][2] << " ";

         pos[i][3]=(double)Nx;
         pos[i][4]=(double)Ny;
         pos[i][5]=(double)Nz;
      
         numSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);
         //cerr<<"Before there were N part = "<< NpartInSquare[numSquare-1] << " in square = "<< numSquare <<endl;
         NpartInSquare[numSquare-1]++;
         //cerr<<"Now there are N part = "<< NpartInSquare[numSquare-1] << " in square = "<< numSquare <<endl;
         if (NpartInSquare[numSquare-1]>MAX_SPHERES){
            cerr<<"Too many spheres in square # "<< numSquare<<endl;
            exit(1);
            }

         PartInSquarePos[numSquare-1][NpartInSquare[numSquare-1]]=i;
         //getchar();
         //cerr<< "PartPos["<<numSquare-1<<"]["<< NpartInSquare[numSquare-1] <<"] = "<<
             //           PartInSquarePos[numSquare-1][NpartInSquare[numSquare-1]]<<endl;

        // myfile<< endl;


   } //end for i-Nsphere

   //myfile.close();
}


void MC_spherocyl_pure(double diam, long int Nsphere, double Lato, double length, double **pos,  double **direct, int * NpartInSquare, int ** PartInSquarePos, int Nsweep_cyl, double &pas, double &ds){

   int Nsquare =(int)(0.99*Lato/(length+diam));
   double Lsquare = Lato/Nsquare;

   cerr<<"****************I'm beginning Monte Carlo*****************"<<endl;
   //getchar();
   sgenrand(time(NULL));

   //ofstream distFile;
   //distFile.open("distFile.dat");
   int flag=0;
   double move[3];
   double wold[3];
   double passo=0.2*diam;
   double dotstep=0.02;
   double count_accept=0;
   
   cerr<<"Passo = "<<pas<<" Dotstep = "<<ds<<endl;
   //int step_count=0;
   
   /*double *timestep=allocate_double_vector(Nsweep_cyl);
   for (int l=0; l<Nsweep_cyl; l++)
            timestep[l]=0; 
   */
   for (int step=1; step<=Nsweep_cyl; step++){
    
      if (pas==0){
      
       //timestep[step-1]=check_order(direct, Nsphere, length);
       if(step%IRATIO == 0 && passo<Lato/2){
          
             double Ratio = (double)count_accept / (Nsphere * IRATIO);
            
              if( Ratio > 0.3 ){
                   passo=passo + passo*10/100;
                   dotstep=dotstep+dotstep*10/100;
                   }       
              else{
                 passo=passo - passo*10/100;
                 dotstep=dotstep-dotstep*10/100;    
                }
              
              count_accept = 0.0;
              cerr<<"Ratio "<<Ratio<<endl;    
        }
      }
      else {
      dotstep=ds;
      passo=pas;
      
      }
      
      
      for ( long int  i=1 ; i <= Nsphere; i++ ){

         long int numSquare=0;

         move[0]=pos[i][0]+ passo*(genrand()-0.5)*2;
         move[1]=pos[i][1]+ passo*(genrand()-0.5)*2;
         move[2]=pos[i][2]+ passo*(genrand()-0.5)*2;
         
         wold[0]=direct[i][0]; 
         wold[1]=direct[i][1]; 
         wold[2]=direct[i][2];
         
         flip(wold, dotstep);
         
         for(int j=0; j<3; j++){

            if(move[j]>Lato)
               move[j]=move[j]-Lato;
            if(move[j]<0.0)
               move[j]=move[j]+Lato;
         }

         int Nx= (int)(move[0]/Lsquare) +1;
         int Ny= (int)(move[1]/Lsquare) +1;
         int Nz= (int)(move[2]/Lsquare) +1;

         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;

         int count =0;


         for (int n1 =Nx-1; n1<=Nx+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
              for (int n2 =Ny-1; n2<=Ny+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =Nz-1; n3<=Nz+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }
                      count++;
                      //find position in numerized squares from 1 to Nsquare^3
                      numSquare=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);
                      //cerr<<"There are "<< NpartInSquare[numSquare-1] <<" neighbors in square = "<< numSquare <<endl;
                  for (long int j=1; j<=NpartInSquare[numSquare-1]; j++){
                     double  delta=0.0;

                     int k=PartInSquarePos[numSquare-1][j];
                     if (k==i){continue;}
                     double centredist=0.0;
                           //find distance of part i from each particle in the around squares
                     double rx=(move[0]-pos[k][0]+Lato*k1);
                     double ry=(move[1]-pos[k][1]+Lato*k2);
                     double rz=(move[2]-pos[k][2]+Lato*k3);
                     
                      centredist=rx*rx+ry*ry+rz*rz;
                      if (centredist<=(length+diam)*(length+diam)){
                     
                          //double vr[]={fabs(rx), fabs(ry), fabs(rz)};
                          double vr[]={rx, ry, rz};
                          double w1[]={wold[0], wold[1], wold[2]};
                          double w2[]={direct[k][0], direct[k][1], direct[k][2]};
                      
                          double dist = dist2_rods(centredist, vr, w1, w2, length/2, length/2);
                           
                          
                          if (dist<diam*diam){
                             flag=1;

                             //cerr<<"Overlap..part i = "<<i<<" !! "<<delta<<" < "<< diam<<endl;
                             break;
                          }
                     }

                  }/*end for j particles*/

               if (flag==1) break;

               }//end Nz for
            if (flag==1) break;
           }//end Ny for
         if (flag==1) break;
        }//end Nx for
         if (flag==1) {flag=0; continue;};

         count_accept=count_accept+1; 
         /*delete old information*/
         int oldSquare = ((int)pos[i][3]+Nsquare*((int)pos[i][4]-1)+Nsquare*Nsquare*((int)pos[i][5]-1));
         int k=0;
         for (long int j=1; j<=NpartInSquare[oldSquare-1]-1; j++){
               k++;
               if (PartInSquarePos[oldSquare-1][j]==i) k++;
               PartInSquarePos[oldSquare-1][j]=PartInSquarePos[oldSquare-1][k];
         }

         PartInSquarePos[oldSquare-1][NpartInSquare[oldSquare-1]]=0;
         NpartInSquare[oldSquare-1]--;

         /*update with new information*/
         pos[i][0]=move[0];
         pos[i][1]=move[1];
         pos[i][2]=move[2];

         pos[i][3]=(double)Nx;
         pos[i][4]=(double)Ny;
         pos[i][5]=(double)Nz;
         
         direct[i][0]=wold[0];
         direct[i][1]=wold[1];
         direct[i][2]=wold[2];
         
         int newSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);
         NpartInSquare[newSquare-1]++;

         if (NpartInSquare[newSquare-1]>MAX_SPHERES){
            cerr<<"Too many spheres ("<<NpartInSquare[newSquare-1]<<") in square # "<< newSquare<<endl;
            exit(1);
         }

         PartInSquarePos[newSquare-1][NpartInSquare[newSquare-1]]=i;

      } /*end for i-Nsphere*/

   }/*end for sweep*/

   /*for (int i=0; i<Nsweep_cyl; i++)
      distFile<<timestep[i]<<" ";
    distFile<<endl;

   distFile.close();*/
   if (pas==0) {
   pas=passo;
   ds=dotstep;
   }
  
   //DestroyDoubleVector(timestep);
}//end function



void MC_spherocyl_SW(double diam, long int Nsphere, double Lato, double length, double **pos, double **direct, int * NpartInSquare, int ** PartInSquarePos, int Nsweep_cyl, double &pas, double &ds, double Tstar, double lambda){

   int Nsquare =(int)(0.99*Lato/(length+diam+2*lambda));
   double Lsquare = Lato/Nsquare;
   
   if (Nsquare<2){
      cerr<<"MC_spherocyl_SW: Too few squares! "<<"Nsquares = "<<Nsquare<< " of length "<<Lsquare<<endl;
      exit(1);
      
   }
   
   
   double cc=0, a, b, c;
   cerr<<"****************I'm beginning Monte Carlo*****************"<<endl;
   //getchar();
   sgenrand(time(NULL));

   //ofstream distFile;
   //distFile.open("distFile.dat");
   
   int flag=0;
   double move[3];
   double wold[3];
   double passo=0.2*diam;
   double dotstep=0.02;
   double count_accept=0;
   
   cerr<<"Passo = "<<pas<<" Dotstep = "<<ds<<endl;
   //int step_count=0;
   
    
   
   for (int step=1; step<=Nsweep_cyl; step++){
    
      if (pas==0){
       if(step%IRATIO == 0 && passo<Lato/2){
          
             double Ratio = (double)count_accept / (Nsphere * IRATIO);
            
              if( Ratio > 0.3 ){
                   passo=passo + passo*10/100;
                   dotstep=dotstep+dotstep*10/100;
                   }       
              else{
                 passo=passo - passo*10/100;
                 dotstep=dotstep-dotstep*10/100;    
                }
              
              count_accept = 0.0;
              cerr<<"Ratio "<<Ratio<<endl;    
        }
      }
      else {
      dotstep=ds;
      passo=pas;
      
      }
      
      
      for ( long int  i=1 ; i <= Nsphere; i++ ){

         long int numSquare=0;
        
         double sumEnergyNew=0;
         double sumEnergyOld=0;
         double diff;
         
         move[0]=pos[i][0]+ passo*(genrand()-0.5)*2;
         move[1]=pos[i][1]+ passo*(genrand()-0.5)*2;
         move[2]=pos[i][2]+ passo*(genrand()-0.5)*2;
         
         wold[0]=direct[i][0]; 
         wold[1]=direct[i][1]; 
         wold[2]=direct[i][2];
         
         flip(wold, dotstep);
         
         for(int j=0; j<3; j++){

            if(move[j]>Lato)
               move[j]=move[j]-Lato;
            if(move[j]<0.0)
               move[j]=move[j]+Lato;
         }

         int Nx= (int)(move[0]/Lsquare) +1;
         int Ny= (int)(move[1]/Lsquare) +1;
         int Nz= (int)(move[2]/Lsquare) +1;

         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;
         
      
   

         int count =0;
         
         for (int n1 =Nx-1; n1<=Nx+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
              for (int n2 =Ny-1; n2<=Ny+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =Nz-1; n3<=Nz+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }
                      count++;
                      //find position in numerized squares from 1 to Nsquare^3
                      numSquare=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);
                      //cerr<<"There are "<< NpartInSquare[numSquare-1] <<" neighbors in square = "<< numSquare <<endl;
                  for (long int j=1; j<=NpartInSquare[numSquare-1]; j++){
                     double  delta=0.0;

                     int k=PartInSquarePos[numSquare-1][j];
                     if (k==i){continue;}
                     double centredist=0.0;
                           //find distance of part i from each particle in the around squares
                     double rx=(move[0]-pos[k][0]+Lato*k1);
                     double ry=(move[1]-pos[k][1]+Lato*k2);
                     double rz=(move[2]-pos[k][2]+Lato*k3);
                     
                      centredist=rx*rx+ry*ry+rz*rz;
                      
                      if (centredist<=(length+diam+2*lambda)*(length+diam+2*lambda)){
                     
                          //double vr[]={fabs(rx), fabs(ry), fabs(rz)};
                          double vr[]={rx, ry, rz};
                          double w1[]={wold[0], wold[1], wold[2]};
                          double w2[]={direct[k][0], direct[k][1], direct[k][2]};
                     
                          double dist = dist2_rods(centredist, vr, w1, w2, length/2, length/2);
                           
                          
                           if (dist<diam*diam){
                              flag=1;

                             //cerr<<"Overlap..part i = "<<i<<" !! "<<delta<<" < "<< diam<<endl;
                             break;
                           }
                           else if (dist>=diam*diam && dist<= (diam+lambda)*(diam+lambda)){
                                      
                              sumEnergyNew=sumEnergyNew-1/Tstar;
                           }
                     
                     }
                     
                     
                  }/*end for j particles*/

               if (flag==1) break;

               }//end Nz for
            if (flag==1) break;
           }//end Ny for
         if (flag==1) break;
        }//end Nx for
         if (flag==1) {flag=0; continue;};
         
         //cerr << "-------------- No overlap " <<endl;
         // Calculate the energy relative to the old position
         int Nxold= (int)(pos[i][0]/Lsquare) +1;
         int Nyold= (int)(pos[i][1]/Lsquare) +1;
         int Nzold= (int)(pos[i][2]/Lsquare) +1;
    
         if (Nxold==Nsquare+1) Nxold=Nsquare;
         if (Nyold==Nsquare+1) Nyold=Nsquare;
         if (Nzold==Nsquare+1) Nzold=Nsquare;

      

         for (int n1 =Nxold-1; n1<=Nxold+1; n1++){
             int m1=n1;
             int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
             }
             for (int n2 =Nyold-1; n2<=Nyold+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                  for (int n3 =Nzold-1; n3<=Nzold+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquare=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);

                      //if (n1==Nx && n2==Ny && n3==Nz || NpartInSquare[numSquare]==0)
                      if ( NpartInSquare[numSquare-1]!=0){
 
                        for (long int j=1; j<=NpartInSquare[numSquare-1]; j++){
                          int k=PartInSquarePos[numSquare-1][j];
                          double centredist=0.0;
                          //find distance of part i from each particle in the around squares
                          double rx=(pos[i][0]-pos[k][0]+Lato*k1);
                          double ry=(pos[i][1]-pos[k][1]+Lato*k2);
                          double rz=(pos[i][2]-pos[k][2]+Lato*k3);
                          centredist=rx*rx+ry*ry+rz*rz;
                          
                          if (centredist<=(length+diam+2*lambda)*(length+diam+2*lambda)){
                              //double vr[]={fabs(rx), fabs(ry), fabs(rz)};
                              double vr[]={rx, ry, rz};
                              double w1[]={direct[i][0], direct[i][1], direct[i][2]};
                              double w2[]={direct[k][0], direct[k][1], direct[k][2]};
                           
                          
                              double dist = dist2_rods(centredist, vr, w1, w2, length/2, length/2);
                                                     
                              if (dist<diam*diam){
                                 cerr<<"Overlap..!! BIG PROBLEM"<<endl;
                                  exit(1);
                              }
                             else if (dist>=diam*diam && dist<= (diam+lambda)*(diam+lambda)){
                                       
                                 sumEnergyOld=sumEnergyOld-1/Tstar;
                              }
                              
                           }//end if centredist
                              
                        }//end for
                      }//end if numSquare
                    

                   }//end Nz for
               
                }//end Ny for
               
             }//end Nx for
        
                  
         diff=sumEnergyNew-sumEnergyOld;
         
         cc=genrand();
         
         if (diff<=0 || cc < exp(-diff)){
             count_accept=count_accept+1; 
             /*delete old information*/
             int oldSquare = ((int)pos[i][3]+Nsquare*((int)pos[i][4]-1)+Nsquare*Nsquare*((int)pos[i][5]-1));
             int k=0;
             for (long int j=1; j<=NpartInSquare[oldSquare-1]-1; j++){
                   k++;
                   if (PartInSquarePos[oldSquare-1][j]==i) k++;
                   PartInSquarePos[oldSquare-1][j]=PartInSquarePos[oldSquare-1][k];
             }

             PartInSquarePos[oldSquare-1][NpartInSquare[oldSquare-1]]=0;
             NpartInSquare[oldSquare-1]--;

             /*update with new information*/
             pos[i][0]=move[0];
             pos[i][1]=move[1];
             pos[i][2]=move[2];

             pos[i][3]=(double)Nx;
             pos[i][4]=(double)Ny;
             pos[i][5]=(double)Nz;
         
             direct[i][0]=wold[0];
             direct[i][1]=wold[1];
             direct[i][2]=wold[2];
         
             int newSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);
             NpartInSquare[newSquare-1]++;

             if (NpartInSquare[newSquare-1]>MAX_SPHERES){
                cerr<<"Too many spheres ("<<NpartInSquare[newSquare-1]<<") in square # "<< newSquare<<endl;
                exit(1);
             }

             PartInSquarePos[newSquare-1][NpartInSquare[newSquare-1]]=i;

        }//end if sumEnergy 

      } /*end for i-Nsphere*/

   }/*end for sweep*/

    
   if (pas==0) {
   pas=passo;
   ds=dotstep;
   }

}//end function

double check_order(double ** direct, int Nsphere, double length){
  
  double Xij, Yij, Zij, dist;
 
  double  ord=0, dot;
  double w1[]={0, 0, 0};
  double w2[]={0, 0, 0};
  
  
  for (int i=1; i<=Nsphere-1; i++){
     w1[0]= direct[i][0];
     w1[1]= direct[i][1];
     w1[2]= direct[i][2];
     for (int j=i+1; j<=Nsphere; j++){
     
         w2[0]= direct[j][0];
         w2[1]= direct[j][1];
         w2[2]= direct[j][2];
      
        dot=VecDot(w1,w2);
       
        ord = ord + (3*dot*dot-1)/2; 
             
     }
  }
  
  return ord=ord/(Nsphere*(Nsphere-1));
  
 
}




double dist2_rods(double normr12, double *r12, double *w1, double *w2,double lh1, double lh2)
{
 
 
 double xla,xmu;
 double rr=normr12;
  //rr= VecDot(r12),
 double rw1= VecDot(r12,w1);
 double rw2= VecDot(r12,w2);
 double w1w2= VecDot(w1,w2);
 double cc= 1-PW2(w1w2);

// Checking whether the rods are or not parallel:
// The original code is modified to have symmetry:

 if(cc<1e-6) {
  if(rw1 && rw2) {
   xla= rw1/2;
   xmu= -rw2/2;
  }
  else return rr;
 }

 else {

// Step 1

  xla= (rw1-w1w2*rw2)/cc;
  xmu= (-rw2+w1w2*rw1)/cc;
 }

// Step 2

if( fabs(xla)>lh1 || fabs(xmu)>lh2 ) {

// Step 3 - 7

  if(fabs(xla)-lh1>fabs(xmu)-lh2) {
   xla= SIGN(lh1,xla);
   xmu= xla*w1w2-rw2;
   if( fabs(xmu)>lh2 ) xmu= SIGN(lh2,xmu);
  }
  else {
   xmu= SIGN(lh2,xmu);
   xla= xmu*w1w2+rw1;
   if( fabs(xla)>lh1 ) xla= SIGN(lh1,xla);
  }
 }

// Step 8

 return rr+PW2(xla)+PW2(xmu) + 2*(xmu*rw2 -xla*(rw1+xmu*w1w2));
}

 void flip ( double *wold,  double dotstep)
 {

        double dot = 0.0;
        double xisq, xi1, xi2, xi;
        double wnew[]={0, 0, 0};
    while ((1.0 - dot) >= dotstep ){ 

           xisq = 1.0;
           while ( xisq >= 1.0 ){ 

              xi1  = (genrand())*2-1;
              xi2  = (genrand())*2-1;
              xisq = xi1 * xi1 + xi2 * xi2;
           }

           xi = sqrt(1.0-xisq);
           wnew[0] = 2.0 * xi1 * xi;
           wnew[1] = 2.0 * xi2 * xi;
           wnew[2] = 1.0 - 2.0 * xisq;
           dot= VecDot(wold,wnew); 

    }
    wold[0]=wnew[0];  
    wold[1]=wnew[1];
    wold[2]=wnew[2];
}


/*Find Radial Distribution Function g(r)*/

void update_hist_gr(double ** pos, int Nsphere, double *hist, double Lato, double deltar, int Nbin){

  int bin;
  double Xij, Yij, Zij, dist;


  for (int i=1; i<=Nsphere-1; i++){
     for (int j=i+1; j<=Nsphere; j++){
        Xij=pos[i][0]-pos[j][0];
        Xij-= Lato*rint(Xij/Lato);

        Yij=pos[i][1]-pos[j][1];
       Yij-=Lato*rint(Yij/Lato);

        Zij=pos[i][2]-pos[j][2];
        Zij-=Lato*rint(Zij/Lato);

        dist=Xij*Xij+Yij*Yij+Zij*Zij;
        dist=sqrt(dist);

        bin = (int) ( dist/deltar );
        if(bin<Nbin) hist[bin] = hist[bin] + 2;
     }
  }
  cerr<<"Update G(r) hist"<<endl;
 /* for (int i=0; i<Nbin; i++){
   cerr<<hist[i]<< " ";
  }
  cerr<<endl;*/
  //getchar();
}


void normal_write_gr(int Nsphere, double *hist, double Lato, double deltar, int Nbin, int realiz){

   double  vbin, nid;
   double rho=Nsphere/(Lato*Lato*Lato);
   ofstream grFile;
   grFile.open("grHomo.dat", fstream::app);
   grFile<<realiz<<" "<<Nsphere<<" "<<Lato<<" "<<Nbin<<" "<<deltar<<" "; 
  /* Normalize  g(r) */
  for (int i=0;i<Nbin;i++) {
    
    vbin=((i+1)*(i+1)*(i+1)-i*i*i)*deltar*deltar*deltar;
    nid=4*PI*vbin*rho/3;

   //xaxis = i*deltar,
   hist[i]=hist[i]/(Nsphere*nid*realiz);
   grFile <<hist[i]<<" "; 
  } grFile<<endl;

   grFile.close();
}

void normal_write_cc(int Nsphere, double length, double *hist, double Lato, double deltar, int Nbin, int realiz){

   double  vbin, nid;
   double rho=Nsphere/(Lato*Lato*Lato);
   ofstream grFile;
   grFile.open("grHomo_cc.dat", fstream::app);
   grFile<<realiz<<" "<<Nsphere<<" "<<length<<" "<<Lato<<" "<<Nbin<<" "<<deltar<<" "; 
  /* Normalize  g(r) */
  for (int i=0;i<Nbin;i++) {
    
    vbin=((i+1)*(i+1)*(i+1)-i*i*i)*deltar*deltar*deltar;
    nid=4*PI*vbin*rho/3;

   //xaxis = i*deltar,
   hist[i]=hist[i]/(Nsphere*nid*realiz);
   grFile <<hist[i]<<" "; 
  } grFile<<endl;

   grFile.close();
}
void update_hist_gr_angle(double ** direct, double **pos, int Nsphere, double *hist, double Lato, double deltar, int Nbin){
  
  double Xij, Yij, Zij, dist;
  int bin;
  double  ord, dot;
  double w1[]={0, 0, 0};
  double w2[]={0, 0, 0};
  int *count;
  double *num;
  count=allocate_int_vector(Nbin);
  num=allocate_double_vector(Nbin);
  for (int l=0; l<Nbin; l++){
            count[l]=0;
            num[l]=0;
  }
  
  for (int i=1; i<=Nsphere-1; i++){
     w1[0]= direct[i][0];
     w1[1]= direct[i][1];
     w1[2]= direct[i][2];
     for (int j=i+1; j<=Nsphere; j++){
     
         w2[0]= direct[j][0];
         w2[1]= direct[j][1];
         w2[2]= direct[j][2];

        Xij=pos[i][0]-pos[j][0];
        Xij-= Lato*rint(Xij/Lato);

        Yij=pos[i][1]-pos[j][1];
        Yij-=Lato*rint(Yij/Lato);

        Zij=pos[i][2]-pos[j][2];
        Zij-=Lato*rint(Zij/Lato);

        dist=Xij*Xij+Yij*Yij+Zij*Zij;
        dist=sqrt(dist);

        bin = (int) ( dist/deltar );
        
        dot=VecDot(w1,w2);
        //cerr<<dot*dot<<endl; getchar();
        ord =  (3*dot*dot-1)/2; 
        
        if(bin<Nbin){ 
           count[bin]=count[bin]+1; 
     
           num[bin] = num[bin]+1*ord;
        
        }
        
        
     }
  }
  cerr<<"Update G(r) hist"<<endl;
  
  for (int i=0; i<Nbin; i++){
   if (count[i]!=0)
        hist[i] = hist[i]+(double)num[i]/(count[i]);
   //cerr<<hist[i]<< " ";
  }
  //cerr<<endl;
  
 
  //cerr<<endl;
  //getchar();
  DestroyIntVector(count);
  DestroyDoubleVector(num);
}


void normal_write_gr_angle(int Nsphere, double length, double *hist, double Lato, double deltar, int Nbin, int realiz){

   double  vbin, nid;
   double rho=Nsphere/(Lato*Lato*Lato);
   ofstream grFile;
   grFile.open("grHomo_angle.dat", fstream::app);
   grFile<<realiz<<" "<<Nsphere<<" "<<length<<" "<<Lato<<" "<<Nbin<<" "<<deltar<<" "; 
  /* Normalize  g(r) */
  for (int i=0;i<Nbin;i++) {
    
    //vbin=((i+1)*(i+1)*(i+1)-i*i*i)*deltar*deltar*deltar;
    //nid=4*PI*rho*vbin/3;

   //xaxis = i*deltar,
   hist[i]=hist[i]/(realiz);
   grFile <<hist[i]<<" "; 
  } grFile<<endl;

   grFile.close();
}


void update_hist_gr_mindist(double ** direct, double **pos, int Nsphere, double length, double *hist, double Lato, double deltar, int Nbin){
  
  double Xij, Yij, Zij, dist;
  int bin;
  double  ord, dot;
  double w1[]={0, 0, 0};
  double w2[]={0, 0, 0};
  double vr[]={0, 0, 0};
  
  
  for (int i=1; i<=Nsphere-1; i++){
     w1[0]= direct[i][0];
     w1[1]= direct[i][1];
     w1[2]= direct[i][2];
     for (int j=i+1; j<=Nsphere; j++){
     
         w2[0]= direct[j][0];
         w2[1]= direct[j][1];
         w2[2]= direct[j][2];
         
        Xij=pos[i][0]-pos[j][0];
        Xij-= Lato*rint(Xij/Lato);

        Yij=pos[i][1]-pos[j][1];
        Yij-=Lato*rint(Yij/Lato);

        Zij=pos[i][2]-pos[j][2];
        Zij-=Lato*rint(Zij/Lato);
        vr[0]= Xij;
        vr[1]= Yij; 
        vr[2]= Zij;
        
        dist=Xij*Xij+Yij*Yij+Zij*Zij;
        
        double distTheta = dist2_rods(dist, vr, w1, w2, length/2, length/2);
        distTheta=sqrt(distTheta);
        
        
        bin = (int) ( distTheta/deltar );
        
        if(bin<Nbin) hist[bin] = hist[bin] + 2;
        
     }
  }
  cerr<<"Update G(r) hist"<<endl;
  
}


void normal_write_gr_mindist(int Nsphere, double length, double *hist, double Lato, double deltar, int Nbin, int realiz, string suffix){

   double  vbin2, vbin3, nid;
   double rho=Nsphere/(Lato*Lato*Lato);
   ofstream grFile;
   
   string file ("connectFile.dat");
   grFile.open((file+suffix).c_str(), fstream::app);
   
   
   grFile<<realiz<<" "<<Nsphere<<" "<<length<<" "<<Lato<<" "<<Nbin<<" "<<deltar<<" "; 
  /* Normalize  g(r) */
  for (int i=0;i<Nbin;i++) {
    
    vbin3=((i+1)*(i+1)*(i+1)-i*i*i)*deltar*deltar*deltar;
    vbin2=((i+1)*(i+1)-i*i)*deltar*deltar;
    nid=4*PI*rho*vbin3/3+ rho*PI*length*vbin2;

   //xaxis = i*deltar,
   hist[i]=hist[i]/(Nsphere*realiz);
   grFile <<hist[i]<<" "; 
  } grFile<<endl;

   grFile.close();
}


/*Nearest neighbor Distribution Function Hp(r)*/

void update_hist_Hp(double **pos, int ** AdjList, int Nsphere, double Lato, double *hist, double deltar, int Nbin){

  int bin;
  double Xij, Yij, Zij, dist=0;
  int count =0, ns=0, index;

  for (int i=1; i<=Nsphere; i++){
        ns++;
        double min=9282726;
        dist=0;
        int k=0;

        while(AdjList[i][k]!=0){

           Xij=pos[i][0]- pos[AdjList[i][k]][0];
           //Xij-= Lato*rint(Xij/Lato);

           Yij=pos[i][1]-pos[AdjList[i][k]][1];
           //Yij-=Lato*rint(Yij/Lato);

           Zij=pos[i][2]-pos[AdjList[i][k]][2];
           //Yij-=Lato*rint(Yij/Lato);

           dist=Xij*Xij+Yij*Yij+Zij*Zij;
           dist=sqrt(dist);


           if(dist<min ){
              min=dist;
              index=k; }
           k++;
        }

        bin = (int) ( min/deltar );
        if(bin<Nbin) hist[bin] = hist[bin] + 1;

  }

  //cerr<<"print Hp(r) hist"<<endl;
  cerr<<"print Hp(r) hist"<<endl;
  for (int i=0; i<Nbin; i++){
   cerr<<hist[i]<< " ";
  }
  cerr<<endl;
  //cerr<<endl <<"# part without neigh = "<<count<<" #part = "<<ns <<endl;
  //getchar();
}


void normal_write_Hp( int Nsphere, double *hist, double deltar, int Nbin, int realiz){


   ofstream HpFile;
   HpFile.open("Hp3d_impenetrable.dat");

  // Normalize  Hp(r)
  for (int i=0;i<Nbin;i++) {
     hist[i]=hist[i]/(Nsphere*realiz*deltar);
     HpFile << i*deltar<<" "<<hist[i]<<endl;
  }

   HpFile.close();
}


int  percolation (int Nsphere, int ** connect, int *clusterHK, double *count_perc, int i, int realiz){

   int max_el=0;
   for(int z=1; z<=Nsphere; z++){
      if (clusterHK[z]>max_el)
        max_el=clusterHK[z];
   }

   //cerr<<"max el="<<max_el<<endl;

   int cl=0, ind_cluster;
   for (int z=1; z<=max_el; z++){
      int count_el=0;
      for (int w=1; w<=Nsphere; w++){
         if (clusterHK[w]==z)
            count_el++;
      }
      if (count_el>cl){
            cl=count_el;
            ind_cluster=z;
       }

   }
   //cerr<<"ci sono ="<<cl<<" occorrenze di "<< ind_cluster<< endl;
   int perc_up=0, perc_down=0;
   for(int z=1; z<=Nsphere; z++){
      if(clusterHK[z]==ind_cluster){
         if(connect[z-1][0]==1)
            perc_up=1;
            //cerr<<"********* PERCOLATION SOPRA **********"<<endl;
         if(connect[z-1][1]==1)
            perc_down=1;
            //cerr<<"********* PERCOLATION SOTTO **********"<<endl;
         }
   }

   if (perc_up==1 && perc_down==1) {
      //count_perc[i]=count_perc[i]+1;
      return ind_cluster;
      }
   else return 0;

}


void percolation_wrap (int Nsphere, int ** connect, int *clusterHK, int *count_perc, int realiz, double** posCond, double Lato, double diam, double shell){

   int max_el=0;
   for(int z=1; z<=Nsphere; z++){
      if (clusterHK[z]>max_el)
        max_el=clusterHK[z];
   }

   //cerr<<"max el="<<max_el<<endl;

   int cl=0, ind_cluster;
   for (int z=1; z<=max_el; z++){
      int count_el=0;
      for (int w=1; w<=Nsphere; w++){
         if (clusterHK[w]==z)
            count_el++;
      }
      if (count_el>cl){
            cl=count_el;
            ind_cluster=z;
       }

   }
   //cerr<<"ci sono ="<<cl<<" occorrenze di "<< ind_cluster<< endl;
   int perc_up=0, perc_down=0;
   for(int z=1; z<=Nsphere; z++){
      perc_up=0;
      perc_down=0;

      if(clusterHK[z]==ind_cluster && connect[z-1][0]==1){
            perc_up=1;
            for(int j=1; j<=Nsphere; j++){
                if(connect[j-1][1]==1 && clusterHK[j]==ind_cluster){
                    double dist=pow((posCond[j][0]-posCond[z][0]-Lato),2)+pow((posCond[j][1]-posCond[z][1]),2)
                            +pow((posCond[j][2]-posCond[z][2]),2);
                    if (dist<diam+2*shell){
                         *count_perc=*count_perc+1;
                         return;
                    }
                }

            }


      }
   }



}



double adaptMC(double shell, double diamCond, int NsphereCond, double Lato, double **posCond, int * NpartInSquareCond, int ** PartInSquarePosCond, std::ofstream& percFile){

   int NsquareCond =(int)(0.98 *Lato/(diamCond));
   double LsquareCond = Lato/NsquareCond;
   double sumdist=0;
   double coord_num=0;

for ( long int  i=1 ; i <= NsphereCond; i++ ){
      double mindist=DBL_MAX;
      int NeighSquare3d=0;
      int countPart=0;
      long int numSquareCond=0;

      int NxCond= (int)(posCond[i][0]/LsquareCond)+1;
      int NyCond= (int)(posCond[i][1]/LsquareCond)+1;
      int NzCond= (int)(posCond[i][2]/LsquareCond)+1;

      if (NxCond==NsquareCond+1) NxCond=NsquareCond;
      if (NyCond==NsquareCond+1) NyCond=NsquareCond;
      if (NzCond==NsquareCond+1) NzCond=NsquareCond;

      for (int n1 =NxCond-1; n1<=NxCond+1; n1++){
               int m1=n1;
               int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=NsquareCond;
                k1=1;
               }
             if (n1>NsquareCond){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =NyCond-1; n2<=NyCond+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=NsquareCond;
                     k2=1;
                   }
                  if (n2>NsquareCond){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =NzCond-1; n3<=NzCond+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=NsquareCond;
                        k3=1;
                       }
                      if (n3>NsquareCond){
                        m3=1;
                        k3=-1;
                        }
                      double mindistquad=DBL_MAX;
                      //find position in numerized squares from 1 to Nsquare^3
                      numSquareCond=m1+NsquareCond*(m2-1)+NsquareCond*NsquareCond*(m3-1);

                      if ( NpartInSquareCond[numSquareCond-1]!=0){

                      for (long int j=1; j<=NpartInSquareCond[numSquareCond-1]; j++){
                           double dist;
                           int k=PartInSquarePosCond[numSquareCond-1][j];
                          //cerr<<"++++++Vicino  "<<k <<endl;
                           if (k==i) continue;

                            dist=pow((posCond[i][0]-posCond[k][0]+Lato*k1),2)+pow((posCond[i][1]-posCond[k][1]+Lato*k2),2)
                            +pow((posCond[i][2]-posCond[k][2]+Lato*k3),2);

                            dist=sqrt(dist);
                            if (dist>diamCond && dist<diamCond+2*shell)
                              countPart++;       ////////////////////DEFINIRE countPart ?????????*(*(**(*=*(Ã§=

                            if (dist<mindistquad)
                               mindistquad=dist;

                     }/*end for j particles*/
                      //cerr<<"----------scelgo "<<mindistquad<<endl;
                    }/*end if NpartInSquare*/


                    if (mindistquad<mindist)
                         mindist=mindistquad;


         }/*end Ny for*/
       }/*end Nx for*/
      }/*end Nz for*/
     //cerr<<"per la part i="<<i<<" "<<mindist<<endl;
     if (mindist==DBL_MAX){

         mindist=0;
         }
     sumdist=sumdist+mindist;

     coord_num=coord_num+countPart;

   } //end for i-Nsphere
   //cerr<<"***********dist finale= "<<sumdist/NsphereCond <<endl;
    //getchar();
   percFile<<sumdist/NsphereCond<<" "<<coord_num/NsphereCond<< endl;
   return sumdist/NsphereCond;

}

int **  createAdj(double shell, double diamCond, int NsphereCond, double Lato, double **posCond, int **connect){


   //int Nsquare =(int)(0.98*Lato/(diamCond+2));
   int Nsquare =(int)(0.98*Lato/(diamCond+2*shell));
   double Lsquare = Lato/Nsquare;

   int max_spheres=(int)((diamCond+2*shell)*(diamCond+2*shell)*(diamCond+2*shell)*27*0.58*6/PI);
   
   
   
   int ** AdjList;
   AdjList=allocate_int_matrix(NsphereCond+2, max_spheres);
   for (int l=0; l<NsphereCond+2; l++)
         for (int k=0; k<max_spheres; k++)
            AdjList[l][k]=0;
   
   
   int *NpartInSquare;
   NpartInSquare=allocate_int_vector(Nsquare*Nsquare*Nsquare);
   
   for (int l=0; l<Nsquare*Nsquare*Nsquare; l++)
            NpartInSquare[l]=0;

   int **PartInSquarePos;
   PartInSquarePos=allocate_int_matrix(Nsquare*Nsquare*Nsquare,max_spheres);
   

   for (int l=0; l<Nsquare*Nsquare*Nsquare; l++)
         for (int k=0; k<max_spheres; k++)
            PartInSquarePos[l][k]=0;
   


   for (int l=0; l<NsphereCond; l++)
         for (int k=0; k<2; k++)
            connect[l][k]=0;

 


    /*###############Update All Neighbors-Structures####################*/

    int count1=0;
    int count2=0;
    for ( long int  i=1 ; i <= NsphereCond; i++ ){
         long int numSquare=0;


         int Nx= (int)(posCond[i][0]/Lsquare) +1;
         int Ny= (int)(posCond[i][1]/Lsquare) +1;
         int Nz= (int)(posCond[i][2]/Lsquare) +1;


         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;

       /* cerr<<"Part "<<i<<endl;
        cerr<<"pos "<<posCond[i][0]<<endl;
        cerr<<(int)(posCond[i][0]/Lsquare)<<endl;
        cerr<<(int)(posCond[i][0]/Lsquare)+1<<endl;
        cerr<<"Nx = "<<Nx<<endl;
        getchar(); */
         numSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);

         NpartInSquare[numSquare-1]++;

         if (NpartInSquare[numSquare-1]>max_spheres){
            cerr<<"Too many spheres in square # "<< numSquare<<endl;
            cerr<<NpartInSquare[numSquare-1]<<endl;
            exit(1);
            }

         PartInSquarePos[numSquare-1][NpartInSquare[numSquare-1]]=i;
         //Update connect
         if (posCond[i][0]<=diamCond/2+shell){ connect[i-1][0]=1; count1++;};
         if (Lato-posCond[i][0]<=diamCond/2+shell ){ connect[i-1][1]=1; count2++;};

   } //end for i-Nsphere

    //cerr<<"trovati #contatti bordi "<< count1<<" & "<<count2<<endl;

   
   /*###############Create Adjacency List##########################*/
   for ( long int  i=1 ; i <= NsphereCond; i++ ){

      int NeighSquare3d=0;
      long int numSquareCond=0;
      int countPart=0;
      int countPart2=0;
      int NxCond= (int)(posCond[i][0]/Lsquare)+1;
      int NyCond= (int)(posCond[i][1]/Lsquare)+1;
      int NzCond= (int)(posCond[i][2]/Lsquare)+1;

      if (NxCond==Nsquare+1) NxCond=Nsquare;
      if (NyCond==Nsquare+1) NyCond=Nsquare;
      if (NzCond==Nsquare+1) NzCond=Nsquare;



      for (int n1 =NxCond-1; n1<=NxCond+1; n1++){
               int m1=n1;
               int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =NyCond-1; n2<=NyCond+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =NzCond-1; n3<=NzCond+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquareCond=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);

                      if ( NpartInSquare[numSquareCond-1]!=0){

                      for (long int j=1; j<=NpartInSquare[numSquareCond-1]; j++){
                          double dist=0.0;
                          int k=PartInSquarePos[numSquareCond-1][j];
                            if (k==i) continue;
                            
                            dist=pow((posCond[i][0]-posCond[k][0]),2)+pow((posCond[i][1]-posCond[k][1]+Lato*k2),2)
                            +pow((posCond[i][2]-posCond[k][2]+Lato*k3),2);
                            dist=sqrt(dist);
                           //dataFile<<dist<<" ";
                           if (dist>=diamCond && dist<=diamCond+2*shell){
                           countPart++;

                            //cerr<<k<<endl;
                              for (int index=0;index<max_spheres; index++){
                                   if (AdjList[i][index]==0){
                                       AdjList[i][index]=k;

                                       break;
                                       }
                              }

                      }

                 }/*end for j particles*/
              }/*end if NpartInSquare*/

              NeighSquare3d++;

         }/*end Ny for*/
       }/*end Nx for*/
      }/*end Nz for*/


   } //end for i-Nsphere

    DestroyIntMatrix(PartInSquarePos, Nsquare*Nsquare*Nsquare);
    DestroyIntVector(NpartInSquare);

  return AdjList;

}

void createAdj2(double shell, double diamCond, int NsphereCond, double Lato, double **posCond, int **connect, int **AdjList){


   int Nsquare =(int)(0.98*Lato/(diamCond+2));
   //int Nsquare =(int)(0.98*Lato/(diamCond+2*shell));
   double Lsquare = Lato/Nsquare;

   cerr << "Nsquare = "<< Nsquare << " Lsquare = " << Lsquare <<endl;
   cerr << "Nsquare*Lsquare = " <<  Nsquare*Lsquare << endl;

   int *NpartInSquare;
   NpartInSquare=allocate_int_vector(Nsquare*Nsquare*Nsquare);

   for (int l=0; l<Nsquare*Nsquare*Nsquare; l++)
            NpartInSquare[l]=0;

   int **PartInSquarePos;
   PartInSquarePos=allocate_int_matrix(Nsquare*Nsquare*Nsquare,MAX_SPHERES);
   cerr << "Nsquare = "<< Nsquare << " Lsquare = " << Lsquare <<endl;

   for (int l=0; l<Nsquare*Nsquare*Nsquare; l++)
         for (int k=0; k<MAX_SPHERES; k++)
            PartInSquarePos[l][k]=0;
   cerr << "Nsquare = "<< Nsquare << " Lsquare = " << Lsquare <<endl;


   for (int l=0; l<NsphereCond; l++)
         for (int k=0; k<2; k++)
            connect[l][k]=0;

   for (int l=0; l<NsphereCond+2; l++)
         for (int k=0; k<MAX_SPHERES; k++)
            AdjList[l][k]=0;



    /*###############Update All Neighbors-Structures####################*/

    int count1=0;
    int count2=0;
    for ( long int  i=1 ; i <= NsphereCond; i++ ){
         long int numSquare=0;


         int Nx= (int)(posCond[i][0]/Lsquare) +1;
         int Ny= (int)(posCond[i][1]/Lsquare) +1;
         int Nz= (int)(posCond[i][2]/Lsquare) +1;


         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;


         numSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);

         NpartInSquare[numSquare-1]++;

         if (NpartInSquare[numSquare-1]>MAX_SPHERES){
            cerr<<"Too many spheres in square # "<< numSquare<<endl;
            exit(1);
            }

         PartInSquarePos[numSquare-1][NpartInSquare[numSquare-1]]=i;
         //Update connect
         if (posCond[i][0]<diamCond/2+shell){ connect[i-1][0]=1; count1++;};
         if (Lato-posCond[i][0]<diamCond/2+shell ){ connect[i-1][1]=1; count2++;};

   } //end for i-Nsphere

    //cerr<<"trovati #contatti bordi "<< count1<<" & "<<count2<<endl;


   /*###############Create Adjacency List##########################*/
   for ( long int  i=1 ; i <= NsphereCond; i++ ){

      int NeighSquare3d=0;
      long int numSquareCond=0;

      int NxCond= (int)(posCond[i][0]/Lsquare)+1;
      int NyCond= (int)(posCond[i][1]/Lsquare)+1;
      int NzCond= (int)(posCond[i][2]/Lsquare)+1;



      for (int n1 =NxCond-1; n1<=NxCond+1; n1++){
               int m1=n1;
               int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =NyCond-1; n2<=NyCond+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =NzCond-1; n3<=NzCond+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquareCond=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);

                      if ( NpartInSquare[numSquareCond-1]!=0){

                      for (long int j=1; j<=NpartInSquare[numSquareCond-1]; j++){
                          double dist=0.0;
                          int k=PartInSquarePos[numSquareCond-1][j];
                            if (k==i) continue;
                            dist=pow((posCond[i][0]-posCond[k][0]),2)+pow((posCond[i][1]-posCond[k][1]),2)
                            +pow((posCond[i][2]-posCond[k][2]),2);
                            //dist=pow((pos[i][0]-pos[k][0]),2)+pow((pos[i][1]-pos[k][1]),2)
                            //   +pow((pos[i][2]-pos[k][2]),2);
                            dist=sqrt(dist);
                           //dataFile<<dist<<" ";
                           if (dist>diamCond && dist<diamCond+2*shell){
                                for (int index=0;index<MAX_SPHERES; index++){
                                   if (AdjList[i][index]==0){
                                       AdjList[i][index]=k;
                                       //cerr<<"AdjList["<<i<<"]["<<index<<"] = "<< AdjList[i][index]<<endl;
                                       break;
                                       }
                                 }
                            }
                 }/*end for j particles*/
              }/*end if NpartInSquare*/

              NeighSquare3d++;

         }/*end Ny for*/
       }/*end Nx for*/
      }/*end Nz for*/



   } //end for i-Nsphere

    DestroyIntMatrix(PartInSquarePos, Nsquare*Nsquare*Nsquare);
    DestroyIntVector(NpartInSquare);


}

int **  createAdjcylinder(double shell, double diamCond, double length, int NsphereCond, double Lato, double **posCond, double **direct, int **connect){


   
   int Nsquare =(int)(0.99*Lato/(length+2*shell+diamCond));
   double Lsquare = Lato/Nsquare;
    
    
   if (Nsquare<2){
   
      cerr<<"CreateAdjCyl: Too few squares! "<<"Nsquares = "<<Nsquare<< " of length "<<Lsquare<<endl;
      exit(1);
      
   }
   
   
   int max_spheres=MAX_SPHERES;
   
   
   
   int ** AdjList;
   AdjList=allocate_int_matrix(NsphereCond+2, max_spheres);
   for (int l=0; l<NsphereCond+2; l++)
         for (int k=0; k<max_spheres; k++)
            AdjList[l][k]=0;
   
   
   int *NpartInSquare;
   NpartInSquare=allocate_int_vector(Nsquare*Nsquare*Nsquare);
   
   for (int l=0; l<Nsquare*Nsquare*Nsquare; l++)
            NpartInSquare[l]=0;

   int **PartInSquarePos;
   PartInSquarePos=allocate_int_matrix(Nsquare*Nsquare*Nsquare,max_spheres);
   

   for (int l=0; l<Nsquare*Nsquare*Nsquare; l++)
         for (int k=0; k<max_spheres; k++)
            PartInSquarePos[l][k]=0;
   


   for (int l=0; l<NsphereCond; l++)
         for (int k=0; k<2; k++)
            connect[l][k]=0;

 


    /*###############Update All Neighbors-Structures####################*/

    int count1=0;
    int count2=0;
    for ( long int  i=1 ; i <= NsphereCond; i++ ){
         long int numSquare=0;


         int Nx= (int)(posCond[i][0]/Lsquare) +1;
         int Ny= (int)(posCond[i][1]/Lsquare) +1;
         int Nz= (int)(posCond[i][2]/Lsquare) +1;


         if (Nx==Nsquare+1) Nx=Nsquare;
         if (Ny==Nsquare+1) Ny=Nsquare;
         if (Nz==Nsquare+1) Nz=Nsquare;

       /* cerr<<"Part "<<i<<endl;
        cerr<<"pos "<<posCond[i][0]<<endl;
        cerr<<(int)(posCond[i][0]/Lsquare)<<endl;
        cerr<<(int)(posCond[i][0]/Lsquare)+1<<endl;
        cerr<<"Nx = "<<Nx<<endl;
        getchar(); */
         numSquare=Nx+Nsquare*(Ny-1)+Nsquare*Nsquare*(Nz-1);

         NpartInSquare[numSquare-1]++;

         if (NpartInSquare[numSquare-1]>max_spheres){
            cerr<<"Too many spheres in square # "<< numSquare<<endl;
            cerr<<NpartInSquare[numSquare-1]<<endl;
            exit(1);
            }

         PartInSquarePos[numSquare-1][NpartInSquare[numSquare-1]]=i;
         //Update connect
         if (posCond[i][0]+length/2*direct[i][0]<=diamCond/2+shell || posCond[i][0]-length/2*direct[i][0]<=diamCond/2+shell)
         { 
       //cerr<<"contact 0 Nodo "<<i<<" dist1: "<<posCond[i][0]+length/2*direct[i][0]<<" - dist2: " <<posCond[i][0]-length/2*direct[i][0]<<endl;      
            connect[i-1][0]=1; count1++;};
         if (Lato-(posCond[i][0]+length/2*direct[i][0])<=diamCond/2+shell || Lato-(posCond[i][0]-length/2*direct[i][0])<=diamCond/2+shell)
         { connect[i-1][1]=1; count2++;};

   } //end for i-Nsphere

    //cerr<<"trovati #contatti bordi "<< count1<<" & "<<count2<<endl;

   
   /*###############Create Adjacency List##########################*/
   for ( long int  i=1 ; i <= NsphereCond; i++ ){
   
      //if (i==1 || i==2|| i==3) cerr<<endl<<endl<<"node "<<i<<"-----> ";

      int NeighSquare3d=0;
      long int numSquareCond=0;
      int countPart=0;
      int countPart2=0;
      int NxCond= (int)(posCond[i][0]/Lsquare)+1;
      int NyCond= (int)(posCond[i][1]/Lsquare)+1;
      int NzCond= (int)(posCond[i][2]/Lsquare)+1;

      if (NxCond==Nsquare+1) NxCond=Nsquare;
      if (NyCond==Nsquare+1) NyCond=Nsquare;
      if (NzCond==Nsquare+1) NzCond=Nsquare;



      for (int n1 =NxCond-1; n1<=NxCond+1; n1++){
               int m1=n1;
               int k1=0;
             if (n1<1){ //i.e. Nx==1
                m1=Nsquare;
                k1=1;
               }
             if (n1>Nsquare){ //i.e. Nx==Nsquare
                m1=1;
                k1=-1;
               }
             for (int n2 =NyCond-1; n2<=NyCond+1; n2++){
                  int m2=n2;
                  int k2=0;
                  if (n2<1){
                     m2=Nsquare;
                     k2=1;
                   }
                  if (n2>Nsquare){
                     m2=1;
                     k2=-1;
                   }
                   for (int n3 =NzCond-1; n3<=NzCond+1; n3++){
                      int m3=n3;
                      int k3=0;
                      if (n3<1){
                        m3=Nsquare;
                        k3=1;
                       }
                      if (n3>Nsquare){
                        m3=1;
                        k3=-1;
                        }

                      //find position in numerized squares from 1 to Nsquare^3
                      numSquareCond=m1+Nsquare*(m2-1)+Nsquare*Nsquare*(m3-1);

                      if ( NpartInSquare[numSquareCond-1]!=0){

                      for (long int j=1; j<=NpartInSquare[numSquareCond-1]; j++){
                          
                          int k=PartInSquarePos[numSquareCond-1][j];
                            if (k==i) continue;
                           
                           double centredist=0.0;
                           //find distance of part i from each particle in the around squares
                           double rx=(posCond[i][0]-posCond[k][0]);
                           double ry=(posCond[i][1]-posCond[k][1]+Lato*k2);
                           double rz=(posCond[i][2]-posCond[k][2]+Lato*k3);
                           
                           centredist=rx*rx+ry*ry+rz*rz;
                           if (centredist<=(length+diamCond+2*shell)*(length+diamCond+2*shell)){
                            
                               //double vr[]={fabs(rx), fabs(ry), fabs(rz)};
                               double vr[]={rx, ry, rz};
                               double w1[]={direct[i][0], direct[i][1], direct[i][2]};
                               double w2[]={direct[k][0], direct[k][1], direct[k][2]};
                           
                               double dist = dist2_rods(centredist, vr, w1, w2, length/2, length/2);
                           
                               if (dist>=diamCond*diamCond && dist<=(diamCond+2*shell)*(diamCond+2*shell)){
                                  //if (dist>=1) {cerr<<sqrt(dist)<<" "<<endl; getchar();}
                                 //if (i==1 || i==2 || i==3 ) {
                                 //cerr<<endl<<" node "<<k<<"(x,y,z) = ("<<rx<<","<<ry<<","<<rz<<")"<<"DIST = "<<sqrt(dist);
                                 //cerr<<endl<<"DIST = "<<sqrt(dist);
                                 //}
                               
                                countPart++;

                                //cerr<<k<<endl;
                                  for (int index=0;index<max_spheres; index++){
                                       if (AdjList[i][index]==0){
                                           AdjList[i][index]=k;

                                           break;
                                       }
                                  }

                               }
                           }//end if centredist    

                 }/*end for j particles*/
              }/*end if NpartInSquare*/

              NeighSquare3d++;

         }/*end Ny for*/
       }/*end Nx for*/
      }/*end Nz for*/


   } //end for i-Nsphere

    DestroyIntMatrix(PartInSquarePos, Nsquare*Nsquare*Nsquare);
    DestroyIntVector(NpartInSquare);

  return AdjList;

}



