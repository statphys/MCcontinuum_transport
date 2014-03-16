#include "mylib.h"

void sort(int *array, int begin, int end, int *index) {
   /*** Use of static here will reduce memory footprint, but will make it thread-unsafe ***/
   static int pivot;
   static int t;     /* temporary variable for swap */

   if (end > begin) {
      int l = begin + 1;
      int r = end;

      swap(array[begin], array[pivot_index()], t); /*** choose arbitrary pivot ***/
      swap(index[begin], index[pivot_index()], t);

      pivot = array[begin];
      while(l < r) {
         if (array[l] <= pivot) {
            l++;
         } else {
            while(l < --r && array[r] >= pivot) /*** skip superfluous swaps ***/
               ;
            swap(array[l], array[r], t);
            swap(index[l], index[r], t);
         }
      }
      l--;
      swap(array[begin], array[l], t);
      swap(index[begin], index[l], t);

      sort(array, begin, l, index);
      sort(array, r, end, index);
   }
}

void sortDeque(deque <int> &array, int begin, int end, deque <double> &index) {
   /*** Use of static here will reduce memory footprint, but will make it thread-unsafe ***/
   static int pivot;
   static int t;     /* temporary variable for swap */
   static double tt;

   if (end > begin) {
      int l = begin + 1;
      int r = end;

      swap(array[begin], array[pivot_index()], t); /*** choose arbitrary pivot ***/
      swapDouble(&index[begin], &index[pivot_index()], tt);

      pivot = array[begin];
      while(l < r) {
         if (array[l] <= pivot) {
            l++;
         } else {
            while(l < --r && array[r] >= pivot) /*** skip superfluous swaps ***/
               ;
            swap(array[l], array[r], t);
            swapDouble(&index[l], &index[r], tt);
         }
      }
      l--;
      swap(array[begin], array[l], t);
      swapDouble(&index[begin], &index[l], tt);

      sortDeque(array, begin, l, index);
      sortDeque(array, r, end, index);
   }
}


void swapDouble(double * a, double * b, double tt){

  tt=*a; 
  *a=*b; 
  *b=tt;

}


 int ricercaSequenziale(int *array, int x, int n) {
     int i=0;

     while (i<n && x!=array[i])
          i++;
     if (i==n)
          return -1;
     return i;
 }

int ricercaSequenzialeDeque(deque <int> array, int x, int n) {
     int i=0;

     while (i<n && x!=array[i])
          i++;
     if (i==n)
          return -1;
     return i;
 }


void createCondStruct(double diamCond, int Nsphere, double Lato, double **pos, int *cluster, int **AdjList, int **connect, deque <deque <int> > &Adj, deque < deque <double> > &CondList, int * NumNeigh, int span_index, int *countPart){

 double Xij, Yij, Zij, dist=0;
 int countElectrodeNeighUp=0, countElectrodeNeighDown=0;

 for (int i=1; i<=Nsphere; i++)
 {

    int k=0;
    if (cluster[i]==span_index){


       while (AdjList[i][k]){

           Adj[i].push_back(AdjList[i][k]);

           Xij=fabs(pos[i][0]-pos[AdjList[i][k]][0]);
           Yij=fabs(pos[i][1]-pos[AdjList[i][k]][1]);
           Zij=fabs(pos[i][2]-pos[AdjList[i][k]][2]);

           Yij=min(Yij,Lato-Yij);
           Zij=min(Zij,Lato-Zij);

           dist=Xij*Xij+Yij*Yij+Zij*Zij;
           dist=sqrt(dist);

           CondList[i].push_back(exp(-2*(dist-diamCond)/XI));

           k++;

       }

    NumNeigh[i]=k;


    if (connect[i-1][0]==1){
      Adj[0].push_back(i);
      Adj[i].push_back(0);
      NumNeigh[i]=NumNeigh[i]+1;
      double espo=exp(-2*(pos[i][0]-diamCond/2)/XI);
      if (espo<1){
          CondList[0].push_back(espo);
          CondList[i].push_back(espo);
      }
      else{
          CondList[0].push_back(1);
          CondList[i].push_back(1);
      }
      countElectrodeNeighUp++;
      }


    if (connect[i-1][1]==1){
      Adj[Nsphere+1].push_back(i);
      Adj[i].push_back(Nsphere+1);
      NumNeigh[i]=NumNeigh[i]+1;
      double espo=exp(-2*(Lato-pos[i][0]-diamCond/2)/XI);
      if (espo<1){
         CondList[Nsphere+1].push_back(espo);
         CondList[i].push_back(espo);
      }
      else{
         CondList[Nsphere+1].push_back(1);
         CondList[i].push_back(1);
      }


      countElectrodeNeighDown++;
      }

     *countPart=*countPart+1;
    }//end if span_index



 }//end for Nsphere

    NumNeigh[0]=countElectrodeNeighUp;
    NumNeigh[Nsphere+1]=countElectrodeNeighDown;

    cluster[0]=span_index;
    cluster[Nsphere+1]=span_index;

}

void createCondStructCyl(double diamCond, double length, int Nsphere, double Lato, double **pos, double **direct, int *cluster, int **AdjList, int **connect, deque <deque <int> > &Adj, deque < deque <double> > &CondList, int * NumNeigh, int span_index, int *countPart){

 double Xij, Yij, Zij, dist=0;
 int countElectrodeNeighUp=0, countElectrodeNeighDown=0;

 for (int i=1; i<=Nsphere; i++)
 {

    int k=0;
    if (cluster[i]==span_index){

       //if (i==1 || i==2) cerr<<endl<<"************* "<<i<<" --> ";
       while (AdjList[i][k]){

           Adj[i].push_back(AdjList[i][k]);

           Xij=(pos[i][0]-pos[AdjList[i][k]][0]);
           Yij=(pos[i][1]-pos[AdjList[i][k]][1]);
           Zij=(pos[i][2]-pos[AdjList[i][k]][2]);

            Yij=MIN_RIJ(Yij, Lato);
            Zij=MIN_RIJ(Zij, Lato);
           //Yij=min(Yij-Lato,Yij+Lato);
           //Zij=min(Zij-Lato,Zij+Lato);

           double vr[]={Xij, Yij, Zij};
           
           double w1[]={direct[i][0], direct[i][1], direct[i][2]};
           
           double w2[]={direct[AdjList[i][k]][0], direct[AdjList[i][k]][1], direct[AdjList[i][k]][2]};
           
           

           dist=Xij*Xij+Yij*Yij+Zij*Zij;
           //dist=sqrt(dist);
          
           double distTheta = dist2_rods(dist, vr, w1, w2, length/2, length/2);
           
          
           
           distTheta=sqrt(distTheta);
           
      //if (i==1 || i==2 || i==3) cerr<<endl<<endl<<AdjList[i][k]<<", "<<"(x,y,z) = ("<<Xij<<","<<Yij<<","<<Zij<<")"<<" / "<< "DIST = "<<distTheta;
            //if (i==1 || i==2 || i==3 || i==4) {cerr<<endl<<distTheta<<" "; getchar();}
           CondList[i].push_back(exp(-2*(distTheta-diamCond)/XI));

           k++;

       } //if (i==1 || i==2) cerr<<endl;

    NumNeigh[i]=k;

    
         
  
    if (connect[i-1][0]==1){
      Adj[0].push_back(i);
      Adj[i].push_back(0);
      NumNeigh[i]=NumNeigh[i]+1;
      
      double a = pos[i][0]+length/2*direct[i][0]; 
      double b = pos[i][0]-length/2*direct[i][0];
       
      
      double dist_edge=min(a,b);

      double espo=exp(-2*(dist_edge-diamCond/2)/XI);
                    
      if (espo<1){
          CondList[0].push_back(espo);
          CondList[i].push_back(espo);
      }
      else{
          //espo=1;
          CondList[0].push_back(1);
          CondList[i].push_back(1);
      }
      
      //cerr<<"##### Nodo "<<i<<" dist min: "<<dist_edge<<" exp ="<<espo<<endl;   
      countElectrodeNeighUp++;
      }


    if (connect[i-1][1]==1){
      Adj[Nsphere+1].push_back(i);
      Adj[i].push_back(Nsphere+1);
      NumNeigh[i]=NumNeigh[i]+1;
      double a =  Lato-(pos[i][0]+length/2*direct[i][0]); 
      double b =  Lato-(pos[i][0]-length/2*direct[i][0]);
      
      
      
      double dist_edge=min(a,b);
      
      
      double espo=exp(-2*(dist_edge-diamCond/2)/XI);
      if (espo<1){
         CondList[Nsphere+1].push_back(espo);
         CondList[i].push_back(espo);
      }
      else{
         CondList[Nsphere+1].push_back(1);
         CondList[i].push_back(1);
      }


      countElectrodeNeighDown++;
      }

     *countPart=*countPart+1;
    }//end if span_index
    
    


 }//end for Nsphere

    NumNeigh[0]=countElectrodeNeighUp;
    NumNeigh[Nsphere+1]=countElectrodeNeighDown;

    cluster[0]=span_index;
    cluster[Nsphere+1]=span_index;

}

double decimation(double diamCond, int Nsphere, double Lato, deque <deque <int> > &Adj, deque < deque <double> > &CondList,  int * NumNeigh, int *NumNeighIndex, int span_index, int countPart, int max_decimated, double min_cond){

    /*-----------Decimation of all nodes without the electrodes -----------*/
    /*-----------Decimation starts from nodes with low connectivity--------*/
    time_t start=0, end=0;
    double dif=0;

   for (int i=Nsphere+1-countPart; i<=Nsphere; i++){

           if (dif>max_decimated) {
            return 0;}

          if (NumNeigh[i]==0) continue;
          else if (NumNeigh[i]==1){

              for (int j=0; j<Adj[NumNeighIndex[i]].size(); j++){

                 deque <int>::iterator it;
                 it=find(Adj[Adj[NumNeighIndex[i]][j]].begin(), Adj[Adj[NumNeighIndex[i]][j]].end(), NumNeighIndex[i]);

                 int position;
                 position = distance(Adj[Adj[NumNeighIndex[i]][j]].begin(),it);

                 /*Delete i from i's neighbor*/

                 Adj[Adj[NumNeighIndex[i]][j]].erase (it);
                 CondList[Adj[NumNeighIndex[i]][j]].erase(CondList[Adj[NumNeighIndex[i]][j]].begin()+position );

                 /*Delete node i, i.e. its lonely neighbors*/
                 Adj[NumNeighIndex[i]].erase(Adj[NumNeighIndex[i]].begin()+j);
                 CondList[NumNeighIndex[i]].erase(CondList[NumNeighIndex[i]].begin()+j);

                 /*Update the connectivity of i and of the neighbor of i*/
                 NumNeigh[i]=NumNeigh[i]-1;
                 int a=ricercaSequenziale(NumNeighIndex, Adj[NumNeighIndex[i]][j], Nsphere+2);
                 NumNeigh[a]=NumNeigh[a]-1;

              }//end for Adj[numNeigh]

           }//end else if NumNeigh[i]==1
           else {
                 time(&start);
                 double Gsum=0;
                 double gi_neighi=0;
                 double g_new=0;

                 for (int j=0; j<CondList[NumNeighIndex[i]].size(); j++){

                   Gsum=Gsum+CondList[NumNeighIndex[i]][j];
                 }

                 //check if each i's neighbors is also a neighbors of neighbors

                 for (int j=0; j<Adj[NumNeighIndex[i]].size(); j++){

                   gi_neighi=CondList[NumNeighIndex[i]][j];
                   for (int k=0; k<Adj[NumNeighIndex[i]].size(); k++){

                      if(Adj[NumNeighIndex[i]][j]!=Adj[NumNeighIndex[i]][k] ){

                        int a = ricercaSequenzialeDeque(Adj[Adj[NumNeighIndex[i]][k]], Adj[NumNeighIndex[i]][j],
                                +Adj[Adj[NumNeighIndex[i]][k]].size());

                        if (a==-1 && Adj[Adj[NumNeighIndex[i]][k]].size()!=0){
                           g_new=gi_neighi*CondList[NumNeighIndex[i]][k]/Gsum; 
                           if (g_new>0){
                              
                              Adj[Adj[NumNeighIndex[i]][k]].push_back(Adj[NumNeighIndex[i]][j]);
                              
                              CondList[Adj[NumNeighIndex[i]][k]].push_back(g_new);

                              int b=ricercaSequenziale(NumNeighIndex, Adj[NumNeighIndex[i]][k], Nsphere+2);
                              NumNeigh[b]=NumNeigh[b]+1;
                             }

                        }
                        else if (a!=-1 && Adj[Adj[NumNeighIndex[i]][k]].size()!=0){

                        g_new=gi_neighi*CondList[NumNeighIndex[i]][k]/Gsum;
                        CondList[Adj[NumNeighIndex[i]][k]][a]=CondList[Adj[NumNeighIndex[i]][k]][a]+g_new;

                        }



                     } //end if
                   }//end for k

                 }//end for j

               /*Free node i */
               for (int j=0; j<Adj[NumNeighIndex[i]].size(); j++){
                 int a =ricercaSequenzialeDeque(Adj[Adj[NumNeighIndex[i]][j]], NumNeighIndex[i], Adj[Adj[NumNeighIndex[i]][j]].size());

                   Adj[Adj[NumNeighIndex[i]][j]].erase(Adj[Adj[NumNeighIndex[i]][j]].begin()+a);
                   CondList[Adj[NumNeighIndex[i]][j]].erase(CondList[Adj[NumNeighIndex[i]][j]].begin()+a);
               }

               Adj[NumNeighIndex[i]].erase(Adj[NumNeighIndex[i]].begin(), Adj[NumNeighIndex[i]].end());
               CondList[NumNeighIndex[i]].erase(CondList[NumNeighIndex[i]].begin(), CondList[NumNeighIndex[i]].end());
               NumNeigh[i]=0;
               time(&end);
               dif = difftime (end,start);
                
              //cerr<<"Decimato nodo "<< i <<"time == "<<dif<<" seconds"<<endl; 

          }//end else
          

   }//End for Nsphere



   cerr<<"-------------------------Fine decimazione"<<endl;
    double conductance;
    
   
    int a = ricercaSequenzialeDeque(Adj[0], Nsphere+1, Adj[0].size());
    conductance=CondList[0][a];
    cerr<<"Conductance = "<<conductance<<endl;
    //getchar();
    return log10(conductance);


}



double createCGMstruct(double diamCond, int Nsphere, double Lato, deque <deque <int> > &Adj, deque < deque <double> > &CondList, double toll){
   

    cerr<<endl<<"Sono dentro createCGM"<<endl;
    cerr<<"toll = "<<toll<<endl;
    

   //int *IA=allocate_int_vector(Nsphere+2+1);
  // ofstream myfile;
   //myfile.open ("err.dat");
   //myfile.setf(ios::scientific,ios::floatfield);
   //myfile.precision(50);

   /****create vectors ia (rows) ja (coloumns) and sysmat (matrix) from the adjacency matrix CondList and Adj***/

   vector <int> IA(Nsphere+2+1);
   for (int l=0; l<=Nsphere+2; l++)
       IA[l]=0;
   
   int countLink=0;
   for (int i=0; i<=Nsphere+1; i++){
    
      IA[i]=countLink+i;
       
       for (int j=0; j<Adj[i].size(); j++){
          
          if (Adj[i][j]>i){   
             countLink=countLink+1;
          }
           
       }
   }
      
   IA[Nsphere+2]=countLink+Nsphere+2;
   
  
   
    
   vector <int> JA(countLink+Nsphere+2);
   for (int l=0; l<countLink+Nsphere+2; l++)
            JA[l]=0;
  
   vector <double> SYSMAT(countLink+Nsphere+2);
   for (int l=0; l<countLink+Nsphere+2; l++)
            SYSMAT[l]=0;


   for (int i=0; i<=Nsphere+1; i++){
       sortDeque(Adj[i], 0, Adj[i].size(), CondList[i]);
       int k=IA[i];
       JA[k]=i;
       int m=0;
       double sumG=0;
       for (int j=0; j<Adj[i].size(); j++){
            sumG=sumG+CondList[i][j];
            if (Adj[i][j]>i){
                JA[k+1+m]=Adj[i][j];
                SYSMAT[k+1+m]=-CondList[i][j];
                m=m+1;
            }
         SYSMAT[k]=sumG;
       }
   }
   
    
     
   /******adjust vectors ia ja and sys by removing the decimated nodes**********/
   
   

   for(int i=0; i<SYSMAT.size(); i++){
                   
           if (SYSMAT[i]==0){
               

               int a = JA[i];
              
               for(int j=0; j<JA.size(); j++){
                       if (JA[j]>a)
                           JA[j]--;
               }
             
               
               for(int k=a+1; k<IA.size(); k++){
                       IA[k]--;
               }
               
               
               IA.erase(IA.begin()+a);

               SYSMAT.erase(SYSMAT.begin()+i);
               JA.erase(JA.begin()+i);
               i--;
           }
   }
   
   int numNode = IA.size()-1;
   int numLink= JA.size();
   

   
   /*boundary condition: add g0A=1 to the first summation SYS[0][0] and gNB to the last SYS[N][N] */
   /*****We have introduced two fictious nodes A and B linked respectively to node O and N+1  *****/
   SYSMAT[0]=SYSMAT[0]+1;
   SYSMAT[numLink-1]=SYSMAT[numLink-1]+1;
   
  
  /*********create cholesky matrix and store it in vector chol***********/

   double *chol=allocate_double_vector(numLink);
   for (int l=0; l<numLink; l++)
            chol[l]=0;
    
    
   kersh(numNode, numLink, IA,JA, SYSMAT, chol);
  
   
   /*declare and initialize with boundary conditions the vector b (Ax=b): just set b[0]= g0a*Va = 1 */
   double *Bvec=allocate_double_vector(numNode);
   for (int l=0; l<numNode; l++)
            Bvec[l]=0;
   
   Bvec[0]=1;
   
   /**** residual vector init Ro=B-A*Vo****/
   double *Rvec=allocate_double_vector(numNode);
   for (int l=0; l<numNode; l++)
            Rvec[l]=0;

   /**** F=A*V ****/
   double *Fvec=allocate_double_vector(numNode);
   for (int l=0; l<numNode; l++)
            Fvec[l]=0;
   
   /**** S=clol(-1)*R into the while ****/
   double *Svec=allocate_double_vector(numNode);
   for (int l=0; l<numNode; l++)
            Svec[l]=0;
   
   /**** D = chol^(-1)*R ****/
   double *Dvec=allocate_double_vector(numNode);
   for (int l=0; l<numNode; l++)
            Dvec[l]=0;

   /**** Solution vector containing the potentials Vo V1 V2 ... Vn ****/
   double *Vvec=allocate_double_vector(numNode);
      
   for (int l=0; l<numNode; l++)
            Vvec[l]=0;

   

   /**************************/
   /*       start CGM        */
   /**************************/


   /**** initialize CGM: first iteration****/
  
   double eps2= toll*toll;
   //double eps2= 1e-10;
   double sum0=0, sum_new=0, sum_r=0, sum_old=0;

   /**** F = A*Vo ****/
   MatVecMult(numNode,numLink, IA , JA, SYSMAT, Vvec, Fvec);


   /**** Ro=Bo-AVo****/
   
   for (int k =0; k<numNode; k++){
       Rvec[k]=Bvec[k]-Fvec[k];
   }

  
   /**** d = chol(-1)*Ro ****/
   lsolve(numNode, numLink, IA,JA, chol, Rvec, Dvec);
   
   
   /****  R*D   ****/
   for (int k =0; k<numNode; k++){

       sum_new=sum_new+Dvec[k]*Rvec[k];
   }
   sum0=sum_new;   
   //cerr<<"sum0 = "<<sum0<<endl;
   
    double sum_nocg=0;
    for (int k =0; k<numNode; k++){

       sum_nocg=sum_nocg+Rvec[k]*Rvec[k];
   }

   


  double v0=0, vn=0;
   int l=0;

 
   /**********Iterate CGM untill precision is reached (EPS_CGM=1e-3) ***********/
   while (sum_new/sum0 > eps2 && l < ITER_MAX){
  
           
        
       //getchar();
       
       //cerr<<" sum "<<sum_new/sum0<<endl;
       MatVecMult(numNode,numLink, IA , JA, SYSMAT, Dvec, Fvec);

       
       sum_r=0;
       for (int k =0; k<numNode; k++){

           sum_r=sum_r+Dvec[k]*Fvec[k];
       }
       

       double alpha= sum_new/sum_r;
       

       for (int k =0; k<numNode; k++){

           Vvec[k]=Vvec[k]+alpha*Dvec[k];
           Rvec[k]=Rvec[k]-alpha*Fvec[k];
       }
       
       if (l!=0 && l%20==0){
            MatVecMult(numNode,numLink, IA , JA, SYSMAT, Vvec, Fvec);
            for (int k =0; k<numNode; k++){
                Rvec[k]=Bvec[k]-Fvec[k];
             }
 
            }


        v0=Vvec[0];
        vn=Vvec[numNode-1];
        

       lsolve(numNode, numLink, IA,JA, chol, Rvec, Svec);
              
       sum_old=sum_new;

       sum_new=0;
       for (int k =0; k<numNode; k++){

           sum_new=sum_new+Rvec[k]*Svec[k];
       }
     
       double beta=sum_new/sum_old;
        
     

       for (int k =0; k<numNode; k++){

           Dvec[k]=Svec[k]+beta*Dvec[k];
       }
      
       l++;
       

   }//end while

   /**** extract conductance from bond linked to the ending node ****/

   //myfile.close();

   cout.precision(5);
   double G=0;
   cerr<<endl;
   cerr<<"Numero iterazioni = "<<l<<endl;

   double F=0;
   for(int j=1; j<CondList[0].size()+1; j++){
       F=F+(Vvec[0]-Vvec[JA[j]])*CondList[0][j];
   }
   cerr<<endl;

   cerr<<" FFF = "<<F<<endl;

   cerr<<"V0N = "<<Vvec[0]<<" - "<<Vvec[numNode-1]<<endl;
   G=0.5*(1/(Vvec[0]-Vvec[numNode-1])-1);
   cerr<<" G= 0.5(1/Von-1) = "<<G<<endl;

   G=(Vvec[numNode-1]/(1-2*Vvec[numNode-1]));
   cerr<<" G = Vn/(1-2Vn) = "<<G<<endl;

   G=(Vvec[numNode-1]/(Vvec[0]-Vvec[numNode-1]));
   cerr<<" G = Vn/Von = "<<G<<endl;

   
   cerr<<" logG ="<<log10(fabs(G))<<endl;
   
   //getchar();
   DestroyDoubleVector(Vvec);  
   DestroyDoubleVector(Svec);
   DestroyDoubleVector(Fvec); 
   DestroyDoubleVector(Dvec); 
   DestroyDoubleVector(Rvec); 
   DestroyDoubleVector(Bvec);    
   DestroyDoubleVector(chol);
    //getchar();
    if(G!=0) 

      return log10(fabs(G));
     else return 0;


}  



//soubrutine che calcola la decomposto incompleta di Cholesky secondo Kershaw
void kersh(int nequ,int nterm, vector <int> ia , vector <int> ja,  vector <double> sysmat, double *chol)
{

/*
neq = n
nterm = nt
chol: vettore per la memorizzazione compatta della decomposta incompleta
*/

int i,j,k,kk,k1,k2,i1,j1;
long double a;


//initialize the output vector chol wich conzains the triangular matrix in compact form
for (k=0;k<nterm; k++)
{
    chol[k]=0;
}


//calculate the cholesky matrix with incomplete factor: if Aij=0 => CHOLij=0
for (kk=0;kk<=(nequ-2);kk++)
  {
                            
    k=ia[kk];
    
    a=sysmat[k]-chol[k];
    
    
    if (a<=0)
    {
     
    a=chol[ia[kk-1]]*chol[ia[kk-1]];
  
    }
    
    chol[k]=sqrt(a);
    
    i=ia[kk]+1;
    j=ia[kk+1]-1;

    for (k1=i;k1<=j;k1++)
    {   
        
        chol[k1]=(sysmat[k1]-chol[k1])/chol[k];
        
    }
    
    for (k2=i;k2<=j-1;k2++)
    {
        j1=ia[ja[k2]];
        chol[j1]=chol[j1]+(chol[k2]*chol[k2]);
        i1=k2+1;
        j1++;
        while ((j1<ia[ja[k2]+1])&&(i1<=j))
        {
              if (ja[j1]==ja[i1])
              {
                 chol[j1]=chol[j1]+chol[k2]*chol[i1];
                 i1++;
                 j1++;
              }
              if (ja[j1]<ja[i1])
              {
                 j1++;
              }
              if (ja[j1]>ja[i1])
              {
                 i1++;
              }
        }
    }
       
   
    if (j>=i)
    {
       chol[ia[ja[j]]]=chol[ia[ja[j]]]+(chol[j]*chol[j]);
    }
}

k=ia[nequ-1];

a=sysmat[k]-chol[k];


if (a<=0)
{
   a=(chol[ia[nequ-2]]*chol[ia[nequ-2]]);
}

chol[k]=sqrt(a);


return;
}

void lsolve(int nequ,int nterm, vector <int> ia , vector <int> ja,  double *lmat, double *Bvec, double *Vvec){

//   Vvec:=(LL^T)^-1*Bvec (solves the system LL^T*Vvec=Bvec)

     int i=0,j=0, n1=0, mm=0;
     double a=0;

      for( int k=0; k<nequ; k++)
         Vvec[k] = 0;
      
      for( int k=0; k<nequ; k++){
         i = ia[k];
         j = ia[k+1]-1;
         Vvec[k] = (Bvec[k] - Vvec[k])/lmat[i];
         for( int m = i+1; m<=j; m++){
            Vvec[ja[m]] = Vvec[ja[m]] + lmat[m]*Vvec[k];
         }
      }
      for (int k=0; k<nequ; k++){ 
         n1 = nequ-k-1;
         a = 0;
         i = ia[n1];
         j = ia[n1+1] - 1;
         for (int m=i; m<=j-1; m++){
            mm = j-m+i;
            a = a + lmat[mm]*Vvec[ja[mm]];
         }
         Vvec[n1] = (Vvec[n1] - a)/lmat[i];
      }
return;
}


void lsolve2(int n,int nt, vector <int> ia , vector <int> ja, double *chol, double *v, double *w)
{
int k,i,j,k1;
long double s[n+1],z[n+1];

//algoritmo per il calcolo vettore z

//put to 0 the output vector w and the intermediary vector s
for (j=0;j<=n;++j)
{
    w[j]=0;
    s[j]=0;
}

//calulate z
z[0]=v[0]/chol[0];
for (i=1;i<=n-1;i++)
{
    k=ia[i];
    for (k1=ia[i-1]+1;k1<=ia[i]-1;k1++)
    {
        j=ja[k1];
        s[j]=s[j]+(chol[k1]*z[i-1]);
    }
    z[i]=(v[i]-s[i])/chol[k];
}

//calculate w

w[n-1]=z[n-1]/chol[nt];
for (i=n-2;i>=0;i--)
{
//initialize the vector s to 0
for (j=0;j<=n;j++)
{
    s[j]=0;
}
    k=ia[i];
    for (k1=ia[i]+1;k1<=ia[i+1]-1;++k1)
    {
        j=ja[k1];
        s[i]=s[i]+(chol[k1]*w[j]);
    }
    w[i]=(z[i]-s[i])/chol[k];
}

return;
}

void MatVecMult(int n,int nterm, vector <int> ia , vector <int> ja,  vector <double> sys, double *Bvec, double *Vvec){

   for(int i=0; i<n; i++){
      Vvec[i]=0;
   }
  

   for(int i=0; i<n; i++){
      int m=ia[i];
      Vvec[i]=Vvec[i]+sys[m]*Bvec[i];
     
      for (m =ia[i]+1; m<=ia[i+1]-1; m++){
          int j=ja[m];   
          Vvec[i]=Vvec[i]+sys[m]*Bvec[j];
          Vvec[j]=Vvec[j]+sys[m]*Bvec[i];        
      }
   
   }

}

void calcConductivity(double diamCond, int Nsphere, double Lato, double **pos, int *cluster, int **AdjList, int **connect, int span_index, double *conductance, double toll, double min_cond){

   int countPart=0;
   int max_decimated=MAX_DECIMATED;
   //int max_decimated=0;
   //if ((PI/6)*Nsphere/(Lato*Lato*Lato) <0.1)
      // max_decimated=0;
    
   
   
   int *NumNeigh;
   NumNeigh=allocate_int_vector(Nsphere+2);
   for (int l=0; l<Nsphere+2; l++)
            NumNeigh[l]=0;

   int *NumNeighIndex;
   NumNeighIndex=allocate_int_vector(Nsphere+2);
   for (int l=0; l<Nsphere+2; l++)
            NumNeighIndex[l]=l;


   deque< deque<int> > Adj(Nsphere+2, deque<int>(0,-1));
   deque< deque<double> > CondList(Nsphere+2, deque<double>(0,-1));
   createCondStruct(diamCond,Nsphere, Lato, pos, cluster, AdjList, connect, Adj, CondList, NumNeigh, span_index, &countPart);

  
   sort(NumNeigh, 1, Nsphere+1, NumNeighIndex);
  
  
        //cerr<<"Sono qui "<<endl; getchar();
      *conductance=*conductance+decimation(diamCond, Nsphere, Lato, Adj, CondList, NumNeigh, NumNeighIndex, span_index, countPart, max_decimated, min_cond);
       
      if (*conductance==0){
          cerr<<"Entro in CGM "<<endl;
          *conductance=*conductance+createCGMstruct(diamCond, Nsphere, Lato, Adj, CondList,toll);
      }
      //getchar();
 


   DestroyIntVector(NumNeighIndex);
   DestroyIntVector(NumNeigh);
  
}

void calcConductivityCyl(double diamCond, double length, int Nsphere, double Lato, double **pos, double **direct, int *cluster, int **AdjList, int **connect, int span_index, double *conductance, double toll, double min_cond){

   int countPart=0;
   int max_decimated=MAX_DECIMATED;
   //int max_decimated=0;
   //if ((PI/6)*Nsphere/(Lato*Lato*Lato) <0.1)
      // max_decimated=0;
   if (min_cond<1e-10){
       
      max_decimated=10;
   }
   
   int *NumNeigh;
   NumNeigh=allocate_int_vector(Nsphere+2);
   for (int l=0; l<Nsphere+2; l++)
            NumNeigh[l]=0;

   int *NumNeighIndex;
   NumNeighIndex=allocate_int_vector(Nsphere+2);
   for (int l=0; l<Nsphere+2; l++)
            NumNeighIndex[l]=l;


   deque< deque<int> > Adj(Nsphere+2, deque<int>(0,-1));
   deque< deque<double> > CondList(Nsphere+2, deque<double>(0,-1));
   createCondStructCyl(diamCond, length, Nsphere, Lato, pos, direct, cluster, AdjList, connect, Adj, CondList, NumNeigh, span_index, &countPart);

  
   sort(NumNeigh, 1, Nsphere+1, NumNeighIndex);
  
  
        //cerr<<"Sono qui "<<endl; getchar();
      *conductance=*conductance+decimation(diamCond, Nsphere, Lato, Adj, CondList, NumNeigh, NumNeighIndex, span_index, countPart, max_decimated, min_cond);
       
      if (*conductance==0){
          cerr<<"Entro in CGM "<<endl;
          *conductance=*conductance+createCGMstruct(diamCond, Nsphere, Lato, Adj, CondList,toll);
      }
      //getchar();
 


   DestroyIntVector(NumNeighIndex);
   DestroyIntVector(NumNeigh);
  
}


