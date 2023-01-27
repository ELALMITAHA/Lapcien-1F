#include "laplacien.hpp"



//--------------------------------------------------------------------//
//***************Constructeur ****************************************//
//Argument d'entrer :xmin borne inf  du domaine 
//                   xmax borne sup du domaine 
//	             N    Taille du maillage 
//                   h    pas en espace (constant car maillage uniforme)
//Argument de sortie:void 
//********************************************************************//

laplacien1d::laplacien1d(double xmin,double xmax,int monN){

	 xmin=0;
	 xmax=1;
	 N=monN;
	 h=(xmax-xmin)/(N+1);
}


//----------------------------------------------------------------------//
//**************************** Matrice du laplacien*********************//
//Argument d'entrer   :void                                             //
//Argument de sortie  :void                                             //
//----------------------------------------------------------------------//

void laplacien1d::Matlaplacien(){
	A.resize(N-1,N-1);
	for (int i(0);i<N-1;++i)
	{
		for (int j(0);j<N-1;++j)
		{
			if (i==j)
				A(i,j)=2;
			else if ((i==j+1)||(i+1==j))
				A(i,j)=-1;
			else
				A(i,j)=0;
		}
	}
	A=(1/pow(h,2))*A;
}



//-----------------------------------------------------------------------//
//******************Second Membre ***************************************//
//Argument d'entrer  :Void 
//Argument de sortie :Void 
//-----------------------------------------------------------------------//

void   laplacien1d::Termesource(){
	
	X.resize(N+1);
	F.resize(N-1);
	for(int i(0);i<N+1;++i)
		{X(i)=xmin+i*h;

		}
	for (int i=0;i<N-1;++i){
		F(i)=(-pow(X(i+1),2)+5*X(i+1)-4)*exp(-X(i+1));
	}
}
//------------------------------------------------------------------------//
//****************Méthode du gradient conjugué****************************//
//Argument d'entrer  :Void 
//Argument de sortie :Void 
//------------------------------------------------------------------------//

/*void laplacien1d::gradientconj(){

	U.resize(N-1);
	for(int i(0);i<N-1;++i)
		U(i)=0;
	double eps(0.001);
	int kmax(1000);
	Eigen::VectorXd d,r,w,r1;
	d.resize(N-1);r.resize(N-1);w.resize(N-1);r1.resize(N-1);
	double alpha,beta;
	int l(0);
	r=A*U-F;
	d=r;
	
	while(l<=kmax||r.norm()>eps){
		w=A*d;
		alpha=(d.dot(r))/(d.dot(w));
		U=U-alpha *d;
		r1=r-alpha*w;
		beta=(r1.dot(r1))/(r.dot(r));
		d=r1+beta*d;		l=l+1;
		r=r1;
	}
	cout<<U<<endl;


}
*/

//******************Solveur directe pour AU=F ****************************//
//Arguement d'entrer : A Matrice du problème (A doit etre inversible )
//                   : F vecteur second membre             
//Argument de sortie : U Vecteur Soltion du systeme lineaire 
//La résolution se fera à l'aide la méthode Lu 
//------------------------------------------------------------------------//

void  laplacien1d::solveur(){
	
	U.resize(N-1);
	int k;
	cout<<"\n";
	cout<<" Ce programme offre deux possibilites pour resoudre le systeme lineaire AU=F."<<endl;
	cout<<"\n";
	cout <<" Une methode directe (LU) et une methode iterative (Jacobie) :"<<endl ;
	cout<<"\n";
	cout<<" Tapez 0 pour la methode directe  "<< ""<<"1 pour la methode iterative " <<endl;
	cout<<"\n";
	cin>>k;
	switch(k){
		case 0:
		U=A.lu().solve(F);
		break;
		case 1:
		U=A.jacobiSvd(ComputeThinU | ComputeThinV).solve(F);
		break;
		default :		cout<<"Je ne comprends pas votre choix  "<<endl;
	}
}


//---------------------------------------------------------------------------//
//*********************Sauvgarde de la soltion*******************************///
//Arguemnt d'entrer :Le nom de fuchier pour sauvgarder la solution           //
//


void laplacien1d:: Save(std::string solution ){
	ofstream fichier("solution.txt", ios::out | ios::trunc);  //déclaration du flux et ouverture du fichier
        
        if(fichier)  // si l'ouverture a réussi
        {
             for(int i(0); i<N-1;i++){
             	fichier<<X(i+1)<<" " <<U(i)<< " "<<(X(i+1)-1)*X(i+1)*exp(-X(i+1))<<endl;
             }
                fichier.close();  // on referme le fichier
        }
        else  // sinon
                cerr << "Erreur à l'ouverture !" << endl;
  }



  //-----------------------------------------------------------------------//
  //***************Calcule de l'erreur de convergence**********************//
  //-----------------------------------------------------------------------//
  //Argument d'entrer  : void    
  //Argument de sortie : Erreur=max_i |U(i)-u(X(i))| 
  //------------------------------------------------------------------------// 

  double laplacien1d:: CalcErreur(){
  	Eigen::VectorXd temp(N-1);

  	for(int i(0);i<N-1;++i){
  		temp(i)=std::abs(U(i)-(X(i+1)-1)*X(i+1)*exp(-X(i+1)));
  	}
  	return temp.maxCoeff();
  }
//-------------------------------------------------------------------------//
//*************Calcule le temps d'execution de notre programme*************//
//-------------------------------------------------------------------------//
//Argument d'entrer : instant initial
//Argument de sortie: instant final 
//la fonction affiche le temps que s'est écoulé entre les deux          
//------------------------------------------------------------------------//



void laplacien1d::duree(time_t _begin, time_t _end) 
{
  double temp; 
  double hours=0, min=0, sec=0; 
  double dureeCalc = difftime(_end, _begin);
  temp = modf(dureeCalc/3600., &hours); 
  temp = modf(temp*60., &min); 
  temp = modf(temp*60., &sec); 
  cout<<"Duree du calcul : "<<hours<<" h "<<min<<" min "<<sec<<" sec"<<endl; 
}

//-------------------------------------------------------------------------//
