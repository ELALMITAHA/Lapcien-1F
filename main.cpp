#include "laplacien.hpp" 




int main(){
	
	double t;
	int Nombre_point;
	std::string fichier;
	



	laplacien1d L(0,1,500);
	L.Matlaplacien();
	L.Termesource();
	//L.gradientconj();
	time_t begin=time(NULL);
	L.solveur();
	time_t end=time(NULL); 
    L.duree(begin,end); 
	L.Save(fichier);
	
	

	cout<<L.CalcErreur()<<endl;

	
	

	return 0;
}