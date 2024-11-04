#include "TimeScheme.h"

int main(int argc, char** argv)
{

   if (argc < 2)
   {
      printf("Please, enter the name of your data file.\n");
      exit(0);
   }
   const std::string data_file_name = argv[1];

   DataFile* df = new DataFile(data_file_name);

   Function* fct = new Function(df);

   Laplacian* lap = new Laplacian(fct, df);

   TimeScheme* ts = new ImplicitScheme(df, lap);

   std::vector<double> U((df->Get_Nx()-1)*(df->Get_Ny()-1),0.0);

   lap->InitialCondition(U);
   double t(0.0);
   int it(0);
   ts->SaveSol(U,"sol",it);
   // ts->SaveSol(lap->ExactSol(t),"exact",it);
   ++it;

   while (t < df->Get_Tf()){
      ts->Integrate(t, U);
      ts->SaveSol(U,"sol",it);
      // ts->SaveSol(lap->ExactSol(t),"exact",it);
      std::cout << "t = " << t << std::endl;
      ++it;
   }

   // ------------------- Pour valider de l'ordre du schéma ----------------

   double Erreur(0.0), Normalise(0.0), ErreurNorm(0.0);
   std::vector<double> ExacteSol(U.size());

   ExacteSol = lap->ExactSol(t);

   for(int k = 0; k < (df->Get_Nx()-1)*(df->Get_Ny()-1); ++k){
      Erreur += (U[k] - ExacteSol[k])*(U[k] - ExacteSol[k]);
      Normalise += ExacteSol[k]*ExacteSol[k];
   }

   ErreurNorm = sqrt(Erreur/Normalise);

   std::cout << "log(Erreur) = " << log(ErreurNorm) << " " << "log(dx) = " << log(df->Get_dx()) << std::endl;

   // ----------------------------------------------------------------------

   delete df, delete fct, delete lap, delete ts;

   return 0;
}
