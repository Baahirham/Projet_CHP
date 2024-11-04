#ifndef _FUNCTION_H

#include "DataFile.h"

class Function {
private:
   // Pointeur de la classe DataFile pour récupérer toutes les
   // valeurs de paramètres
   const DataFile* _df;
   // Diffusion coefficient

   public: // Méthodes et opérateurs de la classe
   Function(DataFile* df);
   double Exact_solution(const double x, const double y, const double t) const;
   double Initial_condition(const double x, const double y) const;
   double Source(const double x, const double y, const double t) const;
   double Dirichlet_Gamma_0(const double x, const double y, const double t) const;
   double Dirichlet_Gamma_1(const double x, const double y, const double t) const;
};

#define _FUNCTION_H
#endif
