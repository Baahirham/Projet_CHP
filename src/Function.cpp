#ifndef _FUNCTION_CPP

#include "Function.h"
#include <cmath>

Function::Function(DataFile* df) :
_df(df)
{
   
}

double Function::Initial_condition(const double x, const double y) const
{
   if (this->_df->Get_cas() == 0)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 1)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 2)
   {  
      return 0.0;
   }
   else if (this->_df->Get_cas() == 3)
   {  
      return 0.0;
   }
   else
   {
      return 0.0;
   }
}

double Function::Exact_solution(const double x, const double y, const double t) const
{
   if (_df->Get_cas() == 0)
   {
      double pi(std::acos(-1.0));
      double xmax(_df->Get_xmax()), ymax(_df->Get_ymax()), xmin(_df->Get_xmin()), ymin(_df->Get_ymin()); 
      return t*sin(2.0*pi*x/(_df->Get_xmax()-_df->Get_xmin()))*sin(2.0*pi*y/(_df->Get_ymax()-_df->Get_ymin()));
   }
   else if (this->_df->Get_cas() == 1)
   {
      return x*(1-x)*y*(1-y);
   }
   else if (this->_df->Get_cas() == 2)
   {  
      return sin(x)+cos(y);
   }
   else if (this->_df->Get_cas() == 3)
   {  
      return 0.0;
   }
   else
   {
      return 0.0;
   }
}

double Function::Source(const double x, const double y, const double t) const
{
   double pi(std::acos(-1.0));
   if (_df->Get_cas() == 0)
   {
      double xmax(_df->Get_xmax()), ymax(_df->Get_ymax()), xmin(_df->Get_xmin()), ymin(_df->Get_ymin()); 
      return sin(2.0*pi*x/(_df->Get_xmax()-_df->Get_xmin()))*sin(2.0*pi*y/(_df->Get_ymax()-_df->Get_ymin()))*(1.0 + t*_df->Get_D()*4.0*pi*pi*(1.0/((xmax - xmin)*(xmax - xmin)) + 1.0/((ymax - ymin)*(ymax - ymin))));
   }
   else if (this->_df->Get_cas() == 1)
   {
      return 2.0*(x - x*x + y - y*y);
   }
   else if (this->_df->Get_cas() == 2)
   {  
      return sin(x)+cos(y);
   }
   else if (this->_df->Get_cas() == 3)
   {  
      double Lx(_df->Get_xmax()-_df->Get_xmin()), Ly(_df->Get_ymax()-_df->Get_ymin());
      return exp(-(x-Lx/2.0)*(x-Lx/2.0))*exp(-(y-Ly/2.0)*(y-Ly/2.0))*cos((pi/2.0)*t);
   }
   else
   {
      return 0.0;
   }
}

double Function::Dirichlet_Gamma_0(const double x, const double y, const double t) const
{
   if (_df->Get_cas() == 0)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 1)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 2)
   {  
      return sin(x)+cos(y);
   }
   else if (this->_df->Get_cas() == 3)
   {  
      return 0.0;
   }
   else
   {
      return 0.0;
   }
}

double Function::Dirichlet_Gamma_1(const double x, const double y, const double t) const
{
   if (_df->Get_cas() == 0)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 1)
   {
      return 0.0;
   }
   else if (this->_df->Get_cas() == 2)
   {  
      return sin(x)+cos(y);
   }
   else if (this->_df->Get_cas() == 3)
   {  
      return 1.0;
   }
   else
   {
      return 0.0;
   }
}


#define _FUNCTION_CPP
#endif
