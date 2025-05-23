#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <stdlib.h>
// Définition de la classe

class DataFile {

   private:
      const std::string _file_name;
      int _cas, _Nx, _Ny; 
      double _xmin, _xmax, _ymin, _ymax, _D, _Tf;
      double _dx, _dy, _dt;
      std::string _solver;


   public: // Méthodes et opérateurs de la classe
   DataFile(std::string file_name);

   const int Get_cas() const {return _cas;};
   const int Get_Nx() const {return _Nx;};
   const int Get_Ny() const {return _Ny;};
   const double Get_xmin() const {return _xmin;};
   const double Get_xmax() const {return _xmax;};
   const double Get_ymin() const {return _ymin;};
   const double Get_ymax() const {return _ymax;};
   const double Get_D() const {return _D;};
   const double Get_Tf() const {return _Tf;};
   const double Get_dx() const {return (_xmax - _xmin)/(double(_Nx));};
   const double Get_dy() const {return (_ymax - _ymin)/(double(_Ny));};
   const double Get_dt() const {return _dt;};
   const std::string Get_Solver() const {return _solver;};
   
};

#define _DATA_FILE_H
#endif
