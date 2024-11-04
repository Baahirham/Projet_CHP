#ifndef _LAPLACIAN_CPP

#include "Laplacian.h"

using namespace std;

// Constructeur
Laplacian::Laplacian(Function* fct, DataFile* df) :
_fct(fct), _df(df)
{
    
}

void Laplacian::InitialCondition(std::vector<double> &U){
    int k(0);
        for (int j = 1; j<_df->Get_Ny(); ++j){
            for (int i = 1; i<_df->Get_Nx(); ++i){
                U[k] = _fct->Initial_condition(_df->Get_xmin() + i*_df->Get_dx(),_df->Get_ymin() + j*_df->Get_dy());
                ++k;
            }
        }
}

std::vector<double> Laplacian::MatVecProd(const std::vector<double> &U){

    std::vector<double> X(U.size());
    int k(0);
    int Nx(_df->Get_Nx()), Ny(_df->Get_Ny());
    double dx(_df->Get_dx()), dy(_df->Get_dy()), xmin(_df->Get_xmin()), ymin(_df->Get_ymin()); 
    double dt(_df->Get_dt()), D(_df->Get_D());
    double a(1.0 + 2.0*dt*D*(1.0/(dx*dx) + 1.0/(dy*dy))), b(-D*dt/(dx*dx)), c(-D*dt/(dy*dy));
    int i, j;
    for(int k = 0; k < (Nx-1)*(Ny-1); ++k)
    {
        i = k%(Nx-1);
        j = k/(Nx-1);
        X[k] = a * U[k];
        if (i > 0) X[k] += b * U[k - 1];
        if (i < Nx-2) X[k] += b * U[k + 1];
        if (j > 0) X[k] += c * U[k - Nx+1];
        if (j < Ny-2) X[k] += c * U[k + Nx-1];
    }
    return X;
}

std::vector<double> Laplacian::RHS(const double t){
    const int Nx(_df->Get_Nx()), Ny(_df->Get_Ny());
    std::vector<double> rhs((Nx-1)*(Ny-1));
    int k(0);
    double dx(_df->Get_dx()), dy(_df->Get_dy());
    double x_i, xmin(_df->Get_xmin());
    double y_j, ymin(_df->Get_ymin());
    for(int j = 1; j < Ny; ++j){
        for(int i = 1; i < Nx; ++i){
            x_i = xmin + i*dx;
            y_j = ymin + j*dy;
            rhs[k] = _fct->Source(x_i,y_j,t); 
            if(i == 1){
                rhs[k] += _df->Get_D()*_fct->Dirichlet_Gamma_1(x_i - dx,y_j,t)/(dx*dx);
            }        
            if(i == Nx-1){
                rhs[k] += _df->Get_D()*_fct->Dirichlet_Gamma_1(x_i + dx,y_j,t)/(dx*dx);
            } 
            if(j == 1){
                rhs[k] += _df->Get_D()*_fct->Dirichlet_Gamma_0(x_i,y_j - dy,t)/(dy*dy);
            } 
            if(j == Ny-1){
                rhs[k] += _df->Get_D()*_fct->Dirichlet_Gamma_0(x_i,y_j + dy,t)/(dy*dy);
            } 
            ++k;
        }
    }
    return rhs;
}

std::vector<double> Laplacian::ExactSol(const double t){
    std::vector<double> ExactSol(_df->Get_Nx()*_df->Get_Ny());
        int k(0);
        for (int j = 1; j<_df->Get_Ny(); ++j){
            for (int i = 1; i<_df->Get_Nx(); ++i){
                ExactSol[k] = _fct->Exact_solution(_df->Get_xmin() + i*_df->Get_dx(),_df->Get_ymin() + j*_df->Get_dy(),t);
                ++k;
            }
        }
        return ExactSol;
}

#define _LAPLACIAN_CPP
#endif
