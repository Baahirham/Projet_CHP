#ifndef _TIME_SCHEME_CPP

#include "TimeScheme.h"

TimeScheme::TimeScheme(DataFile* df, Laplacian* lap) :
_df(df), _lap(lap)
{

}

void TimeScheme::SaveSol(const std::vector<double> &U, std::string n_sol, int n){
    std::string n_file = "../res/" + n_sol + "." + std::to_string(n) + ".dat";
    std::ofstream monflux;
    monflux.open(n_file, std::ios::out);
    int k(0);
    for(int j = 1; j < _df->Get_Ny(); j++){
        for(int i = 1; i < _df->Get_Nx(); i++){
            monflux << _df->Get_xmin() + i*_df->Get_dx() << " " << _df->Get_ymin() + j*_df->Get_dy() << " " << U[k] << std::endl; 
            k++;
        }
    }
    monflux.close();
}

ImplicitScheme::ImplicitScheme(DataFile* df, Laplacian* lap) : 
TimeScheme(df, lap)
{

}

std::vector<double> ImplicitScheme::Jacobi(const std::vector<double> &U, const std::vector<double> &F){
    std::vector<double> x(U.size());
    //  A faire
    return x;
}

std::vector<double> ImplicitScheme::CG(const std::vector<double> &U, const std::vector<double> &F){
    std::vector<double> x(U.size());
    //  A faire
    return x;
}

std::vector<double> ImplicitScheme::BiCGstab(std::vector<double> &U, const std::vector<double> &F){
    std::vector<double> r(U.size()), r_old(U.size()), r_tilde(U.size()), p(U.size()), nu(U.size()), h(U.size()), s(U.size());
    std::vector<double> t1(U.size()), x(U.size()), p_old(U.size());
    double rho(0.0), rho_old(0.0), alpha(0.0), omega(0.0), beta(0.0);
    int Nmax(10000), k(0);
    double dt(_df->Get_dt());

    x = U;
    r = SubVector(AddVector(U,MultiplyBy(F,dt)),_lap->MatVecProd(x));
    r_tilde = r;
    rho = DotProduct(r_tilde,r);
    p = r;
    while (k <= Nmax){
        nu = _lap->MatVecProd(p);
        alpha = rho/DotProduct(r_tilde,nu);
        h = AddVector(x,MultiplyBy(p,alpha));
        s = SubVector(r,MultiplyBy(nu,alpha));
        if (sqrt(DotProduct(s,s)) <= 1.e-12){
            x = h;
            // std::cout << "it = " << k << std::endl;
            return x;
        }
        t1 = _lap->MatVecProd(s);
        omega = DotProduct(t1,s)/DotProduct(t1,t1);
        x = AddVector(h,MultiplyBy(s,omega));
        r_old = r;
        r = SubVector(s,MultiplyBy(t1,omega));
        if (sqrt(DotProduct(r,r)) <= 1.e-12) {
            // std::cout << "it = " << k << std::endl;
            return x;
        }
        rho_old = rho;
        rho = DotProduct(r_tilde,r);
        beta = (rho/rho_old)*(alpha/omega);
        p_old = p;
        p = AddVector(r,MultiplyBy(SubVector(p_old,MultiplyBy(nu,omega)),beta));
        k++;
    }

    // std::cout << "it = " << k << std::endl;
    if (k >= Nmax){
        std::cout << "Pas de convergence" << std::endl;
    }
    return x;
}

void ImplicitScheme::Integrate(double &t, std::vector<double> &U){
    t += _df->Get_dt();
    U = BiCGstab(U, _lap->RHS(t));
}

#define _TIME_SCHEME_CPP
#endif