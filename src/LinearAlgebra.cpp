#ifndef _LINEAR_ALGEBRA_CPP

#include "LinearAlgebra.h"

double DotProduct(const std::vector<double> a, const std::vector<double> b){
    double result(0.0);
    for(int i = 0; i < a.size(); i++){
        result += a[i]*b[i];
    }
    return result;
}

std::vector<double> MultiplyBy(const std::vector<double> a, const double b){
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = b*a[i];
    }
    return result;
}

std::vector<double> AddVector(const std::vector<double> a, const std::vector<double> b){
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] + b[i];
    }
    return result;
}

std::vector<double> SubVector(const std::vector<double> a, const std::vector<double> b){
    std::vector<double> result(a.size());
    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] - b[i];
    }
    return result;
}

#define _LINEAR_ALGEBRA_CPP
#endif