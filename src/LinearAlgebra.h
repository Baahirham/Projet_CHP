#ifndef _LINEAR_ALGEBRA_H

#include <vector>

double DotProduct(const std::vector<double> a, const std::vector<double> b);
std::vector<double> MultiplyBy(const std::vector<double> a, const double b);
std::vector<double> AddVector(const std::vector<double> a, const std::vector<double> b);
std::vector<double> SubVector(const std::vector<double> a, const std::vector<double> b);

#define _LINEAR_ALGEBRA_H
#endif
