#ifndef intracellularFBA_H
#define intracellularFBA_H
#include <vector>

void buildModelFromMat(const char *file_path);
void changeBound(int col, double lb, double ub);
void optimizeModel();
std::vector<double> getFluxes();
#endif
