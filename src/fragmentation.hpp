#ifndef FRAGMENTATION_HPP
#define FRAGMENTATION_HPP

#include <string>

double BCFY_DP(double z, double r);

double BCFY_DV(double z, double r);

double Dc_to_D0(double z, double r);

double D_kniehl_kramer(double z, double N, double eps);

// Configure the LHAPDF scale once from the main setup.
void SetLhapdfScale(double Q);

// LHAPDF-based fragmentation function (central member), D(z,Q) = xfxQ(pid, z, Q) / z
// The scale is taken from main.cpp.
double D_lhapdf(double z, const std::string& setname, int pid);

#endif