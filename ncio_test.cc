// ncio_test.cc

#include <iostream>
#include <string>
#include <complex>
#include "ncio.h"



int main()
{
    const double Dim1 [4] = {1., 2., 3., 4.};
    const double Dim2_1 [6] = {10., 11., 12., 13., 14., 15.};
    const double Dim2_2 [6] = {20., 21., 22., 23., 24., 25.};
    
    const std::string GlobalAttr = "The global attributes";
    
    const size_t dims_num = 2;
    const std::string dim_names [dims_num] = {"Dim_1", "Dim_2"};
    const size_t dim_lengths [dims_num] = {4, 6};
    
    const size_t vars_num = 3;
    const std::string var_names [vars_num] = {"Var_1", "Var_2", "Var_3"};
    const bool var_complex [vars_num] = {false, false, true};
    
    const std::string path="";
    
    newDS_t newDS(dims_num, dim_names, dim_lengths, vars_num, var_names, var_complex, GlobalAttr, path);
    
    newDS.DefCoordVar(0, "Dim1");
    
    const int varid_Dim2_1 = newDS.DefCoordVar_custom(1, "Dim_2_1");
    const int varid_Dim2_2 = newDS.DefCoordVar_custom(1, "Dim_2_2");
    
    newDS.EndDef();
    
    newDS.WriteCoordVar(0, Dim1);
    newDS.WriteCoordVar_custom(varid_Dim2_1, Dim2_1);
    newDS.WriteCoordVar_custom(varid_Dim2_2, Dim2_2);
    
    const double Var1 [4*6] = {11, 12, 13, 14, 15, 16,
                               21, 22, 23, 24, 25, 26,
                               31, 32, 33, 34, 35, 36,
                               41, 42, 43, 44, 45, 46};
    const double Var2 [4*6] = {11, 12, 13, 14, 15, 16,
                               21, 22, 23, 24, 25, 26,
                               31, 32, 33, 34, 35, 36,
                               41, 42, 43, 44, 45, 46};
    std::complex<double> Var3 [4*6] = {{0.,0.},{0.,1.},{0.,2.},{0.,3.},{0.,4.},{0.,5.},
                                       {1.,0.},{1.,1.},{1.,2.},{1.,3.},{1.,4.},{1.,5.},
                                       {2.,0.},{2.,1.},{2.,2.},{2.,3.},{2.,4.},{2.,5.},
                                       {3.,0.},{3.,1.},{3.,2.},{3.,3.},{3.,4.},{3.,5.}};
    
    
    const double*const vars [vars_num] = {Var1, Var2, reinterpret_cast<const double*const>(Var3)};
    newDS.WriteVars(vars);
    
    const double energ [4*6] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
                                1.1, 1.2, 1.3, 1.4, 1.5, 1.6,
                                2.1, 2.2, 2.3, 2.4, 2.5, 2.6,
                                3.1, 3.2, 3.3, 3.4, 3.5, 3.6};
    
    newDS.WriteEnergy(energ);
    
    const double mu [4*6] = {0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
                                0.11, 0.12, 0.13, 0.14, 0.15, 0.16,
                                0.21, 0.22, 0.23, 0.24, 0.25, 0.26,
                                0.31, 0.32, 0.33, 0.34, 0.35, 0.36};
    
    newDS.Writemu(mu);
    
    const int loops [4*6] = {1, 2, 3, 4, 5, 6,
                             11, 12, 13, 14, 15, 16,
                             21, 22, 23, 24, 25, 26,
                             31, 32, 33, 34, 35, 36};
    
    newDS.WriteLoops(loops);
    
    return 0;
}
