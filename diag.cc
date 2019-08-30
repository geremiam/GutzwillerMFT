// diag.cc
/* Routines for finding eigenvalues and eigenvectors of matrices. */
 
#include <complex> // Complex numbers
/* Workaround to get LAPACKE to take std::complex types... Redefines the former as the 
latter. Has to come *before* including lapacke.h and *after* including complex. */
#define lapack_complex_double std::complex<double>
#define MKL_Complex16 std::complex<double>
//#include <lapacke.h> // Linear algebra routines (here: diagonalization)
#include <mkl.h>
#include <iostream> // Input from/output to command line
#include "diag.h" // Include for consistency check with header

void ErrorMessage(const int info);

int simple_zheev(const int MatOrder, std::complex<double>*const a, double*const w, 
                 const bool OutputEvecs, std::complex<double>*const z)
{
    /* https://software.intel.com/en-us/mkl-developer-reference-c-heev
    
    Finds the eigenvalues and optionally eigenvectors of a Hermitian matrix in full 
    storage layout. The first argument is the order (size) of the Hermitian matrix. The 
    second argument is the matrix in row-major full storage layout (i.e. a 
    MatOrder*MatOrder array), though only the lower-triangle of the matrix is referenced 
    (could have chosen the upper triangle). The evals are written to the array w (of 
    length MatOrder) in ascending order. If OutputEvecs is true, the orthonormal evecs 
    are written to z as column vectors (in row-major storage of course), in the order 
    corresponding to the evals in w; otherwise, z is not referenced. */
    
    /* We presume we are working with full arrays and not subarrays, so leading dimension 
    is the same as MatOrder. */
    const int lda=MatOrder;
    // Initialize info to unlikely value.
    int info=-99;
    // We use the lower triangle of the Hermitian matrix.
    const char uplo = 'L';
    
    // Choose the job type
    if (!OutputEvecs)
        info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'N', uplo, MatOrder, a, lda, w);
    else
    {
        info = LAPACKE_zheev(LAPACK_ROW_MAJOR, 'V', uplo, MatOrder, a, lda, w);
        for (int i=0; i<MatOrder*MatOrder; ++i)
            z[i] = a[i];
    }
    
    // Give error message if 'info' is nonzero.
    ErrorMessage(info);
    return info;
}

int simple_zheevr(const int MatOrder, std::complex<double>*const a, double*const w, 
                 const bool OutputEvecs, std::complex<double>*const z)
{
    /* https://software.intel.com/en-us/mkl-developer-reference-c-heevr 
    
    Still not sure what the output isuppz does... seems to give correct answers.
    
    Finds the eigenvalues and optionally eigenvectors of a Hermitian matrix in full 
    storage layout. The first argument is the order (size) of the Hermitian matrix. The 
    second argument is the matrix in row-major full storage layout (i.e. a 
    MatOrder*MatOrder array), though only the lower-triangle of the matrix is referenced 
    (could have chosen the upper triangle). The evals are written to the array w (of 
    length MatOrder) in ascending order. If OutputEvecs is true, the orthonormal evecs 
    are written to z as column vectors (in row-major storage of course), in the order 
    corresponding to the evals in w; otherwise, z is not referenced. */
    
    /* We presume we are working with full arrays and not subarrays, so leading dimension 
    is the same as MatOrder. */
    const int lda = MatOrder, ldz = MatOrder;
    // Initialize info to unlikely value.
    int info=-99;
    // We use the lower triangle of the Hermitian matrix.
    const char uplo = 'L';
    
    const double vl=0., vu=0.; // Not referenced
    const int il=0, iu=0; // Not referenced
    
    // 'abstol' sets the abs. tol. This choice gives the highest precision.
    const double abstol = LAPACKE_dlamch('S');
    
    
    // The number of evals returned is assigned to 'm' (should always be MatOrder here).
    int64_t m=-99;
    
    int64_t isuppz [2*MatOrder]; // Support of evecs in z
    
    // Choose the job type
    if (!OutputEvecs)
        info = LAPACKE_zheevr(LAPACK_ROW_MAJOR, 'N', 'A', uplo, MatOrder, a, lda, vl, vu, 
                              il, iu, abstol, &m, w, z, ldz, isuppz);
    else
        info = LAPACKE_zheevr(LAPACK_ROW_MAJOR, 'V', 'A', uplo, MatOrder, a, lda, vl, vu, 
                              il, iu, abstol, &m, w, z, ldz, isuppz);
    
    // Checks:
    if (m!=MatOrder) std::cout << "simple_zheevr: something went wrong. A" << std::endl;
    
    if (OutputEvecs)
      for (int i=0; i<MatOrder; ++i)
        if ( (isuppz[2*i]!=1) || (isuppz[2*i+1]!=MatOrder) )
          std::cout << "simple_zheevr: something went wrong. B" << std::endl;
    
    // Give error message if 'info' is nonzero.
    ErrorMessage(info);
    return info;
}

void ErrorMessage(const int info)
{
    if (info!=0)
        std::cout << "**** Warning -- LAPACKE error\
                      \nSee documentation for error code meaning.\
                      \ninfo = " << info << std::endl;
}
