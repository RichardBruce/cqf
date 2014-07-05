#include <iostream>

#include "matrix.h"

using namespace std;


int main()
{
    cout << "Test Started" << endl;
    float data_a[16] = { 1, 3, 3, 4, 5, 6, 7, 8, 9, 10, 2, 12, 13, 14, 15, 16  };
    float data_b[9] = { 2, 1, 1, 2, 9, 10, 7, 8, 34 };
    float data_c[25] = { 1, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 4, 5, 6, 0, 0, 0, 7, 8, 9, 0, 0, 0, 0, 1 };
    float data_d[9] = { 2, 2, 4, 1, 3, -2, 3, 1, 3 };
    float data_e[100] = { 1,  2,  3, 4, 5, 6,  7, 8,  9, 10,
                          2,  1,  6, 7, 8, 9, 10, 9,  8,  7,
                          3,  6,  1, 6, 5, 4,  3, 5,  6,  7,
                          4,  7,  6, 1, 9, 8,  7, 3,  7,  5,
                          5,  8,  5, 9, 1, 4,  3, 8,  2,  6, 
                          6,  9,  4, 8, 4, 1,  2, 5,  6,  7,
                          7, 10,  3, 7, 3, 2,  1, 8, 11,  5,
                          8,  9,  5, 3, 8, 5,  8, 1,  5,  5,
                          9,  8,  6, 7, 2, 6, 11, 5,  1,  3,
                          10, 7,  7, 5, 6, 7,  5, 5,  3,  1 };
    float data_f[50] = { 1, 2, 3, 4, 5,
                         6, 7, 8, 9, 0,
                         5, 4, 3, 2, 1,
                         0, 9, 8, 7, 6,
                         3, 4, 5, 6, 7,
                         7, 6, 5, 4, 3,
                         8, 9, 0, 1, 2,
                         3, 2, 1, 0, 9, 
                         4, 7, 9, 4, 3,
                         4, 7, 3, 6, 5 };
    matrix<float> a(&data_a[0], 4, 4);
    matrix<float> b(&data_b[0], 3, 3);
    matrix<float> d(&data_c[0], 5, 5);
    matrix<float> e(&data_c[0], 5, 5);
    matrix<float> g(&data_e[0], 10, 10);
    matrix<float> h(&data_f[0], 5, 5);
    
//    a.dump(cout);
//    cout << endl;
//    b.dump(cout);
//    cout << endl;
    d.dump(cout);
    cout << endl;

    float sol[5];
    float eq[5] = { 1, 2, 3, 4, 5 };
    d.gauss_solve(&eq[0], &sol[0], Tridiagonal);
    d.dump(cout);
    cout << endl;
    for (int i = 0; i < 5; i++)
    {
        cout << sol[i] << ", ";
    }
    cout << endl;
//    matrix<float> soln(&sol[0], 5, 1);
//    matrix<float> res = e *  soln;
//    cout << endl << "Multiplied : " << endl;
//    res.dump(cout);
//    cout << endl;
//    matrix<float> c(a);
//    
//    cout << endl << "Householder Solver" << endl;
//    matrix<float> f(&data_d[0], 3, 3);
//    float hh_eq[3] = { 18, 1, 14 };
//    f.householder_solve(&hh_eq[0], &sol[0]);
//    f.dump(cout);
//    cout << endl;
//    for (int i = 0; i < 3; i++)
//    {
//        cout << sol[i] << ", ";
//    }
//    cout << endl;

    cout << endl << "SoR Solver" << endl;
    matrix<float> sor_m(&data_c[0], 5, 5);
    sor_m.dump(cout);
    cout << endl;
    
    float sor_eq[5] = { 1, 2, 3, 4, 5 };
    sor_m.sor_solve(&sor_eq[0], &sol[0], Tridiagonal);

    sor_m.dump(cout);
    cout << endl;
    for (int i = 0; i < 5; i++)
    {
        cout << sol[i] << ", ";
    }
    cout << endl;
    
//    cout << endl << "Householder tridiagonal" << endl;
//    float hh_e[10];
//    float hh_d[10];
//    matrix<float> hh_a(&data_e[0], 10, 10);
//    cout << "Before: " << endl;
//    hh_a.dump(cout);
//    cout << endl;
//    
//    hh_a.eigen_system(&hh_d[0], Symetric);
//    //cout << "Tridiagonal: " << endl;
//    //hh_a.dump(cout);
//    //cout << endl;
//    
//    cout << "a: " << endl;
//    hh_a.dump(cout);
//    cout << endl;
//    
//    cout << "d: " << endl;
//    for (int i = 0; i < 10; i++)
//    {
//        cout << hh_d[i] << ", ";
//    }
//    cout << endl;
//    
//    
//    cout << "Jacobi eigen vector" << endl;
//    cout << "Before:" << endl;
//    g.dump(cout);
////    cout << "Sum of upper squares: " << g.sum_of_upper_triangular_squares() << endl;
//    
//    matrix<float> ev = g;
//    matrix<float> *eigen_vecs = g.eigen_vectors(0.0001f, Symetric);
//    
//    cout << endl << "Eigen vectors: " << endl;
//    eigen_vecs->dump(cout);    
//    delete eigen_vecs;
//    
//    float *eigen_vals = ev.eigen_values(0.0001f, Symetric);
//    cout << endl << "Eigen values: " << endl;
//    for (int i = 0; i < 10; i++)
//    {
//        cout << eigen_vals[i] << ", ";
//    }
//    cout << endl;
//    delete [] eigen_vals;
//    
//    
//    cout << "Single value decomposition" << endl;
//    h.dump(cout);
//    float svd_w[5];
//    matrix<float> svd_v = h.single_value_decomposition(&svd_w[0]);
//    
//    cout << endl << "v matrix:" << endl;
//    svd_v.dump(cout);
//    
//    cout << endl << "u matrix" << endl;
//    h.dump(cout);
//    
//    cout << endl << "w matrix" << endl;
//    matrix<float> svd_wm = matrix<float>::diagonal_matrix(&svd_w[0], 5, 5);
//    svd_wm.dump(cout);
//    cout << endl;
//    
//    
//    cout << "after: " << endl;
//    h = h * (svd_wm * svd_v.transpose());
//    h.dump(cout);
//    cout << endl;
    
//    matrix<float>* identity = matrix<float>::identity_matrix(10);
//    identity->dump(cout);
//    delete identity;
//    
//    g.dump(cout);
//    cout << endl;
//    float diagonal[5];
//    float off_diag[5];
//    g.symetric_to_tridiagonal(&diagonal[0], &off_diag[0]);
//    g.dump(cout);
//    cout << endl;
//    for (int i = 0; i < 5; i++)
//    {
//        cout << diagonal[i] << ", ";
//    }
//    cout << endl;
//    for (int i = 0; i < 5; i++)
//    {
//        cout << off_diag[i] << ", ";
//    }
//    cout << endl;
    
    
//    a.transpose();
    //a.dump(cout);
    //cout << endl;
//    a.transpose();
    //a.dump(cout);
    //cout << endl;
    
//    cout << a.determinant() << endl << endl;
//    a.invert();
//    a.dump(cout);
//    cout << endl;
//    matrix<float> *lower;
//    b.lu_decomposition(&lower);
//    cout << "Upper: " << endl;
//    b.dump(cout);
//    cout << endl;
//    cout << "Lower: " << endl;
//    lower->dump(cout);
//    cout << endl;
//    cout << "Recomposition: " << endl;
//    matrix<float> recomp = (*lower) * b;
//    recomp.dump(cout);
//    cout << endl;
//    cout << "Cholesky Decomposition: " << endl;
//    b = recomp;
//    b.cholesky_decomposition();
//    b.dump(cout);
//    cout << endl;
//    cout << "Transpose: " << endl;
//    matrix<float> b_tran = b;
//    b_tran.transpose();
//    b_tran.dump(cout);
//    cout << endl;
//    cout << "Recomposition: " << endl;
//    b *= b_tran;
//    b.dump(cout);
//    cout << endl;
    
    //c.dump(cout);
    //cout << endl;
//    c = a + b;
    //c.dump(cout);
    //cout << endl;
//    c += c;
    //c.dump(cout);
    //cout << endl;
//    c -= b;
    //c.dump(cout);
    //cout << endl;
//    matrix<int> d = a * b;
    //d.dump(cout);
    //cout << endl;
//    a *= b;
    //a.dump(cout);
    //cout << endl;

    return 1;
}
