#include "matrix.h"
#include <fstream>

void Matrix::fill_arr()
{
    for (int i = 0; i < m; ++i){
        for (int j = 0; j < n; ++j){
        //    std::cin >> A[i][j];
    }
  }
}

void Matrix::show_arr()
{
    for (int i = 0; i < m; ++i){
        //if(i!=0)std::cout << '\n';
        for (int j = 0; j < n; ++j){
            //std::cout << A[i][j] << " ";
        }
    }
}

void Matrix::mult_arr(Matrix a, Matrix b)
{
    for (int i = 0; i < a.m; i++){
        for (int j = 0; j < b.n; j++){
            A[i][j] = 0;
            for (int k = 0; k < a.n; k++)
                A[i][j] += a.A[i][k] * b.A[k][j];
        }
    }
}

void Matrix::add_arr(Matrix a, Matrix b)
{
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
            A[i][j] = a.A[i][j] + b.A[i][j];
        }
    }
}

void Matrix::transp(){
    double temp{};
    for (int i = 0; i < m; ++i){
        for (int j = 0; j < n; ++j){
            if(i!=j && i < j){
            temp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = temp;
            }
        }
    }
}

int Matrix::Get_m()
{
    return m;
};

int Matrix::Get_n()
{
    return n;
};

void Matrix::Del()
{
   for (int i = 0; i < m; i++)
   delete[] A[i];
   delete[] A;
}

double Matrix::shpur()
{
    double sh{};
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            if(i==j)sh+=A[i][j];
    return sh;
}

void Matrix::merge(Matrix a, Matrix b)
{
    Del();
    n = a.n + b.n;
    m = a.m;
    Create();
    for (int i = 0; i < m; i++){
        for (int j = 0; j < n; j++){
        if(j<a.n) A[i][j] = a.A[i][j];
        else A[i][j] = b.A[i][j-a.n];
        }
    }
}

double Matrix::Get_el(int i, int j)
{
    return A[i][j];
}

void Matrix::Set_el(int i, int j, double el)
{
    A[i][j] = el;
}

void Matrix::unit()
{
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if(i==j) A[i][j] = 1;
            else A[i][j] = 0;
        }
    }
}

void Matrix::Add_to_el(int i, int j, double el)
{
    A[i][j]+=el;
}

void Matrix::reverse()
{
    double **I;
    I = new double*[m];
    for (int i = 0; i < m; ++i)
    I[i] = new double[m];
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if(i==j) I[i][j] = 1;
            else I[i][j] = 0;
        }
    }

    double tem{};
    tem = A[0][0];
    for (int kk = 0; kk < m; kk++) {
        A[0][kk]/=tem;
        I[0][kk]/=tem;
    }

    int k{};
    double temp{};

    // Прямой ход

    for (int s = 1; s < m; s++) {
        if(A[s][s]!=0){
        for (int i = s; i < m; i++) {
            ++k;
            temp = A[i][s-1];
            for (int j = 0; j < m; j++) {
                A[i][j] = A[i][j] - A[i-k][j]*temp;
                I[i][j] = I[i][j] - I[i-k][j]*temp;
            }
        }
        if(A[s][s]!=0){
            tem = A[s][s];
            for (int kk = 0; kk < m; kk++) {
                A[s][kk]/=tem;
                I[s][kk]/=tem;
            }
        }
        k=0;
        }
    }

    // Обратный ход

    for (int s = m-2; s >= 0 ; s--) {
        for (int i = s; i >= 0; i--) {
            ++k;
            temp = A[i][s+1];
            for (int j = m-1; j >= 0; j--) {
                A[i][j] = A[i][j] - A[i+k][j]*temp;
                I[i][j] = I[i][j] - I[i+k][j]*temp;
            }
        }
        k=0;
      }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            A[i][j] = I[i][j];

        }
    }
}

void Matrix::LU(char x)
{
    double **L, **U;
    L = new double*[n];
    U = new double*[n];
    for (int i = 0; i < n; i++){
    L[i] = new double[n];
    U[i] = new double[n];
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            U[i][j] = 0;
            L[i][j] = 0;
        }
    }
    for (int j = 0; j < n; j++) {
        U[0][j] = A[0][j];
        if(U[0][0]!=0) L[j][0] = A[j][0]/U[0][0];
        //else qDebug() << "Division on zero" << '\n'; //std::cout << "Division on zero" << std::endl;
    }
    double pr_summ1{}, pr_summ2{};
    int i{};
    for (i = 1; i < n; i++) {
        for (int j = i; j < n; j++) {
        for (int k = 0; k <= i-1; k++) {
            pr_summ1+=L[i][k]*U[k][j];
            pr_summ2+=L[j][k]*U[k][i];
            }
        U[i][j] = A[i][j] - pr_summ1;
        if(U[i][i]!=0)L[j][i] = (A[j][i] - pr_summ2)/U[i][i];
        pr_summ1 = pr_summ2 = 0;
        //std::cout << U[i][j] << "   " << L[i][j] << std::endl;
        //std::cout << pr_summ1 << "   " << pr_summ2 << std::endl;
        }
    }
    if(x=='U')
    {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = U[i][j];
        }
    }
    }
    if(x=='L')
    {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = L[i][j];
        }
    }
    }
}

double Matrix::det()
{
    double determ = 1;
    this->LU('U');
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if(i==j)determ*=A[i][j];
        }
    }
    return determ;
}
