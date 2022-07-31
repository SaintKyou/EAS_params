#ifndef METROPOLIS_H
#define METROPOLIS_H
#include <random>
#include <vector>
#include <fstream>

class Metropolis        // гиперпараметры для поиска максимума правдоподобия - число итераций (N) и значение шага по параметрам (eps)
{
public:
    Metropolis(std::vector<double> &x1, std::vector<double> &y1,  std::vector<int> &xt1,std::vector<double> &yt) : x(x1), y(y1), t(xt1), y_t(yt)
    {
        std::ifstream conf("config_rect.txt"); //config_rect config
        for(int i = 0; i < Num_det; ++i) conf >> x_det[i] >> y_det[i];
    }
    ~Metropolis()
    {
    x.clear();
    y.clear();
    t.clear();
    y_t.clear();
    }

    std::vector<double> &x;
    std::vector<double> &y;     // TRUE values
    std::vector<int> &t;
    std::vector<double> &y_t;
    void start_init();
    void find_min(int N);
    void arr_dir();
    double calc_del_r();
    double calc_psi();
    double nkg(int &i);
    double ne_calc();
    double L();
    void print_params();
    std::vector<double> par;
    std::vector<double> angl;
private:
    double S = 0.36;
    int Num_det = 16;
    int threshold = 7;                      // порог сработавших кластеров
    double pi = 3.14159265358979323846;
    double c = 0.2998; // m/nsec
    std::random_device rd;
    double *x_det = new double[Num_det];
    double *y_det = new double[Num_det];
    std::vector<double> par_prom;
    std::vector<double> eps;
};

#endif // METROPOLIS_H
