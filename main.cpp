#include "functions.h"
#include "metropolis.h"
#include <random>
#include <iostream>
#include <future>
#include <thread>
#include <ctime>
#include <mutex>
#include <sstream>

const int Num_det = 16;
const int number_of_threads = 8;

std::vector <double> del_R;
std::vector <double> del_psi;
std::vector<std::thread> thread_pool;

std::mutex mtx;
std::ofstream out("pred.txt");
std::string s;

void loop3(int j, std::vector<double> &xx, std::vector<double> &yy, std::vector<int> &xx_t ,std::vector<double> &yy_t)
{
    while (j < xx.size()/Num_det){
    std::vector<double> x1(&xx[Num_det*j], &xx[Num_det*(j+1)]), y1(&yy[2*j], &yy[2*(j+1)]), yt1(&yy_t[2*j], &yy_t[2*(j+1)]);
    std::vector<int> xt1(&xx_t[Num_det*j], &xx_t[Num_det*(j+1)]);
    Metropolis M(x1, y1, xt1 ,yt1, s);
    M.arr_dir();
    M.start_init();
    M.find_min(10000);
    if(!(M.par[0] > 9.5 or M.par[0] < -9.5 or M.par[1] > 9.5 or M.par[1] < -9.5)){      // отбор внутри границ установки
    del_R.push_back(M.calc_del_r());
    if(M.calc_psi()!=-1) del_psi.push_back(M.calc_psi());
    mtx.lock();
    out << M.par[0] << ','<< M.par[1] << ','<< M.par[2] << ',' << M.par[3] << ',' << M.ne_calc() << '\n';
    mtx.unlock();
    }
    std::cout << j << ' ' << M.calc_del_r() << ' ' << M.calc_psi()  << '\n';
    j+=number_of_threads;
    }
}

int main(int argc, char **argv)
{
    if(argc != 3){
    std::cout << "Please choose config and data file" << '\n';
    return 0;
    }
    std::fstream f1(std::string(argv[1]).c_str());
    s = std::string(argv[2]);

    std::vector<double> xx, yy, yy_t;
    std::vector<int> xx_t;
    std::string line;
    std::vector<std::string> parsedRow2;
    int count{};
    while(std::getline(f1,line))
    {
        std::stringstream lineStream(line);
        std::string cell;
        std::vector<std::string> parsedRow;
        while(std::getline(lineStream,cell, ',') ) parsedRow.push_back(cell);
        for(int i = 0; i < 2*Num_det + 4; ++i)
        parsedRow2.push_back(parsedRow[i]);
        ++count;
    }
    std::cout << "Number of events: " << count << '\n';

    for(int i = 0; i < parsedRow2.size(); ++i){
        if(i%(2*Num_det + 4) ==0){
            for(int k = 0; k < 2*Num_det + 4; ++k){
                if(k < Num_det) xx.push_back(stod(parsedRow2[k+i]));
                else if(k < 2*Num_det and k >= Num_det)  xx_t.push_back(stoi(parsedRow2[k+i]));
                else if (k >= 2*Num_det and k < 2*Num_det+2 ) yy.push_back(stoi(parsedRow2[k+i]));
                else yy_t.push_back(stoi(parsedRow2[k+i]));
            }
        }
    }

    for(int i = 0; i < number_of_threads; ++i) thread_pool.emplace_back(loop3, i, std::ref(xx), std::ref(yy), std::ref(xx_t), std::ref(yy_t));
    for(int i = 0; i < number_of_threads; ++i){
        std::thread& entry = thread_pool[i];
        entry.join();
    }

    std::cout << "Events in border" << del_R.size() << '\n';
    std::cout << "Distance error: " << '\n';
    std::cout << "Mean "  << MX(del_R) << ' ' << " Mean Square Deviation " << std::sqrt(DX(del_R)) << '\n';
    std::cout << "0.72 quantile: " << quantile(del_R, 0.72) << '\n';

    print(del_R, "del_R.txt");

    std::cout << "psi error: " << '\n';
    std::cout << "Mean "  << MX(del_psi) << ' ' << " Mean Square Deviation " << std::sqrt(DX(del_psi)) << '\n';
    std::cout << "0.72 quantile: " << quantile(del_psi, 0.72) << '\n';

    std::vector<double> den;
    den = fun_dence(del_R);
    print(den, "dence.txt");
    return 0;
}
