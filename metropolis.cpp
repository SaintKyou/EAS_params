#include "metropolis.h"
#include "functions.h"
#include "matrix.h"
#include <iostream>

double Metropolis::nkg(int &i){
    double rm = 80.0;
    double delt = 0.1;
    double r = std::sqrt(std::pow(par[0]-x_det[i], 2) + std::pow(par[1]-y_det[i], 2));
    //if(r < 0.1) r = 0.1;
//    return (1/rm/rm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-par[2]))
//           *(std::pow((r/rm),(par[2]-2.)))*std::pow((1+r/rm),(par[2]-4.5));
    return (1/rm/rm)*std::tgamma(4.5-par[2])/(2*pi*tgamma(par[2])*tgamma(4.5-par[2]))
           *(std::pow(((r+delt)/rm),(par[2]-2.)))*std::pow((1+(r+delt)/rm),(par[2]-4.5));
}

double Metropolis::L(){
    double F{}, N_i{};
    for(int i = 0; i < Num_det; ++i){
        N_i = par[3]*nkg(i)*S;
        if(x[i]< 15) F += x[i]* std::log(N_i) - N_i;
        else F += std::log(1/std::sqrt(2*pi*N_i)) - std::pow((x[i]-N_i)/(2*std::sqrt(N_i)),2);
    }
    return F;
}

void Metropolis::start_init(){      // Some zero initialization
        std::mt19937 mt(rd());
        double Sum_x{}, Sum_y{}, Q_all{}, Ne_0{};
        for(int i = 0; i < Num_det; ++i){
            Sum_x += x[i]*x_det[i];
            Sum_y += x[i]*y_det[i];
            Q_all += x[i];              // x - n
        }
        par.push_back(Sum_x/Q_all);
        par.push_back(Sum_y/Q_all);
        double s = 1.0;
        par.push_back(s);
        Sum_x = 0;
        for(int i = 0; i < Num_det; ++i)
            Sum_x+=S*nkg(i);
        Ne_0 = Q_all/Sum_x;
        par.push_back(Ne_0);

        eps.push_back(10);
        eps.push_back(10);
        eps.push_back(0.1);
        eps.push_back(1000);
}

double Metropolis::ne_calc(){
    double Sum_x{},Q_all{};
    for(int i = 0; i < Num_det; ++i) Q_all += x[i];
    for(int i = 0; i < Num_det; ++i)
        Sum_x+=S*nkg(i);
    return Q_all/Sum_x;
}

double Metropolis::calc_del_r(){
    double loss{};
    loss = std::sqrt(std::pow(par[0] - y[0],2)+std::pow(par[1] - y[1],2));
    return loss;
}

void Metropolis::find_min(int N){
    int j{};
    std::mt19937 mt(rd());
    int num_par{};
    std::uniform_real_distribution<double> distribution{0, 1};
    std::uniform_int_distribution<> it(0, 3);
    double gamma{};
    double step{};
    par_prom = par;
    double loss{}, loss2{};
    while (j < N){
        loss = L();
        num_par = it(mt);
        gamma = distribution(mt);
        step = (-1+2*gamma)*eps[num_par];
        par[num_par] +=step;
        loss2 = L();
        if (loss2 <= loss) par[num_par] -=step;
//        if(std::abs(loss2-loss) < 1e-3){
//            for(int i = 0; i < 4; ++i) eps[i]*=0.1;
//        }
        //if(j%(N/10)==0) std::cout << j << ' ' << par[0] << ' ' << par[1] << ' ' << par[2] << ' ' << loss << '\n';
        //if(loss <= 1e-3) break;
        j+=1;
    }
}

void Metropolis::arr_dir(){
    int count{}, u{}, min{};
    double tetha{}, phi{};
    double alpha{}, betta{}, C{};
    int* tim = &t[0];
    Matrix M(3,3), Omega(3,1), X(3,1);

    // Нормировочка
    for (int i = 0; i < Num_det; ++i) {
        if(tim[i]!=-1) {
            min = tim[i];
            break;
        }
    }
    for (int i = 0; i < Num_det; ++i) {
        if(tim[i] < min && tim[i]!=-1) min = tim[i];
    }
    for (int i = 0; i < Num_det; ++i) {
        if(tim[i]!=-1) tim[i]-=min;
    }

    //Устранение выбросов
    for (int j = 0; j < Num_det; ++j ) {
        if(tim[j] > 1000) tim[j] = -1;
        if(tim[j]>=0) count+=1;
    }

    if(count > threshold){
        for(int j = 0; j < Num_det; ++j){
            if(tim[j]!=-1){
                ++u;
                M.Add_to_el(0,0,x_det[j]*x_det[j]);
                M.Add_to_el(0,1,x_det[j]*y_det[j]);
                M.Add_to_el(0,2,x_det[j]);
                M.Add_to_el(1,1,y_det[j]*y_det[j]);
                M.Add_to_el(1,2,y_det[j]);
                M.Add_to_el(1,0,x_det[j]*y_det[j]);
                M.Add_to_el(2,0,x_det[j]);
                M.Add_to_el(2,1,y_det[j]);
                M.Set_el(2,2,u);
                Omega.Add_to_el(0,0,c*tim[j]*x_det[j]);
                Omega.Add_to_el(1,0,c*tim[j]*y_det[j]);
                Omega.Add_to_el(2,0,c*tim[j]);
            }
        }
        M.reverse();
        X.mult_arr(M,Omega);

        alpha = X.Get_el(0,0);
        betta = X.Get_el(1,0);
        C = sqrt(1-alpha*alpha-betta*betta);

        tetha = acos(C)*180.0/M_PI;
        phi = atan2(betta, alpha)*180.0/M_PI;
        if(phi < 0) phi+=360.0;

        if(tetha>7) angl.insert(angl.end(), {tetha, phi});//out << tetha << "," << phi << "," << count << std::endl;
        else angl.insert(angl.end(), {-1, -1});//out << tetha << "," << -1 << "," << count << std::endl;
        //std::cout << tetha << "," << phi << "," << count  << std::endl;
        count = u  = 0;

        for (int i = 0; i < 3; ++i) {
            Omega.Set_el(i,0,0);
            for (int j = 0; j < 3; ++j) {
                M.Set_el(i,j,0);
            }
        }
    }
    else{
        angl.insert(angl.end(), {-1, -1}); //out << -1 << "," << -1 << "," << count << std::endl;
        count = 0;
    }
}

void Metropolis::print_params(){
    if(par.size()==4)  std::cout << "EAS Parameters: " << par[0] << ' ' << par[1] << ' '<< par[2] << ' '<< par[3] << '\n';
    if(angl.size()==2) std::cout << "Arrival Direction: " << angl[0] << ' ' << angl[1] << ' ' << "Real: " << y_t[0] << ' ' << y_t[1] << '\n';
}

double Metropolis::calc_psi(){
    if(angl[0]!=-1 and angl[1]!=-1){
    double t_v = angl[0]*pi/180., p_v = angl[1]*pi/180., t = y_t[0]*pi/180., p = y_t[1]*pi/180.;
    return std::acos(std::sin(t_v)*std::cos(p_v)*std::sin(t)*std::cos(p)+
    std::sin(t_v)*std::sin(p_v)*std::sin(t)*std::sin(p)+
    std::cos(t_v)*std::cos(t))*180./pi;
    }
    else return -1;
}
