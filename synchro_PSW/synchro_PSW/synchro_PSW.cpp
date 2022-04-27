#include <iostream>
#include <math.h>
#include <vector>
#include <random>
#include <functional>
#include<fstream>

using namespace std;

//Константы
const double b = 4.15;    //Скользящие движения появлются примерно на 4,15
const double alpha = 2.0;
const double v = 0.65;
const double lambda = 0.294;
const double w = 2.0;
const double delta = 0.588;
const double k = 1.5;

template<typename T>
int sign(T var)
{
    if (var < 0)
    {
        return -1;
    }
    else { return 1; }
}

std::vector<double> lorenzModel_PWS(vector<double> point, vector<double> parameters) {
    //Обозначения:
    //point={x,y,z}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={N(0),b(1),alpha(2),beta(3),lambda(4),nu(5),omega(6)}

    //вектор производных {dx/dt, dy/dt, dz/dt}
    vector<double> dPoint(9);
    double sum = 0;
    for(int i = 0;i < 3; i++) {
        sum += point[3*i];
    }
    sum = k * sum;
    for (int i = 0; i < 3; i++) {
        if ((abs(point[3 * i]) < 1) && (point[3 * i + 2] < parameters[1]))
        {
            //система Аs
            dPoint[3 * i] = point[3*i] + sum - k * parameters[0]*point[3*i];
            dPoint[3 * i + 1] = -parameters[2] * point[3 * i + 1];
            dPoint[3 * i + 2] = -parameters[5] * point[3 * i + 2];
        }
        else
        {
            if (((point[3 * i] < -1) && (point[3 * i + 2] <= parameters[1]))
                || ((point[3 * i] < -sign(point[3 * i + 1])) && (point[3 * i + 2] > parameters[1])))
            {
                //система Al
                dPoint[3 * i] = -parameters[4] * (point[3 * i] + 1) + parameters[6] * (point[3 * i + 2] - parameters[1]) + sum - k * parameters[0] * point[3 * i];
                dPoint[3 * i + 1] = -parameters[3] * (point[3 * i + 1] + 1);
                dPoint[3 * i + 2] = -parameters[6] * (point[3 * i] + 1) - parameters[4] * (point[3 * i + 2] - parameters[1]);
            }
            else
            {
                //система Ar
                dPoint[3 * i] = -parameters[4] * (point[3 * i] - 1) - parameters[6] * (point[3 * i + 2] - parameters[1]) + sum - k * parameters[0] * point[3 * i];
                dPoint[3 * i + 1] = -parameters[3] * (point[3 * i + 1] - 1);
                dPoint[3 * i + 2] = parameters[6] * (point[3 * i] - 1) - parameters[4] * (point[3 * i + 2] - parameters[1]);
            }
        }
    }

    return dPoint;
}

std::vector<double> lorenzModel(vector<double> point, vector<double> parameters) {
    //Обозначения:
    //point={x,y,z}
    //dPoint={dx/dt, dy/dt, dz/dt}
    //parameters={N(0),sigma(1), r(2), beta(3)}

    //вектор производных {dx/dt, dy/dt, dz/dt}
    vector<double> dPoint(9);
    double sum = 0;
    for (int i = 0; i < 3; i++) {
        sum += point[3 * i];
    }
    sum = k * sum;
    for (int i = 0; i < 3; i++) {
        dPoint[3 * i] = parameters[1] * (point[3 * i + 1] - point[3 * i]) + sum - k * parameters[0] * point[3 * i]; //parameters[1] * (point[3 * i + 1] - point[3 * i])
        dPoint[3 * i + 1] = point[3 * i] * (parameters[2] - point[3 * i + 2]) - point[3 * i + 1]; //point[3 * i] * (parameters[2] - perem[3 * i + 2]) - perem[3 * i + 1]
        dPoint[3 * i + 2] = point[3 * i] * point[3 * i + 1] - parameters[3] * point[3 * i + 2]; // point[3 * i] * point[3 * i + 1] - parameters[3] * point[3 * i + 2]
    }

    return dPoint;
}

//интегратор Рунге-Кутты 4-го порядка с фиксированным шагом
std::vector<double> odeRK4(vector<double>(*fun)(vector<double>, vector<double>),
    vector<double> initCon, vector<double> parameters, double t_step)
{
    //Векторы "точки" длинной с вектор начальных условий
    vector<double> point(initCon.size());
    vector<double> buf_point(initCon.size());

    point = initCon;

    //Вектор производных
    vector<double> dpoint(initCon.size());

    //Двумерный массив коэффициентов интегрирования
    double** K = new double* [initCon.size()];
    for (int i = 0; i < point.size(); i++)
    {
        K[i] = new double[4];
    }

    /////////////////////////////////////////
    //K1
    for (int i = 0; i < point.size(); i++)
    {
        K[i][0] = t_step * fun(point, parameters)[i];
    }

    //приращение buf_point для K2
    for (int i = 0; i < point.size(); i++)
    {
        buf_point[i] = point[i] + 0.5 * K[i][0];
    }

    //K2
    for (int i = 0; i < point.size(); i++)
    {
        K[i][1] = t_step * fun(buf_point, parameters)[i];
    }

    //приращение buf_point для K3
    for (int i = 0; i < point.size(); i++)
    {
        buf_point[i] = point[i] + 0.5 * K[i][1];
    }

    //K3
    for (int i = 0; i < point.size(); i++)
    {
        K[i][2] = t_step * fun(buf_point, parameters)[i];
    }

    //приращение buf_point для K4
    for (int i = 0; i < point.size(); i++)
    {
        buf_point[i] = point[i] + K[i][2];
    }

    //K4
    for (int i = 0; i < point.size(); i++)
    {
        K[i][3] = t_step * fun(buf_point, parameters)[i];
    }

    //интегрирование
    for (int i = 0; i < point.size(); i++)
    {
        point[i] += (K[i][0] + 2 * K[i][1] + 2 * K[i][2] + K[i][3]) / 6;
    }
    /////////////////////////////////////////

    //очистка памяти
    for (int i = 0; i < point.size(); i++)
    {
        delete[] K[i];
    }
    delete[] K;

    return point;
}


int main()
{
    double dt = 0.001;
    int time = 100000;
    vector<double> point(9);
    vector<double> parameters(7);
    vector<double> param(4);

    ofstream f, f1, f2, f3, f4;
    //f.open("res.txt", ios::out);
    //f1.open("res_diff.txt", ios::out);
    f2.open("res_diff_orig.txt", ios::out);
    f3.open("res_orig.txt", ios::out);
    f4.open("res_gradient.txt", ios_base::out);
    f.open("point.txt", ios::out);

    //начальные значения
    point[0] = 0.1;
    point[1] = 0.1;
    point[2] = 0.1;

    point[3] = 0.2;
    point[4] = 0.5;
    point[5] = 0.8;

    point[6] = 0.1;
    point[7] = 1;
    point[8] = 1;

    parameters[0] = 3;  // размер сети N
    parameters[1] = b;
    parameters[2] = alpha;
    parameters[3] = delta;
    parameters[4] = lambda;
    parameters[5] = v;
    parameters[6] = w;

    param[0] = 3;
    param[1] = 10.0;
    param[2] = 28.0;
    param[3] = 8.0 / 3.0;

    int N = 10; // число экспериментов
    int n = 0; // число систем с синхронизацией


    /*for (int i = 0; i < time; i++)
    {
        point = odeRK4(lorenzModel, point, parameters, dt);
        //cout << point[0] << ", " << point[1] << ", " << point[2] << endl;
        f << i << " " << point[0] << " " << point[3] << " " << point[6] << " " << endl;
    }*/

    /*for (int i = 0; i < time; i++)
    {
        point = odeRK4(lorenzModel_PWS, point, parameters, dt);
        //cout << point[0] << ", " << point[1] << ", " << point[2] << endl;
        f1 << i << " " << point[0] - point[6] << " " << point[3] - point[6] << endl;
    }*/

   /* for (int i = 0; i < time; i++) {
       point = odeRK4(lorenzModel, point, param, dt);
       //cout << point[0] << ", " << point[1] << ", " << point[2] << endl;
       f2 << i << " " << point[0] - point[6] << " " << point[3] - point[6] << endl;
    }

    for (int i = 0; i < time; i++)
    {
        point = odeRK4(lorenzModel, point, param, dt);
        f3 << i << " " << point[0] << " " << point[3] << " " << point[6] << " " << endl;
    }*/

    for (double eps = 0; eps < 1; eps += 0.1) {
        for (int k = 0; k < 10; k++) {
            for (int j = 0; j < 10; j++) {
                point[0] = (double)(rand()) / RAND_MAX;
                point[1] = (double)(rand()) / RAND_MAX;
                point[2] = (double)(rand()) / RAND_MAX;
                point[3] = (double)(rand()) / RAND_MAX;
                point[4] = (double)(rand()) / RAND_MAX;
                point[5] = (double)(rand()) / RAND_MAX;
                point[6] = (double)(rand()) / RAND_MAX;
                point[7] = (double)(rand()) / RAND_MAX;
                point[8] = (double)(rand()) / RAND_MAX;
                if ((abs(point[0] - point[6]) < eps) & ((point[3] - point[6]) < eps)) {
                    n++;
                }
                f << eps << " " << k << " " << j << " " << n << " " << point[0] << " " << point[3] << " " << point[6] << " " << endl;
                //вопрос ,как определить есть ли синхронизация. Возможно нужна проверка 
                //принадлежит ли i-я координата х диапазону [-eps;eps]. 
            }
            if (n < 10) {
                double p = n / N;
                f4 << eps << " " << k << " " << n << endl;
            }
            else {
                n = 0;
            }
        }
    }
}
