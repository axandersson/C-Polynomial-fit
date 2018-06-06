#include <vector>
#include <iostream>
#include "PolyFit.hpp"
using namespace std;
int main() {
    const int order = 3;
    vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    vector<double> y = {-7.5740, 5.2055, 11.3065, 9.4567, 7.0056, 1.2991, 20.3567, 52.3649, 105.1265, 196.0492};
    vector<double> polCoef;
    
    const int error = polyFit(polCoef,x,y,order);

    if(!error)
    {
        cout << "Pol coeficient (c_0 + c_1 * x + c_2 x^2 + ... + c_n * x^n" << endl;
        cout << "c = [  ";
        for(const auto& c : polCoef)
            cout << c << " ";
        cout << " ]" << endl;;
    }



    return 0;
}