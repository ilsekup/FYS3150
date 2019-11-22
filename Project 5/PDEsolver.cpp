#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <armadillo>


using namespace std;
using namespace arma;


//explicit
vec explicitsch(int n, int t_steps)
{
double dx = 1.0 / (n+1);
double dt = 1.0 /(t_steps);
double alpha = dt / ( dx * dx) ;
if( alpha < 0.5)
  {
    cout << "Warning alpha is bigger than 0.5" << alpha << endl;
  }

vec u_xx = zeros<vec>(n+1);
vec u_t = zeros<vec>(n+1);

u_xx(0) = u_t(0) = 0.0;
u_xx(n) = u_t(n) = 1.0;

for (int j = 1; j < t_steps; j++)
  {
    for (int i = 1; i < n; i++)
    {
    u_t(i) = alpha*( u_xx(i-1) + u_xx(i+1) ) + (1 - 2*alpha) * u_xx(i);
    u_xx(i) = u_t(i);
    }
    cout << u_t << endl;
  }
return u_t;
}

int main(int argc, char* argv[])
{
int n, t;
cout << "n = " << endl;
cin >> n;
cout << "t = " << endl;
cin >> t;
vec u_t = explicitsch(n,t);
return 0;
}
