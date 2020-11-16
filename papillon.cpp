#include <iostream>
#include <fstream>
#include <math.h>

#include "rk4.cpp"

#define X   0
#define Y   1
#define Z   2

using namespace std;

double rho = 28.0;
double sigma = 10.0;
double beta = 8.0 / 3.0;

void deriv(int n, double t, double y[], double dy[])
{
	dy[X] = sigma * (y[Y] - y[X]);
	dy[Y] = y[X] * (rho - y[Z]) - y[Y];
	dy[Z] = y[X] * y[Y] - beta * y[Z];
}

int main()
{
	ofstream fd("papillon.res");

	int i, n = 3;
	double Tmax = 300, dt = 0.01, t = 0;
	double y[n];

	y[X] = 1.0;
	y[Y] = 1.0;
	y[Z] = 1.0;

	while (t < Tmax)
	{
		rk4(n, t, y, dt, deriv);
		t += dt;

		fd << t;
		for (i = 0; i < n; i++)
		{
			fd << " " << y[i];
		}
		fd << endl;
	}

	fd.close();
	return 0;
}
