#include <gmp.h>
#include <gmpxx.h>
#include <time.h>
#include <fstream>
#include <string.h>
#include <iostream>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include "gaussian_integers.h"
using namespace std;
int main()
{
	//mpz p0 = 2, ps = 5, pt = 11, qs = 7, qt = 19; //GCD(g-i,N)=Norm(p);
	mpz p0 = 2, ps = 5, pt = 47, qs = 7, qt = 19;
	int d = 3;
	mpz p0d = pow(p0, d);
	mpz p = 2 * p0d * ps * pt + 1;
	mpz q = 2 * p0d * qs * qt + 1;
	mpz N = p * q;
	cout << "The binary representation of ps is " << decToBin(ps) << endl;
	cout << "The binary representation of qs is " << decToBin(qs) << endl;
	cout << "The binary representation of pt is " << decToBin(pt) << endl;
	cout << "The binary representation of qt is " << decToBin(qt) << endl;
	mpz phiN = Eulerphi(N);
	cout << "N= ";
	factor(N);
	cout << "p= ";
	factor(p);
	cout << "q= ";
	factor(q);
	cout << "N = ";
	printGIproduct(GIprimefactor({ N,0 }));
	cout << "Eulerphi(N) = " << phiN.get_str() << endl;
	mpz g = 1;
	for (; g < N; g++)
	{
		if (isorder(g, p0d, N))
			break;
	}
	cout << "g = " << g.get_str() << endl;
	mpz h = modularExponentiation(g, p0d / 4, N);
	cout << "h = " << h.get_str() << endl;
	GI gamma = { h, -1 };
	cout << "pho = GCD(" << N.get_str() << ", " << gamma << ") = ";
	GI beta = GIgcd(gamma, { N, 0 });
	cout << beta << " = ";
	printGIproduct(GIprimefactor(beta));
	for (int i = 200003; i < N; i++)
	{
		GI alpha{ i, 0 };
		if (JacobiSymbol(i, p) == 1 && JacobiSymbol(i, q) == 1)
		{
			cout << "x = " << i << endl;
			for (int j = 1; j < N; j++)
			{
				if (i == modularExponentiation(j, 2, N))
				{
					cout << "sqrt(x) = " << j << endl;
					break;
				}
			}
			for (int j = 1;; j++)
			{
				if (modularExponentiation(i, j, N) == 1)
				{
					cout << "The order of x in Z_N^* is " << j << endl;
					break;
				}
			}
			/*
			for (int j = 1; j < N; j++)
			{
				if (modularExponentiation(j, p0d * pt * qt, N) == i)
				{
					cout << "x is of the form y^{p_0^dp_tq_t}" << j << endl;
					break;
				}
			}
			*/
			QuarticResSymbol(alpha, beta, 1);
			break;
		}
	}
	return 0;
}
