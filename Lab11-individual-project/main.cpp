#include <cstdlib>
#include <iostream>
#include <cmath>
#include <stdio.h>
#include <cstdio>
#include <fstream>
#include <iomanip>

#include "math.h"
#include "iostream"
//#include "calerf.h"

using namespace std;
double** get_condition_matrix(int n, int m, double dt);
double** get_analitical_values(int n, int m, double dt, double h);
void print_matrix(double** A, int n, int m);
double** alloc_0_matrix(int n, int m);
double* alloc_0_vector(int n);
void dealloc(double **M, int n);

void laasonen_lu(double **U, int t, int x, double h);
void laasonen_thomas(double **U, int t, int x, double h);
void thomas_algorithm(double *l, double *d, double *u, double *b, double *x, int m);

void classical(double **A, int n, int m, double h);


// LU
void gauss(double **U, double **L, double** T, int N);
int partial_pivot(double** M, double** L, int row, double** T, int N);
double* gauss_resoult(double **U, double **L, double **P, double *b, int N);


namespace calerf
{
long double CALERFL(const long double arg, const int jint);
long double ERFL(const long double x);
long double ERFCL(const long double x);
long double EREXL(const long double x);
};

//deklaracja stalych i zmiennych projektowych
static double lambda_indirect = 1.0, D = 1.0;
static double lambda_direct = 0.4;
static int r = 1, a = 10, tmax = 2;
static int GAUSS = 0, THOMAS = 1;

// UK£AD MACIERZY U
//		x0 x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11
//	t0
//	t1
//	t2
//  t3
//  t4
//  t5

int main()
{
	double dt;
	int x_max;
	int t_max;
	double error;
	double error_lu;
	double error_thomas;
	double temp;

	double** U_analitical;
	double** U_laasonen_lu;
	double** U_laasonen_thomas;
	double** U_direct;


	/////////////////////////////////////////////////////////////////////////// (1) ////////////////////////////////////////////////////////////////////////////////////

	
	ofstream o1("1_laasonen_gnu.txt");
	ofstream o11("1_laasonen.txt");
	o11 << setprecision(16) <<"B³êdy dla Lassonen\n" << setw(24) << "log10(h)" << setw(24) << "log10(blad_lu)" << setw(24) << "log10(blad_thomas)\n\n";
	o1 << setprecision(16) << endl;
	for(double h=0.05; h<= 0.50; h+=0.01){
		error_lu=0;
		error_thomas=0;
		dt = (lambda_indirect * h*h) / D;
		x_max = (int)(a / h) + 1;
		t_max = (int)(tmax / dt) + 1;
	
		U_analitical = get_analitical_values(t_max, x_max, dt, h);
		U_laasonen_lu = get_condition_matrix(t_max, x_max, dt);
		U_laasonen_thomas = get_condition_matrix(t_max, x_max, dt);
		
		laasonen_lu(U_laasonen_lu, t_max, x_max, h);
		laasonen_thomas(U_laasonen_thomas, t_max, x_max, h);
		
		for(int i=0; i<x_max; i++){
			temp = fabs(U_analitical[t_max-1][i] - U_laasonen_lu[t_max-1][i]);
			if (temp > error_lu)
				error_lu = temp;
			temp = fabs(U_analitical[t_max-1][i] - U_laasonen_thomas[t_max-1][i]);
			if (temp > error_thomas)
				error_thomas = temp;
		}
		dealloc(U_analitical, t_max);
		dealloc(U_laasonen_lu, t_max);
		dealloc(U_laasonen_thomas, t_max);
		cout << h << endl;
		o1 << log10(h) << " " << log10(error_lu) << " " << log10(error_thomas)<< endl;
		o11 << setw(24) << log10(h) << setw(24) << log10(error_lu) << setw(24) << log10(error_thomas) << endl;	
	}
	o1.close();
	o11.close();
	

	ofstream o2("1_direct_gnu.txt");
	ofstream o22("1_direct.txt");
	o22 << setprecision(16) <<"B³êdy dla Metody klasycznej bezpoœredniej\n" << setw(24) << "log10(h)" << setw(24) << "log10(|blad|)\n\n";
	o2 << setprecision(16) << endl;
	for(double h=0.05; h<= 0.51; h+=0.01){
		error=0;
		dt = (lambda_direct * h*h) / D;
		x_max = (int)(a / h) + 1;	
		t_max = (int)(tmax / dt) + 1;	
	
		U_analitical = get_analitical_values(t_max, x_max, dt, h);
		U_direct = get_condition_matrix(t_max, x_max, dt);	
		classical(U_direct, t_max, x_max, h);

		
		for(int i=0; i<x_max; i++){
			temp = fabs(U_analitical[t_max-1][i] - U_direct[t_max-1][i]);
			if (temp > error)
				error = temp;
		}
		dealloc(U_analitical, t_max);
		dealloc(U_direct, t_max);
		o2 << log10(h) << " " << log10(error) << endl;
		o22 << setw(24) << log10(h) << setw(24) << log10(error) << endl;	
	}
	o2.close();
	o22.close();
	
	
		/////////////////////////////////////////////////////////////////////////// (2) ////////////////////////////////////////////////////////////////////////////////////
	
	
	
	double h=0.04;
	ofstream o3("2_laasonen_lu_gnu.txt");
	o3 << setprecision(16);
	ofstream o33("2_laasonen_lu.txt");
	o33 << setprecision(16) <<"Rozwi¹zania anlityczne i numeryczne dla ró¿nych te[0,2] \n" << setw(20) << "x" << setw(20) << "(anal)t=0.2s" << setw(20) << "(anal)t=0.4s";
	o33 << setw(20) << "(anal)t=1s" << setw(20) << "(anal)t=1.6s" << setw(20) << "(anal)t=2s" << setw(20) << "(numer)t=0.2s" << setw(20) << "(numer)t=0.4s";
	o33 << setw(20) << "(numer)t=1s" << setw(20) << "(numer)t=1.6s" << setw(20) << "(numer)t=2s\n\n"; 
	ofstream o4("2_laasonen_thomas_gnu.txt");
	o4 << setprecision(16);
	ofstream o44("2_laasonen_thomas.txt");
	o44 << setprecision(16) <<"Rozwi¹zania anlityczne i numeryczne dla ró¿nych te[0,2] \n" << setw(20) << "x" << setw(20) << "(anal)t=0.2s" << setw(20) << "(anal)t=0.4s";
	o44 << setw(20) << "(anal)t=1s" << setw(20) << "(anal)t=1.6s" << setw(20) << "(anal)t=2s" << setw(20) << "(numer)t=0.2s" << setw(20) << "(numer)t=0.4s";
	o44 << setw(20) << "(numer)t=1s" << setw(20) << "(numer)t=1.6s" << setw(20) << "(numer)t=2s\n\n";
	dt = (lambda_indirect * h*h) / D;
	x_max = (int)(a / h) + 1;
	t_max = (int)(tmax / dt) + 1;

	U_analitical = get_analitical_values(t_max, x_max, dt, h);
	U_laasonen_lu = get_condition_matrix(t_max, x_max, dt);
	U_laasonen_thomas = get_condition_matrix(t_max, x_max, dt);	
	laasonen_lu(U_laasonen_lu, t_max, x_max, h);
	laasonen_thomas(U_laasonen_thomas, t_max, x_max, h);

	double x=r;
	for(int i=0; i<x_max;i++){
//		dt =0.0004,   [125] = 0.2s   [250] = 0.4s   [625] = 1s   [1000] = 1.6s    [t_max-1] = 2s
		o3 << x << " " << U_analitical[125][i] << " " << U_analitical[250][i] << " " << U_analitical[625][i] << " " << U_analitical[1000][i] << " " << U_analitical[t_max-1][i] << " ";
		o3 << U_laasonen_lu[125][i] << " " << U_laasonen_lu[250][i] << " " << U_laasonen_lu[625][i] << " " << U_laasonen_lu[1000][i] << " " << U_laasonen_lu[t_max-1][i] << endl;
		o4 << x << " " << U_analitical[125][i] << " " << U_analitical[250][i] << " " << U_analitical[625][i] << " " << U_analitical[1000][i] << " " << U_analitical[t_max-1][i] << " ";
		o4 << U_laasonen_thomas[125][i] << " " << U_laasonen_thomas[250][i] << " " << U_laasonen_thomas[625][i] << " " << U_laasonen_thomas[1000][i] << " " << U_laasonen_thomas[t_max-1][i] << endl;
		
		o33 << setw(20) << x << setw(20) << U_analitical[125][i] << setw(20) << U_analitical[250][i] << setw(20) << U_analitical[625][i] << setw(20) << U_analitical[1000][i] << setw(20) <<  U_analitical[t_max-1][i];
		o33 << setw(20) << U_laasonen_lu[125][i] << setw(20) << U_laasonen_lu[250][i] << setw(20) << U_laasonen_lu[625][i] << setw(20) << U_laasonen_lu[1000][i] << setw(20) << U_laasonen_lu[t_max-1][i] << endl;
		o44 << setw(20) << x << setw(20) << U_analitical[125][i] << setw(20) << U_analitical[250][i] << setw(20) << U_analitical[625][i] << setw(20) << U_analitical[1000][i] << setw(20) << U_analitical[t_max-1][i];
		o44 << setw(20) << U_laasonen_thomas[125][i] << setw(20) << U_laasonen_thomas[250][i] << setw(20) << U_laasonen_thomas[625][i] << setw(20) << U_laasonen_thomas[1000][i] << setw(20) << U_laasonen_thomas[t_max-1][i] << endl;
		x+=h;
	}
	o3.close();
	o33.close();
	o4.close();
	o44.close();
	
	
	
	
	
	double h=0.04;
	ofstream o5("2_direct_gnu.txt");
	o5 << setprecision(16);
	ofstream o55("2_direct.txt");
	o55 << setprecision(16) <<"Rozwi¹zania anlityczne i numeryczne dla ró¿nych te[0,2] \n" << setw(20) << "x" << setw(20) << "(anal)t=0.2s" << setw(20) << "(anal)t=0.4s";
	o55 << setw(20) << "(anal)t=1s" << setw(20) << "(anal)t=1.6s" << setw(20) << "(anal)t=2s" << setw(20) << "(numer)t=0.2s" << setw(20) << "(numer)t=0.4s";
	o55 << setw(20) << "(numer)t=1s" << setw(20) << "(numer)t=1.6s" << setw(20) << "(numer)t=2s\n\n"; 
	dt = (lambda_direct * h*h) / D;
	x_max = (int)(a / h) + 1;	
	t_max = (int)(tmax / dt) + 1;	
	cout << t_max << endl;
	cout << dt << endl;
	
	U_analitical = get_analitical_values(t_max, x_max, dt, h);
	U_direct = get_condition_matrix(t_max, x_max, dt);	
	classical(U_direct, t_max, x_max, h);

	
	double x=r;
	for(int i=0; i<x_max;i++){
//		dt =0.0004,   [313] = 0.2s   [625] = 0.4s   [1562] = 1s   [2500] = 1.6s    [t_max-1] = 2s
		o5 << x << " " << U_analitical[313][i] << " " << U_analitical[625][i] << " " << U_analitical[1562][i] << " " << U_analitical[2500][i] << " " << U_analitical[t_max-1][i] << " ";
		o5 << U_direct[313][i] << " " << U_direct[625][i] << " " << U_direct[1562][i] << " " << U_direct[2500][i] << " " << U_direct[t_max-1][i] << endl;
		
		o55 << setw(20) << x << setw(20) << U_analitical[313][i] << setw(20) << U_analitical[625][i] << setw(20) << U_analitical[1562][i] << setw(20) << U_analitical[2500][i] << setw(20) << U_analitical[t_max-1][i];
		o55 << setw(20) << U_direct[313][i] << setw(20) << U_direct[625][i] << setw(20) << U_direct[1562][i] << setw(20) << U_direct[2500][i] << setw(20) << U_direct[t_max-1][i] << endl;
		x+=h;
	}
	o5.close();
	o55.close();
	
	
	/////////////////////////////////////////////////////////////////////////// (3) ////////////////////////////////////////////////////////////////////////////////////

	
	double h=0.04;
	dt = (lambda_indirect * h*h) / D;
	x_max = (int)(a / h) + 1;
	t_max = (int)(tmax / dt) + 1;
	U_analitical = get_analitical_values(t_max, x_max, dt, h);
	U_laasonen_lu = get_condition_matrix(t_max, x_max, dt);
	U_laasonen_thomas = get_condition_matrix(t_max, x_max, dt);	
	laasonen_lu(U_laasonen_lu, t_max, x_max, h);
	laasonen_thomas(U_laasonen_thomas, t_max, x_max, h);
		
	
	ofstream o6("3_laasonen_gnu.txt");
	ofstream o66("3_laasonen.txt");
	o66 << setprecision(16) <<"B³êdy dla Lassonen\n" << setw(24) << "t" << setw(24) << "blad_lu" << setw(24) << "blad_thomas\n\n";
	o6 << setprecision(16) << endl;
	
	for(int t=5; t <t_max; t+=5){
		error_lu=0;
		error_thomas=0;
		
		for(int i=0; i<x_max; i++){
			temp = fabs(U_analitical[t][i] - U_laasonen_lu[t][i]);
			if (temp > error_lu)
				error_lu = temp;
			temp = fabs(U_analitical[t][i] - U_laasonen_thomas[t][i]);
			if (temp > error_thomas)
				error_thomas = temp;
		}
		o6 << t*dt << " " << error_lu << " " << error_thomas << endl;
		o66 << setw(24) << t*dt << setw(24) << error_lu << setw(24) << error_thomas << endl;	
	
	o6.close();
	o66.close();
	*/
		
	
	double h = 0.04;
	dt = (lambda_direct * h*h) / D;
	x_max = (int)(a / h) + 1;	
	t_max = (int)(tmax / dt) + 1;	
	
	U_analitical = get_analitical_values(t_max, x_max, dt, h);
	U_direct = get_condition_matrix(t_max, x_max, dt);	
	classical(U_direct, t_max, x_max, h);
		
	ofstream o7("3_direct_gnu.txt");
	ofstream o77("3_direct.txt");
	o77 << setprecision(16) <<"B³êdy dla Lassonen\n" << setw(24) << "t" << setw(24) << "blad\n\n";
	o7 << setprecision(16) << endl;
	for(int t=5; t <t_max; t+=1){
		error=0;
		
		for(int i=0; i<x_max; i++){
			temp = fabs(U_analitical[t][i] - U_direct[t][i]);
			if (temp > error)
				error = temp;
		}

		o7 << t*dt << " " << error << endl;
		o77 << setw(24) << t*dt << setw(24) << error << endl;	
	}
	o7.close();
	o77.close();
	
	
}

void dealloc(double **M, int n){
	for(int i; i<n; i++)
		delete M[i];
	delete M;
}

void classical(double **U, int n, int m, double h) {
	double x;
	for (int k = 0; k < n - 1; k++) {	// t
		x = r + h;
		for (int i = 1; i < m - 1 ; i++) { // x
			U[k + 1][i] = lambda_direct * (U[k][i - 1] - 2.0 * U[k][i] + U[k][i + 1]) + (2.0 * lambda_direct * h / x) * (U[k][i + 1] - U[k][i]) + U[k][i];
			//U[k+1][i]= lambda_direct*U[k][i-1] + (-2.0*lambda_direct*(1.0+h/x)+1.0)*U[k][i] + lambda_direct*(1+2*h/x)*U[k][i+1]; przekszta³cone - taki sam wynik
			x += h;
		}

	}
}

void laasonen_lu(double **U, int t, int x, double h) {

	// tworze macierz i wektory dla równania A*y=b
	double **A = alloc_0_matrix(x, x);
	double *b = alloc_0_vector(x);
	double *y = alloc_0_vector(x);
	// dodatkowe macierze dla metody LU
	double **L = alloc_0_matrix(x, x);
	double **P = alloc_0_matrix(x, x);
	
	// ZMIENNA x zale¿na od iteracji
	double xx;

	for (int k = 1; k < t; k++) {
		// zaczynamy od k=1 wiêc wartoœæ pocz¹tkowa r + 1*h
		xx = r + h;
		A[0][0] = 1.0;
		A[x - 1][x - 1] = 1.0;
		for (int i = 1; i < x - 1; i++) {
			A[i][i - 1] = lambda_indirect * (1.0 - h/xx);
			A[i][i] = -(2.0*lambda_indirect + 1.0);
			A[i][i + 1] = lambda_indirect * (1.0 + h / xx);
			b[i] = -U[k-1][i];
			xx += h;
		}
		b[0] = U[k][0];
		b[x - 1] = U[k][x-1];
	
		// rozwi¹zuje uk³ad przy pomocy dekompozycji LU
		gauss(A, L, P, x);
		y = gauss_resoult(A, L, P, b, x);
			
		// przypisanie wyników
		for (int i = 1; i < x - 1; i++)
			U[k][i] = y[i];
	}
	for (int i = 0; i < x; i++)
		delete A[i];
	delete A;
	delete b;
	delete y;
}

void laasonen_thomas(double **U, int t, int x, double h) {

	double *l = alloc_0_vector(x);
	double *d = alloc_0_vector(x);
	double *u = alloc_0_vector(x);
	double *b = alloc_0_vector(x);
	double *y = alloc_0_vector(x);

	double xx;
	for (int k = 1; k < t; k++) {
		d[0] = 1.0;
		d[x - 1] = 1.0;

		xx = r + h;
		for (int i = 1; i < x - 1; i++) {
			l[i] = lambda_indirect * (1.0 - h / xx);
			d[i]= -(2.0*lambda_indirect + 1.0);
			u[i] = lambda_indirect * (1.0 + h / xx);
			b[i] = -U[k - 1][i];
			xx += h;
		}
		b[0] = U[k][0];
		b[x - 1] = U[k][x - 1];

		thomas_algorithm(l, d, u, b, y, x);
		for (int i = 1; i < x - 1; i++)
			U[k][i] = y[i];
	}
}

void thomas_algorithm(double *l, double *d, double *u, double *b, double *x, int N) {
	for (int i = 1; i < N; i++)
		d[i] -= l[i] / d[i - 1] * u[i - 1];
	
	for (int i = 1; i < N; i++)
		b[i] -= l[i] / d[i - 1] * b[i - 1];

	x[N - 1] = b[N - 1] / d[N - 1];
	for (int i = N - 2; i >= 0; i--)
		x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
}

void gauss(double **U, double **L, double** P, int N) {

	// inicjalizuje Macierz przekszta³ceñ
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) P[i][j] = 1.0;
			else P[i][j] = 0.0;
		}
	}

	// Wyznaczam macierz LU za pomoc¹ eliminacji Gaussa
	// obliczam temp ( gdy¿ licznik zerujemy w pierwszym przejœciu pêtli i w kolejnych kolumnach program nic by nie zmieni³)
	for (int k = 0; k < N; k++) {
		// if el. podstawowy == 0 --> pivoting (czêœciowy wybór elementu podstawowego).
		if (U[k][k] == 0.0 && partial_pivot(U, L, k, P, N) == 0) {
			// nie znaleziono elementu podstawowego. b³¹d.
			return;
		}

		// zerowanie kolejnej kolumny U i uzupe³nianie L
		for (int i = 1 + k; i < N; i++) {
			double temp = U[i][k] / U[k][k];
			for (int j = k; j < N; j++) {
				U[i][j] -= U[k][j] * temp;
				L[i][j] = temp;
			}
		}
	}

	// z Macierzy L tworze macierz trojk¹tn¹ doln¹ (obliczone wsp. poni¿ej przek¹tnej (diagonali) pozostaja bez zmian)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			L[i][i] = 1;
			if (j > i)
				L[i][j] = 0.0;
		}
	}
}

int partial_pivot(double** U, double **L, int x, double** P, int N) {

	// Szukam max elementu w kolumnie pod M[x][x]
	int max_index = x;
	for (int i = x + 1; i < N; i++) {
		if (fabs(U[i][x]) > U[x][x])
			max_index = i;
	}

	if (max_index == x) {
		// nie w³aœciwa macierz 
		return 0;
	}

	else {
		double temp;
		// zamieniam wiersze w macierzy U i L i T
		for (int i = 0; i < N; i++) {
			temp = U[x][i];
			U[x][i] = U[max_index][i];
			U[max_index][i] = temp;

			temp = L[x][i];
			L[x][i] = L[max_index][i];
			L[max_index][i] = temp;

			temp = P[x][i];
			P[x][i] = P[max_index][i];
			P[max_index][i] = temp;
		}
		return 1;
	}
}

double* gauss_resoult(double **U, double **L, double **P, double *b, int N) {
	double *y = new double[N];
	double *x = new double[N];

	// przekszta³cam macierz wyników b,  zgodnie z macierz¹ przekszta³ceñ P
	double *temp = new double[N];
	for (int i = 0; i < N; i++) {
		temp[i] = 0.0;
		for (int j = 0; j < N; j++) {
			temp[i] += P[i][j] * b[j];
		}
	}
	b = temp;

	// równanie Ly = b 
	y[0] = b[0] / L[0][0];
	for (int i = 1; i < N; i++) {
		y[i] = b[i];
		double summ = 0;
		for (int j = 0; j < i; j++) {
			summ += L[i][j] * y[j];
		}
		summ = summ / L[i][i];
		y[i] -= summ;
	}

	// równanie Ux = y
	x[N - 1] = y[N - 1] / U[N - 1][N - 1];
	for (int i = N - 2; i >= 0; i--) {
		double value = y[i];
		for (int j = i + 1; j < N; j++) {
			value -= U[i][j] * x[j];
		}
		value = value / U[i][i];
		x[i] = value;
	}

	delete[] y;
	delete[] b;
	return x;
}


double** get_condition_matrix(int n, int m, double dt) {

	// tworze macierz n x m 
	double** A = new double *[n];
	for (int i = 0; i < n; i++)
		A[i] = new double[m];

	// wype³niam macie¿ 0
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			A[i][j] = 0.0;

	//warunek pocz¹tkowy U(x,0) = 1
	for (int j = 0; j < m; j++)
		A[0][j] = 1;

	//warunek brzegowy U(r,t) = 0
	for (int i = 1; i < n; i++)
		A[i][0] = 0;

	//warunek brzegowy U(r+a,t) = 1-r(r+a)erfc(a/(2*sqrt(D*t)))
	double t = dt;
	for (int i = 1; i < n; i++) {
		A[i][m - 1] = 1.0 - (((double)r / (double)(r + a)*erfc(a / (2.0*sqrt(D*t)))));
		t += dt;
	}
	return A;
}

double** get_analitical_values(int n, int m, double dt, double h) {
	double** A = get_condition_matrix(n, m, dt);
	double t = 0.0;
	double x;
	for (int i = 1; i < n; i++) { // t
		x = r + h;
		t += dt;
		for (int j = 1; j < m - 1; j++) { // x
			A[i][j] = 1.0 - r / x * calerf::ERFCL((x - r) / (2.0 * sqrt(D*t)));
			x += h;
		}
	}
	return A;
}

void print_matrix(double** A, int n, int m) {
	cout << "===========================================" << endl << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			cout << setw(10) << A[i][j] << " ";
		}
		cout << endl;
	}
}

double** alloc_0_matrix(int n, int m) {
	double **M = NULL;
	M = new double*[n];
	for (int i = 0; i < n; i++)
		M[i] = new double[m];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			M[i][j] = 0.0;
	return M;
}

double* alloc_0_vector(int n) {
	double *V = NULL;
	V = new double[n];
	for (int i = 0; i < n; i++)
		V[i] = 0.0;
	return V;
}



const long double aintl_(const long double x)
{ // This function has been added by L. K. Bieniasz
if(x >= 0.0L)return floorl(x);
else return ceill(x);
}





////////////////////////
// long double versions
////////////////////////


long double calerf::CALERFL(const long double arg, const int jint)
{
//-----------------------------------------------------------------
//
//  This packet evaluates ERF(x), ERFC(x), and exp(x*x)*ERFC(x)
//  for a real argument x.  It contains four long double functions:
//  ERFL(), ERFCL(), and erexl(), and CALERFL().
//
//  The calling statements for the primary entries are:
//
//         y = ERFL(x),
//
//         y = ERFCL(x),
//  and
//         y = erexl(x).
//
//  The function CALERFL() is intended for internal packet use only,
//  all computations within the packet being concentrated in this
//  routine. The functions invoke CALERFL() with the
//  statement
//
//         CALERFL(arg,jint)
//
//  where the parameter usage is as follows
//
//      Function                     Parameters for CALERFL()
//       call       l       arg            jint
//
//     ERFL(arg)     ANY REAL ARGUMENT      0
//     ERFCL(arg)    fabs(arg)  < XBIG      1
//     EREXL(arg)    XNEG < arg < XMAX      2
//
//  The main computation evaluates near-minimax approximations
//  from "Rational Chebyshev approximations for the error function"
//  by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
//  transportable program uses rational functions that theoretically
//  approximate  ERFL(x) and  ERFCL(x) to at least 18 significant
//  decimal digits. The accuracy achieved depends on the arithmetic
//  system, the compiler, the intrinsic functions, and proper
//  selection of the machine-dependent constants.
//
//******************************************************************
//
//  Explanation of machine-dependent constants
//
//   XMIN   = the smallest positive floating-point number.
//   XINF   = the largest positive finite floating-point number.
//   XNEG   = the largest negative argument acceptable to EREXL();
//            the negative of the solution to the equation
//            2*exp(x*x) = XINF.
//   XSMALL = argument below which ERFL(x) may be represented by
//            2*x/sqrt(pi) and above which  x*x  will not underflow.
//            A conservative value is the largest machine number x
//            such that   1.0 + x = 1.0   to machine precision.
//   XBIG   = largest argument acceptable to ERFCL();  solution to
//            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
//            W(x) = exp(-x*x)/[x*sqrt(pi)].
//   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
//            machine precision.  A conservative value is
//            1/[2*sqrt(XSMALL)]
//   XMAX   = largest acceptable argument to EREXL(); the minimum
//            of XINF and 1/[sqrt(pi)*XMIN].
//
//  Approximate values for some important machines are:
//
//                          XMIN       XINF        XNEG     XSMALL
//
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
//  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
//  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
//  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
//  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
//
//
//                          XBIG       XHUGE       XMAX
//
//  IEEE (IBM/XT,
//    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
//  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
//  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
//  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
//  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
//
//******************************************************************
//
//  Error returns
//
//  The program returns  erfc = 0      for  arg >= XBIG;
//
//                       erex = XINF   for  arg <  XNEG;
//      and
//                       erex = 0      for  arg >= XMAX.
//
//
//  Intrinsic functions required are:
//
//     fabsl(), expl()
//
//-----------------------------------------------------------------
//  Based on the netlib FORTRAN package by W. J. Cody,
//  Mathematics and Computer Science Division
//  Argonne National Laboratory
//  Argonne, IL 60439
//
//  Latest modification of the above package: March 19, 1990
//-----------------------------------------------------------------

//-----------------------------------------------------------------
//  Mathematical constants
//-----------------------------------------------------------------
const long double ZERO    =  0.0e0L;
const long double HALF    =  0.5e0L;
const long double ONE     =  1.0e0L;
const long double TWO     =  2.0e0L;
const long double FOUR    =  4.0e0L;
const long double SIXTEEN = 16.0e0L;

static const long double SQRPI  = 5.6418958354775628695e-1L;
static const long double THRESH = 0.46875e0L;

//-----------------------------------------------------------------
//  Machine-dependent constants (for IBM PC)
//-----------------------------------------------------------------
static const long double XINF   =    1.79e308L;
static const long double XNEG   = -26.628e0L;
static const long double XSMALL =  1.11e-16L; 
static const long double XBIG   =  26.543e0L;
static const long double XHUGE  =  1.0e10L;    // 6.71e7L;    // Modified by L. K. Bieniasz
static const long double XMAX   =  0.5e2466L;  // 2.53e307L;  // Modified by L.K. Bieniasz

//-----------------------------------------------------------------
//  Coefficients for approximation to  erf  in first interval
//-----------------------------------------------------------------
static const long double A[5] = {
                                3.16112374387056560e00L,
                                1.13864154151050156e02L,
                                3.77485237685302021e02L,
                                3.20937758913846947e03L,
                                1.85777706184603153e-1L
                                };

static const long double B[4] = {
                                2.36012909523441209e01L,
                                2.44024637934444173e02L,
                                1.28261652607737228e03L,
                                2.84423683343917062e03L
                                };
//-----------------------------------------------------------------
//  Coefficients for approximation to  erfc  in second interval
//-----------------------------------------------------------------
static const long double C[9] = {
                                5.64188496988670089e-1L,
                                8.88314979438837594e0L,
                                6.61191906371416295e01L,
                                2.98635138197400131e02L,
                                8.81952221241769090e02L,
                                1.71204761263407058e03L,
                                2.05107837782607147e03L,
                                1.23033935479799725e03L,
                                2.15311535474403846e-8L
                                };

static const long double D[8] = {
                                1.57449261107098347e01L,
                                1.17693950891312499e02L,
                                5.37181101862009858e02L,
                                1.62138957456669019e03L,
                                3.29079923573345963e03L,
                                4.36261909014324716e03L,
                                3.43936767414372164e03L,
                                1.23033935480374942e03L
                                };
//-----------------------------------------------------------------
//  Coefficients for approximation to  erfc  in third interval
//-----------------------------------------------------------------
static const long double P[6] = {
                                3.05326634961232344e-1L,
                                3.60344899949804439e-1L,
                                1.25781726111229246e-1L,
                                1.60837851487422766e-2L,
                                6.58749161529837803e-4L,
                                1.63153871373020978e-2L
                                };

static const long double Q[5] = {
                                2.56852019228982242e00L,
                                1.87295284992346047e00L,
                                5.27905102951428412e-1L,
                                6.05183413124413191e-2L,
                                2.33520497626869185e-3L
                                };
//-----------------------------------------------------------------

register int i;
long double del,x,xden,xnum,y,ysq;
long double result;

x = arg;
y = fabsl(x);
if(y <= THRESH)
  {
  //------------------------------------
  //  Evaluate  erf  for  |x| <= 0.46875
  //------------------------------------
  ysq = ZERO;
  if(y > XSMALL)ysq = y * y;
  xnum = A[4]*ysq;
  xden = ysq;
  for(i=0; i<3; i++)
     {
     xnum = (xnum + A[i])*ysq;
     xden = (xden + B[i])*ysq;
     }
  result = x * (xnum + A[3])/(xden + B[3]);
  if(jint != 0)result = ONE - result;
  if(jint == 2)result = expl(ysq)*result;

  return result;
  }
//-------------------------------------------
//  Evaluate  erfc  for 0.46875 <= |x| <= 4.0
//-------------------------------------------
else
if(y <= FOUR)
  {
  xnum = C[8]*y;
  xden = y;
  for(i=0; i<7; i++)
     {
     xnum = (xnum + C[i])*y;
     xden = (xden + D[i])*y;
     }
  result = (xnum + C[7])/(xden + D[7]);
  if(jint != 2)
    {   
    ysq = aintl_(y*SIXTEEN)/SIXTEEN; 
    
    del = (y-ysq)*(y+ysq);
    result = expl(-ysq*ysq) * expl(-del)*result;
    }
  }
//-------------------------------
//  Evaluate  erfc  for |x| > 4.0
//-------------------------------
else{
    result = ZERO;
    if(y >= XBIG)
      {
      if((jint != 2) || (y >= XMAX))goto LAB300;
      if(y >= XHUGE)
        {
        result = SQRPI/y;
        goto LAB300;
        }
      }
    ysq = ONE/(y * y);
    xnum = P[5]*ysq;
    xden = ysq;
    for(i=0; i<4; i++)
       {
       xnum = (xnum + P[i])*ysq;
       xden = (xden + Q[i])*ysq;
       }
    result = ysq *(xnum + P[4])/(xden + Q[4]);
    result = (SQRPI - result)/y;
    if(jint != 2)
      {
      ysq = aintl_(y*SIXTEEN)/SIXTEEN; 
      del = (y-ysq)*(y+ysq);
      result = expl(-ysq*ysq) * expl(-del)*result;
      }
    }

//-----------------------------------------
//  Fix up for negative argument, erf, etc.
//-----------------------------------------
LAB300:
if(jint == 0)
  {
  result = (HALF - result)+HALF;
  if(x < ZERO)result = -result;
  }
else
if(jint == 1)
  {
  if(x < ZERO)result = TWO - result;
  }
else{ // jint == 2

    if(x < ZERO)
      {
      if(x < XNEG)result = XINF;
      else{
          ysq = aintl_(x*SIXTEEN)/SIXTEEN;
          
          
          del = (x-ysq)*(x+ysq);
          y = expl(ysq*ysq) * expl(del);
          result = (y+y)-result;
          }
      }
    }

return result;
}





long double calerf::ERFL(const long double x)
{
//-------------------------------------------------------------------
//  This function computes approximate values for ERF(x).
//  (see comments heading CALERF()).
//
//  Based on the netlib package by W. J. Cody, January 8, 1985
//-------------------------------------------------------------------
return calerf::CALERFL(x,0);
}






long double calerf::ERFCL(const long double x)
{
//-------------------------------------------------------------------
//  This function computes approximate values for ERFC(x).
//  (see comments heading CALERF()).
//
//  Based on the netlib package by W. J. Cody, January 8, 1985
//-------------------------------------------------------------------
return calerf::CALERFL(x,1);
}





long double calerf::EREXL(const long double x)
{
//-----------------------------------------------------------------
//  This function computes approximate values for
//  exp(x*x) * ERFC(x).
//  (see comments heading CALERF()).
//
//  Based on the netlib package by W. J. Cody, March 30, 1987
//-----------------------------------------------------------------
return calerf::CALERFL(x,2);
}
