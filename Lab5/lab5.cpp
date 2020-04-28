#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#define N 4

using namespace std;
string fname = "data.txt";

// funkcje pomocnicze
void print_vector(double *V, string what);
void print_matrix(double **M, string what);
double** alloc_matrix();
double** task_matrix();

// główne funkcje rozwiązujące równanie Ax = b metodą LU(z gaussem) z częsciowym wyborem el. podstawowego 
void gauss(double **U, double **L, double** T);
void resoult(double **U, double **L, double **T, double* x);
int partial_pivot(double** M, double** L, int row, double** T);

int main()
{
	cout << fixed;

	// U - upper triangular, L - lower triangular, P - macierz przekształceń ( PA = LU )
	double** U = task_matrix();
	double** L = alloc_matrix();
	double** P = alloc_matrix();
	double* x = new double[N];
	print_matrix(U, "A");

	gauss(U, L, P);
	resoult(U, L, P, x);


	print_vector(x, "x");
}



void gauss(double **U, double **L, double** P) {
	fstream o;
	o.open(fname, ios::out | ios::app);
	o << "Eliminacja gaussa z częsciowym wyborem el. podstawowego: " << endl;

	// inicjalizuje Macierz przekształceń
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (i == j) P[i][j] = 1.0;
			else P[i][j] = 0.0;
		}
	}


	// Wyznaczam macierz LU za pomocą eliminacji Gaussa
	// obliczam temp ( gdyż licznik zerujemy w pierwszym przejściu pętli i w kolejnych kolumnach program nic by nie zmienił)
	for (int k = 0; k < N; k++) {
		print_matrix(U," Etap: " + to_string(k));
		// if el. podstawowy == 0 --> pivoting (częściowy wybór elementu podstawowego).
		if (U[k][k] == 0.0 && partial_pivot(U, L, k, P) == 0) {
				// nie znaleziono elementu podstawowego. błąd.
				o << "error" << endl;
				return;
		}

		// zerowanie kolejnej kolumny U i uzupełnianie L
		for (int i = 1 + k; i < N; i++) {
			double temp = U[i][k] / U[k][k];
			for (int j = k; j < N; j++) {
				U[i][j] -= U[k][j] * temp;
				L[i][j] = temp;
			}
		}
	}

	// z Macierzy L tworze macierz trojkątną dolną (obliczone wsp. poniżej przekątnej (diagonali) pozostaja bez zmian)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			L[i][i] = 1;
			if (j > i)
				L[i][j] = 0.0;
		}
	}

	print_matrix(U, "U");
	print_matrix(L, "L");
	print_matrix(P, "Macierz przekształceń wierszy");

	o << "\n\n\n\n WYNIK:\n";
}

int partial_pivot(double** U, double **L, int x, double** P) {

	// Szukam max elementu w kolumnie pod M[x][x]
	int max_index = x;
	for (int i = x + 1; i < N; i++) {
		if (fabs(U[i][x]) > U[x][x])
			max_index = i;
	}

	if (max_index == x) {
		// nie właściwa macierz 
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

void resoult(double **U, double **L, double **P, double* x) {
	double *y = new double[N];
	double *b = new double[N];

	b[0] = 35.0;
	b[1] = 104.0;
	b[2] = -366.0;
	b[3] = -354.0;

	// przekształcam macierz wyników b,  zgodnie z macierzą przekształceń P
	double *temp = new double[N];
	for (int i = 0; i < N; i++) {
		temp[i] = 0.0;
		for (int j = 0; j < N; j++) {
			temp[i] += P[i][j] * b[j];
		}
	}
	b = temp;

	print_vector(b, "b (przekształcone)");

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

	print_vector(y, "y");

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
}


void print_matrix(double **M, string what) {
	fstream o;
	o.open("data.txt", ios::out | ios::app);
	o << "------------------- " << what << " -------------------------" << endl;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			o << setw(11) << M[i][j];
		}
		o << endl;
	}
	o << "\n\n\n";
}

void print_vector(double *V, string what) {
	fstream o;
	o.open("data.txt", ios::out | ios::app);
	o << "------------------- " << what << " -------------------------" << endl;

	for (int i = 0; i < N; i++)
		o << V[i] << endl;
	o << "\n\n\n";
}


double** alloc_matrix() {
	double **M = NULL;
	M = new double*[N];
	for (int i = 0; i < N; i++)
		M[i] = new double[N];
	return M;
}


double** task_matrix() {
	double** A = alloc_matrix();
	A[0][0] = 1.0;
	A[0][1] = -20.0;
	A[0][2] = 30.0;
	A[0][3] = -4.0;
	A[1][0] = 2.0;
	A[1][1] = -40.0;
	A[1][2] = -6.0;
	A[1][3] = 50.0;
	A[2][0] = 9.0;
	A[2][1] = -180.0;
	A[2][2] = 11.0;
	A[2][3] = -12.0;
	A[3][0] = -16.0;
	A[3][1] = 15.0;
	A[3][2] = -140.0;
	A[3][3] = 13.0;
	return A;
}


