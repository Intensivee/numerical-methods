#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#define N 4
#define param 0.5

using namespace std;

void print_vector(double *V, string what);
void print_matrix(double **M, string what);

double** alloc_matrix();
double** alloc_A_matrix();
double* alloc_b_vector();
double* alloc_x_vector();
double* vector_plus_vector(double* v, double* u);
double** matrix_plus_matrix(double** A, double** B);
double** matrix_multiply_matrix(double** A, double** B);
double* matrix_multiply_vector(double** A, double* v);
double** matrix_negation(double** A);
double get_estimator(double* Old, double* New);
double get_residuum(double* x);

void Jacobiego(int n_max, double tol_x, double tol_f);
void Gaussa_Seidela(int n_max, double tol_x, double tol_f);
void SOR(int n_max, double tol_x, double tol_f);

int main()
{
	Jacobiego(15, 0.000000000000001, 0.000000000000001);
	Gaussa_Seidela(15, 0.00000000000001, 0.00000000000001);
	SOR(60, 0.00000000000001, 0.00000000000001);
}

void Jacobiego(int n_max, double tol_x, double tol_f) {
	fstream o;
	o.open("data.txt", ios::out | ios::app);

	//A = L_U + D, gdzie: L_U = L + U
	double **L_U = alloc_A_matrix();
	double **D_inv = alloc_A_matrix();
	double *b = alloc_b_vector();
	double *x = alloc_x_vector();

	// tworze macierze L_U = L+U oraz odwróconą macierz D * (-1)
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {

			// wszystko poza przekątną
			if (j != i)
				D_inv[i][j] = 0;

			// przekątna
			else {
				if (D_inv[i][j] != 0)
					D_inv[i][j] = 1.0 / D_inv[i][j];
				else {
					o << " Nie istnieje macierz odwrotna do macierzy D." << endl;
					return;
				}
				L_U[i][j] = 0;
			}
		}
	}


	// wyznaczam macierz M = -D^-1 * (L+U) oraz c = D^-1 * b
	double* c = matrix_multiply_vector(D_inv, b);
	double** M = matrix_negation(matrix_multiply_matrix(D_inv, L_U));
	double* new_x = new double[N];


	o << "-----------------------------------------------------   Jacobiego   -------------------------------------------------------" << endl;
	o << fixed << setw(7) << "  n |" << setw(25) << "          x1            |" << setw(25) << "          x2            |";
	o << setw(25) << "          x3            |" << setw(25) << "          x4            |" << setw(25) << "          e_n           |";
	o << setw(25) << "       residuum         |" << endl;
	o << setprecision(16) << setw(5) << 0 << " |" << setw(23) << x[0] << " |" << setw(23) << x[1] << " |";
	o << setw(23) << x[2] << " |" << setw(23) << x[3] << " |" << setw(23) << "------------------------|";
	o << setw(23) << "------------------------|" << endl;

	double e_n;
	double residuum;

	// głowna pętla iterująca
	for (int i = 1; i <= n_max; i++) {
		new_x = vector_plus_vector(matrix_multiply_vector(M, x), c);
		e_n = get_estimator(new_x, x);
		residuum = get_residuum(new_x);
		x = new_x;


		// zapisywanie do pliku
		o << setprecision(16) << setw(5) << i << " |" << setw(23) << x[0] << " |" << setw(23) << x[1] << " |";
		o << setw(23) << x[2] << " |" << setw(23) << x[3] << " |" << setw(23) << e_n << " |";
		o << setw(23) << residuum << " |" << endl;


		// warunek tol_x i residuum
		if (e_n < tol_x && residuum < tol_f) {
			o << "--------------------Spełniono Warunek dokladnosci X* oraz warunek residuum-------------------\n" << endl;
			print_vector(x, "WYNIK");
			return;
		}
	}

	o << "--------------------Spełniono warunek ilości iteracji-------------------\n" << endl;
	print_vector(x, "WYNIK");


}

void Gaussa_Seidela(int n_max, double tol_x, double tol_f) {
	fstream o;
	o.open("data.txt", ios::out | ios::app);

	//A = L_D + U, gdzie: L_D = L + D
	double **L_D_inv = alloc_matrix();
	double **U = alloc_A_matrix();
	double *b = alloc_b_vector();
	double *x = alloc_x_vector();

	// Macierz L_D_inv = (L+D)^-1 
	double temp[4][4] = { {1.0 / 100.0, 0.0, 0.0, 0.0},
					{-1.0 / 20000.0, 1.0 / 200.0, 0.0, 0.0 },
					{ 101.0 / 1500000.0, -1.0 / 15000.0, 1.0 / 300.0, 0.0},
					{-15327.0 / 200000000.0, 127.0 / 2000000.0,-1.0 / 20000.0, 1.0 / 400.0} };

	// setting temp into L_D_inv and creating Upper from A
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			L_D_inv[i][j] = temp[i][j];
			if (j <= i) U[i][j] = 0.0;
		}
	}

	// wyznaczam macierz M = -D^-1 * (L+U) oraz c = D^-1 * b
	double* c = matrix_multiply_vector(L_D_inv, b);
	double** M = matrix_negation(matrix_multiply_matrix(L_D_inv, U));



	o << "-----------------------------------------------------   Gaussa-Seidela   -------------------------------------------------------" << endl;
	o << fixed << setw(7) << "  n |" << setw(25) << "          x1            |" << setw(25) << "          x2            |";
	o << setw(25) << "          x3            |" << setw(25) << "          x4            |" << setw(25) << "          e_n           |";
	o << setw(25) << "       residuum         |" << endl;
	o << setprecision(16) << setw(5) << 0 << " |" << setw(23) << x[0] << " |" << setw(23) << x[1] << " |";
	o << setw(23) << x[2] << " |" << setw(23) << x[3] << " |" << setw(23) << "------------------------|";
	o << setw(23) << "------------------------|" << endl;

	double e_n;
	double residuum;
	double* new_x = new double[N];

	// głowna pętla iterująca
	for (int i = 1; i <= n_max; i++) {
		new_x = vector_plus_vector(matrix_multiply_vector(M, x), c);
		e_n = get_estimator(new_x, x);
		residuum = get_residuum(new_x);
		x = new_x;


		// zapisywanie do pliku
		o << setprecision(16) << setw(5) << i << " |" << setw(23) << x[0] << " |" << setw(23) << x[1] << " |";
		o << setw(23) << x[2] << " |" << setw(23) << x[3] << " |" << setw(23) << e_n << " |";
		o << setw(23) << residuum << " |" << endl;


		// warunek tol_x i residuum
		if (e_n < tol_x && residuum < tol_f) {
			o << "--------------------Spełniono Warunek dokladnosci X* oraz warunek residuum-------------------\n" << endl;
			print_vector(x, "WYNIK");
			return;
		}
	}

	o << "--------------------Spełniono warunek ilości iteracji-------------------\n" << endl;
	print_vector(x, "WYNIK");


}

void SOR(int n_max, double tol_x, double tol_f) {
	fstream o;
	o.open("data.txt", ios::out | ios::app);

	//A = L_D + U, gdzie: L_D = L + D
	double **L_param_D_inv = alloc_matrix();
	double **D = alloc_A_matrix();
	double **U = alloc_A_matrix();
	double *b = alloc_b_vector();
	double *x = alloc_x_vector();

	// przekształcam D w (1-1/param) D oraz ustawiam U z macierzy A
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (j <= i) U[i][j] = 0.0;
			if (j != i)	D[i][j] = 0.0;
			if (j == i) D[i][j] = (1.0 - 1.0 / param) * D[i][j];
		}
	}

	// Macierz L_param_D_inv = (L+(1/param)*D)^-1 
	double temp[4][4] = { {1.0 / 200.0, 0.0 ,0.0 ,0.0},
					{-1.0 / 80000.0, 1.0 / 400.0, 0.0 , 0.0 },
					{67.0 / 4000000.0, -1.0 / 60000.0, 1.0 / 600.0, 0.0 },
					{-15163.0 / 800000000.0, 63.0 / 4000000.0, -1.0 / 80000.0, 1.0 / 800.0} };
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			L_param_D_inv[i][j] = temp[i][j];
		}
	}


	// wyznaczam macierz M oraz c
	double* c = matrix_multiply_vector(L_param_D_inv, b);
	double** M = matrix_multiply_matrix(L_param_D_inv, matrix_plus_matrix(D, U));
	M = matrix_negation(M);
	double* new_x = new double[N];




	o << "-----------------------------------------------------   SOR   -------------------------------------------------------" << endl;
	o << fixed << setw(7) << "  n |" << setw(25) << "          x1            |" << setw(25) << "          x2            |";
	o << setw(25) << "          x3            |" << setw(25) << "          x4            |" << setw(25) << "          e_n           |";
	o << setw(25) << "       residuum         |" << endl;
	o << setprecision(16) << setw(5) << 0 << " |" << setw(23) << x[0] << " |" << setw(23) << x[1] << " |";
	o << setw(23) << x[2] << " |" << setw(23) << x[3] << " |" << setw(23) << "------------------------|";
	o << setw(23) << "------------------------|" << endl;

	double e_n;
	double residuum;

	// głowna pętla iterująca
	for (int i = 1; i <= n_max; i++) {
		new_x = vector_plus_vector(matrix_multiply_vector(M, x), c);
		e_n = get_estimator(new_x, x);
		residuum = get_residuum(new_x);
		x = new_x;


		// zapisywanie do pliku
		o << setprecision(16) << setw(5) << i << " |" << setw(23) << x[0] << " |" << setw(23) << x[1] << " |";
		o << setw(23) << x[2] << " |" << setw(23) << x[3] << " |" << setw(23) << e_n << " |";
		o << setw(23) << residuum << " |" << endl;


		// warunek tol_x i residuum
		if (e_n < tol_x && residuum < tol_f) {
			o << "--------------------Spełniono Warunek dokladnosci X* oraz warunek residuum-------------------\n" << endl;
			print_vector(x, "WYNIK");
			return;
		}
	}

	o << "--------------------Spełniono warunek ilości iteracji-------------------\n" << endl;
	print_vector(x, "WYNIK");


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
	o << fixed  << setprecision(16);
	o.open("data.txt", ios::out | ios::app);
	o << "------------------- " << what << " -------------------------" << endl;

	for (int i = 0; i < N; i++)
		o << V[i] << endl;
	o << "\n\n\n";
}

double* alloc_b_vector() {
	double *V = NULL;
	V = new double[N];
	V[0] = 116.0;
	V[1] = -226.0;
	V[2] = 912.0;
	V[3] = -1174.0;
	return V;
}

double* alloc_x_vector() {
	double *V = NULL;
	V = new double[N];
	V[0] = 2.0;
	V[1] = 2.0;
	V[2] = 2.0;
	V[3] = 2.0;
	return V;
}

double** alloc_A_matrix() {
	// nie zczytuje z pliku żeby program był bardziej mobilny (mail)
	double **M = NULL;
	M = new double*[N];
	for (int i = 0; i < N; i++)
		M[i] = new double[N];

	M[0][0] = 100.0;
	M[0][1] = -1.0;
	M[0][2] = 2.0;
	M[0][3] = -3.0;
	M[1][0] = 1.0;
	M[1][1] = 200.0;
	M[1][2] = -4.0;
	M[1][3] = 5.0;
	M[2][0] = -2.0;
	M[2][1] = 4.0;
	M[2][2] = 300.0;
	M[2][3] = -6.0;
	M[3][0] = 3.0;
	M[3][1] = -5.0;
	M[3][2] = 6.0;
	M[3][3] = 400.0;
	return M;
}

double** alloc_matrix() {
	double **M = NULL;
	M = new double*[N];
	for (int i = 0; i < N; i++)
		M[i] = new double[N];
	return M;
}

double** matrix_multiply_matrix(double** A, double** B) {
	double **M = alloc_matrix();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			M[i][j] = 0.0;
			for (int k = 0; k < N; k++) {
				M[i][j] += A[i][k] * B[k][j];
			}

		}
	}
	return M;
}

double* matrix_multiply_vector(double** A, double* V) {
	double *U = new double[N];
	for (int i = 0; i < N; i++) {
		U[i] = 0.0;
		for (int j = 0; j < N; j++) {
			U[i] += A[i][j] * V[j];
		}
	}
	return U;
}

double* vector_plus_vector(double* v, double* u) {
	double* x = new double[N];
	for (int i = 0; i < N; i++)
		x[i] = v[i] + u[i];
	return x;
}

double** matrix_negation(double** A) {
	double **M = alloc_matrix();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			M[i][j] = -A[i][j];
		}
	}
	return M;
}

double get_estimator(double* v, double* u) {

	double max = fabs(v[0] - u[0]);
	for (int i = 1; i < N; i++) {
		if (fabs(v[i]-u[i]) > max)
			max = fabs(v[i]-u[i]);
	}
	return max;
}

double get_residuum(double* x) {

	// y = Ax - b = 0
	double * y = matrix_multiply_vector(alloc_A_matrix(), x);
	double * b = alloc_b_vector();
	for (int i = 0; i < N; i++) {
		y[i] -= b[i];
	}

	// znajduje najwieksza wartość z pośród residuum 
	double max = fabs(y[0]);
	for (int i = 1; i < N; i++) {
		if (fabs(y[i]) > max)
			max = fabs(y[i]);
	}
	return max;
}

double** matrix_plus_matrix(double** A, double** B) {
	double** M = alloc_matrix();
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			M[i][j] = A[i][j] + B[i][j];
		}
	}
	return M;
}

