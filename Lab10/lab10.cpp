#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace std;


double analitical_value(double t) {
	return 1 - exp(-10.0 * (t + atan(t)));
}


// bezpośrednia metoda Eulera
void solve(double step, double val) {
	ofstream o("data.txt");

	double BME = 0.0;
	double PME = 0.0;
	double PMT = 0.0;
	double temp;
	int iter = 0;

	for (double t = 0; t < val; t += step, iter++) {

		// BME
		BME = BME - (10.0*t*t + 20.0) / (t*t + 1.0) * (BME - 1.0) * step;

		// PME - po rozwiązaniu równania liniowego 
		temp = (10.0*(t + step)*(t + step) + 20.0) / ((t + step) * (t + step) + 1.0);
		PME = (PME + temp * step) / (1.0 + temp * step);

		//PMT - po rozwiązaniu równania liniowego 
		temp = (10.0 * (t + step) * (t + step) + 20.0) / ((t + step) * (t + step) + 1.0);
		PMT = ((-step / 2.0) * (((10.0 * t * t + 20.0) / (t * t + 1.0)) * (PMT - 1.0) - temp) + PMT) / (1.0 + (step / 2.0) * temp);
		
		// zapisuje tylko co 1000 wartość (żeby wykres był czytelny)
		if (iter % 1000 == 0) {

			// wartości funkcji
			o << fixed << setprecision(16) << t << " " << BME << " " << PME << " ";
			o <<	PMT << " " << analitical_value(t) << endl;
		}
	}
}

// w tym przypadku limitem dla pętli jest ilosc iteracji
// a nie max wartość, bo dla bardzo małych kroków trwało by to za długo
void error_solve(double step, double N) {
	ofstream o2("bledy.txt");

	double BME, PME, PMT;
	double BME_error, PME_error, PMT_error;
	double temp, temp_error, t;

	while (step + 1 > 1) {
		BME = PME = PMT = 0.0;
		BME_error = PME_error = PMT_error = 0.0;

		// BME
		t = step;
		for (int i = 0; i < N; i++) {

			BME = BME - (10.0*t*t + 20.0) / (t*t + 1.0) * (BME - 1.0) * step;

			temp_error = fabs(BME - analitical_value(t));
			if (temp_error > BME_error)
				BME_error = temp_error;

			t += step;
		}

		// PME
		t = step;
		for (int i = 0; i < N; i++) {

			temp = (10.0*(t + step)*(t + step) + 20.0) / ((t + step) * (t + step) + 1.0);
			PME = (PME + temp * step) / (1.0 + temp * step);

			temp_error = fabs(PME - analitical_value(t));
			if (temp_error > PME_error)
				PME_error = temp_error;

			t += step;
		}

		// PMT
		t = step;
		for (int i = 0; i < N; i++) {

			temp = (10.0 * (t + step) * (t + step) + 20.0) / ((t + step) * (t + step) + 1.0);
			PMT = ((-step / 2.0) * (((10.0 * t * t + 20.0) / (t * t + 1.0)) * (PMT - 1.0) - temp) + PMT) / (1.0 + (step / 2.0) * temp);

			temp_error = fabs(PMT - analitical_value(t));
			if (temp_error > PMT_error)
				PMT_error = temp_error;

			t += step;
		}
		o2 << fixed << setprecision(16) << log10(step) << " " << log10(BME_error) << " " << log10(PME_error) << " " << log10(PMT_error) << endl;
		step *= 0.9;
	}
}


void BME_single_solve(double step, double val) {
	ofstream o("BME_unstable.txt");

	double BME = 0.0;
	double temp;
	int iter = 0;

	for (double t = 0; t < val; t += step, iter++) {

		// BME
		BME = BME - (10.0*t*t + 20.0) / (t*t + 1.0) * (BME - 1.0) * step;

		o << fixed << setprecision(16) << t << " " << BME << " " << analitical_value(t) << endl;
	}
}

int main()
{
	//solve(0.00005, 4);
	//BME_single_solve(0.13, 4);
	error_solve(0.2, 10000);
}


