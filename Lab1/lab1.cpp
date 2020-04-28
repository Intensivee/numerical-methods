#include <iostream>

using namespace std;

int main()
{

//--------------DOUBLE--------------------

	int t = 0;					// liczba bitów mantysy, zaczynamy od 0, w while sprawdzamy czy dla t=1 czyli liczby (1.1)bin  -> (1.5)dec, jak tak to inkrementujemy t
	double dx = 0.5;			// epsilon maszynowy, zaczynamy od wartosci 0.5    czyli 2^-t, gdzie t = 1 (sprawdzamy czy dla t=1 zachodzi while);
	while (1.0 + dx > 1.0) {
		t++;				
		dx /= 2.0;				// sprawdzamy dla większej dokładności liczby (f=0.f1f2f3.. ft)bin dla t=2 -> (f=0.01)bin = (f=0.25)dec itd. 
	}
					
	dx *= 2.0;					// cofamy sie do ostatniego epsilon kiedy wszedł do pentli
	cout << "DOUBLE: " << endl;
	cout << "Epsilon maszynowy:      " << dx << endl;
	cout << "DBL_EPSILON:            " << DBL_EPSILON << endl;
	cout << "Liczba bitow mantysy:   " << t << endl;


//-------------FLOAT----------------------
	
	int t2 = 0;
	float dxx = 0.5f;
	while (1.0f + dxx > 1.0f) {
		t2++;
		dxx /= 2.0f;
	}

	dxx *= 2.0f;

	cout << "\nFLOAT: " << endl;
	cout << "Epsilon maszynowy:      " << dxx << endl;
	cout << "FLT_EPSILON:            " << FLT_EPSILON << endl;
	cout << "Liczba bitow mantysy:   " << t2 << endl;




}

