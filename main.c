/***********************************************************************
 *H4.1
 * Kompilieren mit
 *
 *   gcc  -Wall -pedantic main.c -o main -lm
 *
 *
 * Ausführen mit
 *   ./main
 *
 **********************************************************************/


#include <stdio.h>
#include <math.h>

#define N 1000  // Anzahl der Gitterpunkte
#define x0 0.0  // Startpunkt
#define x1 60.0  // Endpunkt
#define u0 0.0  // Randbedingung am Startpunkt
#define u1 0.0  // Randbedingung am Endpunkt

double g(double x) {
    // Funktion g(x) für g(x) ≡ 0
    return 0.0;
}

double eigenvalue(int n, double h, double *u) {
    // Berechnung des Eigenwerts für die gegebene Eigenfunktion u
    // Hier können Sie das Verfahren Ihrer Wahl zur Eigenwertbestimmung verwenden
    // Beispiel: Bisektionsverfahren
    double a = -10.0;  // Untere Grenze des Eigenwertbereichs
    double b = 10.0;  // Obere Grenze des Eigenwertbereichs
    double tol = 1e-6;  // Toleranz für die Genauigkeit
    double lambda = 0.0;  // Initialer Schätzwert für den Eigenwert
    
    while (fabs(b - a) > tol) {
        double c = (a + b) / 2.0;
        double f_a = u[n-1] + c * u[n];
        double f_c = u[1] + c * u[0];
        
        if (f_a * f_c < 0) {
            b = c;
        } else {
            a = c;
        }
        
        lambda = (a + b) / 2.0;
    }
    
    return lambda;
}

void numerov_method() {
    double h = (x1 - x0) / N;  // Schrittweite
    
    double u[N+1];  // Eigenfunktion u(x)
    double v[N+1];  // Funktion v(x) = g(x) * u(x)
    
    // Randbedingungen
    u[0] = u0;
    u[1] = u[0] + h;
    
    for (int i = 1; i <= N; i++) {
        double x = x0 + i * h;
        v[i] = g(x) * u[i];
    }
    
    for (int i = 1; i < N; i++) {
        double Q_i = 1 - (5.0/12.0) * h * h * v[i];
        double S_i = 1 + (1.0/12.0) * h * h * v[i-1];
        
        u[i+1] = (2 * (1 - (1.0/12.0) * h * h * v[i+1]) * u[i] - S_i * u[i-1]) / Q_i;
    }
    
    double eigenvalues[10];
    
    // Eigenwerte berechnen
    for (int i = 0; i < 10; i++) {
        eigenvalues[i] = eigenvalue(N, h, u);
        printf("Eigenwert %d: %f\n", i+1, eigenvalues[i]);
    }
}

int main() {
    numerov_method();
    return 0;
}
