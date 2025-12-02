#ifndef GAUSS_H
#define GAUSS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 100

void gauss_elimination(int n, double a[MAX_SIZE+1][MAX_SIZE+1], double b[MAX_SIZE+1], double X[MAX_SIZE+1]) {
    // Implementación similar a la tuya pero sin errores
    for (int i = 0; i < n; i++) {
        // Pivoteo parcial
        int max_row = i;
        for (int k = i+1; k < n; k++) {
            if (fabs(a[k][i]) > fabs(a[max_row][i])) {
                max_row = k;
            }
        }
        
        // Intercambiar filas
        for (int k = 0; k < n; k++) {
            double temp = a[i][k];
            a[i][k] = a[max_row][k];
            a[max_row][k] = temp;
        }
        double temp = b[i];
        b[i] = b[max_row];
        b[max_row] = temp;
        
        // Eliminación gaussiana
        for (int k = i+1; k < n; k++) {
            double factor = a[k][i] / a[i][i];
            for (int j = i; j < n; j++) {
                a[k][j] -= factor * a[i][j];
            }
            b[k] -= factor * b[i];
        }
    }
    
    // Sustitución hacia atrás
    for (int i = n-1; i >= 0; i--) {
        X[i] = b[i];
        for (int j = i+1; j < n; j++) {
            X[i] -= a[i][j] * X[j];
        }
        X[i] /= a[i][i];
    }
}

#endif
