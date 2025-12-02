#ifndef GAUSS_H
#define GAUSS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 100

// Eliminación gaussiana (corregida)
void gauss_elimination(int n, double a[MAX_SIZE+1][MAX_SIZE+1], 
                       double b[MAX_SIZE+1], double X[MAX_SIZE+1]) {
    int i, j, k;
    double factor, temp;
    
    // Eliminación gaussiana con pivoteo parcial
    for (i = 0; i < n; i++) {
        // Pivoteo parcial: encontrar la fila con el mayor elemento en la columna i
        int max_row = i;
        double max_val = fabs(a[i][i]);
        
        for (j = i + 1; j < n; j++) {
            if (fabs(a[j][i]) > max_val) {
                max_val = fabs(a[j][i]);
                max_row = j;
            }
        }
        
        // Intercambiar filas si es necesario
        if (max_row != i) {
            for (k = 0; k < n; k++) {
                temp = a[i][k];
                a[i][k] = a[max_row][k];
                a[max_row][k] = temp;
            }
            temp = b[i];
            b[i] = b[max_row];
            b[max_row] = temp;
        }
        
        // Verificar que el pivote no sea cero
        if (fabs(a[i][i]) < 1e-10) {
            printf("Error: Matriz singular o casi singular.\n");
            exit(1);
        }
        
        // Eliminación hacia adelante
        for (j = i + 1; j < n; j++) {
            factor = a[j][i] / a[i][i];
            for (k = i; k < n; k++) {
                a[j][k] -= factor * a[i][k];
            }
            b[j] -= factor * b[i];
        }
    }
    
    // Sustitución hacia atrás
    for (i = n - 1; i >= 0; i--) {
        X[i] = b[i];
        for (j = i + 1; j < n; j++) {
            X[i] -= a[i][j] * X[j];
        }
        X[i] /= a[i][i];
    }
}

#endif
