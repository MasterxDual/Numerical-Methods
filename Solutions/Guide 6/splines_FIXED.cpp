#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 50
#define MAX_SIZE 100 
#include "gauss_FIXED.h"

/**
 * Function to read Xi, Yi data pairs from a file
 * Expected format: first line contains number of points, then xi yi pairs
 * @param filename Name of the file to read
 * @param X Array to store X values
 * @param Y Array to store Y values
 * @param n Pointer to store the number of data points read
 * @return 1 if successful, 0 otherwise
 */
int read_data_points(const char* filename, double X[], double Y[], int* n);

/**
 * Function to display the data points
 * @param X Array of X values
 * @param Y Array of Y values
 * @param n Number of data points
 */
void print_data_points(double X[], double Y[], int n);

/**
 * Function to print the cubic spline coefficients
 * @param X Array of X values (data points)
 * @param solution Array containing the solution coefficients
 * @param n Number of data points
 */
void print_cubic_splines(double X[], double solution[], int n);

/**
 * Function to evaluate the cubic spline at a given point x
 * @param X Array of X values (data points)
 * @param solution Array of spline coefficients
 * @param n Number of data points
 * @param x The x value to evaluate
 * @return The interpolated y value
 */
double evaluate_spline(double X[], double solution[], int n, double x);

/**
 * Function to compute and print the linear splines
 * @param X Array of X values (data points)
 * @param Y Array of Y values (data points)
 * @param n Number of data points
 */
void linear_spline(double X[], double Y[], int n);

int main(int argc, char const *argv[]) {
    // Arrays for data points to read from text file
    double X[MAX_POINTS], Y[MAX_POINTS];
    // Arrays for polynomial coefficients calculation
    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
    // Solution of Least Squares polynomial
    double a[MAX_SIZE+1];
    // n = number of data points and degree is the polynomial degree
    int n, option, row1, row2;
    

    // Read data points from file
    if (!read_data_points("data.txt", X, Y, &n)) {
        printf("Failed to read data from file. Exiting.\n");
        return 1;
    }
    
    // Print the data points
    print_data_points(X, Y, n);

    // Initialize A and b to zero
    for (int i = 0; i < 4*(n-1); i++) {
        b[i] = 0.0;
        for (int j = 0; j < 4*(n-1); j++) {
            A[i][j] = 0.0;
        }
    }

    printf("Choose an option for spline interpolation:\n");
    printf("1. Linear Spline\n");
    printf("2. Cubic Spline\n");
    scanf("%d", &option);

    switch(option) {
        case 1:
            linear_spline(X, Y, n);
            // We evalute points x_hat to verify if the spline works correctly
            printf("Do you want to evaluate a point x_hat? (1 for yes, 0 for no): ");
            scanf("%d", &option);
            while(option) {
                double x_hat;
                printf("Enter the value of x_hat: ");
                scanf("%lf", &x_hat);
                // Find the right interval for x_hat
                if(x_hat < X[0] || x_hat > X[n-1]) {
                    printf("x_hat is out of bounds [%lf, %lf]. Please enter a value within the range.\n", X[0], X[n-1]);
                } else {
                    // Locate the interval
                    int k = 0;
                    while(k < n-1 && x_hat > X[k+1]) {
                        k++;
                    }
                    // Linear interpolation formula
                    double y_hat = Y[k] + (Y[k+1] - Y[k]) * (x_hat - X[k]) / (X[k+1] - X[k]);
                    printf("The interpolated value at x_hat = %.3f is y_hat = %.3f\n", x_hat, y_hat);
                }
                printf("Do you want to evaluate another point x_hat? (1 for yes, 0 for no): ");
                scanf("%d", &option);
            }
            break;
        case 2: {
            int num_coeffs = 4 * (n - 1);
            
            // Inicializar matriz A y vector b a cero
            for (int i = 0; i < num_coeffs; i++) {
                b[i] = 0.0;
                for (int j = 0; j < num_coeffs; j++) {
                    A[i][j] = 0.0;
                }
            }
            
            // 1. Condiciones de interpolación (cada spline pasa por sus puntos extremos)
            for (int i = 0; i < n-1; i++) {
                // Condición en el punto izquierdo X[i]
                int row = 2 * i;
                A[row][4*i] = pow(X[i], 3);
                A[row][4*i+1] = pow(X[i], 2);
                A[row][4*i+2] = X[i];
                A[row][4*i+3] = 1.0;
                b[row] = Y[i];
                
                // Condición en el punto derecho X[i+1]
                row = 2 * i + 1;
                A[row][4*i] = pow(X[i+1], 3);
                A[row][4*i+1] = pow(X[i+1], 2);
                A[row][4*i+2] = X[i+1];
                A[row][4*i+3] = 1.0;
                b[row] = Y[i+1];
            }
            
            // 2. Continuidad de primera derivada en nodos interiores
            for (int i = 1; i < n-1; i++) {
                int row = 2 * (n-1) + (i-1);
                
                // Primera derivada del spline i-1 en X[i]
                A[row][4*(i-1)] = 3 * pow(X[i], 2);    // 3a*x²
                A[row][4*(i-1)+1] = 2 * X[i];          // 2b*x
                A[row][4*(i-1)+2] = 1.0;               // c
                
                // Primera derivada del spline i en X[i] (negativa para igualar)
                A[row][4*i] = -3 * pow(X[i], 2);
                A[row][4*i+1] = -2 * X[i];
                A[row][4*i+2] = -1.0;
                
                b[row] = 0.0;
            }
            
            // 3. Continuidad de segunda derivada en nodos interiores
            for (int i = 1; i < n-1; i++) {
                int row = 2 * (n-1) + (n-2) + (i-1);
                
                // Segunda derivada del spline i-1 en X[i]
                A[row][4*(i-1)] = 6 * X[i];    // 6a*x
                A[row][4*(i-1)+1] = 2.0;       // 2b
                
                // Segunda derivada del spline i en X[i] (negativa para igualar)
                A[row][4*i] = -6 * X[i];
                A[row][4*i+1] = -2.0;
                
                b[row] = 0.0;
            }
            
            // 4. Condiciones de frontera naturales (segunda derivada = 0 en extremos)
            // En x = X[0] (primer spline)
            int row1 = num_coeffs - 2;
            A[row1][0] = 6 * X[0];
            A[row1][1] = 2.0;
            b[row1] = 0.0;
            
            // En x = X[n-1] (último spline)
            int row2 = num_coeffs - 1;
            A[row2][4*(n-2)] = 6 * X[n-1];
            A[row2][4*(n-2)+1] = 2.0;
            b[row2] = 0.0;
            
            // Resolver el sistema
            printf("\nSolving system of %d equations...\n", num_coeffs);
            gauss_elimination(num_coeffs, A, b, solution);
            
            // Verificar la solución evaluando en los puntos originales
            printf("\nVerification at original points:\n");
            printf("X\t\tY actual\tY spline\tError\n");
            for (int i = 0; i < n; i++) {
                double y_spline = evaluate_spline(X, solution, n, X[i]);
                double error = fabs(Y[i] - y_spline);
                printf("%.4f\t%.6f\t%.6f\t%.6e\n", X[i], Y[i], y_spline, error);
            }
            
            print_cubic_splines(X, solution, n);
            
            // Evaluar puntos adicionales
            printf("Do you want to evaluate a point x_hat? (1 for yes, 0 for no): ");
            scanf("%d", &option);
            while(option) {
                double x_hat;
                printf("Enter the value of x_hat: ");
                scanf("%lf", &x_hat);
                double y_hat = evaluate_spline(X, solution, n, x_hat);
                printf("The interpolated value at x_hat = %.3f is y_hat = %.6f\n", x_hat, y_hat);
                printf("Do you want to evaluate another point x_hat? (1 for yes, 0 for no): ");
                scanf("%d", &option);
            }
            break;
        }
        default:
            printf("Invalid option. Please try again.\n");
            break;
    }


    return 0;
}


int read_data_points(const char* filename, double X[], double Y[], int* n) {
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: Cannot open file '%s'\n", filename);
        return 0;
    }
    
    // Leer número de puntos
    if (fscanf(fp, "%d", n) != 1) {
        printf("Error: Cannot read number of data points\n");
        fclose(fp);
        return 0;
    }
    
    if (*n <= 0 || *n > MAX_POINTS) {
        printf("Error: Invalid number of points (%d)\n", *n);
        fclose(fp);
        return 0;
    }
    
    // Leer puntos
    for (int i = 0; i < *n; i++) {
        if (fscanf(fp, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: Cannot read data point %d\n", i + 1);
            fclose(fp);
            return 0;
        }
    }
    
    fclose(fp);
    printf("Successfully read %d data points\n\n", *n);
    return 1;
}

void print_data_points(double X[], double Y[], int n) {
    printf("Data Points:\n");
    printf("=============\n");
    printf("   i  |      Xi      |      Yi      \n");
    printf("------|--------------|-------------\n");
    for (int i = 0; i < n; i++) {
        printf("%4d  | %12.6f | %12.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
}

void print_cubic_splines(double X[], double solution[], int n) {
    printf("\n------------------CUBIC SPLINE COEFFICIENTS------------------\n");
    for(int k = 0; k < n-1; k++) {
        printf("Spline %d (interval [%.3f, %.3f]):\n", k+1, X[k], X[k+1]);
        printf("  S%d(x) = %.6fx³ + %.6fx² + %.6fx + %.6f\n", 
               k+1, solution[4*k], solution[4*k+1], solution[4*k+2], solution[4*k+3]);
        printf("  Valid for x in [%.3f, %.3f]\n\n", X[k], X[k+1]);
    }
}

double evaluate_spline(double X[], double solution[], int n, double x) {
    int k;
    
    // Encontrar el intervalo correcto
    if (x <= X[0]) {
        k = 0;
    } else if (x >= X[n-1]) {
        k = n - 2;
    } else {
        for (k = 0; k < n-1; k++) {
            if (x >= X[k] && x <= X[k+1]) {
                break;
            }
        }
    }
    
    // Evaluar polinomio cúbico: a*x³ + b*x² + c*x + d
    double x2 = x * x;
    double x3 = x2 * x;
    double y = solution[4*k] * x3 + 
               solution[4*k+1] * x2 + 
               solution[4*k+2] * x + 
               solution[4*k+3];
    
    return y;
}

/**
 * Function to compute and print the linear splines
 * @param X Array of X values (data points)
 * @param Y Array of Y values (data points)
 * @param n Number of data points
 */
void linear_spline(double X[], double Y[], int n) {
    printf("------------------LINEAR SPLINES------------------\n");
    for (int k = 0; k < n - 1; k++) {
        double mk = (Y[k+1] - Y[k]) / (X[k+1] - X[k]);
        printf("Spline %d (from X[%d]=%.3f to X[%d]=%.3f):\n", k+1, k, X[k], k+1, X[k+1]);
        printf("  f_%d(X) = %.6f + %.6f * (X - %.6f)\n\n", k+1, Y[k], mk, X[k]);
    }
}
