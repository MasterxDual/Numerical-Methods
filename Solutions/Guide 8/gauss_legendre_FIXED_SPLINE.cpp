#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 20
#define MAX_SIZE 100

#include "gauss_FIXED.h"

/**
 * Function to define function f(x)
 * @param x The point at which to evaluate the function (in our case, X̂)
 * @return The value of the function at x
 */
double f(double x);

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
 * Function to evaluate the cubic spline at a given point x
 * @param X Array of X values (data points)
 * @param solution Array of spline coefficients
 * @param n Number of data points
 * @param x The x value to evaluate
 * @return The interpolated y value
 */
double evaluate_spline(double X[], double solution[], int n, double x);

/* 
Summary:
    2 points:
        Exact for: Polynomials of degree ≤ 3
        Ideal functions:
        f(x) = ax³ + bx² + cx + d
        Very smooth and monotonic functions
        Simple exponentials in small intervals
        When it fails: Oscillating functions, abrupt changes
        Example: f(x) = x³ + 2x² → Exact result
    3 points:
        Exact for: Polynomials of degree ≤ 5
        Ideal functions:
        Moderate oscillations (like your sin(2x)*e^(-x))
        Functions with 1-2 local extrema
        Exponentials with complex behavior
        Optimal trade-off: Accuracy vs. efficiency
        Example: Your function → Very good result
    4 points:
        Exact for: Polynomials of degree ≤ 7
        Ideal functions:
        Greater curvature and complexity
        Trigonometric functions with medium frequency
        Combinations of exponential and trigonometric functions
        When to use: You need high precision without being excessive
    5 points:
        Exact for: Polynomials of degree ≤ 9
        Ideal functions:
        Highly oscillatory (e.g., sin(10x))
        Multiple local extrema
        Functions with smooth singularities
        Cost-benefit: Begins to be expensive
    6 points:
        Exact for: Polynomials of degree ≤ 11
        Ideal functions:
        Extremely complex
        Many oscillations in the interval
        When 5 points do not converge
        
*/

// Global variables to store spline data for integration
double spline_X[MAX_POINTS];
double spline_solution[MAX_SIZE + 1];
int spline_n = 0;

int main(int argc, char const *argv[]) {
    // Weighting factor
    double c0, c1, c2, c3, c4, c5;
    // Function arguments
    double x0, x1, x2, x3, x4, x5;
    // Number of points
    int number_of_points;
    // Result of the integral
    double integral;

    // Arrays for polynomial coefficients calculation to use for Spline
    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
    // Data points to read from text file
    double X[MAX_POINTS], Y[MAX_POINTS];
    // Data points to calculate the integral
    double new_X[MAX_POINTS], new_Y[MAX_POINTS];
    // Number of points
    int n;
    // Limits of integration
    double a_limit, b_limit;
    // Distance between two consecutive points
    double h;
    // Choice between function or data table
    int choice;

    printf("¿Do you have a function or Do you have a data table?\n");
    printf("1. I have a data table\n");
    printf("2. I have a function\n");
    scanf("%d", &choice);

    if(choice == 1) {
        // Read data points from file
        if (!read_data_points("data.txt", X, Y, &n)) {
            printf("Failed to read data from file. Exiting.\n");
            return 1;
        }
        // Print the data points
        print_data_points(X, Y, n);

        int num_coeffs = 4 * (n - 1);

        for (int i = 0; i < num_coeffs; i++) {
            b[i] = 0.0;
            for (int j = 0; j < num_coeffs; j++) {
                A[i][j] = 0.0;
            }
        }

        // 1. Condiciones de interpolación (2 por intervalo)
        for (int i = 0; i < n-1; i++) {
            int interval = i;

            // Condición en el extremo izquierdo del intervalo
            int row = 2*i;
            A[row][4*interval] = pow(X[i], 3);
            A[row][4*interval+1] = pow(X[i], 2);
            A[row][4*interval+2] = X[i];
            A[row][4*interval+3] = 1.0;
            b[row] = Y[i];

            // Condición en el extremo derecho del intervalo
            row = 2*i + 1;
            A[row][4*interval] = pow(X[i+1], 3);
            A[row][4*interval+1] = pow(X[i+1], 2);
            A[row][4*interval+2] = X[i+1];
            A[row][4*interval+3] = 1.0;
            b[row] = Y[i+1];
        }

        // 2. Continuidad de primera derivada en nodos interiores
        for (int i = 1; i < n-1; i++) {
            int row = 2*(n-1) + (i-1);

            // Derivada del polinomio anterior en X[i]
            A[row][4*(i-1)] = 3 * pow(X[i], 2);    // 3a*x²
            A[row][4*(i-1)+1] = 2 * X[i];          // 2b*x
            A[row][4*(i-1)+2] = 1.0;               // c

            // Derivada del polinomio actual en X[i] (negativo)
            A[row][4*i] = -3 * pow(X[i], 2);
            A[row][4*i+1] = -2 * X[i];
            A[row][4*i+2] = -1.0;

            b[row] = 0.0;
        }

        // 3. Continuidad de segunda derivada en nodos interiores
        for (int i = 1; i < n-1; i++) {
            int row = 2*(n-1) + (n-2) + (i-1);

            // Segunda derivada del polinomio anterior en X[i]
            A[row][4*(i-1)] = 6 * X[i];    // 6a*x
            A[row][4*(i-1)+1] = 2.0;       // 2b

            // Segunda derivada del polinomio actual en X[i] (negativo)
            A[row][4*i] = -6 * X[i];
            A[row][4*i+1] = -2.0;

            b[row] = 0.0;
        }

        // 4. Condiciones de frontera naturales (segunda derivada = 0 en extremos)
        // En x = X[0]
        int row = num_coeffs - 2;
        A[row][0] = 6 * X[0];
        A[row][1] = 2.0;
        b[row] = 0.0;

        // En x = X[n-1]
        row = num_coeffs - 1;
        A[row][4*(n-2)] = 6 * X[n-1];
        A[row][4*(n-2)+1] = 2.0;
        b[row] = 0.0;

        // Resolver el sistema
        gauss_elimination(num_coeffs, A, b, solution);

        // It's not necessary to use equally spaced points for the integral
        // Equally spaced points are only necessary to graph or show the interpolation
        // Divide the interval [X[0], X[n-1]] into n-1 subintervals of equal length
        // h = (X[n - 1] - X[0]) / (n - 1);
        /* double start_x = 1.0;
        double end_x = 2.0;
        h = 0.1;
        int num_points = (int)((end_x - start_x) / h) + 1; */

        
        // Divide the interval [X[0], X[n-1]] into n-1 subintervals of equal length
        /* h = (X[n - 1] - X[0]) / (n - 1);

        // Calculate X[n] and Y[n] equally spaced
        for(int i = 0; i < n; i++) {
            new_X[i] = X[0] + i * h;
            new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
        } */

        // We save the spline coefficients and points for f(x) and then to integrate
        for (int i = 0; i < n; i++) {
            spline_X[i] = X[i];
        }
        for (int i = 0; i < 4*(n-1); i++) {
            spline_solution[i] = solution[i];
        }
        spline_n = n;
    }

    // Agregar verificación en los puntos originales
    printf("\nVerificación del spline en puntos originales:\n");
    printf("X\tY dato\tY spline\tError\n");
    for(int i = 0; i < n; i++) {
        double y_spline = evaluate_spline(X, solution, n, X[i]);
        double error = fabs(Y[i] - y_spline);
        printf("%.1f\t%.4f\t%.4f\t%.6f\n", X[i], Y[i], y_spline, error);
    }

    printf("Insert the limits of integration:\n");
    scanf("%lf %lf", &a_limit, &b_limit);

    printf("Insert the number of points (between 2 and 6)\n");
    scanf("%d", &number_of_points);

    switch(number_of_points) {
        case 2: 
            // Implement 2-point Gauss-Legendre quadrature
            c0 = 1.0;
            c1 = 1.0;
            x0 = -0.577350269;
            x1 = 0.577350269;
            // integral = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)));
            integral = ((b_limit-a_limit)/2) * (
            c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2) + 
            c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)
            );
            break;
        case 3:
            // Implement 3-point Gauss-Legendre quadrature
            c0 = 0.5555556;
            c1 = 0.8888889;
            c2 = 0.5555556;
            x0 = -0.774596669;
            x1 = 0.0;
            x2 = 0.774596669;
            // integral = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)) + (c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2)));
            integral = ((b_limit-a_limit)/2) * (
            c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2) + 
            c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2) + 
            c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2)
            );
            break;
        case 4:
            // Implement 4-point Gauss-Legendre quadrature
            c0 = 0.3478548;
            c1 = 0.6521452;
            c2 = 0.6521452;
            c3 = 0.3478548;
            x0 = -0.861136312;
            x1 = -0.339981044;
            x2 = 0.339981044;
            x3 = 0.861136312;
            // integral = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)) + (c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2)) + (c3 * f((((b_limit-a_limit)*x3) + (b_limit+a_limit))/2)));
            integral = ((b_limit-a_limit)/2) * (
            c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2) + 
            c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2) + 
            c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2) + 
            c3 * f((((b_limit-a_limit)*x3) + (b_limit+a_limit))/2)
            );
            break;
        case 5:
            // Implement 5-point Gauss-Legendre quadrature
            c0 = 0.2369269;
            c1 = 0.4786287;
            c2 = 0.5688889;
            c3 = 0.4786287;
            c4 = 0.2369269;
            x0 = -0.906179846;
            x1 = -0.538469310;
            x2 = 0.0;
            x3 = 0.538469310;
            x4 = 0.906179846;
            // integral = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)) + (c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2)) + (c3 * f((((b_limit-a_limit)*x3) + (b_limit+a_limit))/2)) + (c4 * f((((b_limit-a_limit)*x4) + (b_limit+a_limit))/2)));
            integral = ((b_limit-a_limit)/2) * (
            c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2) + 
            c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2) + 
            c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2) + 
            c3 * f((((b_limit-a_limit)*x3) + (b_limit+a_limit))/2) + 
            c4 * f((((b_limit-a_limit)*x4) + (b_limit+a_limit))/2)
            );
            break;
        case 6:
            // Implement 6-point Gauss-Legendre quadrature
            c0 = 0.1713245;
            c1 = 0.3607616;
            c2 = 0.4679139;
            c3 = 0.4679139;
            c4 = 0.3607616;
            c5 = 0.1713245;
            x0 = -0.932469514;
            x1 = -0.661209386;
            x2 = -0.238619186;
            x3 = 0.238619186;
            x4 = 0.661209386;
            x5 = 0.932469514;
            // integral = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)) + (c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2)) + (c3 * f((((b_limit-a_limit)*x3) + (b_limit+a_limit))/2)) + (c4 * f((((b_limit-a_limit)*x4) + (b_limit+a_limit))/2)) + (c5 * f((((b_limit-a_limit)*x5) + (b_limit+a_limit))/2)));
            integral = ((b_limit-a_limit)/2) * (
            c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2) + 
            c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2) + 
            c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2) + 
            c3 * f((((b_limit-a_limit)*x3) + (b_limit+a_limit))/2) + 
            c4 * f((((b_limit-a_limit)*x4) + (b_limit+a_limit))/2) + 
            c5 * f((((b_limit-a_limit)*x5) + (b_limit+a_limit))/2)
            );
            break;
        default:
            printf("Error: Number of points must be between 2 and 6\n");
            exit(0);
    }

    printf("The integral is: %lf\n", integral);
    return 0;
}


double f(double x) {
    if(spline_n == 0) {
        return exp(2 * sin(x)) * (1 + log(x));
    }
    return evaluate_spline(spline_X, spline_solution, spline_n, x);
}

/* int read_data_points(const char* filename, double X[], double Y[], int* n) {
    FILE *fp;
    
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: Cannot open file '%s'\n", filename);
        return 0;
    }
    
    printf("File '%s' opened successfully\n", filename);
    
    // Read number of data points
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
    
    // Read data points
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
} */

int read_data_points(const char* filename, double X[], double Y[], int* n) {
    FILE *fp;
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: Cannot open file '%s'\n", filename);
        return 0;
    }
    
    *n = 0;
    while (*n < MAX_POINTS && fscanf(fp, "%lf %lf", &X[*n], &Y[*n]) == 2) {
        (*n)++;
    }
    
    fclose(fp);
    return 1;
}

/* void print_data_points(double X[], double Y[], int n) {
    printf("Data Points:\n");
    printf("=============\n");
    printf("   i  |      Xi      |      Yi      \n");
    printf("------|--------------|-------------\n");
    for (int i = 0; i < n; i++) {
        printf("%4d  | %12.6f | %12.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
} */

void print_data_points(double X[], double Y[], int n) {
    printf("Data Points:\n");
    for (int i = 0; i < n; i++) {
        printf("X[%d] = %lf, Y[%d] = %lf\n", i, X[i], i, Y[i]);
    }
    printf("\n");
}


/**
 * Function to evaluate the cubic spline at a given point x
 * @param X Array of X values (data points)
 * @param solution Array of spline coefficients
 * @param n Number of data points
 * @param x The x value to evaluate
 * @return The interpolated y value
 */
/* double evaluate_spline(double X[], double solution[], int n, double x) {
    int k;

    // Find the correct spline interval
    if (x <= X[0]) {
        k = 0;  // Use first spline for x smaller than first point
    } else if (x >= X[n-1]) {
        k = n - 2; // Use last spline for x larger than last point
    } else {
        for (k = 0; k < n-1; k++) {
            if (x >= X[k] && x <= X[k+1]) {
                break;
            }
        }
    }

    // Evaluate S_k(x) = a*x^3 + b*x^2 + c*x + d
    double y = solution[4*k] * pow(x, 3) +
               solution[4*k+1] * pow(x, 2) +
               solution[4*k+2] * x +
               solution[4*k+3];

    return y;
} */

 double evaluate_spline(double X[], double coeffs[], int n, double x) {
    // Encontrar el intervalo correcto
    int interval = -1;
    for (int i = 0; i < n-1; i++) {
        if (x >= X[i] && x <= X[i+1]) {
            interval = i;
            break;
        }
    }
    
    // Si x está fuera del rango, usar el primer o último intervalo
    if (interval == -1) {
        if (x < X[0]) interval = 0;
        else interval = n-2;
    }
    
    // Calcular coeficientes para este intervalo
    double a = coeffs[4*interval];
    double b = coeffs[4*interval+1];
    double c = coeffs[4*interval+2];
    double d = coeffs[4*interval+3];
    
    // Evaluar polinomio cúbico: a*x³ + b*x² + c*x + d
    return a*pow(x, 3) + b*pow(x, 2) + c*x + d;
}

