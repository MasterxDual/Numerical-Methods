#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 20
#define MAX_SIZE 100

#include "gauss.h"

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

/* We know the formula of height of center of pressure like:
    h_cp = (∫(from 0 to H) h * p(h) dh) / (∫(from 0 to H) p(h) dh)
    where p(h) = pressure in function of height h
    H = 112 meters. Total height
 */

// Global variables to store spline data for integration
double spline_X[MAX_POINTS];
double spline_solution[MAX_SIZE + 1];
int spline_n = 0;
int mode = 0; // 0 -> p(h), 1 -> h*p(h)

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
        if (!read_data_points("data1.txt", X, Y, &n)) {
            printf("Failed to read data from file. Exiting.\n");
            return 1;
        }
        // Print the data points
        print_data_points(X, Y, n);

        // 2. Inicializar A y b
        for (int i = 0; i < 4*(n-1); i++) {
            b[i] = 0.0;
            for (int j = 0; j < 4*(n-1); j++) {
                A[i][j] = 0.0;
            }
        }

        // We calculate A[4(n-1)][4(n-1)] and b[4(n-1)]
        // First 2(n-1) equations: Each spline passes through its two endpoints
        for(int k = 0; k < n-1; k++) {
            // Spline k passes through point (X[k], Y[k])
            for(int j = 0; j <= 3; j++) {
                A[2*k][4*k+j] = pow(X[k], 3-j);
            }
            b[2*k] = Y[k];

            // Spline k passes through point (X[k+1], Y[k+1])
            for(int j = 0; j <= 3; j++) {
                A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
            }
            b[2*k+1] = Y[k+1];
        }


        // Next n-2 equations: Continuity of first derivatives
        for(int k = 0; k < n-2; k++) {
            int row = 2*(n-1) + k;
            // First derivative of spline k at X[k+1]
            for(int j = 0; j <= 2; j++) {
                A[row][4*k+j] = (3-j) * pow(X[k+1], 2-j);
            }
            // First derivative of spline k+1 at X[k+1]
            for(int j = 0; j <= 2; j++) {
                A[row][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
            }
            b[row] = 0.0;
        }

        // Next n-2 equations: Continuity of second derivatives
        for(int k = 0; k < n-2; k++) {
            int row = 2*(n-1) + (n-2) + k;
            // Second derivative of spline k at X[k+1]
            A[row][4*k] = 6 * X[k+1];
            A[row][4*k+1] = 2;
            // Second derivative of spline k+1 at X[k+1]
            A[row][4*(k+1)] = -6 * X[k+1];
            A[row][4*(k+1)+1] = -2;
            b[row] = 0.0;
        }

        // Two boundary conditions: Natural spline (second derivatives = 0 at endpoints)
        int row1 = 4*(n-1) - 2;
        int row2 = 4*(n-1) - 1;

        // Second derivative = 0 at X[0] (first spline)
        A[row1][0] = 6 * X[0];
        A[row1][1] = 2;
        b[row1] = 0.0;

        // Second derivative = 0 at X[n-1] (last spline)
        A[row2][4*(n-2)] = 6 * X[n-1];
        A[row2][4*(n-2)+1] = 2;
        b[row2] = 0.0;

        // We use the function from gauss.h to solve the system with Gaussian elimination
        gauss_elimination(4*(n-1), A, b, solution);

        // It's not necessary to use equally spaced points for the integral
        // Equally spaced points are only necessary to graph or show the interpolation
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
        
        // DEBUG: Test spline interpolation at some points
        printf("\n=== SPLINE VERIFICATION ===\n");
        for(int i = 0; i < n; i++) {
            double test_p = evaluate_spline(spline_X, spline_solution, spline_n, X[i]);
            printf("h=%.1f: Original p=%.1f, Spline p=%.1f, Error=%.6f\n", 
                   X[i], Y[i], test_p, fabs(Y[i] - test_p));
        }
        printf("===========================\n\n");
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
            integral = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)));
            break;
        case 3:
            // Implement 3-point Gauss-Legendre quadrature
            c0 = 0.5555556;
            c1 = 0.8888889;
            c2 = 0.5555556;
            x0 = -0.774596669;
            x1 = 0.0;
            x2 = 0.774596669;
            integral = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)) + (c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2)));
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
            integral = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)) + (c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2)) + (c3 * f((((b_limit-a_limit)*x3) + (b_limit+a_limit))/2)));
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
            integral = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)) + (c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2)) + (c3 * f((((b_limit-a_limit)*x3) + (b_limit+a_limit))/2)) + (c4 * f((((b_limit-a_limit)*x4) + (b_limit+a_limit))/2)));
            break;
        case 6: {
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

            // Denominator: ∫ p(h) dh
            mode = 0;
            double integral_p = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)) + (c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2)) + (c3 * f((((b_limit-a_limit)*x3) + (b_limit+a_limit))/2)) + (c4 * f((((b_limit-a_limit)*x4) + (b_limit+a_limit))/2)) + (c5 * f((((b_limit-a_limit)*x5) + (b_limit+a_limit))/2)));/* Llamás al bloque Gauss-Legendre que ya tenés */
                
            // Numerator: ∫ h·p(h) dh
            mode = 1;
            double integral_hp = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)) + (c2 * f((((b_limit-a_limit)*x2) + (b_limit+a_limit))/2)) + (c3 * f((((b_limit-a_limit)*x3) + (b_limit+a_limit))/2)) + (c4 * f((((b_limit-a_limit)*x4) + (b_limit+a_limit))/2)) + (c5 * f((((b_limit-a_limit)*x5) + (b_limit+a_limit))/2)));/* Volvés a llamar al bloque Gauss-Legendre */
                
            // Debug: Print intermediate results
            printf("DEBUG: ∫ p(h) dh = %.6lf\n", integral_p);
            printf("DEBUG: ∫ h·p(h) dh = %.6lf\n", integral_hp);
            printf("DEBUG: Integration limits: [%.1lf, %.1lf]\n", a_limit, b_limit);
            
            // Center of pressure
            double hc = integral_hp / integral_p;
            printf("Center of pressure: %.6lf m\n", hc);
            break;
        }
        default:
            printf("Error: Number of points must be between 2 and 6\n");
            exit(0);
    }

    // printf("The integral is: %lf\n", integral);
    return 0;
}

double f(double x) {
    double p = evaluate_spline(spline_X, spline_solution, spline_n, x); // interpolated pressure

    if (mode == 0) { // 0 -> p(h), 1 -> h*p(h)
        return p;        // For denominator ∫ p(h) dh
    } else {
        return x * p;    // For numerator ∫ h·p(h) dh
    }
}

int read_data_points(const char* filename, double X[], double Y[], int* n) {
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

/**
 * Function to evaluate the cubic spline at a given point x
 * @param X Array of X values (data points)
 * @param solution Array of spline coefficients
 * @param n Number of data points
 * @param x The x value to evaluate
 * @return The interpolated y value
 */
double evaluate_spline(double X[], double solution[], int n, double x) {
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
}