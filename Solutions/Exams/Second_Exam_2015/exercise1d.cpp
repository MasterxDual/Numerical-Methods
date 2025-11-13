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

double interpolate_linear(double X[], double Y[], int n, double x);

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

        
    }

    printf("Insert the limits of integration:\n");
    scanf("%lf %lf", &a_limit, &b_limit);

    printf("Insert the number of points (between 2 and 6)\n");
    scanf("%d", &number_of_points);

    switch(number_of_points) {
        case 2: 
        {
            // Implement 2-point Gauss-Legendre quadrature
            // Usar esto en caso que tengamos la función
            /* c0 = 1.0;
            c1 = 1.0;
            x0 = -0.577350269;
            x1 = 0.577350269;
            // integral = (((b_limit-a_limit)/2) * (c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2)) + (c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)));
            integral = ((b_limit-a_limit)/2) * (
            c0 * f((((b_limit-a_limit)*x0) + (b_limit+a_limit))/2) + 
            c1 * f((((b_limit-a_limit)*x1) + (b_limit+a_limit))/2)
            ); */

            c0 = 1.0;
            c1 = 1.0;
            x0 = -0.577350269;
            x1 = 0.577350269;
            
            // Transformar nodos de [-1,1] a [a_limit,b_limit]
            double t0_2 = ((b_limit - a_limit) * x0 + (b_limit + a_limit)) / 2.0;
            double t1_2 = ((b_limit - a_limit) * x1 + (b_limit + a_limit)) / 2.0;
            
            // Usar interpolación lineal en lugar de spline
            double f0_2 = interpolate_linear(X, Y, n, t0_2);
            double f1_2 = interpolate_linear(X, Y, n, t1_2);
            
            integral = ((b_limit - a_limit) / 2.0) * (c0 * f0_2 + c1 * f1_2);
            break;
        }
        case 3:
        {
            // Implement 3-point Gauss-Legendre quadrature
            /* c0 = 0.5555556;
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
            ); */
            c0 = 0.5555556;
            c1 = 0.8888889;
            c2 = 0.5555556;
            x0 = -0.774596669;
            x1 = 0.0;
            x2 = 0.774596669;
            
            double t0_3 = ((b_limit - a_limit) * x0 + (b_limit + a_limit)) / 2.0;
            double t1_3 = ((b_limit - a_limit) * x1 + (b_limit + a_limit)) / 2.0;
            double t2_3 = ((b_limit - a_limit) * x2 + (b_limit + a_limit)) / 2.0;
            
            double f0_3 = interpolate_linear(X, Y, n, t0_3);
            double f1_3 = interpolate_linear(X, Y, n, t1_3);
            double f2_3 = interpolate_linear(X, Y, n, t2_3);
            
            integral = ((b_limit - a_limit) / 2.0) * (c0 * f0_3 + c1 * f1_3 + c2 * f2_3);
            break;
        }
        case 4:
        {
            // Implement 4-point Gauss-Legendre quadrature
            /* c0 = 0.3478548;
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
            ); */
            c0 = 0.3478548;
            c1 = 0.6521452;
            c2 = 0.6521452;
            c3 = 0.3478548;
            x0 = -0.861136312;
            x1 = -0.339981044;
            x2 = 0.339981044;
            x3 = 0.861136312;

            double t0_4 = ((b_limit - a_limit) * x0 + (b_limit + a_limit)) / 2.0;
            double t1_4 = ((b_limit - a_limit) * x1 + (b_limit + a_limit)) / 2.0;
            double t2_4 = ((b_limit - a_limit) * x2 + (b_limit + a_limit)) / 2.0;
            double t3_4 = ((b_limit - a_limit) * x3 + (b_limit + a_limit)) / 2.0;

            double f0_4 = interpolate_linear(X, Y, n, t0_4);
            double f1_4 = interpolate_linear(X, Y, n, t1_4);
            double f2_4 = interpolate_linear(X, Y, n, t2_4);
            double f3_4 = interpolate_linear(X, Y, n, t3_4);

            integral = ((b_limit - a_limit) / 2.0) * (c0 * f0_4 + c1 * f1_4 + c2 * f2_4 + c3 * f3_4);
            break;
        }
        case 5:
        {
            // Implement 5-point Gauss-Legendre quadrature
            /* c0 = 0.2369269;
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
            ); */
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

            double t0_5 = ((b_limit - a_limit) * x0 + (b_limit + a_limit)) / 2.0;
            double t1_5 = ((b_limit - a_limit) * x1 + (b_limit + a_limit)) / 2.0;
            double t2_5 = ((b_limit - a_limit) * x2 + (b_limit + a_limit)) / 2.0;
            double t3_5 = ((b_limit - a_limit) * x3 + (b_limit + a_limit)) / 2.0;
            double t4_5 = ((b_limit - a_limit) * x4 + (b_limit + a_limit)) / 2.0;

            double f0_5 = interpolate_linear(X, Y, n, t0_5);
            double f1_5 = interpolate_linear(X, Y, n, t1_5);
            double f2_5 = interpolate_linear(X, Y, n, t2_5);
            double f3_5 = interpolate_linear(X, Y, n, t3_5);
            double f4_5 = interpolate_linear(X, Y, n, t4_5);

            integral = ((b_limit - a_limit) / 2.0) * (c0 * f0_5 + c1 * f1_5 + c2 * f2_5 + c3 * f3_5 + c4 * f4_5);
            break;
        }
        case 6:
        {
            // Implement 6-point Gauss-Legendre quadrature
            /* c0 = 0.1713245;
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
            ); */
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

            double t0_6 = ((b_limit - a_limit) * x0 + (b_limit + a_limit)) / 2.0;
            double t1_6 = ((b_limit - a_limit) * x1 + (b_limit + a_limit)) / 2.0;
            double t2_6 = ((b_limit - a_limit) * x2 + (b_limit + a_limit)) / 2.0;
            double t3_6 = ((b_limit - a_limit) * x3 + (b_limit + a_limit)) / 2.0;
            double t4_6 = ((b_limit - a_limit) * x4 + (b_limit + a_limit)) / 2.0;
            double t5_6 = ((b_limit - a_limit) * x5 + (b_limit + a_limit)) / 2.0;

            double f0_6 = interpolate_linear(X, Y, n, t0_6);
            double f1_6 = interpolate_linear(X, Y, n, t1_6);
            double f2_6 = interpolate_linear(X, Y, n, t2_6);
            double f3_6 = interpolate_linear(X, Y, n, t3_6);
            double f4_6 = interpolate_linear(X, Y, n, t4_6);
            double f5_6 = interpolate_linear(X, Y, n, t5_6);

            integral = ((b_limit - a_limit) / 2.0) * (c0 * f0_6 + c1 * f1_6 + c2 * f2_6 + c3 * f3_6 + c4 * f4_6 + c5 * f5_6);
            break;
        }
        default:
            printf("Error: Number of points must be between 2 and 6\n");
            exit(0);
    }

    printf("The integral is: %lf\n", integral);
    return 0;
}


double f(double x) {
    return (((2.0*x) / (pow(x, 2) + 1.0)) - (cos(x)));
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

double interpolate_linear(double X[], double Y[], int n, double x) {
    if (x <= X[0]) return Y[0];
    if (x >= X[n-1]) return Y[n-1];
    
    for (int i = 0; i < n-1; i++) {
        if (x >= X[i] && x <= X[i+1]) {
            double t = (x - X[i]) / (X[i+1] - X[i]);
            return Y[i] + t * (Y[i+1] - Y[i]);
        }
    }
    return Y[n-1];
}