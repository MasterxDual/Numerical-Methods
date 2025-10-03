#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 20
#define MAX_SIZE 100

#include "gauss.h"

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
 * Function to define function f(x)
 * @param x The point at which to evaluate the function (in our case, X̂)
 * @return The value of the function at x
 */
double f(double x);

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
 * Second derivative using aproximation
 * @param func Pointer to the function
 * @param x The point at which to evaluate the second derivative
 * @param h A small value for the finite difference approximation (default is 1e-5)
 * @return The second derivative of the function at point x
 *  */ 
double second_derivative(double (*func)(double), double x, double h = 1e-5);

/** 
 * Fourth derivative using aproximation
 * @param func Pointer to the function
 * @param x The point at which to evaluate the fourth derivative
 * @param h A small value for the finite difference approximation (default is 1e-5)
 * @return The fourth derivative of the function at point x
 *  */
double fourth_derivative(double (*func)(double), double x, double h = 0.01);

/* 
In this exercise we saw that using composition simpson with cubic spline interpolation gave us an error because the determinant was zero.
So, we decided to use linear interpolation which is more stable than cubic splines.
We recommend to use linear interpolation when:
--> Experimental data (like yours)
--> Very different scales (X: 0.001, Y: 800+)
--> Few points (< 15)
--> Large jumps in the data
--> You need numerical stability
--> Robust numerical integration
We recommend to use the original method with cubic spline when:
--> Perfect theoretical/mathematical data
--> Similar scales (X: 1-10, Y: 1-10)
--> Many points (15+) well distributed
--> Smooth data without noise
--> You need continuous derivatives
--> Visualization/animation

We used Composition Trapeze because we dont need to use equally spaced points.
So we eliminate code that equally spaces the points and we add the change of velocity calculation.
The result of V = 44.5127 m/s with equally spaced points using 12 equally spaced points
The result of V = 44.5421 m/s without equally spaced points 
Note: The implementation of equally spaced points is commented out

*/

int main(int argc, char const *argv[]) {
    int choice, subintervals;
    // Limits of integration
    double a, b;
    // Integral of Composition Trapeze
    double sum = 0.0;
    // Values to divide the interval of integration [a,b]
    double x = 0.0;
    // Distance between two consecutive points
    double h = 0.0;

    printf("Choose an option (1.Composition Trapeze 2.Simple Trapeze 3.Composition Simpson 4.Simpson 1/3)\n");
    scanf("%d", &choice);
    if (choice == 1) {
        printf("Do you have a function or Do you have a data table?\n");
        printf("1. I have a function\n");
        printf("2. I have a data table\n");
        scanf("%d", &choice);
        switch (choice) {
            case 1:
                printf("Insert the limits of integration:\n");
                scanf("%lf %lf", &a, &b);
                printf("Please enter the number of subintervals:\n");
                scanf("%d", &subintervals);

                // Calculate I
                sum = f(a) + f(b);
                h = (b - a) / subintervals;        
                for(int i = 1; i <= subintervals-1; i++) {
                    x = a + i * h;
                    sum += 2 * f(x);
                }
                sum = ((b-a)/(2*subintervals)) * sum;
                // Print integral using Trapeze
                printf("The integral is: %lf\n", sum);
                break;

            case 2: {
                // Data points to read from text file
                double X[MAX_POINTS], Y[MAX_POINTS];
                // Data points to calculate the integral using linear interpolation
                double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                // Number of points
                int n;
                
                // Read data points from file
                if (!read_data_points("data.txt", X, Y, &n)) {
                    printf("Failed to read data from file. Exiting.\n");
                    return 1;
                }
                // Print the data points
                print_data_points(X, Y, n);

                printf("Using LINEAR INTERPOLATION (more stable than cubic splines)\n");
                
                // Test with different number of points
                printf("Original data points: %d\n", n);
                printf("Choose number of interpolated points (recommended: 10-50): ");
                int num_points;
                scanf("%d", &num_points);
                
                // Generate equally spaced points using linear interpolation
                /* h = (X[n-1] - X[0]) / (num_points - 1);
                
                for(int i = 0; i < num_points; i++) {
                    new_X[i] = X[0] + i * h;
                    
                    // Find the correct interval for linear interpolation
                    int k = 0;
                    while(k < n-1 && new_X[i] > X[k+1]) {
                        k++;
                    }
                    
                    // Linear interpolation: y = y0 + (y1-y0)*(x-x0)/(x1-x0)
                    if(k < n-1) {
                        new_Y[i] = Y[k] + (Y[k+1] - Y[k]) * (new_X[i] - X[k]) / (X[k+1] - X[k]);
                    } else {
                        new_Y[i] = Y[n-1]; // Last point
                    }
                }

                // Apply trapezoidal rule to interpolated points
                sum = new_Y[0] + new_Y[num_points - 1];
                for(int i = 1; i < num_points - 1; i++) {
                    sum += 2 * new_Y[i];
                }
                sum = ((new_X[num_points - 1] - new_X[0]) / (2 * (num_points - 1))) * sum; */
                // Aplicar regla del trapecio directamente a los puntos desiguales
                double sum = 0.0;
                for(int i = 0; i < n-1; i++) {
                    sum += 0.5 * (Y[i] + Y[i+1]) * (X[i+1] - X[i]);
                }

                double V = 0.0;

                // In this part of the code we modify it like this:
                V= sqrt((2 * sum) / 0.075); // We consider m = 0.075 kg

                // Print integral using Trapeze with linear interpolation
                printf("The integral using linear interpolation is: %lf\n", sum);
                printf("The velocity is: %lf m/s\n", V);
                printf("Note: Used %d equally spaced points for integration\n", num_points);
                break;
            }
            default:
                printf("Invalid choice\n");
                break;
        }
    } else if(choice == 2) {
        // Point to calculate the error
        double c;
        // Integral of Simple Trapeze
        double Iaprox = 0.0;
        // Error of Simple Trapeze
        double aprox_error = 0.0;
        // Exactly error and porcentual error
        double exact_error = 0.0, porcentual_error = 0.0;
        // Exact Integral (needs to be calculated manually)
        double Iexact = 0.0;

        printf("Insert the limits of integration:\n");
        scanf("%lf %lf", &a, &b);
        printf("Insert a value between the interval [a,b] to calculate the error:\n");
        scanf("%lf", &c);
        printf("Insert the exact value of the integral to calculate the exact error:\n");
        scanf("%lf", &Iexact);

        // Calculate I
        Iaprox = (b-a) * ((f(b) + f(a)) / 2.0);
        aprox_error = fabs(-(1.0/12.0) * second_derivative(f, c) * pow((b - a), 3));
        exact_error = fabs(Iexact - Iaprox);
        porcentual_error = (fabs(Iexact - Iaprox) / fabs(Iexact)) * 100.0;

        // Print results
        printf("The aproximated integral is: %lf\n", Iaprox);
        printf("The aproximated error is: %lf\n", aprox_error);
        printf("The exact error is: %lf\n", exact_error);
        printf("The porcentual error is: %lf%%\n", porcentual_error);

    } else if(choice == 3) {
        printf("Do you have a function or Do you have a data table?\n");
        printf("1. I have a function\n");
        printf("2. I have a data table\n");
        scanf("%d", &choice);
        switch(choice) {
            case 1:
                // Values to divide the interval of integration [a,b]
                double x[MAX_POINTS];
                printf("Insert the limits of integration:\n");
                scanf("%lf %lf", &a, &b);
                printf("Please enter the number of subintervals (it needs to be an even number):\n");
                scanf("%d", &subintervals);
                if(subintervals % 2 != 0) {
                    printf("I said that it needs to be an even number, not an odd number");
                    exit(0);
                }
                
                // Calculate I
                sum = f(a) + f(b);
                h = (b - a) / subintervals;        
                for(int i = 1; i < subintervals; i++) {
                    x[i] = a + (i * h);
                    if (i % 2 == 0) {
                        sum += 2 * f(x[i]);  // índices pares
                    } else {
                        sum += 4 * f(x[i]);  // índices impares
                    }
                }
                sum = (h/3) * sum;

                // Print integral using Composition Simpson
                printf("The integral is: %lf\n", sum);
                break;
            case 2:
                // Data points to read from text file
                double X[MAX_POINTS], Y[MAX_POINTS];
                // Data points to calculate the integral using linear interpolation
                double new_X[MAX_POINTS], new_Y[MAX_POINTS];
                // Number of points
                int n;

                // Read data points from file
                if (!read_data_points("data.txt", X, Y, &n)) {
                    printf("Failed to read data from file. Exiting.\n");
                    return 1;
                }
                // Print the data points
                print_data_points(X, Y, n);

                printf("Using LINEAR INTERPOLATION (more stable than cubic splines)\n");
                
                // For Simpson, we need an even number of subintervals
                printf("Original data points: %d\n", n);
                printf("Choose number of subintervals (must be EVEN, recommended: 10-50): ");
                scanf("%d", &subintervals);
                if(subintervals % 2 != 0) {
                    printf("Making it even by adding 1...\n");
                    subintervals++;
                }
                int num_points = subintervals + 1; // Number of points = subintervals + 1
                
                // Generate equally spaced points using linear interpolation
                h = (X[n-1] - X[0]) / subintervals;
                
                for(int i = 0; i <= subintervals; i++) {
                    new_X[i] = X[0] + i * h;
                    
                    // Find the correct interval for linear interpolation
                    int k = 0;
                    while(k < n-1 && new_X[i] > X[k+1]) {
                        k++;
                    }
                    
                    // Linear interpolation: y = y0 + (y1-y0)*(x-x0)/(x1-x0)
                    if(k < n-1) {
                        new_Y[i] = Y[k] + (Y[k+1] - Y[k]) * (new_X[i] - X[k]) / (X[k+1] - X[k]);
                    } else {
                        new_Y[i] = Y[n-1]; // Last point
                    }
                }

                // Apply Simpson's rule to interpolated points
                sum = new_Y[0] + new_Y[subintervals];
                for(int i = 1; i < subintervals; i++) {
                    if(i % 2 == 0) {
                        sum += 2 * new_Y[i];  // Even indices
                    } else {
                        sum += 4 * new_Y[i];  // Odd indices
                    }
                }
                sum = (h/3.0) * sum;

                // Print integral using Simpson with linear interpolation
                printf("The integral using linear interpolation is: %lf\n", sum);
                printf("Note: Used %d equally spaced points (%d subintervals) for Simpson integration\n", 
                       num_points, subintervals);
                break;
        }
    } else if(choice == 4) {
        // This method works perfectly with grade 3 polynomials or less
        // Point to calculate the error
        double c;
        // Integral of Simpson 1/3
        double Iaprox = 0.0;
        // Error of Simpson 1/3
        double aprox_error = 0.0;
        // Exactly error and porcentual error
        double exact_error = 0.0, porcentual_error = 0.0;
        // Exact Integral (needs to be calculated manually)
        double Iexact = 0.0;

        printf("Insert the limits of integration:\n");
        scanf("%lf %lf", &a, &b);
        printf("Insert a value between the interval [a,b] to calculate the error:\n");
        scanf("%lf", &c);
        printf("Insert the exact value of the integral to calculate the exact error:\n");
        scanf("%lf", &Iexact);

        // Calculate I
        Iaprox = ((b-a)/6.0) * (f(a) + 4.0*f((a+b)/2.0) + f(b));
        aprox_error = fabs(-(1.0/2880.0) * pow(b-a, 5) * fourth_derivative(f, c));
        exact_error = fabs(Iexact - Iaprox);
        porcentual_error = (fabs(Iexact - Iaprox) / fabs(Iexact)) * 100.0;

        // Print results
        printf("The aproximated integral is: %lf\n", Iaprox);
        printf("The aproximated error is: %lf\n", aprox_error);
        printf("The exact error is: %lf\n", exact_error);
        printf("The porcentual error is: %lf%%\n", porcentual_error);
    } else {
        printf("You should insert a number between the interval [1, 4]");
    }
    
    return 0;
}

double f(double x) {
    return (sin(2*x) * exp(-x));
}

double second_derivative(double (*func)(double), double x, double h) {
    return (func(x + h) - 2 * func(x) + func(x - h)) / (h * h);
}

double fourth_derivative(double (*func)(double), double x, double h) {
    return (func(x - 2*h) - 4*func(x - h) + 6*func(x) - 4*func(x + h) + func(x + 2*h)) / (pow(h, 4));
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