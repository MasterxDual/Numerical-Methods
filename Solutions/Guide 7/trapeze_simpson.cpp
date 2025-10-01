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

int main(int argc, char const *argv[]) {
    int choice, subintervals;
    // Limits of integration
    double a, b;
    // Integral of Trapeze
    double sum = 0.0;
    // Values to divide the interval of integration [a,b]
    double x = 0.0;
    // Distance between two consecutive points
    double h = 0.0;

    

    printf("Trapeze or Simpson 1/3?(1.Trapeze 2.Simpson 1/3 3.Simpson Composition)\n");
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
                // Arrays for polynomial coefficients calculation to use for Spline
                double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
                // Data points to read from text file
                double X[MAX_POINTS], Y[MAX_POINTS];
                // Data points to calculate the integral
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

                // Divide the interval [X[0], X[n-1]] into n-1 subintervals of equal length
                h = (X[n - 1] - X[0]) / (n - 1);

                // Calculate X[n] and Y[n] equally spaced
                for(int i = 0; i < n; i++) {
                    new_X[i] = X[0] + i * h;
                    new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
                }

                // Calculate I
                sum = new_Y[0] + new_Y[n - 1];
                for(int i = 1; i <= n - 2; i++) {
                    sum += 2 * new_Y[i];
                }
                sum = ((new_X[n - 1] - new_X[0]) / (2 * (n - 1))) * sum;

                // Print integral using Trapeze with made equally spaced points
                printf("The integral is: %lf\n", sum);
                break;
            }
            default:
                printf("Invalid choice\n");
                break;
        }
    } else if(choice == 2) {
        
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
                for(int i = 1; i <= subintervals-3; i+=2) {
                    x[i] = a + (i * h);
                    sum += (4 * f(x[i])) + (2 * f(x[i]+h));
                }
                sum = (h/3) * (sum + 4*f(a + ((subintervals-1)*h)));

                // Print integral using Composition Simpson
                printf("The integral is: %lf\n", sum);
                break;
            case 2:
                // Arrays for polynomial coefficients calculation to use for Spline
                double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
                // Data points to read from text file
                double X[MAX_POINTS], Y[MAX_POINTS];
                // Data points to calculate the integral
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

                // 2. Inicializar A y b
                for (int i = 0; i < 4*(n-1); i++) {
                    b[i] = 0.0;
                    for (int j = 0; j < 4*(n-1); j++) {
                        A[i][j] = 0.0;
                    }
                }

                // If number of points minus 1 is an even
                if((n-1) % 2 == 0) {
                    subintervals = n-1;
                } else {
                    // If number of points minus 1 is an odd
                    subintervals = n;
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


                // Divide the interval [X[0], X[n-1]] into n-1 subintervals of equal length
                h = (X[n - 1] - X[0]) / subintervals;

                // Calculate X[n] and Y[n] equally spaced
                for(int i = 0; i <= subintervals; i++) {
                    new_X[i] = X[0] + i * h;
                    new_Y[i] = evaluate_spline(X, solution, n, new_X[i]);
                }

                // Calculate I applying Simpson rule
                sum = new_Y[0] + new_Y[subintervals];
                for(int i = 1; i < subintervals; i++) {
                    if(i % 2 == 0) {
                        sum += 2 * new_Y[i];
                    } else {
                        sum += 4 * new_Y[i];
                    }
                }
                sum = (h/3) * sum;

                // Print integral using Composition Simpson with made equally spaced points
                printf("The integral is: %lf\n", sum);
                break;
        }
    } else {
        printf("You should insert a number between the interval [1, 3]");
    }
    
    return 0;
}

double f(double x) {
    return x * x;
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

