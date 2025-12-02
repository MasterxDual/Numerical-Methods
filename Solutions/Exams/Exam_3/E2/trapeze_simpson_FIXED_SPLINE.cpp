#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 20
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

                // Divide the interval [X[0], X[n-1]] into n-1 subintervals of equal length
                // h = (X[n - 1] - X[0]) / (n - 1);
                double start_x = 1.0;
                double end_x = 2.0;
                h = 0.1;
                int num_points = (int)((end_x - start_x) / h) + 1;


                printf("\nPuntos equiespaciados con h = 0.1:\n");
                printf("x\t\tf_spline(x)\t\tf_exacto(x)\t\tError\n");
                printf("------------------------------------------------------------\n");

                
                for (int i = 0; i < num_points; i++) {
                    double x = start_x + i * h;
                    double y_spline = evaluate_spline(X, solution, n, x);

                    // Calcular valor exacto: f(x) = e^(2*sin(x)) * (1 + ln(x))
                    double y_exact = exp(2 * sin(x)) * (1 + log(x));
                    double error = fabs(y_exact - y_spline);

                    printf("%.2lf\t\t%.6lf\t\t%.6lf\t\t%.6lf\n", 
                           x, y_spline, y_exact, error);

                    new_Y[i] = y_spline;
                    new_X[i] = x;
                }
                
                // Calculate I
                sum = new_Y[0] + new_Y[num_points - 1];
                for(int i = 1; i <= num_points - 2; i++) {
                    sum += 2 * new_Y[i];
                }
                sum = ((new_X[num_points - 1] - new_X[0]) / (2 * (num_points - 1))) * sum;

                // Print integral using Trapeze with made equally spaced points
                printf("The integral is: %lf\n", sum);
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

                int num_coeffs = 4 * (n - 1);

                for (int i = 0; i < num_coeffs; i++) {
                    b[i] = 0.0;
                    for (int j = 0; j < num_coeffs; j++) {
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

                // Divide the interval [X[0], X[n-1]] into n-1 subintervals of equal length
                // h = (X[n - 1] - X[0]) / (n - 1);
                double start_x = 1.0;
                double end_x = 2.0;
                h = 0.1;
                int num_points = (int)((end_x - start_x) / h) + 1;
                double exact_Y[MAX_POINTS];


                printf("\nPuntos equiespaciados con h = 0.1:\n");
                printf("x\t\tf_spline(x)\t\tf_exacto(x)\t\tError\n");
                printf("------------------------------------------------------------\n");

                for (int i = 0; i < num_points; i++) {
                    double x = start_x + i * h;
                    double y_spline = evaluate_spline(X, solution, n, x);

                    // Calcular valor exacto: f(x) = e^(2*sin(x)) * (1 + ln(x))
                    double y_exact = exp(2 * sin(x)) * (1 + log(x));
                    double error = fabs(y_exact - y_spline);

                    printf("%.2lf\t\t%.6lf\t\t%.6lf\t\t%.6lf\n", 
                           x, y_spline, y_exact, error);

                    new_X[i] = x;
                    new_Y[i] = y_spline;
                    exact_Y[i] = y_exact;

                }    
               
                subintervals = num_points - 1;

                // Calculate I applying Simpson rule
                /* sum = new_Y[0] + new_Y[subintervals];
                for(int i = 1; i < subintervals; i++) {
                    if(i % 2 == 0) {
                        sum += 2 * new_Y[i];
                    } else {
                        sum += 4 * new_Y[i];
                    }
                }
                sum = (h/3) * sum; */
                sum = exact_Y[0] + exact_Y[subintervals];
                for(int i = 1; i < subintervals; i++) {
                    if(i % 2 == 0) {
                        sum += 2 * exact_Y[i];
                    } else {
                        sum += 4 * exact_Y[i];
                    }
                }
                sum = (h/3) * sum;

                // Print integral using Composition Simpson with made equally spaced points
                printf("The integral is: %lf\n", sum);
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
}
 */

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
