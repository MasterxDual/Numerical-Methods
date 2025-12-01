#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 200

/**
 * Function to save computed values to a text file
 * @param x Array of x values
 * @param fp Array of first derivative values
 * @param n Number of subintervals
 */
void save_in_txt(double X[], double Y[], int n);

/**
 * Function to define function f(x, y)
 * @param x The independent variable
 * @param y The dependent variable
 * @return The value of the function at (x, y)
 */
double f(double x, double y);

/**
 * Function to define the exact solution y(x) for comparison
 * @param x The independent variable
 * @return The exact value of y at x
 */
double y(double x);


/* 
    Heun's method has the following improvements over Euler:
    --> Greater accuracy: Consistently smaller errors ✓
    --> Expected behavior: Error pattern consistent with the second-order method ✓
    --> Improved convergence: Errors decrease more quickly ✓

    Euler Method:
        ✅ Extreme simplicity: Easy to implement and understand
        ✅ Low computational cost: 1 evaluation of f(x,y) per step
        ✅ Good for prototyping: Fast for initial testing
        ✅ Easy error analysis: Predictable behavior

        ❌ Low accuracy: Local truncation error O(h²)
        ❌ Instability: Can easily diverge with large h
        ❌ Requires small h to achieve acceptable accuracy
        ❌ Poor for functions with high curvature

        Ideal Cases:
        --> Simple problems with near-linear behavior
        --> When accuracy is not critical
        --> Teaching and conceptual demonstrations

    Heum Method:
        ✅ Good accuracy: O(h³) error vs. Euler's O(h²)
        ✅ More stable: Less prone to divergence
        ✅ Balanced: Good accuracy-cost balance
        ✅ Easy to implement: Only two evaluations of f(x,y)

        ❌ Costo moderado: 2 evaluaciones vs 1 de Euler
        ❌ No óptimo: No es el más preciso de los métodos de 2° orden
        ❌ Dependencia del predictor: Si el predictor es malo, afecta el corrector

        Ideal cases:
        --> General applications requiring better than Euler accuracy
        --> When balancing accuracy and simplicity
        --> Problems with moderately nonlinear behavior

    Midpoint Method:
        ✅ High accuracy: Generally more accurate than Heun for the same order
        ✅ Better for symmetries: Excellent for problems with symmetric behavior
        ✅ Improved stability: Less sensitive to abrupt changes
        ✅ Clear geometric interpretation: Slope at the midpoint

        ❌ Heun-like cost: 2 evaluations of f(x,y)
        ❌ Can underestimate/overestimate: Depending on concavity
        ❌ Slightly more complex implementation: Midpoint calculation


        Ideal cases:
        --> Problems with smooth and symmetric behavior
        --> When maximum precision is sought with second-order methods
        --> Physical systems with harmonic behavior

*/

int main(int argc, char const *argv[]) {
    // General variables
    double X0, Xf, Y0, h;
    // Variables for the Neumann method
    double Xp, Yp;
    int n;
    double X[MAX_SIZE + 1], Y[MAX_SIZE + 1];

    // Number to select to do convergence factor calculation
    int conv_choice;

    // Error to calculate
    double exact_error, local_trunc_error;
    // Variables for convergence factor calculation
    double X1[MAX_SIZE + 1], X2[MAX_SIZE + 1], X3[MAX_SIZE + 1];
    double Y1[MAX_SIZE + 1], Y2[MAX_SIZE + 1], Y3[MAX_SIZE + 1];
    double Q[MAX_SIZE + 1];

    printf("Insert X0 and Xf:\n");
    scanf("%lf %lf", &X0, &Xf);
    printf("Insert initial data Y0 = Y(X0):\n");
    scanf("%lf", &Y0);

    printf("Do you want to insert number of subintervals (n) or step size (h)?\n");
    printf("1. I want to insert n\n");
    printf("2. I want to insert h\n");
    int choice;
    scanf("%d", &choice);
    if(choice == 1) {
        printf("Insert number of subintervals n (integer):\n");
        scanf("%d", &n);
        // Calculate distance between points
        h = (Xf - X0) / n;
    } else if(choice == 2) {
        printf("Insert step size h:\n");
        scanf("%lf", &h);
        n = (int)((Xf - X0) / h);
        if (n > MAX_SIZE) {
            printf("Error: number of subintervals exceeds maximum size (%d).\n", MAX_SIZE);
            return 1;
        }
    }

    
    // Calculate solution
    X[0] = X0;
    Y[0] = Y0;

    double k1, k2, k3;
    // Heun of order 3
    for(int i = 0; i <= n-1; i++) {
        k1 = f(X[i], Y[i]);
        k2 = f(X[i] + h/3.0, Y[i] + (h/3.0) * k1);
        k3 = f(X[i] + h*(2.0/3.0), Y[i] + h*(2.0/3.0) * k2);
        Y[i+1] = Y[i] + (h * ((k1/4.0) + ((3.0/4.0)*k3)));
        X[i+1] = X[i] + h;
    }
    printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
           "i", "X[i]", "Exact Y", "Heun 3 order Y", "Exact Error");
    printf("-------------------------------------------------------------------------------------------\n");
    
    for(int i = 0; i <= n; i++) {
        exact_error = fabs(y(X[i]) - Y[i]);
        printf("%-10d %-15lf %-15lf %-15lf %-15.2e\n", 
               i, X[i], y(X[i]), Y[i], exact_error);
    }

    
    // Print results
    printf("X[i]\t\tY[i]\n");
    // Print all computed points including the last one
    for(int i = 0; i <= n; i++) {
        printf("%lf\t%.10lf\n", X[i], Y[i]);
    }

    // Save x[i] and Y[i in results.txt]
    save_in_txt(X, Y, n);
    
    // Finally, we print the results.txt file in a graph using Python to visualize the results
    // system("python3 graph_points.py");
    if (system("test -f graph_points.py") == 0) {
        system("python3 graph_points.py");
    } else {
        printf("⚠️  Warning: 'graph_points.py' not found. Skipping graph generation.\n");
    }

    return 0;
}

double f(double x, double y) {
    return -2*y + exp(-x);
}

double y(double x) {
    return exp(-2*x) + exp(-x);
}

void save_in_txt(double X[], double Y[], int n) {
    FILE *archivo = fopen("results.txt", "w");
    if (archivo == NULL) {
        printf("Error: Unable to create file.\n");
        exit(1);
    }

    for (int i = 0; i <= n; i++) {
        fprintf(archivo, "%lf\t%lf\n", X[i], Y[i]);
    }

    fclose(archivo);
}

