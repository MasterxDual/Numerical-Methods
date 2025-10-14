#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 100

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

/**
 * Function to define the total derivative of f with respect to x
 * @param x The independent variable
 * @param y The dependent variable
 * @return The value of the total derivative at (x, y)
 */
double fprima(double x, double y);  // Nueva función: derivada total de f respecto a x

/**
 * Function to define the third derivative of y with respect to x
 * @param x The independent variable
 * @param y The dependent variable
 * @return The value of the third derivative at (x, y)
 */
double y3prima(double x, double y);

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
    // Error to calculate
    double exact_error, local_trunc_error;

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
    }

    
    // Calculate solution
    X[0] = X0;
    Y[0] = Y0;

    printf("Insert the method to use: 1. Euler's method 2. Neum Method 3. Midpoint method\n");
    scanf("%d", &choice);

    switch(choice) {
        case 1:
            for(int i = 1; i <= n; i++) {
                // Compute next X incrementally to avoid rounding surprises
                X[i] = X[i-1] + h; // X[i] = X0 + i*h; equivalently
                // Euler method
                Y[i] = Y[i-1] + h * f(X[i-1], Y[i-1]); 
            }
            for(int i = 0; i <= n; i++) {
                exact_error = fabs(y(X[i]) - Y[i]);
                printf("At X = %lf, Exact Y = %lf, Euler Y = %lf, Exact Error (e%d) = %lf\n", X[i], y(X[i]), Y[i], i, exact_error);
            }
        
            printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                   "i", "X[i]", "Exact Y", "Euler Y", "Exact Error", "Local Trunc. Err");
            printf("-------------------------------------------------------------------------------------------\n");
        
            // Euler's method approximates the curve by means of straight line segments tangent to each point.
            // If these segments lie above the curve → Euler overestimates.
            // If they lie below it → Euler underestimates.
            // If local_trunc_error < 0 => The value calculated by Euler is less than the exact one. Euler underestimates.
            // If local_trunc_error > 0 => The value calculated by Euler is greater than the exact one. Euler overestimates.
            for(int i = 0; i <= n; i++) {
                exact_error = fabs(y(X[i]) - Y[i]);
                local_trunc_error = (h * h / 2.0) * fprima(X[i], y(X[i])); 
                printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n", 
                       i, X[i], y(X[i]), Y[i], exact_error, local_trunc_error);
            }
            break;
        case 2: 
            for(int i = 1; i <= n; i++) {
                X[i] = X0 + (i*h);
                // Neum method
                Xp = X[i] + h;
                Yp = Y[i-1] + h * f(X[i-1], Y[i-1]);
                Y[i] = Y[i-1] + (h/2.0) * (f(X[i-1], Y[i-1]) + f(Xp, Yp));
            }
            printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                    "i", "X[i]", "Exact Y", "Heun Y", "Exact Error", "Local Trunc. Err");
            printf("-------------------------------------------------------------------------------------------\n");

            for(int i = 0; i <= n; i++) {
                exact_error = fabs(y(X[i]) - Y[i]);
                local_trunc_error = (pow(h, 3) / 12.0) * y3prima(X[i], y(X[i]));
                printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n", 
                       i, X[i], y(X[i]), Y[i], exact_error, local_trunc_error);
            }

            break;
        case 3:
            // Midpoint method that gave us the teacher in class
            /* for(int i = 1; i <= n; i++) {
                X[i] = X0 + (i*h/2.0);
                Xp = X[i] + h;
                Yp = Y[i-1] + ((h/2.0) * f(X[i-1], Y[i-1]));
                Y[i] = Y[i-1] + (h * f(Xp, Yp));
            } */

            // Midpoint method that gave me ChatGPT
            for(int i = 1; i <= n; i++) {
                X[i] = X0 + (i * h);
                // Slope (pendiente) at the start of the subinterval
                double k1 = f(X[i-1], Y[i-1]);
                // Slope at midpoint using Euler predictor
                double k2 = f(X[i-1] + h/2.0, Y[i-1] + (h/2.0) * k1);
                // Use the slope at the midpoint to move forward
                Y[i] = Y[i-1] + h * k2;
            }
        
            printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                   "i", "X[i]", "Exact Y", "Midpoint Y", "Exact Error", "Local Trunc. Err");
            printf("-------------------------------------------------------------------------------------------\n");
            
            for(int i = 0; i <= n; i++) {
                exact_error = fabs(y(X[i]) - Y[i]);
                local_trunc_error = (pow(h, 3) / 24.0) * y3prima(X[i], y(X[i]));
                printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n", 
                       i, X[i], y(X[i]), Y[i], exact_error, local_trunc_error);
            }
            break;
    }

    
    // Print results
    printf("X[i]\t\tY[i]\n");
    // Print all computed points including the last one
    for(int i = 0; i <= n; i++) {
        printf("%lf\t%lf\n", X[i], Y[i]);
    }

    // Save x[i] and Y[i in results.txt]
    save_in_txt(X, Y, n);
    
    // Finally, we print the results.txt file in a graph using Python to visualize the results
    system("python3 graph_points.py");

    
    return 0;
}

double f(double x, double y) {
    return (-2 * x * y);
}

double y(double x) {
    return exp(-x*x);
}

double y3prima(double x, double y) {
    return 4 * x * y * (3 - 2 * x * x);
}

double fprima(double x, double y) {
    double fx = -2 * y;
    double fy = -2 * x;
    return fx + fy * f(x, y);
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

