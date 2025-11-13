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
 * Function to save computed exact solution values to a text file
 * @param x Array of x values
 * @param Yfunc Function pointer to exact solution y(x)
 * @param n Number of subintervals
 *  */ 
void save_in_txt2(double X[], double (*Yfunc)(double), int n);

/**
 * Function to define function f(x, y)
 * @param x The independent variable
 * @param y The dependent variable
 * @return The value of the function at (x, y)
 */
double f(double x, double y);


/**
 * Function to calculate convergence factor for Runge-Kutta 4th order method
 *  */ 
void check_and_run_graph();

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

/**
 * Function to perform a single Runge-Kutta 4th order step
 * @param x Current x value
 * @param y Current y value
 * @param step Step size
 * @return The new y value after the RK4 step
 */ 
double rk4_step(double x, double y, double step);

/**
 * Function to calculate convergence factor for Euler's method
 * @param n Number of subintervals
 * @param h Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y0 Initial y value
 */ 
void convergence_factor_euler(int n, double h, double X0, double Xf, double Y0);

/**
 * Function to calculate convergence factor for Heun's method
 * @param n Number of subintervals
 * @param h Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y0 Initial y value
 */
void convergence_factor_heun(int n, double h, double X0, double Xf, double Y0);

/**
 * Function to calculate convergence factor for Midpoint method
 * @param n Number of subintervals
 * @param h Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y0 Initial y value
 */
void convergence_factor_midpoint(int n, double h, double X0, double Xf, double Y0);


/**
 * Function to calculate convergence factor for Runge-Kutta 4th order method
 * @param n Number of subintervals
 * @param h Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y0 Initial y value
 *  */ 
void convergence_factor_rk4(int n, double h, double X0, double Xf, double Y0);


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

    printf("Insert the method to use: 4. Runge-Kutta of order 4 5. Method of Exercise\n");
    scanf("%d", &choice);

    switch(choice) {
        case 4:
            printf("Do you want to calculate convergence factor for Runge Kutta of order 4? (1.Yes 2.No)\n");
            scanf("%d", &conv_choice);

            double k1, k2, k3, k4;
            // Runge-Kutta of order 4
            for(int i = 0; i <= n-1; i++) {
                k1 = f(X[i], Y[i]);
                k2 = f(X[i] + h/2.0, Y[i] + (h/2.0) * k1);
                k3 = f(X[i] + h/2.0, Y[i] + (h/2.0) * k2);
                k4 = f(X[i] + h, Y[i] + h * k3);
                Y[i+1] = Y[i] + (h/6.0) * (k1 + 2*k2 + 2*k3 + k4);
                X[i+1] = X[i] + h;
            }

            printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                   "i", "X[i]", "Exact Y", "RK4 Y", "Exact Error", "Local Trunc. Err");
            printf("-------------------------------------------------------------------------------------------\n");
            
            for(int i = 0; i <= n; i++) {
                exact_error = fabs(y(X[i]) - Y[i]);

                // Here we estimate the local truncation error (LTE) for RK4
                // Practical estimation of LTE (per step) by comparing h vs h/2
                double Y_full = Y[i];  // current value with step h
                double Y_half;
                if (i == 0) {
                    Y_half = Y[0]; // not applicable for the first point
                    local_trunc_error = 0.0;
                } else {
                    // Recalculate from X[i-1], Y[i-1] with two steps of h/2
                    double mid = rk4_step(X[i-1], Y[i-1], h/2.0);
                    Y_half = rk4_step(X[i-1] + h/2.0, mid, h/2.0);

                    // Local truncation error estimator: (Y_half - Y_full) / (2^{p+1}-1) = /31
                    local_trunc_error = fabs((Y_half - Y_full) / 15.0);
                }

                printf("%-10d %-15lf %-15lf %-15lf %-15.2e %-15.2e\n", 
                       i, X[i], y(X[i]), Y[i], exact_error, local_trunc_error);

                // Save x[i] and Y[i] in runge_kutta_4_results.txt]
                save_in_txt(X, Y, n);
                rename("results.txt", "runge_kutta_4_results.txt");
            }
            break;
        case 5:
            X[1] = X0 + h;
            Y[1] = Y[0] + h * f(X[0], Y[0]);

            for(int i = 1; i < n; i++) {
                X[i+1] = X0 + (i+1)*h;
                Y[i+1] = Y[i] + h * (2 * f(X[i], Y[i]) - f(X[i-1], Y[i-1]));
            }
            // Euler's Method
            /* for(int i = 1; i <= n; i++) {
                // Compute next X incrementally to avoid rounding surprises
                X[i] = X[i-1] + h; // X[i] = X0 + i*h; equivalently
                // Euler method
                Y[i] = Y[i-1] + h * f(X[i-1], Y[i-1]); 
            } */
            
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
            // Save x[i] and Y[i] in runge_kutta_4_results.txt]
            save_in_txt(X, Y, n);
            // save_in_txt2(X, y, n);
            rename("results.txt", "method_of_exercise_results.txt");
            // rename("results2.txt", "exact_results.txt");
            break;
    }

    
    // Print results
    printf("X[i]\t\tY[i]\n");
    // Print all computed points including the last one
    for(int i = 0; i <= n; i++) {
        printf("%lf\t%lf\n", X[i], Y[i]);
    }

    // Save x[i] and Y[i in results.txt]
    // save_in_txt(X, Y, n);
    
    // Finally, we print the results.txt file in a graph using Python to visualize the results
    // system("python3 graph_points.py");
    /* if (system("test -f graph_points.py") == 0) {
        system("python3 graph_points.py");
    } else {
        printf("⚠️  Warning: 'graph_points.py' not found. Skipping graph generation.\n");
    } */

    check_and_run_graph();

    return 0;
}

double f(double x, double y) {
    return ((x * exp(x*x)) / y);
}

double y(double x) {
    return sqrt(exp(x*x));
}

/**
 * @brief Verifica si existen los cuatro archivos de resultados
 *        y ejecuta graph_points.py si todos están presentes.
 */
void check_and_run_graph() {
    int all_exist = 1;
    const char* filenames[] = {
        "method_of_exercise_results.txt",
        //"exact_results.txt",
        "runge_kutta_4_results.txt"
    };

    printf("\nVerifying results of files...\n");

    // Verifying existence of each file
    for (int i = 0; i < 2; i++) {
        FILE *file = fopen(filenames[i], "r");
        if (file == NULL) {
            all_exist = 0;
            printf("Missing file: %s\n", filenames[i]);
        } else {
            fclose(file);
            printf("File found: %s\n", filenames[i]);
        }
    }

    // If all files exist, run the Python script
    if (all_exist) {
        printf("\n✅ All files are present. Running combined graph...\n");
        if (system("test -f graph_points_all.py") == 0) {
            system("python3 graph_points_all.py");
        } else {
            printf("'graph_points_all.py' not found. Skipping graph generation.\n");
        }
    } else {
        printf("\nMissing result files. Please run the missing methods before graphing.\n");
    }
}


/**
 * Third derivative of y(x) = e^{-x^2}
 * Computed symbolically as y''' = 4xy(3 - 2x^2)
 */
double y3prima(double x, double y) {
    return 4 * x * y * (3 - 2 * x * x);
}

double fprima(double x, double y) {
    double fx = -2 * y;
    double fy = -2 * x;
    return fx + fy * f(x, y);
}


double rk4_step(double x, double y, double step) {
    double k1 = f(x, y);
    double k2 = f(x + step/2.0, y + (step/2.0) * k1);
    double k3 = f(x + step/2.0, y + (step/2.0) * k2);
    double k4 = f(x + step,     y + step * k3);
    return y + (step/6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
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

// y(X[i])
void save_in_txt2(double X[], double (*Yfunc)(double), int n) {
    FILE *archivo = fopen("results2.txt", "w");
    if (archivo == NULL) {
        printf("Error: Unable to create file.\n");
        exit(1);
    }

    for (int i = 0; i <= n; i++) {
        fprintf(archivo, "%lf\t%lf\n", X[i], Yfunc(X[i]));
    }

    fclose(archivo);
}


void convergence_factor_euler(int n1, double h1, double X0, double Xf, double Y0) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    double Yh[MAX_SIZE + 1], Yh2[MAX_SIZE*2 + 1], Yh4[MAX_SIZE*4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE*2 + 1], Xh4[MAX_SIZE*4 + 1];
    double Q[MAX_SIZE + 1];

    // Inicializaciones
    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Yh[0] = Yh2[0] = Yh4[0] = Y0;

    // Euler con paso h
    for (int i = 0; i < n1; i++) {
        Xh[i+1] = Xh[i] + h1;
        Yh[i+1] = Yh[i] + h1 * f(Xh[i], Yh[i]);
    }

    // Euler con paso h/2
    for (int i = 0; i < 2*n1; i++) {
        Xh2[i+1] = Xh2[i] + h2;
        Yh2[i+1] = Yh2[i] + h2 * f(Xh2[i], Yh2[i]);
    }

    // Euler con paso h/4
    for (int i = 0; i < 4*n1; i++) {
        Xh4[i+1] = Xh4[i] + h3;
        Yh4[i+1] = Yh4[i] + h3 * f(Xh4[i], Yh4[i]);
    }

    printf("\n%-10s %-15s %-15s\n", "i", "x_i", "Q_i");
    printf("------------------------------------------\n");

    // Cálculo del factor de convergencia en los mismos puntos X
    // Note: Q[0] is not defined because at the initial point there is no error
    // (it's the exact initial condition), so we start from i=1
    Q[0] = 0.0;  // Not defined, set to 0 by convention
    
    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;   // posición equivalente para h/2
        int idx4 = 4*i;   // posición equivalente para h/4
        double num = fabs(Yh[i] - Yh2[idx2]);
        double den = fabs(Yh2[idx2] - Yh4[idx4]);

        if (den > 1e-12) {
            Q[i] = log(num / den) / log(2.0);
            printf("%-10d %-15lf %-15lf\n", i, Xh[i], Q[i]);
        } else {
            Q[i] = 0.0;
            printf("%-10d %-15lf %-15s\n", i, Xh[i], "N/A (no error)");
        }
    }

    save_in_txt(Xh, Q, n1);
    rename("results.txt", "convergence_euler.txt");

    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    } else {
        printf("⚠️  Warning: 'graph_convergence.py' not found. Skipping graph generation.\n");
    }
}

void convergence_factor_heun(int n1, double h1, double X0, double Xf, double Y0) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    double Yh[MAX_SIZE + 1], Yh2[MAX_SIZE*2 + 1], Yh4[MAX_SIZE*4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE*2 + 1], Xh4[MAX_SIZE*4 + 1];
    double Q[MAX_SIZE + 1];

    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Yh[0] = Yh2[0] = Yh4[0] = Y0;

    // Heun con paso h
    for (int i = 0; i < n1; i++) {
        double predictor = Yh[i] + h1 * f(Xh[i], Yh[i]);
        Yh[i+1] = Yh[i] + (h1/2.0)*(f(Xh[i], Yh[i]) + f(Xh[i]+h1, predictor));
        Xh[i+1] = Xh[i] + h1;
    }

    // Heun con paso h/2
    for (int i = 0; i < 2*n1; i++) {
        double predictor = Yh2[i] + h2 * f(Xh2[i], Yh2[i]);
        Yh2[i+1] = Yh2[i] + (h2/2.0)*(f(Xh2[i], Yh2[i]) + f(Xh2[i]+h2, predictor));
        Xh2[i+1] = Xh2[i] + h2;
    }

    // Heun con paso h/4
    for (int i = 0; i < 4*n1; i++) {
        double predictor = Yh4[i] + h3 * f(Xh4[i], Yh4[i]);
        Yh4[i+1] = Yh4[i] + (h3/2.0)*(f(Xh4[i], Yh4[i]) + f(Xh4[i]+h3, predictor));
        Xh4[i+1] = Xh4[i] + h3;
    }

    printf("\n%-10s %-15s %-15s\n", "i", "x_i", "Q_i");
    printf("------------------------------------------\n");
    Q[0] = 0.0;

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;
        int idx4 = 4*i;
        double num = fabs(Yh[i] - Yh2[idx2]);
        double den = fabs(Yh2[idx2] - Yh4[idx4]);

        if (den > 1e-12) {
            Q[i] = log(num / den) / log(2.0);
            printf("%-10d %-15lf %-15lf\n", i, Xh[i], Q[i]);
        } else {
            Q[i] = 0.0;
            printf("%-10d %-15lf %-15s\n", i, Xh[i], "N/A");
        }
    }

    save_in_txt(Xh, Q, n1);
    rename("results.txt", "convergence_heun.txt");

    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    }
}

void convergence_factor_midpoint(int n1, double h1, double X0, double Xf, double Y0) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    double Yh[MAX_SIZE + 1], Yh2[MAX_SIZE*2 + 1], Yh4[MAX_SIZE*4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE*2 + 1], Xh4[MAX_SIZE*4 + 1];
    double Q[MAX_SIZE + 1];

    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Yh[0] = Yh2[0] = Yh4[0] = Y0;

    // Midpoint paso h
    for (int i = 0; i < n1; i++) {
        double k1 = f(Xh[i], Yh[i]);
        double k2 = f(Xh[i] + h1/2.0, Yh[i] + (h1/2.0)*k1);
        Yh[i+1] = Yh[i] + h1*k2;
        Xh[i+1] = Xh[i] + h1;
    }

    // Midpoint paso h/2
    for (int i = 0; i < 2*n1; i++) {
        double k1 = f(Xh2[i], Yh2[i]);
        double k2 = f(Xh2[i] + h2/2.0, Yh2[i] + (h2/2.0)*k1);
        Yh2[i+1] = Yh2[i] + h2*k2;
        Xh2[i+1] = Xh2[i] + h2;
    }

    // Midpoint paso h/4
    for (int i = 0; i < 4*n1; i++) {
        double k1 = f(Xh4[i], Yh4[i]);
        double k2 = f(Xh4[i] + h3/2.0, Yh4[i] + (h3/2.0)*k1);
        Yh4[i+1] = Yh4[i] + h3*k2;
        Xh4[i+1] = Xh4[i] + h3;
    }

    printf("\n%-10s %-15s %-15s\n", "i", "x_i", "Q_i");
    printf("------------------------------------------\n");
    Q[0] = 0.0;

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;
        int idx4 = 4*i;
        double num = fabs(Yh[i] - Yh2[idx2]);
        double den = fabs(Yh2[idx2] - Yh4[idx4]);

        if (den > 1e-12) {
            Q[i] = log(num / den) / log(2.0);
            printf("%-10d %-15lf %-15lf\n", i, Xh[i], Q[i]);
        } else {
            Q[i] = 0.0;
            printf("%-10d %-15lf %-15s\n", i, Xh[i], "N/A");
        }
    }

    save_in_txt(Xh, Q, n1);
    rename("results.txt", "convergence_midpoint.txt");

    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    }
}

void convergence_factor_rk4(int n1, double h1, double X0, double Xf, double Y0) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    double Yh[MAX_SIZE + 1], Yh2[MAX_SIZE*2 + 1], Yh4[MAX_SIZE*4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE*2 + 1], Xh4[MAX_SIZE*4 + 1];
    double Q[MAX_SIZE + 1];

    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Yh[0] = Yh2[0] = Yh4[0] = Y0;

    // RK4 paso h
    for (int i = 0; i < n1; i++) {
        Yh[i+1] = rk4_step(Xh[i], Yh[i], h1);
        Xh[i+1] = Xh[i] + h1;
    }

    // RK4 paso h/2
    for (int i = 0; i < 2*n1; i++) {
        Yh2[i+1] = rk4_step(Xh2[i], Yh2[i], h2);
        Xh2[i+1] = Xh2[i] + h2;
    }

    // RK4 paso h/4
    for (int i = 0; i < 4*n1; i++) {
        Yh4[i+1] = rk4_step(Xh4[i], Yh4[i], h3);
        Xh4[i+1] = Xh4[i] + h3;
    }

    printf("\n%-10s %-15s %-15s\n", "i", "x_i", "Q_i");
    printf("------------------------------------------\n");
    Q[0] = 0.0;

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;
        int idx4 = 4*i;
        double num = fabs(Yh[i] - Yh2[idx2]);
        double den = fabs(Yh2[idx2] - Yh4[idx4]);

        if (den > 1e-12) {
            Q[i] = log(num / den) / log(2.0);
            printf("%-10d %-15lf %-15lf\n", i, Xh[i], Q[i]);
        } else {
            Q[i] = 0.0;
            printf("%-10d %-15lf %-15s\n", i, Xh[i], "N/A");
        }
    }

    save_in_txt(Xh, Q, n1);
    rename("results.txt", "convergence_rk4.txt");

    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    }
}


