#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 200

/**
 * Function to save computed values to a text file
 * @param filename Name of the output file
 * @param X Array of x values
 * @param Y Array of y values
 * @param n Number of subintervals
 */
void save_in_txt(const char *filename, double X[], double Y[], int n);

/**
 * Function to define function f1(x, y1, y2)
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @return The value of the function at (x, y1, y2)
 */
double f12(double x, double y1, double y2);

/**
 * Function to define function f2(x, y1, y2)
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @return The value of the function at (x, y1, y2)
 */
double f22(double x, double y1, double y2);

/**
 * Function to define the exact solution y1(x) for comparison
 * @param x The independent variable
 * @return The exact value of y1 at x
 */
double y1(double x);

/**
 * Function to define the exact solution y2(x) for comparison
 * @param x The independent variable
 * @return The exact value of y2 at x
 */
double y2(double x);

/**
 * Function to define the total derivative of f with respect to x
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @param f Pointer to the function f(x, y1, y2)
 * @return The value of the total derivative at (x, y1, y2)
 */
double fprima2(double x, double y1, double y2, double (*f)(double, double, double));

/**
 * Function to define the total derivative of f with respect to x
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @param f Pointer to the function f(x, y1, y2)
 * @return The value of the total derivative at (x, y1, y2)
 */
double fprima3(double x, double y1, double y2, double y3, double (*f)(double, double, double, double));

/**
 * Function to define the third derivative of y with respect to x
 * @param x The independent variable
 * @param y The dependent variable
 * @return The value of the third derivative at (x, y)
 */
double y3prima(double x, double y);

/**
 * Function to perform a single Runge-Kutta 4th order step for a system of 2 EDOs
 * @param x Current x value
 * @param y1 Pointer to the first dependent variable
 * @param y2 Pointer to the second dependent variable
 * @param h Step size
 * @param f1 Pointer to the first function f1(x, y1, y2)
 * @param f2 Pointer to the second function f2(x, y1, y2)
 * @return void
 */
void rk4_step2(double x, double *y1, double *y2, double h, double (*f1)(double, double, double), double (*f2)(double, double, double));

/**
 * 
 * Function to calculate local truncation error for RK4 method for a system of 2 EDOs
 * @param x Current x value
 * @param y1 Current value of the first dependent variable
 * @param y2 Current value of the second dependent variable
 * @param h Step size
 * @param f1 Pointer to the first function f1(x, y1, y2)
 * @param f2 Pointer to the second function f2(x, y1, y2)
 * @param lte1 Pointer to store the local truncation error for the first variable
 * @param lte2 Pointer to store the local truncation error for the second variable
 * @return void
 */ 
void local_trunc_error_rk4_2(double x, double y1, double y2, double h, double (*f1)(double, double, double), double (*f2)(double, double, double), double *lte1, double *lte2);


/**
 * Function to calculate convergence factor for Euler's method for a system of 2 EDOs
 * @param n1 Number of subintervals
 * @param h1 Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y10 Initial y1 value
 * @param Y20 Initial y2 value
 * @param f1 Pointer to the first function f1(x, y1, y2)
 * @param f2 Pointer to the second function f2(x, y1, y2)
 * @param filename Name of the output file for the first variable
 * @param filename2 Name of the output file for the second variable
 *  */ 
void convergence_factor_euler_2(int n1, double h1, double X0, double Xf,
                                double Y10, double Y20,
                                double (*f1)(double, double, double),
                                double (*f2)(double, double, double),
                                const char *filename,
                                const char *filename2);
                                   
/**
 * Function to calculate convergence factor for RK4 method for a system of 2 EDOs
 * @param n1 Number of subintervals
 * @param h1 Step size
 * @param X0 Initial x value
 * @param Y10 Initial y1 value
 * @param Y20 Initial y2 value
 * @param f1 Pointer to the first function f1(x, y1, y2)
 * @param f2 Pointer to the second function f2(x, y1, y2)
 * @param out1 Name of the output file for the first variable
 * @param out2 Name of the output file for the second variable
 */
void convergence_factor_rk4_2(int n1, double h1, double X0, double Y10, double Y20,
                              double (*f1)(double,double,double),
                              double (*f2)(double,double,double), 
                              const char *out1, const char *out2);


int main(int argc, char const *argv[]) {
    // General variables
    double X0, Xf, Y0, h;
    // Variables for the Neumann method
    double Xp, Yp;
    int n;
    // Number of EDOs in the system
    int edo_count;  
    double X[MAX_SIZE + 1], Y[MAX_SIZE + 1], Y1[MAX_SIZE + 1], Y2[MAX_SIZE + 1], Y3[MAX_SIZE + 1];

    // Number to select to do convergence factor calculation
    int conv_choice;

    // Error to calculate
    double exact_error, local_trunc_error;

    double exact_error1[MAX_SIZE + 1];

    printf("Insert X0 and Xf:\n");
    scanf("%lf %lf", &X0, &Xf);
    

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
            n = MAX_SIZE;
            // printf("Error: number of subintervals exceeds maximum size (%d).\n", MAX_SIZE);
            // return 1;
        }
    }

    printf("Insert the method to use: 1. Euler's method (R.K. 1) 2. Heun method (R.K. 2) 3. Midpoint method (R.K. 2) 4. Runge-Kutta of order 4\n");
    scanf("%d", &choice);

    switch(choice) {
        case 1:
            printf("How many EDO's does it have you system? (2 or 3)\n");
            scanf("%d", &edo_count);
            if(edo_count == 2) {
                X[0] = X0;
                printf("Insert initial data Y01 = Y1(X0):\n");
                scanf("%lf", &Y1[0]);
                printf("Insert initial data Y02 = Y2(X0):\n");
                scanf("%lf", &Y2[0]);

                printf("Do you want to calculate convergence factor for Euler's method? (1.Yes 2.No)\n");
                scanf("%d", &conv_choice);
                
                if(conv_choice == 1) {
                   convergence_factor_euler_2(n, h, X0, Xf, Y1[0], Y2[0], f12, f22, "convergence_euler_edo2.txt", "convergence_euler2_edo2.txt");
                }

                for(int i = 0; i <= n-1; i++) {
                    X[i+1] = X[i] + h;
                    Y1[i+1] = Y1[i] + h * f12(X[i], Y1[i], Y2[i]);
                    Y2[i+1] = Y2[i] + h * f22(X[i], Y1[i], Y2[i]);
                }

                double new_X[MAX_SIZE + 1], new_Y1[MAX_SIZE + 1], new_Y2[MAX_SIZE + 1];

                for(int i = 0; i <= n; i++) {
                    new_X[i] = 0.0;
                    new_Y1[i] = 0.0;
                    new_Y2[i] = 0.0;
                }
                int count = 0;

                // Print results for Y1
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y1", "Euler Y1", "Exact Error");
                printf("-------------------------------------------------------------------------------------------\n");
                
                // Euler's method approximates the curve by means of straight line segments tangent to each point.
                // If these segments lie above the curve → Euler overestimates.
                // If they lie below it → Euler underestimates.
                // If local_trunc_error < 0 => The value calculated by Euler is less than the exact one. Euler underestimates.
                // If local_trunc_error > 0 => The value calculated by Euler is greater than the exact one. Euler overestimates.
                for(int i = n; i >= 0; i--) {
                    if(count == 6) break;
                    double exact_error = fabs(y1(X[i]) - Y1[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.2e\n",
                           i, X[i], y1(X[i]), Y1[i], exact_error);
                    count++;
                    new_X[i] = X[i];
                    new_Y1[i] = Y1[i];
                }
                
                count = 0;
                // Imprimir Y2
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y2", "Euler Y2", "Exact Error");
                printf("------------------------------------------------------------\n");
                for(int i = n; i >= 0; i--) {
                    if(count == 6) break;
                    double exact_error = fabs(y2(X[i]) - Y2[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.2e\n",
                           i, X[i], y2(X[i]), Y2[i], exact_error);
                    count++; 
                    new_Y2[i] = Y2[i];
                }

                // Print results
                printf("X[i]\t\tY1[i]\t\tY2[i]\n");
                // Print all computed points including the last one
                for(int i = n-5; i <= n; i++) {
                    printf("%lf\t%lf\t%lf\n", new_X[i], new_Y1[i], new_Y2[i]);
                }
            
                // Save x[i] and Y[i] in results.txt
                save_in_txt("results_Y1.txt", new_X, new_Y1, n);
                save_in_txt("results_Y2.txt", new_X, new_Y2, n);

                // Finally, we print the results.txt file in a graph using Python to visualize the results
                // system("python3 graph_points.py");
                if (system("test -f graph_points_edo2.py") == 0) {
                    system("python3 graph_points_edo2.py");
                } else {
                    printf("⚠️  Warning: 'graph_points_edo2.py' not found. Skipping graph generation.\n");
                }
            }
            break;
        case 2: 
            break;
         case 3:
            break; 
        case 4:
            printf("How many EDO's does it have you system? (2 or 3)\n");
            scanf("%d", &edo_count);
            if(edo_count == 2) {
                X[0] = X0;
                printf("Insert initial data Y01 = Y1(X0):\n");
                scanf("%lf", &Y1[0]);
                printf("Insert initial data Y02 = Y2(X0):\n");
                scanf("%lf", &Y2[0]);

                printf("Do you want to calculate convergence factor for Runge Kutta's 4 method? (1.Yes 2.No)\n");
                scanf("%d", &conv_choice);
                
                if(conv_choice == 1) {
                   convergence_factor_rk4_2(n, h, X0, Y1[0], Y2[0], f12, f22,
                                "convergence_Q1.txt",
                                "convergence_Q2.txt");
                }

                double k11, k12, k21, k22, k31, k32, k41, k42;
                // Runge-Kutta of order 4
                for(int i = 0; i < n; i++) {
                    k11 = f12(X[i], Y1[i], Y2[i]);
                    k12 = f22(X[i], Y1[i], Y2[i]);
                    k21 = f12(X[i] + h/2.0, Y1[i] + (h/2.0) * k11, Y2[i] + (h/2.0) * k12);
                    k22 = f22(X[i] + h/2.0, Y1[i] + (h/2.0) * k11, Y2[i] + (h/2.0) * k12);
                    k31 = f12(X[i] + h/2.0, Y1[i] + (h/2.0) * k21, Y2[i] + (h/2.0) * k22);
                    k32 = f22(X[i] + h/2.0, Y1[i] + (h/2.0) * k21, Y2[i] + (h/2.0) * k22);
                    k41 = f12(X[i] + h, Y1[i] + h * k31, Y2[i] + h * k32);
                    k42 = f22(X[i] + h, Y1[i] + h * k31, Y2[i] + h * k32);
                    Y1[i+1] = Y1[i] + (h/6.0) * (k11 + 2.0 * k21 + 2.0 * k31 + k41);
                    Y2[i+1] = Y2[i] + (h/6.0) * (k12 + 2.0 * k22 + 2.0 * k32 + k42);
                    X[i+1] = X[i] + h;
                }

                // Calcular LTE una vez y guardarlos en arrays
                double lte1_arr[MAX_SIZE+1], lte2_arr[MAX_SIZE+1];
                lte1_arr[0] = lte2_arr[0] = 0.0;
                for(int i = 1; i <= n; i++) {
                    double lte1, lte2;
                    local_trunc_error_rk4_2(X[i-1], Y1[i-1], Y2[i-1], h, f12, f22, &lte1, &lte2);
                    lte1_arr[i] = lte1;
                    lte2_arr[i] = lte2;
                }

                double new_X[MAX_SIZE + 1], new_Y1[MAX_SIZE + 1], new_Y2[MAX_SIZE + 1];

                for(int i = 0; i <= n; i++) {
                    new_X[i] = 0.0;
                    new_Y1[i] = 0.0;
                    new_Y2[i] = 0.0;
                }
                int count = 0;
                // Imprimir Y1
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y1", "RK4 Y1", "Exact Error");
                printf("------------------------------------------------------------\n");
                for(int i = n; i >= 0; i--) {
                    if(count == 6) break;
                    double exact_error = fabs(y1(X[i]) - Y1[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.2e\n",
                           i, X[i], y1(X[i]), Y1[i], exact_error);
                    count++;
                    new_X[i] = X[i];
                    new_Y1[i] = Y1[i];
                }

                count = 0;
                // Imprimir Y2
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y2", "RK4 Y2", "Exact Error");
                printf("------------------------------------------------------------\n");
                for(int i = n; i >= 0; i--) {
                    if(count == 6) break;
                    double exact_error = fabs(y2(X[i]) - Y2[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.2e\n",
                           i, X[i], y2(X[i]), Y2[i], exact_error);
                    count++; 
                    new_Y2[i] = Y2[i];
                }

                // Print results
                printf("X[i]\t\tY1[i]\t\tY2[i]\n");
                // Print all computed points including the last one
                for(int i = n-5; i <= n; i++) {
                    printf("%lf\t%lf\t%lf\n", new_X[i], new_Y1[i], new_Y2[i]);
                }
            
                // Save x[i] and Y[i] in results.txt
                save_in_txt("results_Y1.txt", new_X, new_Y1, n);
                save_in_txt("results_Y2.txt", new_X, new_Y2, n);

                // Finally, we print the results.txt file in a graph using Python to visualize the results
                if (system("test -f graph_points_edo2.py") == 0) {
                    system("python3 graph_points_edo2.py");
                } else {
                    printf("Warning: 'graph_points_edo2.py' not found. Skipping graph generation.\n");
                }
            }
            break; 
    }
    return 0;
}

// Functions for systems of two EDOs
double f12(double X, double Y1, double Y2) {
    // return 3 * X + Y2;
    return Y2;
}

double f22(double X, double Y1, double Y2) {
    return (-(2*Y2 + 5*Y1));
}

double y1(double x) {
    return(exp(-x) * sin(2*x));
}

double y2(double x) {
    return ((2* exp(-x) * cos(2*x)) - (exp(-x) * sin(2*x)));
}

/**
 * Third derivative of y(x) = e^{-x^2}
 * Computed symbolically as y''' = 4xy(3 - 2x^2)
 */
double y3prima(double x, double y) {
    return 4 * x * y * (3 - 2 * x * x);
}

/* double fprima(double x, double y) {
    double fx = -2 * y;
    double fy = -2 * x;
    return fx + fy * f(x, y);
} */

double fprima2(double x, double y1, double y2,
               double (*f)(double, double, double)) {
    double fx = -2 * (y1 + y2);  // ejemplo de ∂f/∂x
    double fy = -2 * x;          // ejemplo de ∂f/∂y
    return fx + fy * f(x, y1, y2);
}

double fprima3(double x, double y1, double y2, double y3,
               double (*f)(double, double, double, double)) {
    double fx = -2 * (y1 + y2 + y3);
    double fy = -2 * x;
    return fx + fy * f(x, y1, y2, y3);
}

void rk4_step2(double x, double *y1, double *y2, double h,
               double (*f1)(double, double, double),
               double (*f2)(double, double, double)) {
    double k1y1 = f1(x, y1[0], y2[0]);
    double k1y2 = f2(x, y1[0], y2[0]);

    double k2y1 = f1(x + h/2.0, y1[0] + h/2.0 * k1y1, y2[0] + h/2.0 * k1y2);
    double k2y2 = f2(x + h/2.0, y1[0] + h/2.0 * k1y1, y2[0] + h/2.0 * k1y2);

    double k3y1 = f1(x + h/2.0, y1[0] + h/2.0 * k2y1, y2[0] + h/2.0 * k2y2);
    double k3y2 = f2(x + h/2.0, y1[0] + h/2.0 * k2y1, y2[0] + h/2.0 * k2y2);

    double k4y1 = f1(x + h, y1[0] + h * k3y1, y2[0] + h * k3y2);
    double k4y2 = f2(x + h, y1[0] + h * k3y1, y2[0] + h * k3y2);

    y1[0] += h/6.0 * (k1y1 + 2*k2y1 + 2*k3y1 + k4y1);
    y2[0] += h/6.0 * (k1y2 + 2*k2y2 + 2*k3y2 + k4y2);
}

void local_trunc_error_rk4_2(double x, double y1, double y2, double h,
                             double (*f1)(double, double, double),
                             double (*f2)(double, double, double),
                             double *lte1, double *lte2) {
    double Y_full[2] = {y1, y2};
    double Y_half[2] = {y1, y2};

    // Paso completo
    rk4_step2(x, Y_full, Y_full+1, h, f1, f2);

    // Dos pasos de h/2
    double mid[2] = {y1, y2};
    rk4_step2(x, mid, mid+1, h/2.0, f1, f2);
    rk4_step2(x + h/2.0, mid, mid+1, h/2.0, f1, f2);

    *lte1 = fabs(mid[0] - Y_full[0]) / 15.0;
    *lte2 = fabs(mid[1] - Y_full[1]) / 15.0;
}


void rk4_step3(double x, double *y1, double *y2, double *y3, double h,
               double (*f1)(double, double, double, double),
               double (*f2)(double, double, double, double),
               double (*f3)(double, double, double, double)) {
    double k1y1 = f1(x, y1[0], y2[0], y3[0]);
    double k1y2 = f2(x, y1[0], y2[0], y3[0]);
    double k1y3 = f3(x, y1[0], y2[0], y3[0]);

    double k2y1 = f1(x + h/2.0, y1[0] + h/2.0*k1y1, y2[0] + h/2.0*k1y2, y3[0] + h/2.0*k1y3);
    double k2y2 = f2(x + h/2.0, y1[0] + h/2.0*k1y1, y2[0] + h/2.0*k1y2, y3[0] + h/2.0*k1y3);
    double k2y3 = f3(x + h/2.0, y1[0] + h/2.0*k1y1, y2[0] + h/2.0*k1y2, y3[0] + h/2.0*k1y3);

    double k3y1 = f1(x + h/2.0, y1[0] + h/2.0*k2y1, y2[0] + h/2.0*k2y2, y3[0] + h/2.0*k2y3);
    double k3y2 = f2(x + h/2.0, y1[0] + h/2.0*k2y1, y2[0] + h/2.0*k2y2, y3[0] + h/2.0*k2y3);
    double k3y3 = f3(x + h/2.0, y1[0] + h/2.0*k2y1, y2[0] + h/2.0*k2y2, y3[0] + h/2.0*k2y3);

    double k4y1 = f1(x + h, y1[0] + h*k3y1, y2[0] + h*k3y2, y3[0] + h*k3y3);
    double k4y2 = f2(x + h, y1[0] + h*k3y1, y2[0] + h*k3y2, y3[0] + h*k3y3);
    double k4y3 = f3(x + h, y1[0] + h*k3y1, y2[0] + h*k3y2, y3[0] + h*k3y3);

    y1[0] += h/6.0 * (k1y1 + 2*k2y1 + 2*k3y1 + k4y1);
    y2[0] += h/6.0 * (k1y2 + 2*k2y2 + 2*k3y2 + k4y2);
    y3[0] += h/6.0 * (k1y3 + 2*k2y3 + 2*k3y3 + k4y3);
}



void save_in_txt(const char *filename, double X[], double Y[], int n) {
    FILE *archivo = fopen(filename, "w");
    if (archivo == NULL) {
        printf("Error: Unable to create file '%s'.\n", filename);
        exit(1);
    }

    for (int i = n-5; i <= n; i++) {
        fprintf(archivo, "%lf\t%lf\n", X[i], Y[i]);
    }

    fclose(archivo);
    printf("File '%s' saved successfully.\n", filename);
}

void convergence_factor_euler_2(int n1, double h1, double X0, double Xf,
                                double Y10, double Y20,
                                double (*f1)(double, double, double),
                                double (*f2)(double, double, double),
                                const char *filename,
                                const char *filename2) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    // Arrays para Y1 e Y2
    double Y1h[MAX_SIZE + 1], Y2h[MAX_SIZE + 1];
    double Y1h2[MAX_SIZE * 2 + 1], Y2h2[MAX_SIZE * 2 + 1];
    double Y1h4[MAX_SIZE * 4 + 1], Y2h4[MAX_SIZE * 4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE * 2 + 1], Xh4[MAX_SIZE * 4 + 1];
    double Q1[MAX_SIZE + 1], Q2[MAX_SIZE + 1];

    // Inicialización
    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Y1h[0] = Y1h2[0] = Y1h4[0] = Y10;
    Y2h[0] = Y2h2[0] = Y2h4[0] = Y20;

    // ---------- Euler con paso h ----------
    for (int i = 0; i < n1; i++) {
        Xh[i+1] = Xh[i] + h1;
        Y1h[i+1] = Y1h[i] + h1 * f1(Xh[i], Y1h[i], Y2h[i]);
        Y2h[i+1] = Y2h[i] + h1 * f2(Xh[i], Y1h[i], Y2h[i]);
    }

    // ---------- Euler con paso h/2 ----------
    for (int i = 0; i < 2*n1; i++) {
        Xh2[i+1] = Xh2[i] + h2;
        Y1h2[i+1] = Y1h2[i] + h2 * f1(Xh2[i], Y1h2[i], Y2h2[i]);
        Y2h2[i+1] = Y2h2[i] + h2 * f2(Xh2[i], Y1h2[i], Y2h2[i]);
    }

    // ---------- Euler con paso h/4 ----------
    for (int i = 0; i < 4*n1; i++) {
        Xh4[i+1] = Xh4[i] + h3;
        Y1h4[i+1] = Y1h4[i] + h3 * f1(Xh4[i], Y1h4[i], Y2h4[i]);
        Y2h4[i+1] = Y2h4[i] + h3 * f2(Xh4[i], Y1h4[i], Y2h4[i]);
    }

    // ---------- Calcular factores Q ----------
    printf("\n%-10s %-15s %-15s %-15s\n", "i", "x_i", "Q1_i", "Q2_i");
    printf("----------------------------------------------------------\n");

    Q1[0] = Q2[0] = 0.0; // no definidos en el punto inicial

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;   // posición equivalente para h/2
        int idx4 = 4*i;   // posición equivalente para h/4

        double num1 = fabs(Y1h[i] - Y1h2[idx2]);
        double den1 = fabs(Y1h2[idx2] - Y1h4[idx4]);
        double num2 = fabs(Y2h[i] - Y2h2[idx2]);
        double den2 = fabs(Y2h2[idx2] - Y2h4[idx4]);

        Q1[i] = (den1 > 1e-12) ? log(num1 / den1) / log(2.0) : 0.0;
        Q2[i] = (den2 > 1e-12) ? log(num2 / den2) / log(2.0) : 0.0;

        printf("%-10d %-15lf %-15lf %-15lf\n", i, Xh[i], Q1[i], Q2[i]);
    }

    // Guardar resultados
    save_in_txt("results_Q1.txt", Xh, Q1, n1);
    save_in_txt("results_Q2.txt", Xh, Q2, n1);
    rename("results_Q1.txt", filename);
    rename("results_Q2.txt", filename2);

    // Generar gráfico si existe script Python
    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    } else {
        printf("⚠️  Warning: 'graph_convergence.py' not found. Skipping graph generation.\n");
    }
}

void convergence_factor_rk4_2(int n1, double h1, double X0, double Y10, double Y20,
                              double (*f1)(double,double,double),
                              double (*f2)(double,double,double),
                              const char *filename1,
                              const char *filename2) {
    double h2 = h1/2.0;
    double h3 = h1/4.0;

    double Y1h[MAX_SIZE+1], Y2h[MAX_SIZE+1];
    double Y1h2[MAX_SIZE*2+1], Y2h2[MAX_SIZE*2+1];
    double Y1h4[MAX_SIZE*4+1], Y2h4[MAX_SIZE*4+1];
    double Xh[MAX_SIZE+1], Xh2[MAX_SIZE*2+1], Xh4[MAX_SIZE*4+1];
    double Q1[MAX_SIZE+1], Q2[MAX_SIZE+1];

    // Inicialización
    Xh[0]=Xh2[0]=Xh4[0]=X0;
    Y1h[0]=Y1h2[0]=Y1h4[0]=Y10;
    Y2h[0]=Y2h2[0]=Y2h4[0]=Y20;

    // RK4 paso h
    for(int i=0;i<n1;i++){
        rk4_step2(Xh[i], &Y1h[i], &Y2h[i], h1, f1, f2);
        Xh[i+1] = Xh[i]+h1;
    }
    // RK4 paso h/2
    for(int i=0;i<2*n1;i++){
        rk4_step2(Xh2[i], &Y1h2[i], &Y2h2[i], h2, f1, f2);
        Xh2[i+1] = Xh2[i]+h2;
    }
    // RK4 paso h/4
    for(int i=0;i<4*n1;i++){
        rk4_step2(Xh4[i], &Y1h4[i], &Y2h4[i], h3, f1, f2);
        Xh4[i+1] = Xh4[i]+h3;
    }

    printf("\n%-10s %-15s %-15s %-15s\n","i","x_i","Q1_i","Q2_i");
    printf("------------------------------------------------\n");

    // Calcular factores de convergencia usando indices calculados a partir de x_i
    for(int i=1;i<=n1;i++){
        double xi = X0 + i*h1;                 // punto exacto x
        int idx2 = (int)((xi - X0)/h2 + 0.5);  // redondeo al índice más cercano en h/2
        int idx4 = (int)((xi - X0)/h3 + 0.5);  // redondeo al índice más cercano en h/4

        double num1 = fabs(Y1h[i]-Y1h2[idx2]);
        double den1 = fabs(Y1h2[idx2]-Y1h4[idx4]);
        double num2 = fabs(Y2h[i]-Y2h2[idx2]);
        double den2 = fabs(Y2h2[idx2]-Y2h4[idx4]);

        Q1[i] = (den1>1e-12)? log(num1/den1)/log(2.0) : 0.0;
        Q2[i] = (den2>1e-12)? log(num2/den2)/log(2.0) : 0.0;

        char s1[64], s2[64];
        if (Q1[i] != 0.0) {
            snprintf(s1, sizeof(s1), "%15.6f", Q1[i]);
        } else {
            snprintf(s1, sizeof(s1), "N/A");
        }
        if (Q2[i] != 0.0) {
            snprintf(s2, sizeof(s2), "%15.6f", Q2[i]);
        } else {
            snprintf(s2, sizeof(s2), "N/A");
        }

        printf("%-10d %-15lf %-15s %-15s\n", i, xi, s1, s2);
    }

    // Guardar resultados en TXT
    save_in_txt(filename1, Xh, Q1, n1);
    save_in_txt(filename2, Xh, Q2, n1);

    // Generar gráfico si existe script Python
    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    } else {
        printf("⚠️  Warning: 'graph_convergence.py' not found. Skipping graph generation.\n");
    }
}




