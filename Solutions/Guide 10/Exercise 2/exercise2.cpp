#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 10000

const double G = 6.672e-11;
const double Me = 5.9742e24;
const double GM = 3.986004418e14;  // GM terrestre en m³/s²
const double R_earth = 6378140.0;  // Radio terrestre en metros

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
 * Function to define function f1(x, y1, y2, y3)
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @param y3 The third dependent variable
 * @return The value of the function at (x, y1, y2, y3)
 */
double f13(double x, double y1, double y2, double y3);

/**
 * Function to define function f2(x, y1, y2, y3)
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @param y3 The third dependent variable
 * @return The value of the function at (x, y1, y2, y3)
 */
double f23(double x, double y1, double y2, double y3);

/**
 * Function to define function f3(x, y1, y2, y3)
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @param y3 The third dependent variable
 * @return The value of the function at (x, y1, y2, y3)
 */
double f33(double x, double y1, double y2, double y3);

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
 * Function to define the exact solution y3(x) for comparison
 * @param x The independent variable
 * @return The exact value of y3 at x
 */
double y3(double x);

/**
 * Function to define function f1(x, y1, y2, y3, y4)
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @param y3 The third dependent variable
 * @param y4 The fourth dependent variable
 * @return The value of the function at (x, y1, y2, y3, y4)
 */
double f14(double x, double y1, double y2, double y3, double y4);

/**
 * Function to define function f2(x, y1, y2, y3, y4)
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @param y3 The third dependent variable
 * @param y4 The fourth dependent variable
 * @return The value of the function at (x, y1, y2, y3, y4)
 */
double f24(double x, double y1, double y2, double y3, double y4);

/**
 * Function to define function f3(x, y1, y2, y3, y4)
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @param y3 The third dependent variable
 * @param y4 The fourth dependent variable
 * @return The value of the function at (x, y1, y2, y3, y4)
 */
double f34(double x, double y1, double y2, double y3, double y4);

/**
 * Function to define function f4(x, y1, y2, y3, y4)
 * @param x The independent variable
 * @param y1 The first dependent variable
 * @param y2 The second dependent variable
 * @param y3 The third dependent variable
 * @param y4 The fourth dependent variable
 * @return The value of the function at (x, y1, y2, y3, y4)
 */
double f44(double x, double y1, double y2, double y3, double y4);

/**
 * Function to perform a single Runge-Kutta 4th order step for a system of 4 EDOs
 * @param x Current x value
 * @param y1 Pointer to the first dependent variable
 * @param y2 Pointer to the second dependent variable
 * @param y3 Pointer to the third dependent variable
 * @param y4 Pointer to the fourth dependent variable
 * @param h Step size
 * @param f1 Pointer to the first function f1(x, y1, y2, y3, y4)
 * @param f2 Pointer to the second function f2(x, y1, y2, y3, y4)
 * @param f3 Pointer to the third function f3(x, y1, y2, y3, y4)
 * @param f4 Pointer to the fourth function f4(x, y1, y2, y3, y4)
 * @return void
 */
void rk4_step4(double x, double *y1, double *y2, double *y3, double *y4, double h,
               double (*f1)(double, double, double, double, double),
               double (*f2)(double, double, double, double, double),
               double (*f3)(double, double, double, double, double),
               double (*f4)(double, double, double, double, double));

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
 * Function to perform a single Runge-Kutta 4th order step for a system of 3 EDOs
 * @param x Current x value
 * @param y1 Pointer to the first dependent variable
 * @param y2 Pointer to the second dependent variable
 * @param y3 Pointer to the third dependent variable
 * @param h Step size
 * @param f1 Pointer to the first function f1(x, y1, y2, y3)
 * @param f2 Pointer to the second function f2(x, y1, y2, y3)
 * @param f3 Pointer to the third function f3(x, y1, y2, y3)
 * @return void
 */
void rk4_step3(double x, double *y1, double *y2, double *y3, double h, 
                double (*f1)(double, double, double, double),
                double (*f2)(double, double, double, double),
                double (*f3)(double, double, double, double));

/** * Function to calculate local truncation error for RK4 method for a system of 4 EDOs
 * @param x Current x value
 * @param y1 Current value of the first dependent variable
 * @param y2 Current value of the second dependent variable
 * @param y3 Current value of the third dependent variable
 * @param y4_val Current value of the fourth dependent variable
 * @param h Step size
 * @param f1 Pointer to the first function f1(x, y1, y2, y3, y4)
 * @param f2 Pointer to the second function f2(x, y1, y2, y3, y4)
 * @param f3 Pointer to the third function f3(x, y1, y2, y3, y4)
 * @param f4 Pointer to the fourth function f4(x, y1, y2, y3, y4)
 * @param lte1 Pointer to store the local truncation error for the first variable
 * @param lte2 Pointer to store the local truncation error for the second variable
 * @param lte3 Pointer to store the local truncation error for the third variable
 * @param lte4 Pointer to store the local truncation error for the fourth variable
 * @return void
 */
void local_trunc_error_rk4_4(double x, double y1, double y2, double y3, double y4_val, double h,
                             double (*f1)(double, double, double, double, double),
                             double (*f2)(double, double, double, double, double),
                             double (*f3)(double, double, double, double, double),
                             double (*f4)(double, double, double, double, double),
                             double *lte1, double *lte2, double *lte3, double *lte4);


/** * Function to calculate convergence factor for RK4 method for a system of 4 EDOs
 * @param n1 Number of subintervals
 * @param h1 Step size
 * @param X0 Initial x value
 * @param Y10 Initial y1 value
 * @param Y20 Initial y2 value
 * @param Y30 Initial y3 value
 * @param Y40 Initial y4 value
 * @param f1 Pointer to the first function f1(x, y1, y2, y3, y4)
 * @param f2 Pointer to the second function f2(x, y1, y2, y3, y4)
 * @param f3 Pointer to the third function f3(x, y1, y2, y3, y4)
 * @param f4 Pointer to the fourth function f4(x, y1, y2, y3, y4)
 * @param out1 Name of the output file for the first variable
 * @param out2 Name of the output file for the second variable
 * @param out3 Name of the output file for the third variable
 * @param out4 Name of the output file for the fourth variable
 */
void convergence_factor_rk4_4(int n1, double h1, double X0,
                              double Y10, double Y20, double Y30, double Y40,
                              double (*f1)(double, double, double, double, double),
                              double (*f2)(double, double, double, double, double),
                              double (*f3)(double, double, double, double, double),
                              double (*f4)(double, double, double, double, double),
                              const char *out1, const char *out2, const char *out3, const char *out4);

/**
 * Function to calculate local truncation error for RK4 method for a system of 3 EDOs
 * @param x Current x value
 * @param y1 Current value of the first dependent variable
 * @param y2 Current value of the second dependent variable
 * @param y3 Current value of the third dependent variable
 * @param h Step size
 * @param f1 Pointer to the first function f1(x, y1, y2, y3)
 * @param f2 Pointer to the second function f2(x, y1, y2, y3)
 * @param f3 Pointer to the third function f3(x, y1, y2, y3)
 * @param lte1 Pointer to store the local truncation error for the first variable
 * @param lte2 Pointer to store the local truncation error for the second variable
 * @param lte3 Pointer to store the local truncation error for the third variable
 * @return void
 */ 
void local_trunc_error_rk4_3(double x, double y1, double y2, double y3, double h,
                             double (*f1)(double, double, double, double),
                             double (*f2)(double, double, double, double),
                             double (*f3)(double, double, double, double), double *lte1, double *lte2, double *lte3);

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
 * Function to calculate convergence factor for Euler's method for a system of 3 EDOs
 * @param n1 Number of subintervals
 * @param h1 Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y10 Initial y1 value
 * @param Y20 Initial y2 value
 * @param Y30 Initial y3 value
 * @param f1 Pointer to the first function f1(x, y1, y2, y3)
 * @param f2 Pointer to the second function f2(x, y1, y2, y3)
 * @param f3 Pointer to the third function f3(x, y1, y2, y3)
 */
void convergence_factor_euler_3(int n1, double h1, double X0, double Xf,
                                double Y10, double Y20, double Y30,
                                double (*f1)(double, double, double, double),
                                double (*f2)(double, double, double, double),
                                double (*f3)(double, double, double, double),
                                const char *filename1,
                                const char *filename2,
                                const char *filename3);


/** * Function to calculate convergence factor for Euler's method for a system of 4 EDOs
 * @param n1 Number of subintervals
 * @param h1 Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y10 Initial y1 value
 * @param Y20 Initial y2 value
 * @param Y30 Initial y3 value
 * @param Y40 Initial y4 value
 * @param f1 Pointer to the first function f1(x, y1, y2, y3, y4)
 * @param f2 Pointer to the second function f2(x, y1, y2, y3, y4)
 * @param f3 Pointer to the third function f3(x, y1, y2, y3, y4)
 * @param f4 Pointer to the fourth function f4(x, y1, y2, y3, y4)
 * @param filename1 Name of the output file for the first variable
 * @param filename2 Name of the output file for the second variable
 * @param filename3 Name of the output file for the third variable
 * @param filename4 Name of the output file for the fourth variable
 */
void convergence_factor_euler_4(int n1, double h1, double X0, double Xf,
                                double Y10, double Y20, double Y30, double Y40,
                                double (*f1)(double, double, double, double, double),
                                double (*f2)(double, double, double, double, double),
                                double (*f3)(double, double, double, double, double),
                                double (*f4)(double, double, double, double, double),
                                const char *filename1, const char *filename2,
                                const char *filename3, const char *filename4);
                                
/**
 * Function to calculate convergence factor for Heun method for a system of 2 EDOs
 * @param n1 Number of subintervals
 * @param h1 Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y10 Initial y1 value
 * @param Y20 Initial y2 value
 * @param f1 Pointer to the first function f1(x, y1, y2)
 * @param f2 Pointer to the second function f2(x, y1, y2)
 * @param filename1 Name of the output file for the first variable
 * @param filename2 Name of the output file for the second variable
 */
void convergence_factor_heun_2(int n1, double h1, double X0, double Xf,
                               double Y10, double Y20,
                               double (*f1)(double, double, double),
                               double (*f2)(double, double, double),
                               const char *filename1,
                               const char *filename2);

/**
 * Function to calculate convergence factor for Heun method for a system of 3 EDOs
 * @param n1 Number of subintervals
 * @param h1 Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y10 Initial y1 value
 * @param Y20 Initial y2 value
 * @param Y30 Initial y3 value
 * @param f1 Pointer to the first function f1(x, y1, y2, y3)
 * @param f2 Pointer to the second function f2(x, y1, y2, y3)
 * @param f3 Pointer to the third function f3(x, y1, y2, y3)
 * @param filename1 Name of the output file for the first variable
 * @param filename2 Name of the output file for the second variable
 * @param filename3 Name of the output file for the third variable
 */
void convergence_factor_heun_3(int n1, double h1, double X0, double Xf,
                               double Y10, double Y20, double Y30,
                               double (*f1)(double, double, double, double),
                               double (*f2)(double, double, double, double),
                               double (*f3)(double, double, double, double),
                               const char *filename1,
                               const char *filename2,
                               const char *filename3);

/**
 * Function to calculate convergence factor for Midpoint method for a system of 2 EDOs
 * @param n1 Number of subintervals
 * @param h1 Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y10 Initial y1 value
 * @param Y20 Initial y2 value
 * @param f1 Pointer to the first function f1(x, y1, y2)
 * @param f2 Pointer to the second function f2(x, y1, y2)
 * @param fileQ1 Name of the output file for the first variable
 * @param fileQ2 Name of the output file for the second variable
 */
void convergence_factor_midpoint_2(int n1, double h1, double X0, double Xf,
                                   double Y10, double Y20,
                                   double (*f1)(double, double, double),
                                   double (*f2)(double, double, double),
                                   const char *fileQ1,
                                   const char *fileQ2);


/**
 * Function to calculate convergence factor for Midpoint method for a system of 3 EDOs
 * @param n1 Number of subintervals
 * @param h1 Step size
 * @param X0 Initial x value
 * @param Xf Final x value
 * @param Y10 Initial y1 value
 * @param Y20 Initial y2 value
 * @param Y30 Initial y3 value
 * @param f1 Pointer to the first function f1(x, y1, y2, y3)
 * @param f2 Pointer to the second function f2(x, y1, y2, y3)
 * @param f3 Pointer to the third function f3(x, y1, y2, y3)
 * @param fileQ1 Name of the output file for the first variable
 * @param fileQ2 Name of the output file for the second variable
 * @param fileQ3 Name of the output file for the third variable
 */
void convergence_factor_midpoint_3(int n1, double h1, double X0, double Xf,
                                   double Y10, double Y20, double Y30,
                                   double (*f1)(double, double, double, double),
                                   double (*f2)(double, double, double, double),
                                   double (*f3)(double, double, double, double),
                                   const char *fileQ1,
                                   const char *fileQ2,
                                   const char *fileQ3);
                                   
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

/**
 * Function to calculate convergence factor for RK4 method for a system of 3 EDOs
 * @param n1 Number of subintervals
 * @param h1 Step size
 * @param X0 Initial x value
 * @param Y10 Initial y1 value
 * @param Y20 Initial y2 value
 * @param Y30 Initial y3 value
 * @param f1 Pointer to the first function f1(x, y1, y2, y3)
 * @param f2 Pointer to the second function f2(x, y1, y2, y3)
 * @param f3 Pointer to the third function f3(x, y1, y2, y3)
 * @param out1 Name of the output file for the first variable
 * @param out2 Name of the output file for the second variable
 * @param out3 Name of the output file for the third variable
 */
void convergence_factor_rk4_3(int n1, double h1, double X0,
                              double Y10, double Y20, double Y30,
                              double (*f1)(double,double,double,double),
                              double (*f2)(double,double,double,double),
                              double (*f3)(double,double,double,double),
                              const char *out1, const char *out2, const char *out3);




int main(int argc, char const *argv[]) {
    // General variables
    double X0, Xf, Y0, h;
    // Variables for the Neumann method
    double Xp, Yp;
    int n;
    // Number of EDOs in the system
    int edo_count;  
    double X[MAX_SIZE + 1], Y[MAX_SIZE + 1], Y1[MAX_SIZE + 1], Y2[MAX_SIZE + 1], Y3[MAX_SIZE + 1], Y4[MAX_SIZE + 1];

    // Number to select to do convergence factor calculation
    int conv_choice;

    // Error to calculate
    double exact_error, local_trunc_error;

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
            printf("Error: number of subintervals exceeds maximum size (%d).\n", MAX_SIZE);
            return 1;
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

                // Print results for Y1
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y1", "Euler Y1", "Exact Error", "Local Trunc. Err");
                printf("-------------------------------------------------------------------------------------------\n");
                
                // Euler's method approximates the curve by means of straight line segments tangent to each point.
                // If these segments lie above the curve → Euler overestimates.
                // If they lie below it → Euler underestimates.
                // If local_trunc_error < 0 => The value calculated by Euler is less than the exact one. Euler underestimates.
                // If local_trunc_error > 0 => The value calculated by Euler is greater than the exact one. Euler overestimates.
                for(int i = 0; i <= n; i++) {
                    exact_error = fabs(y1(X[i]) - Y1[i]);
                    local_trunc_error = (h * h / 2.0) * fprima2(X[i], y1(X[i]), y2(X[i]), f12); 
                    printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n", 
                           i, X[i], y1(X[i]), Y1[i], exact_error, local_trunc_error);
                }

                // Print results for Y2
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y2", "Euler Y2", "Exact Error", "Local Trunc. Err");
                printf("-------------------------------------------------------------------------------------------\n");
                
                for(int i = 0; i <= n; i++) {
                    exact_error = fabs(y2(X[i]) - Y2[i]);
                    local_trunc_error = (h * h / 2.0) * fprima2(X[i], y1(X[i]), y2(X[i]), f22); 
                    printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n", 
                           i, X[i], y2(X[i]), Y2[i], exact_error, local_trunc_error);
                }

                // Print results
                printf("X[i]\t\tY1[i]\t\tY2[i]\n");
                // Print all computed points including the last one
                for(int i = 0; i <= n; i++) {
                    printf("%lf\t%lf\t%lf\n", X[i], Y1[i], Y2[i]);
                }
            
                // Save x[i] and Y[i in results.txt]
                save_in_txt("results_Y1.txt", X, Y1, n);
                save_in_txt("results_Y2.txt", X, Y2, n);

                // Finally, we print the results.txt file in a graph using Python to visualize the results
                // system("python3 graph_points.py");
                if (system("test -f graph_points_edo2.py") == 0) {
                    system("python3 graph_points_edo2.py");
                } else {
                    printf("⚠️  Warning: 'graph_points_edo2.py' not found. Skipping graph generation.\n");
                }
            } else if(edo_count == 3) {
                X[0] = X0;
                printf("Insert initial data Y01 = Y1(X0):\n");
                scanf("%lf", &Y1[0]);
                printf("Insert initial data Y02 = Y2(X0):\n");
                scanf("%lf", &Y2[0]);
                printf("Insert initial data Y03 = Y3(X0):\n");
                scanf("%lf", &Y3[0]);

                printf("Do you want to calculate convergence factor for Euler's method? (1.Yes 2.No)\n");
                scanf("%d", &conv_choice);
                
                if(conv_choice == 1) {
                    convergence_factor_euler_3(n, h, X0, Xf, Y1[0], Y2[0], Y3[0], f13, f23, f33, "convergence_euler_edo3.txt"
                        , "convergence_euler2_edo3.txt",
                        "convergence_euler3_edo3.txt");
                }

                // Calculate Euler method for the system of 3 EDOs
                for(int i = 0; i <= n-1; i++) {
                    X[i+1] = X[i] + h;
                    Y1[i+1] = Y1[i] + h * f13(X[i], Y1[i], Y2[i], Y3[i]);
                    Y2[i+1] = Y2[i] + h * f23(X[i], Y1[i], Y2[i], Y3[i]);
                    Y3[i+1] = Y3[i] + h * f33(X[i], Y1[i], Y2[i], Y3[i]);
                }
                // Print results for Y1
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y1", "Euler Y1", "Exact Error", "Local Trunc. Err");
                printf("-------------------------------------------------------------------------------------------\n");
                
                // Euler's method approximates the curve by means of straight line segments tangent to each point.
                // If these segments lie above the curve → Euler overestimates.
                // If they lie below it → Euler underestimates.
                // If local_trunc_error < 0 => The value calculated by Euler is less than the exact one. Euler underestimates.
                // If local_trunc_error > 0 => The value calculated by Euler is greater than the exact one. Euler overestimates.
                for(int i = 0; i <= n; i++) {
                    exact_error = fabs(y1(X[i]) - Y1[i]);
                    local_trunc_error = (h * h / 2.0) * fprima3(X[i], y1(X[i]), y2(X[i]), y3(X[i]), f13); 
                    printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n", 
                           i, X[i], y1(X[i]), Y1[i], exact_error, local_trunc_error);
                }
                // Print results for Y2
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y2", "Euler Y2", "Exact Error", "Local Trunc. Err");
                printf("-------------------------------------------------------------------------------------------\n");
                
                for(int i = 0; i <= n; i++) {
                    exact_error = fabs(y2(X[i]) - Y2[i]);
                    local_trunc_error = (h * h / 2.0) * fprima3(X[i], y1(X[i]), y2(X[i]), y3(X[i]), f23); 
                    printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n", 
                           i, X[i], y2(X[i]), Y2[i], exact_error, local_trunc_error);
                }
                // Print results for Y3
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y3", "Euler Y3", "Exact Error", "Local Trunc. Err");
                printf("-------------------------------------------------------------------------------------------\n");
                
                for(int i = 0; i <= n; i++) {
                    exact_error = fabs(y3(X[i]) - Y3[i]);
                    local_trunc_error = (h * h / 2.0) * fprima3(X[i], y1(X[i]), y2(X[i]), y3(X[i]), f33); 
                    printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n", 
                           i, X[i], y3(X[i]), Y3[i], exact_error, local_trunc_error);
                }
                // Print results
                printf("X[i]\t\tY1[i]\tY2[i]\tY3[i]\n");
                // Print all computed points including the last one
                for(int i = 0; i <= n; i++) {
                    printf("%lf\t%lf\t%lf\t%lf\n", X[i], Y1[i], Y2[i], Y3[i]);
                }
            
                // Save x[i] and Y[i in results.txt]
                save_in_txt("results_Y1", X, Y1, n);
                save_in_txt("results_Y2", X, Y2, n);
                save_in_txt("results_Y3", X, Y3, n);
                // Finally, we print the results.txt file in a graph using Python to visualize the results
                // system("python3 graph_points.py");
                if (system("test -f graph_points_edo3.py") == 0) {
                    system("python3 graph_points_edo3.py");
                } else {
                    printf("⚠️  Warning: 'graph_points_edo3.py' not found. Skipping graph generation.\n");
                }
            } else if(edo_count == 4) {
                X[0] = X0;
                // CONDICIONES INICIALES CORRECTAS
                double r0 = 6378140.0 + 772000.0;  // R_earth + 772 km
                double vr0 = 0.0;                   // Velocidad radial inicial
                double theta0 = 0.0;                // Ángulo inicial
                double omega0 = 6700.0 / r0;        // Velocidad angular inicial
                        
                printf("Using initial conditions for spacecraft:\n");
                printf("Y01 (r0) = %.1f m\n", r0);
                printf("Y02 (vr0) = %.1f m/s\n", vr0);
                printf("Y03 (theta0) = %.1f rad\n", theta0);
                printf("Y04 (omega0) = %.6f rad/s\n", omega0);
                        
                Y1[0] = r0;
                Y2[0] = vr0;
                Y3[0] = theta0;
                Y4[0] = omega0;
                        
                printf("Do you want to calculate convergence factor for Euler's method? (1.Yes 2.No)\n");
                scanf("%d", &conv_choice);
                        
                if(conv_choice == 1) {
                    convergence_factor_euler_4(n, h, X0, Xf, Y1[0], Y2[0], Y3[0], Y4[0], 
                                             f14, f24, f34, f44,
                                             "convergence_euler_edo4.txt",
                                             "convergence_euler2_edo4.txt",
                                             "convergence_euler3_edo4.txt",
                                             "convergence_euler4_edo4.txt");
                }
            
                int impact_index = n;
                int found_impact = 0;
            
                // Euler con detección de impacto
                for(int i = 0; i < n-1; i++) {
                    X[i+1] = X[i] + h;
                    Y1[i+1] = Y1[i] + h * f14(X[i], Y1[i], Y2[i], Y3[i], Y4[i]);
                    Y2[i+1] = Y2[i] + h * f24(X[i], Y1[i], Y2[i], Y3[i], Y4[i]);
                    Y3[i+1] = Y3[i] + h * f34(X[i], Y1[i], Y2[i], Y3[i], Y4[i]);
                    Y4[i+1] = Y4[i] + h * f44(X[i], Y1[i], Y2[i], Y3[i], Y4[i]);
                    
                    // Verificar conservación de energía (para debug)
                    double energy = 0.5 * (Y2[i+1]*Y2[i+1] + Y1[i+1]*Y1[i+1]*Y4[i+1]*Y4[i+1]) - GM/Y1[i+1];
                    if (i % 100 == 0) {
                        printf("t=%.1f s, r=%.0f m, E=%.6e J/kg\n", X[i+1], Y1[i+1], energy);
                    }
                    
                    // Detección de impacto
                    if (Y1[i+1] <= R_earth && !found_impact) {
                        impact_index = i + 1;
                        found_impact = 1;
                        printf("\n*** IMPACTO DETECTADO ***\n");
                        printf("Tiempo: %.2f s\n", X[i+1]);
                        printf("Ángulo θ: %.6f rad\n", Y3[i+1]);
                        printf("Radio: %.0f m\n", Y1[i+1]);
                        break;
                    }
                    
                    // Seguridad: si el radio se hace muy pequeño, parar
                    if (Y1[i+1] < 0.5 * R_earth) {
                        impact_index = i + 1;
                        found_impact = 1;
                        printf("ERROR: Radio demasiado pequeño\n");
                        break;
                    }
                }
            
                int n_effective = found_impact ? impact_index : n;
            
                // Resultados
                printf("\nRESULTADOS EULER - IMPACTO:\n");
                printf("θ en impacto = %.6f rad (%.2f°)\n", Y3[n_effective], Y3[n_effective] * 180.0/M_PI);
                printf("Tiempo = %.2f s\n", X[n_effective]);
                printf("Radio final = %.0f m\n", Y1[n_effective]);
            }
            break;
        case 2: 
            printf("How many EDO's does it have your system? (2 or 3)\n");
            scanf("%d", &edo_count);

            if (edo_count == 2) {
                X[0] = X0;
                printf("Insert initial data Y01 = Y1(X0):\n");
                scanf("%lf", &Y1[0]);
                printf("Insert initial data Y02 = Y2(X0):\n");
                scanf("%lf", &Y2[0]);

                printf("Do you want to calculate convergence factor for Heun's method? (1.Yes 2.No)\n");
                scanf("%d", &conv_choice);

                if (conv_choice == 1) {
                    convergence_factor_heun_2(
                        n, h, X0, Xf, Y1[0], Y2[0], 
                        f12, f22, 
                        "convergence_heun_edo2.txt", 
                        "convergence_heun2_edo2.txt"
                    );
                }
            
                // Coefficients of prediction and correction
                double k11, k12, k21, k22;

                // Principal for of Heun's method for systems of 2 EDOs
                for (int i = 0; i < n; i++) {
                    double x_i = X[i];
                    double y1_i = Y1[i];
                    double y2_i = Y2[i];

                    // Prediction (Euler)
                    k11 = f12(x_i, y1_i, y2_i);
                    k12 = f22(x_i, y1_i, y2_i);
                    double y1_pred = y1_i + h * k11;
                    double y2_pred = y2_i + h * k12;
                
                    // Correction (Heun)
                    k21 = f12(x_i + h, y1_pred, y2_pred);
                    k22 = f22(x_i + h, y1_pred, y2_pred);
                
                    Y1[i+1] = y1_i + (h / 2.0) * (k11 + k21);
                    Y2[i+1] = y2_i + (h / 2.0) * (k12 + k22);
                    X[i+1] = x_i + h;
                }
            
                // Print results Y1
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y1", "Heun Y1", "Exact Error");
                printf("-------------------------------------------------------------------------------------------\n");
                
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y1(X[i]) - Y1[i]);
                    // double local_trunc_error = (pow(h, 3) / 12.0) * y3prima2(X[i], y1(X[i]), y2(X[i]), f12); 
                    printf("%-10d %-15lf %-15lf %-15lf %-15.3e\n", 
                           i, X[i], y1(X[i]), Y1[i], exact_error);
                }

                // Print results Y2
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y2", "Heun Y2", "Exact Error");
                printf("-------------------------------------------------------------------------------------------\n");
                
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y2(X[i]) - Y2[i]);
                    // double local_trunc_error = (pow(h, 3) / 12.0) * y3prima2(X[i], y1(X[i]), y2(X[i]), f22); 
                    printf("%-10d %-15lf %-15lf %-15lf %-15.3e\n", 
                           i, X[i], y2(X[i]), Y2[i], exact_error);
                }
            
                // Show final results
                printf("\nX[i]\t\tY1[i]\t\tY2[i]\n");
                for (int i = 0; i <= n; i++) {
                    printf("%lf\t%lf\t%lf\n", X[i], Y1[i], Y2[i]);
                }
            
                // Guardar resultados
                save_in_txt("results_Y1.txt", X, Y1, n);
                save_in_txt("results_Y2.txt", X, Y2, n);
            
                // Graficar con Python (opcional)
                if (system("test -f graph_points_edo2.py") == 0) {
                    system("python3 graph_points_edo2.py");
                } else {
                    printf("⚠️  Warning: 'graph_points_edo2.py' not found. Skipping graph generation.\n");
                }
            } else if(edo_count == 3) {
                X[0] = X0;
                printf("Insert initial data Y01 = Y1(X0):\n");
                scanf("%lf", &Y1[0]);
                printf("Insert initial data Y02 = Y2(X0):\n");
                scanf("%lf", &Y2[0]);
                printf("Insert initial data Y03 = Y3(X0):\n");
                scanf("%lf", &Y3[0]);

                printf("Do you want to calculate convergence factor for Heun's method? (1.Yes 2.No)\n");
                scanf("%d", &conv_choice);

                if (conv_choice == 1) {
                    convergence_factor_heun_3(
                        n, h, X0, Xf, Y1[0], Y2[0], Y3[0],
                        f13, f23, f33,
                        "convergence_heun_edo3.txt",
                        "convergence_heun2_edo3.txt",
                        "convergence_heun3_edo3.txt"
                    );
                }

                // Coefficients of prediction and correction
                double k11, k12, k13, k21, k22, k23;

                // Heun method for 3 EDOs
                for (int i = 0; i < n; i++) {
                    double x_i = X[i];
                    double y1_i = Y1[i];
                    double y2_i = Y2[i];
                    double y3_i = Y3[i];
                
                    // Prediction (Euler)
                    k11 = f13(x_i, y1_i, y2_i, y3_i);
                    k12 = f23(x_i, y1_i, y2_i, y3_i);
                    k13 = f33(x_i, y1_i, y2_i, y3_i);
                
                    double y1_pred = y1_i + h * k11;
                    double y2_pred = y2_i + h * k12;
                    double y3_pred = y3_i + h * k13;
                
                    // Correction (Heun)
                    k21 = f13(x_i + h, y1_pred, y2_pred, y3_pred);
                    k22 = f23(x_i + h, y1_pred, y2_pred, y3_pred);
                    k23 = f33(x_i + h, y1_pred, y2_pred, y3_pred);
                
                    Y1[i+1] = y1_i + (h / 2.0) * (k11 + k21);
                    Y2[i+1] = y2_i + (h / 2.0) * (k12 + k22);
                    Y3[i+1] = y3_i + (h / 2.0) * (k13 + k23);
                
                    X[i+1] = x_i + h;
                }

                // ----------- Print results for Y1 -----------
                printf("\n%-10s %-15s %-15s %-15s %-15s\n",
                       "i", "X[i]", "Exact Y1", "Heun Y1", "Exact Error");
                printf("-------------------------------------------------------------------------------------------\n");
                
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y1(X[i]) - Y1[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.3e\n",
                           i, X[i], y1(X[i]), Y1[i], exact_error);
                }

                // ----------- Print results for Y2 -----------
                printf("\n%-10s %-15s %-15s %-15s %-15s\n",
                       "i", "X[i]", "Exact Y2", "Heun Y2", "Exact Error");
                printf("-------------------------------------------------------------------------------------------\n");
                
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y2(X[i]) - Y2[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.3e\n",
                           i, X[i], y2(X[i]), Y2[i], exact_error);
                }

                // ----------- Print results for Y3 -----------
                printf("\n%-10s %-15s %-15s %-15s %-15s\n",
                       "i", "X[i]", "Exact Y3", "Heun Y3", "Exact Error");
                printf("-------------------------------------------------------------------------------------------\n");
                
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y3(X[i]) - Y3[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.3e\n",
                           i, X[i], y3(X[i]), Y3[i], exact_error);
                }

                // Save results
                save_in_txt("results_Y1.txt", X, Y1, n);
                save_in_txt("results_Y2.txt", X, Y2, n);
                save_in_txt("results_Y3.txt", X, Y3, n);

                // Optional graph
                if (system("test -f graph_points_edo3.py") == 0) {
                    system("python3 graph_points_edo3.py");
                } else {
                    printf("⚠️  Warning: 'graph_points_edo3.py' not found. Skipping graph generation.\n");
                }

            }
            break;
         case 3:
            printf("How many EDO's does it have your system? (2 or 3)\n");
            scanf("%d", &edo_count);

            if (edo_count == 2) {
                X[0] = X0;
                printf("Insert initial data Y01 = Y1(X0):\n");
                scanf("%lf", &Y1[0]);
                printf("Insert initial data Y02 = Y2(X0):\n");
                scanf("%lf", &Y2[0]);
    
                printf("Do you want to calculate convergence factor for Heun's method? (1.Yes 2.No)\n");
                scanf("%d", &conv_choice);
    
                if (conv_choice == 1) {
                    convergence_factor_midpoint_2(
                        n, h, X0, Xf, Y1[0], Y2[0],
                        f12, f22,
                        "convergence_midpoint_edo2.txt",
                        "convergence_midpoint2_edo2.txt"
                    );
                }
    
                // Midpoint method for systems of 2 EDOs
                for (int i = 0; i < n; i++) {
                    double x_i = X[i];
                    double y1_i = Y1[i];
                    double y2_i = Y2[i];
                
                    // Paso 1: Calcular k1
                    double k11 = f12(x_i, y1_i, y2_i);
                    double k12 = f22(x_i, y1_i, y2_i);
                
                    // Paso 2: Calcular k2 (en el punto medio)
                    double k21 = f12(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12);
                    double k22 = f22(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12);
                
                    // Paso 3: Actualizar
                    Y1[i+1] = y1_i + h * k21;
                    Y2[i+1] = y2_i + h * k22;
                    X[i+1] = x_i + h;
                }
            
                // Print results Y1
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y1", "Midpoint Y1", "Exact Error");
                printf("-------------------------------------------------------------\n");
                
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y1(X[i]) - Y1[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.3e\n", 
                           i, X[i], y1(X[i]), Y1[i], exact_error);
                }

                // Print results Y2
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y2", "Midpoint Y2", "Exact Error");
                printf("-------------------------------------------------------------\n");
                
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y2(X[i]) - Y2[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.3e\n", 
                           i, X[i], y2(X[i]), Y2[i], exact_error);
                }

                // Save results to files
                save_in_txt("results_Y1.txt", X, Y1, n);
                save_in_txt("results_Y2.txt", X, Y2, n);
    
                // Graph results using Python
                if (system("test -f graph_points_edo2.py") == 0) {
                    system("python3 graph_points_edo2.py");
                } else {
                    printf("⚠️  Warning: 'graph_points_edo2.py' not found. Skipping graph generation.\n");
                }
            } else if(edo_count == 3) {
                X[0] = X0;

                printf("Insert initial data Y01 = Y1(X0):\n");
                scanf("%lf", &Y1[0]);
                printf("Insert initial data Y02 = Y2(X0):\n");
                scanf("%lf", &Y2[0]);
                printf("Insert initial data Y03 = Y3(X0):\n");
                scanf("%lf", &Y3[0]);

                printf("Do you want to calculate convergence factor for Midpoint method? (1.Yes 2.No)\n");
                scanf("%d", &conv_choice);

                if (conv_choice == 1) {
                    convergence_factor_midpoint_3(
                        n, h, X0, Xf,
                        Y1[0], Y2[0], Y3[0],
                        f13, f23, f33,
                        "convergence_midpoint_edo3.txt",
                        "convergence_midpoint2_edo3.txt",
                        "convergence_midpoint3_edo3.txt"
                    );
                }

                // ---------- Midpoint method for systems of 3 EDOs ----------
                for (int i = 0; i < n; i++) {
                    double x_i = X[i];
                    double y1_i = Y1[i];
                    double y2_i = Y2[i];
                    double y3_i = Y3[i];
                
                    // Paso 1: Calcular k1
                    double k11 = f13(x_i, y1_i, y2_i, y3_i);
                    double k12 = f23(x_i, y1_i, y2_i, y3_i);
                    double k13 = f33(x_i, y1_i, y2_i, y3_i);
                
                    // Paso 2: Calcular k2 (en el punto medio)
                    double k21 = f13(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12, y3_i + (h/2.0)*k13);
                    double k22 = f23(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12, y3_i + (h/2.0)*k13);
                    double k23 = f33(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12, y3_i + (h/2.0)*k13);
                
                    // Paso 3: Actualizar
                    Y1[i+1] = y1_i + h * k21;
                    Y2[i+1] = y2_i + h * k22;
                    Y3[i+1] = y3_i + h * k23;
                    X[i+1] = x_i + h;
                }

                // ---------- Print results for Y1 ----------
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y1", "Midpoint Y1", "Exact Error");
                printf("-------------------------------------------------------------\n");
                
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y1(X[i]) - Y1[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.3e\n", 
                           i, X[i], y1(X[i]), Y1[i], exact_error);
                }

                // ---------- Print results for Y2 ----------
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y2", "Midpoint Y2", "Exact Error");
                printf("-------------------------------------------------------------\n");
                
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y2(X[i]) - Y2[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.3e\n", 
                           i, X[i], y2(X[i]), Y2[i], exact_error);
                }

                // ---------- Print results for Y3 ----------
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y3", "Midpoint Y3", "Exact Error");
                printf("-------------------------------------------------------------\n");
                
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y3(X[i]) - Y3[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.3e\n", 
                           i, X[i], y3(X[i]), Y3[i], exact_error);
                }

                // ---------- Save results to files ----------
                save_in_txt("results_Y1.txt", X, Y1, n);
                save_in_txt("results_Y2.txt", X, Y2, n);
                save_in_txt("results_Y3.txt", X, Y3, n);

                // ---------- Graph results using Python ----------
                if (system("test -f graph_points_edo3.py") == 0) {
                    system("python3 graph_points_edo3.py");
                } else {
                    printf("⚠️  Warning: 'graph_points_edo3.py' not found. Skipping graph generation.\n");
                }
            }
            break; 
        case 4:
            printf("How many EDO's does it have you system? (2, 3 or 4)\n");
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

                // Imprimir Y1
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y1", "RK4 Y1", "Exact Error", "Local Trunc. Err");
                printf("------------------------------------------------------------\n");
                for(int i = 0; i <= n; i++) {
                    double exact_error = fabs(y1(X[i]) - Y1[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.2e %-15.2e\n",
                           i, X[i], y1(X[i]), Y1[i], exact_error, lte1_arr[i]);
                }

                // Imprimir Y2
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y2", "RK4 Y2", "Exact Error", "Local Trunc. Err");
                printf("------------------------------------------------------------\n");
                for(int i = 0; i <= n; i++) {
                    double exact_error = fabs(y2(X[i]) - Y2[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.2e %-15.2e\n",
                           i, X[i], y2(X[i]), Y2[i], exact_error, lte2_arr[i]);
                }

                // Print results
                printf("X[i]\t\tY1[i]\t\tY2[i]\n");
                // Print all computed points including the last one
                for(int i = 0; i <= n; i++) {
                    printf("%lf\t%lf\t%lf\n", X[i], Y1[i], Y2[i]);
                }
            
                // Save x[i] and Y[i] in results.txt
                save_in_txt("results_Y1.txt", X, Y1, n);
                save_in_txt("results_Y2.txt", X, Y2, n);

                // Finally, we print the results.txt file in a graph using Python to visualize the results
                if (system("test -f graph_points_edo2.py") == 0) {
                    system("python3 graph_points_edo2.py");
                } else {
                    printf("Warning: 'graph_points_edo2.py' not found. Skipping graph generation.\n");
                }
            } else if(edo_count == 3) {
                X[0] = X0;
                printf("Insert initial data Y01 = Y1(X0):\n");
                scanf("%lf", &Y1[0]);
                printf("Insert initial data Y02 = Y2(X0):\n");
                scanf("%lf", &Y2[0]);
                printf("Insert initial data Y03 = Y3(X0):\n");
                scanf("%lf", &Y3[0]);

                printf("Do you want to calculate convergence factor for Runge-Kutta's 4 method? (1.Yes 2.No)\n");
                scanf("%d", &conv_choice);

                if (conv_choice == 1) {
                    convergence_factor_rk4_3(
                        n, h, X0, Y1[0], Y2[0], Y3[0],
                        f13, f23, f33,
                        "convergence_rk4_edo3.txt",
                        "convergence_rk42_edo3.txt",
                        "convergence_rk43_edo3.txt"
                    );
                }

                // Coefficients of Runge-Kutta 4
                double k11, k12, k13;
                double k21, k22, k23;
                double k31, k32, k33;
                double k41, k42, k43;

                // Principal for of integration
                for (int i = 0; i < n; i++) {
                    double x_i = X[i];
                    double y1_i = Y1[i];
                    double y2_i = Y2[i];
                    double y3_i = Y3[i];
                
                    // --- Step 1 1 ---
                    k11 = f13(x_i, y1_i, y2_i, y3_i);
                    k12 = f23(x_i, y1_i, y2_i, y3_i);
                    k13 = f33(x_i, y1_i, y2_i, y3_i);
                
                    // --- Step 2 ---
                    k21 = f13(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12, y3_i + (h/2.0)*k13);
                    k22 = f23(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12, y3_i + (h/2.0)*k13);
                    k23 = f33(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12, y3_i + (h/2.0)*k13);
                
                    // --- Step 3 ---
                    k31 = f13(x_i + h/2.0, y1_i + (h/2.0)*k21, y2_i + (h/2.0)*k22, y3_i + (h/2.0)*k23);
                    k32 = f23(x_i + h/2.0, y1_i + (h/2.0)*k21, y2_i + (h/2.0)*k22, y3_i + (h/2.0)*k23);
                    k33 = f33(x_i + h/2.0, y1_i + (h/2.0)*k21, y2_i + (h/2.0)*k22, y3_i + (h/2.0)*k23);
                
                    // --- Step 4 ---
                    k41 = f13(x_i + h, y1_i + h*k31, y2_i + h*k32, y3_i + h*k33);
                    k42 = f23(x_i + h, y1_i + h*k31, y2_i + h*k32, y3_i + h*k33);
                    k43 = f33(x_i + h, y1_i + h*k31, y2_i + h*k32, y3_i + h*k33);
                
                    // --- Update ---
                    Y1[i+1] = y1_i + (h/6.0) * (k11 + 2*k21 + 2*k31 + k41);
                    Y2[i+1] = y2_i + (h/6.0) * (k12 + 2*k22 + 2*k32 + k42);
                    Y3[i+1] = y3_i + (h/6.0) * (k13 + 2*k23 + 2*k33 + k43);
                    X[i+1] = x_i + h;
                }

                // ============================================================
                // == Calculus of truncated local error (LTE)
                // ============================================================
                double lte1_arr[MAX_SIZE+1], lte2_arr[MAX_SIZE+1], lte3_arr[MAX_SIZE+1];
                lte1_arr[0] = lte2_arr[0] = lte3_arr[0] = 0.0;

                for (int i = 1; i <= n; i++) {
                    double lte1, lte2, lte3;
                    local_trunc_error_rk4_3(X[i-1], Y1[i-1], Y2[i-1], Y3[i-1], h, f13, f23, f33, &lte1, &lte2, &lte3);
                    lte1_arr[i] = lte1;
                    lte2_arr[i] = lte2;
                    lte3_arr[i] = lte3;
                }

                // ============================================================
                // == Prints of results
                // ============================================================

                // Y1
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y1", "RK4 Y1", "Exact Error", "Local Trunc. Err");
                printf("--------------------------------------------------------------------------\n");
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y1(X[i]) - Y1[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.2e %-15.2e\n",
                           i, X[i], y1(X[i]), Y1[i], exact_error, lte1_arr[i]);
                }

                // Y2
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y2", "RK4 Y2", "Exact Error", "Local Trunc. Err");
                printf("--------------------------------------------------------------------------\n");
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y2(X[i]) - Y2[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.2e %-15.2e\n",
                           i, X[i], y2(X[i]), Y2[i], exact_error, lte2_arr[i]);
                }

                // Y3
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "Exact Y3", "RK4 Y3", "Exact Error", "Local Trunc. Err");
                printf("--------------------------------------------------------------------------\n");
                for (int i = 0; i <= n; i++) {
                    double exact_error = fabs(y3(X[i]) - Y3[i]);
                    printf("%-10d %-15lf %-15lf %-15lf %-15.2e %-15.2e\n",
                           i, X[i], y3(X[i]), Y3[i], exact_error, lte3_arr[i]);
                }

                // ============================================================
                // == Save results and plot
                // ============================================================
                save_in_txt("results_Y1.txt", X, Y1, n);
                save_in_txt("results_Y2.txt", X, Y2, n);
                save_in_txt("results_Y3.txt", X, Y3, n);

                if (system("test -f graph_points_edo3.py") == 0) {
                    system("python3 graph_points_edo3.py");
                } else {
                    printf("Warning: 'graph_points_edo3.py' not found. Skipping graph generation.\n");
                }
            } else if(edo_count == 4) {
                X[0] = X0;
                printf("Insert initial data Y01 = Y1(X0):\n");
                scanf("%lf", &Y1[0]);
                printf("Insert initial data Y02 = Y2(X0):\n");
                scanf("%lf", &Y2[0]);
                printf("Insert initial data Y03 = Y3(X0):\n");
                scanf("%lf", &Y3[0]);
                printf("Insert initial data Y04 = Y4(X0):\n");
                scanf("%lf", &Y4[0]);

                printf("Do you want to calculate convergence factor for Runge-Kutta's 4 method? (1.Yes 2.No)\n");
                scanf("%d", &conv_choice);

                if (conv_choice == 1) {
                    convergence_factor_rk4_4(
                        n, h, X0, Y1[0], Y2[0], Y3[0], Y4[0],
                        f14, f24, f34, f44,
                        "convergence_rk4_edo4.txt",
                        "convergence_rk42_edo4.txt", 
                        "convergence_rk43_edo4.txt",
                        "convergence_rk44_edo4.txt"
                    );
                }
            
                // Coefficients of Runge-Kutta 4 for 4 EDOs
                double k11, k12, k13, k14;
                double k21, k22, k23, k24;
                double k31, k32, k33, k34;
                double k41, k42, k43, k44;

                // ============================================================
                // == MODIFICACIÓN: Integración con criterio de parada
                // ============================================================
                const double R_earth = 6378140.0; // Radio terrestre en metros
                            
                int impact_index = n; // Por defecto, usa todo el intervalo
                int found_impact = 0;
                            
                // Principal for of integration WITH IMPACT DETECTION
                for (int i = 0; i < n; i++) {
                    double x_i = X[i];
                    double y1_i = Y1[i];
                    double y2_i = Y2[i];
                    double y3_i = Y3[i];
                    double y4_i = Y4[i];
                    
                    // --- Steps de RK4 (igual que antes) ---
                    k11 = f14(x_i, y1_i, y2_i, y3_i, y4_i);
                    k12 = f24(x_i, y1_i, y2_i, y3_i, y4_i);
                    k13 = f34(x_i, y1_i, y2_i, y3_i, y4_i);
                    k14 = f44(x_i, y1_i, y2_i, y3_i, y4_i);
                    
                    k21 = f14(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12, y3_i + (h/2.0)*k13, y4_i + (h/2.0)*k14);
                    k22 = f24(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12, y3_i + (h/2.0)*k13, y4_i + (h/2.0)*k14);
                    k23 = f34(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12, y3_i + (h/2.0)*k13, y4_i + (h/2.0)*k14);
                    k24 = f44(x_i + h/2.0, y1_i + (h/2.0)*k11, y2_i + (h/2.0)*k12, y3_i + (h/2.0)*k13, y4_i + (h/2.0)*k14);
                    
                    k31 = f14(x_i + h/2.0, y1_i + (h/2.0)*k21, y2_i + (h/2.0)*k22, y3_i + (h/2.0)*k23, y4_i + (h/2.0)*k24);
                    k32 = f24(x_i + h/2.0, y1_i + (h/2.0)*k21, y2_i + (h/2.0)*k22, y3_i + (h/2.0)*k23, y4_i + (h/2.0)*k24);
                    k33 = f34(x_i + h/2.0, y1_i + (h/2.0)*k21, y2_i + (h/2.0)*k22, y3_i + (h/2.0)*k23, y4_i + (h/2.0)*k24);
                    k34 = f44(x_i + h/2.0, y1_i + (h/2.0)*k21, y2_i + (h/2.0)*k22, y3_i + (h/2.0)*k23, y4_i + (h/2.0)*k24);
                    
                    k41 = f14(x_i + h, y1_i + h*k31, y2_i + h*k32, y3_i + h*k33, y4_i + h*k34);
                    k42 = f24(x_i + h, y1_i + h*k31, y2_i + h*k32, y3_i + h*k33, y4_i + h*k34);
                    k43 = f34(x_i + h, y1_i + h*k31, y2_i + h*k32, y3_i + h*k33, y4_i + h*k34);
                    k44 = f44(x_i + h, y1_i + h*k31, y2_i + h*k32, y3_i + h*k33, y4_i + h*k34);
                    
                    Y1[i+1] = y1_i + (h/6.0) * (k11 + 2*k21 + 2*k31 + k41);
                    Y2[i+1] = y2_i + (h/6.0) * (k12 + 2*k22 + 2*k32 + k42);
                    Y3[i+1] = y3_i + (h/6.0) * (k13 + 2*k23 + 2*k33 + k43);
                    Y4[i+1] = y4_i + (h/6.0) * (k14 + 2*k24 + 2*k34 + k44);
                    X[i+1] = x_i + h;
                    
                    // ===== DETECCIÓN DE IMPACTO =====
                    if (Y1[i+1] <= R_earth && !found_impact) {
                        impact_index = i + 1;
                        found_impact = 1;
                        printf("*** IMPACTO DETECTADO en t = %.2f s, θ = %.6f rad ***\n", X[i+1], Y3[i+1]);
                        break; // Salir del loop - ya impactó
                    }
                }

                // ============================================================
                // == Calculus of truncated local error (LTE)
                // ============================================================
                double lte1_arr[MAX_SIZE+1], lte2_arr[MAX_SIZE+1], lte3_arr[MAX_SIZE+1], lte4_arr[MAX_SIZE+1];
                lte1_arr[0] = lte2_arr[0] = lte3_arr[0] = lte4_arr[0] = 0.0;
            
                for (int i = 1; i <= n; i++) {
                    double lte1, lte2, lte3, lte4;
                    local_trunc_error_rk4_4(X[i-1], Y1[i-1], Y2[i-1], Y3[i-1], Y4[i-1], h, 
                                           f14, f24, f34, f44, &lte1, &lte2, &lte3, &lte4);
                    lte1_arr[i] = lte1;
                    lte2_arr[i] = lte2;
                    lte3_arr[i] = lte3;
                    lte4_arr[i] = lte4;
                }
                
                // ============================================================
                // == IMPRIMIR SOLO HASTA EL IMPACTO
                // ============================================================
                int n_effective = found_impact ? impact_index : n;
                
                // Y1 - Radio
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "Tiempo [s]", "r [m]", "Error Local", "Descripción");
                printf("-------------------------------------------------------------------\n");
                for (int i = 0; i <= n_effective; i++) {
                    printf("%-10d %-15.1f %-15.0f %-15.2e %-15s\n",
                           i, X[i], Y1[i], lte1_arr[i], "Radio");
                }
                
                // Y3 - Ángulo θ (EL QUE TE PIDEN)
                printf("\n%-10s %-15s %-15s %-15s %-15s\n", 
                       "i", "Tiempo [s]", "θ [rad]", "Error Local", "Descripción");
                printf("-------------------------------------------------------------------\n");
                for (int i = 0; i <= n_effective; i++) {
                    printf("%-10d %-15.1f %-15.6f %-15.2e %-15s\n",
                           i, X[i], Y3[i], lte3_arr[i], "Ángulo");
                }
                
                // ============================================================
                // == RESULTADO FINAL - θ EN EL IMPACTO
                // ============================================================
                printf("\n");
                printf("===============================================\n");
                printf("           RESULTADO FINAL - IMPACTO\n");
                printf("===============================================\n");
                printf("Tiempo de impacto:  %.2f s (%.2f min)\n", X[n_effective], X[n_effective]/60.0);
                printf("Ángulo en impacto:  θ = %.6f rad\n", Y3[n_effective]);
                printf("Ángulo en impacto:  θ = %.2f°\n", Y3[n_effective] * 180.0/M_PI);
                printf("Radio en impacto:   r = %.0f m\n", Y1[n_effective]);
                printf("Velocidad radial:   vr = %.2f m/s\n", Y2[n_effective]);
                printf("Velocidad angular:  ω = %.6f rad/s\n", Y4[n_effective]);
                printf("===============================================\n");
                
                // Guardar resultados solo hasta el impacto
                save_in_txt("results_Y1.txt", X, Y1, n_effective);
                save_in_txt("results_Y2.txt", X, Y2, n_effective);
                save_in_txt("results_Y3.txt", X, Y3, n_effective);
                save_in_txt("results_Y4.txt", X, Y4, n_effective);
            
                
            
                /* // ============================================================
                // == Prints of results
                // ============================================================
            
                // Y1
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "RK4 Y1", "Exact Error", "Local Trunc. Err", "Description");
                printf("------------------------------------------------------------------------------------------\n");
                for (int i = 0; i <= n; i++) {
                    // Para el problema de la nave no hay solución exacta conocida
                    double exact_error = 0.0; // fabs(y1(X[i]) - Y1[i]);
                    printf("%-10d %-15lf %-15lf %-15.2e %-15.2e %-15s\n",
                           i, X[i], Y1[i], exact_error, lte1_arr[i], "r (m)");
                }
            
                // Y2
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "RK4 Y2", "Exact Error", "Local Trunc. Err", "Description");
                printf("------------------------------------------------------------------------------------------\n");
                for (int i = 0; i <= n; i++) {
                    double exact_error = 0.0;
                    printf("%-10d %-15lf %-15lf %-15.2e %-15.2e %-15s\n",
                           i, X[i], Y2[i], exact_error, lte2_arr[i], "vr (m/s)");
                }
            
                // Y3
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "RK4 Y3", "Exact Error", "Local Trunc. Err", "Description");
                printf("------------------------------------------------------------------------------------------\n");
                for (int i = 0; i <= n; i++) {
                    double exact_error = 0.0;
                    printf("%-10d %-15lf %-15lf %-15.2e %-15.2e %-15s\n",
                           i, X[i], Y3[i], exact_error, lte3_arr[i], "θ (rad)");
                }
            
                // Y4
                printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", 
                       "i", "X[i]", "RK4 Y4", "Exact Error", "Local Trunc. Err", "Description");
                printf("------------------------------------------------------------------------------------------\n");
                for (int i = 0; i <= n; i++) {
                    double exact_error = 0.0;
                    printf("%-10d %-15lf %-15lf %-15.2e %-15.2e %-15s\n",
                           i, X[i], Y4[i], exact_error, lte4_arr[i], "ω (rad/s)");
                }
            
                // ============================================================
                // == Save results and plot
                // ============================================================
                save_in_txt("results_Y1.txt", X, Y1, n);
                save_in_txt("results_Y2.txt", X, Y2, n);
                save_in_txt("results_Y3.txt", X, Y3, n);
                save_in_txt("results_Y4.txt", X, Y4, n);
            
                printf("\nFinal values:\n");
                printf("r(tf) = %lf m\n", Y1[n]);
                printf("vr(tf) = %lf m/s\n", Y2[n]);
                printf("θ(tf) = %lf rad\n", Y3[n]);
                printf("ω(tf) = %lf rad/s\n", Y4[n]);
            
                if (system("test -f graph_points_edo4.py") == 0) {
                    system("python3 graph_points_edo4.py");
                } else {
                    printf("Warning: 'graph_points_edo4.py' not found. Skipping graph generation.\n");
                }
            } */
            break; 
    }
    return 0;
}
}

// Functions for systems of two EDOs
double f12(double X, double Y1, double Y2) {
    return 3 * X + Y2;
}

double f22(double X, double Y1, double Y2) {
    return pow(X, 2) - Y1 - 1;
}

// Functions for systems of three EDOs
double f13(double X, double Y1, double Y2, double Y3) {
    return X + Y1 + Y2 + Y3;
}

double f23(double X, double Y1, double Y2, double Y3) {
    return X - Y1 + Y2 + Y3;
}

double f33(double X, double Y1, double Y2, double Y3) {
    return X + Y2 - Y3 + Y1;
}

// There's no exact solution
double y1(double x) {
    return 0;
}

double y2(double x) {
    return 0;
}

double y3(double x) {
    return 0;
}

double y4(double x) {
    return 0;
}

double f14(double x, double y1, double y2, double y3, double y4) {
    return y2;  // y1' = y2 (dr/dt = vr)
}

double f24(double x, double y1, double y2, double y3, double y4) {
    return y1 * y4 * y4 - GM / (y1 * y1);  // y2' = r*θ'² - GM/r²
}

double f34(double x, double y1, double y2, double y3, double y4) {
    return y4;  // y3' = y4 (dθ/dt = ω)
}

double f44(double x, double y1, double y2, double y3, double y4) {
    if (fabs(y1) < 1e-9) return 0.0;  // Evitar división por cero
    return -2.0 * y2 * y4 / y1;  // y4' = -2*vr*ω/r
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

void rk4_step4(double x, double *y1, double *y2, double *y3, double *y4, double h,
               double (*f1)(double, double, double, double, double),
               double (*f2)(double, double, double, double, double),
               double (*f3)(double, double, double, double, double),
               double (*f4)(double, double, double, double, double)) {
    
    double k11 = f1(x, *y1, *y2, *y3, *y4);
    double k12 = f2(x, *y1, *y2, *y3, *y4);
    double k13 = f3(x, *y1, *y2, *y3, *y4);
    double k14 = f4(x, *y1, *y2, *y3, *y4);
    
    double k21 = f1(x + h/2.0, *y1 + (h/2.0)*k11, *y2 + (h/2.0)*k12, *y3 + (h/2.0)*k13, *y4 + (h/2.0)*k14);
    double k22 = f2(x + h/2.0, *y1 + (h/2.0)*k11, *y2 + (h/2.0)*k12, *y3 + (h/2.0)*k13, *y4 + (h/2.0)*k14);
    double k23 = f3(x + h/2.0, *y1 + (h/2.0)*k11, *y2 + (h/2.0)*k12, *y3 + (h/2.0)*k13, *y4 + (h/2.0)*k14);
    double k24 = f4(x + h/2.0, *y1 + (h/2.0)*k11, *y2 + (h/2.0)*k12, *y3 + (h/2.0)*k13, *y4 + (h/2.0)*k14);
    
    double k31 = f1(x + h/2.0, *y1 + (h/2.0)*k21, *y2 + (h/2.0)*k22, *y3 + (h/2.0)*k23, *y4 + (h/2.0)*k24);
    double k32 = f2(x + h/2.0, *y1 + (h/2.0)*k21, *y2 + (h/2.0)*k22, *y3 + (h/2.0)*k23, *y4 + (h/2.0)*k24);
    double k33 = f3(x + h/2.0, *y1 + (h/2.0)*k21, *y2 + (h/2.0)*k22, *y3 + (h/2.0)*k23, *y4 + (h/2.0)*k24);
    double k34 = f4(x + h/2.0, *y1 + (h/2.0)*k21, *y2 + (h/2.0)*k22, *y3 + (h/2.0)*k23, *y4 + (h/2.0)*k24);
    
    double k41 = f1(x + h, *y1 + h*k31, *y2 + h*k32, *y3 + h*k33, *y4 + h*k34);
    double k42 = f2(x + h, *y1 + h*k31, *y2 + h*k32, *y3 + h*k33, *y4 + h*k34);
    double k43 = f3(x + h, *y1 + h*k31, *y2 + h*k32, *y3 + h*k33, *y4 + h*k34);
    double k44 = f4(x + h, *y1 + h*k31, *y2 + h*k32, *y3 + h*k33, *y4 + h*k34);
    
    *y1 += h/6.0 * (k11 + 2*k21 + 2*k31 + k41);
    *y2 += h/6.0 * (k12 + 2*k22 + 2*k32 + k42);
    *y3 += h/6.0 * (k13 + 2*k23 + 2*k33 + k43);
    *y4 += h/6.0 * (k14 + 2*k24 + 2*k34 + k44);
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

void local_trunc_error_rk4_4(double x, double y1, double y2, double y3, double y4_val, double h,
                             double (*f1)(double, double, double, double, double),
                             double (*f2)(double, double, double, double, double),
                             double (*f3)(double, double, double, double, double),
                             double (*f4)(double, double, double, double, double),
                             double *lte1, double *lte2, double *lte3, double *lte4) {
    double Y_full[4] = {y1, y2, y3, y4_val};
    double mid[4] = {y1, y2, y3, y4_val};

    rk4_step4(x, Y_full, Y_full+1, Y_full+2, Y_full+3, h, f1, f2, f3, f4);

    double temp[4] = {y1, y2, y3, y4_val};
    rk4_step4(x, temp, temp+1, temp+2, temp+3, h/2.0, f1, f2, f3, f4);
    rk4_step4(x + h/2.0, temp, temp+1, temp+2, temp+3, h/2.0, f1, f2, f3, f4);

    *lte1 = fabs(temp[0] - Y_full[0]) / 15.0;
    *lte2 = fabs(temp[1] - Y_full[1]) / 15.0;
    *lte3 = fabs(temp[2] - Y_full[2]) / 15.0;
    *lte4 = fabs(temp[3] - Y_full[3]) / 15.0;
}

void convergence_factor_rk4_4(int n1, double h1, double X0,
                              double Y10, double Y20, double Y30, double Y40,
                              double (*f1)(double, double, double, double, double),
                              double (*f2)(double, double, double, double, double),
                              double (*f3)(double, double, double, double, double),
                              double (*f4)(double, double, double, double, double),
                              const char *out1, const char *out2, const char *out3, const char *out4) {
    double h2 = h1/2.0;
    double h3 = h1/4.0;

    double Y1h[MAX_SIZE+1], Y2h[MAX_SIZE+1], Y3h[MAX_SIZE+1], Y4h[MAX_SIZE+1];
    double Y1h2[MAX_SIZE*2+1], Y2h2[MAX_SIZE*2+1], Y3h2[MAX_SIZE*2+1], Y4h2[MAX_SIZE*2+1];
    double Y1h4[MAX_SIZE*4+1], Y2h4[MAX_SIZE*4+1], Y3h4[MAX_SIZE*4+1], Y4h4[MAX_SIZE*4+1];
    double Xh[MAX_SIZE+1], Xh2[MAX_SIZE*2+1], Xh4[MAX_SIZE*4+1];
    double Q1[MAX_SIZE+1], Q2[MAX_SIZE+1], Q3[MAX_SIZE+1], Q4[MAX_SIZE+1];

    // Inicialización
    Xh[0]=Xh2[0]=Xh4[0]=X0;
    Y1h[0]=Y1h2[0]=Y1h4[0]=Y10;
    Y2h[0]=Y2h2[0]=Y2h4[0]=Y20;
    Y3h[0]=Y3h2[0]=Y3h4[0]=Y30;
    Y4h[0]=Y4h2[0]=Y4h4[0]=Y40;

    // RK4 paso h
    for(int i=0;i<n1;i++){
        rk4_step4(Xh[i], &Y1h[i], &Y2h[i], &Y3h[i], &Y4h[i], h1, f1, f2, f3, f4);
        Xh[i+1] = Xh[i]+h1;
    }
    // RK4 paso h/2
    for(int i=0;i<2*n1;i++){
        rk4_step4(Xh2[i], &Y1h2[i], &Y2h2[i], &Y3h2[i], &Y4h2[i], h2, f1, f2, f3, f4);
        Xh2[i+1] = Xh2[i]+h2;
    }
    // RK4 paso h/4
    for(int i=0;i<4*n1;i++){
        rk4_step4(Xh4[i], &Y1h4[i], &Y2h4[i], &Y3h4[i], &Y4h4[i], h3, f1, f2, f3, f4);
        Xh4[i+1] = Xh4[i]+h3;
    }

    printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n","i","x_i","Q1_i","Q2_i","Q3_i","Q4_i");
    printf("----------------------------------------------------------------------------\n");

    for(int i=1;i<=n1;i++){
        double xi = X0 + i*h1;
        int idx2 = (int)((xi - X0)/h2 + 0.5);
        int idx4 = (int)((xi - X0)/h3 + 0.5);

        double num1 = fabs(Y1h[i]-Y1h2[idx2]);
        double den1 = fabs(Y1h2[idx2]-Y1h4[idx4]);
        double num2 = fabs(Y2h[i]-Y2h2[idx2]);
        double den2 = fabs(Y2h2[idx2]-Y2h4[idx4]);
        double num3 = fabs(Y3h[i]-Y3h2[idx2]);
        double den3 = fabs(Y3h2[idx2]-Y3h4[idx4]);
        double num4 = fabs(Y4h[i]-Y4h2[idx2]);
        double den4 = fabs(Y4h2[idx2]-Y4h4[idx4]);

        Q1[i] = (den1>1e-12)? log(num1/den1)/log(2.0) : 0.0;
        Q2[i] = (den2>1e-12)? log(num2/den2)/log(2.0) : 0.0;
        Q3[i] = (den3>1e-12)? log(num3/den3)/log(2.0) : 0.0;
        Q4[i] = (den4>1e-12)? log(num4/den4)/log(2.0) : 0.0;

        printf("%-10d %-15lf %-15.6f %-15.6f %-15.6f %-15.6f\n", 
               i, xi, Q1[i], Q2[i], Q3[i], Q4[i]);
    }

    save_in_txt(out1, Xh, Q1, n1);
    save_in_txt(out2, Xh, Q2, n1);
    save_in_txt(out3, Xh, Q3, n1);
    save_in_txt(out4, Xh, Q4, n1);

    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    } else {
        printf("⚠️  Warning: 'graph_convergence.py' not found. Skipping graph generation.\n");
    }
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

void local_trunc_error_rk4_3(double x, double y1, double y2, double y3, double h,
                             double (*f1)(double, double, double, double),
                             double (*f2)(double, double, double, double),
                             double (*f3)(double, double, double, double),
                             double *lte1, double *lte2, double *lte3) {
    double Y_full[3] = {y1, y2, y3};
    double mid[3] = {y1, y2, y3};

    rk4_step3(x, Y_full, Y_full+1, Y_full+2, h, f1, f2, f3);

    double temp[3] = {y1, y2, y3};
    rk4_step3(x, temp, temp+1, temp+2, h/2.0, f1, f2, f3);
    rk4_step3(x + h/2.0, temp, temp+1, temp+2, h/2.0, f1, f2, f3);

    *lte1 = fabs(temp[0] - Y_full[0]) / 15.0;
    *lte2 = fabs(temp[1] - Y_full[1]) / 15.0;
    *lte3 = fabs(temp[2] - Y_full[2]) / 15.0;
}

void save_in_txt(const char *filename, double X[], double Y[], int n) {
    FILE *archivo = fopen(filename, "w");
    if (archivo == NULL) {
        printf("Error: Unable to create file '%s'.\n", filename);
        exit(1);
    }

    for (int i = 0; i <= n; i++) {
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

void convergence_factor_euler_3(int n1, double h1, double X0, double Xf,
                                double Y10, double Y20, double Y30,
                                double (*f1)(double, double, double, double),
                                double (*f2)(double, double, double, double),
                                double (*f3)(double, double, double, double),
                                const char *filename1,
                                const char *filename2,
                                const char *filename3) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    // Arrays
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE * 2 + 1], Xh4[MAX_SIZE * 4 + 1];
    double Y1h[MAX_SIZE + 1], Y2h[MAX_SIZE + 1], Y3h[MAX_SIZE + 1];
    double Y1h2[MAX_SIZE * 2 + 1], Y2h2[MAX_SIZE * 2 + 1], Y3h2[MAX_SIZE * 2 + 1];
    double Y1h4[MAX_SIZE * 4 + 1], Y2h4[MAX_SIZE * 4 + 1], Y3h4[MAX_SIZE * 4 + 1];
    double Q1[MAX_SIZE + 1], Q2[MAX_SIZE + 1], Q3[MAX_SIZE + 1];

    // Inicialización
    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Y1h[0] = Y1h2[0] = Y1h4[0] = Y10;
    Y2h[0] = Y2h2[0] = Y2h4[0] = Y20;
    Y3h[0] = Y3h2[0] = Y3h4[0] = Y30;

    // ---------- Euler con paso h ----------
    for (int i = 0; i < n1; i++) {
        Xh[i+1] = Xh[i] + h1;
        Y1h[i+1] = Y1h[i] + h1 * f1(Xh[i], Y1h[i], Y2h[i], Y3h[i]);
        Y2h[i+1] = Y2h[i] + h1 * f2(Xh[i], Y1h[i], Y2h[i], Y3h[i]);
        Y3h[i+1] = Y3h[i] + h1 * f3(Xh[i], Y1h[i], Y2h[i], Y3h[i]);
    }

    // ---------- Euler con paso h/2 ----------
    for (int i = 0; i < 2*n1; i++) {
        Xh2[i+1] = Xh2[i] + h2;
        Y1h2[i+1] = Y1h2[i] + h2 * f1(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i]);
        Y2h2[i+1] = Y2h2[i] + h2 * f2(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i]);
        Y3h2[i+1] = Y3h2[i] + h2 * f3(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i]);
    }

    // ---------- Euler con paso h/4 ----------
    for (int i = 0; i < 4*n1; i++) {
        Xh4[i+1] = Xh4[i] + h3;
        Y1h4[i+1] = Y1h4[i] + h3 * f1(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i]);
        Y2h4[i+1] = Y2h4[i] + h3 * f2(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i]);
        Y3h4[i+1] = Y3h4[i] + h3 * f3(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i]);
    }

    // ---------- Calcular factores Q ----------
    printf("\n%-10s %-15s %-15s %-15s %-15s\n", "i", "x_i", "Q1_i", "Q2_i", "Q3_i");
    printf("--------------------------------------------------------------------------\n");

    Q1[0] = Q2[0] = Q3[0] = 0.0;

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;
        int idx4 = 4*i;

        double num1 = fabs(Y1h[i] - Y1h2[idx2]);
        double den1 = fabs(Y1h2[idx2] - Y1h4[idx4]);
        double num2 = fabs(Y2h[i] - Y2h2[idx2]);
        double den2 = fabs(Y2h2[idx2] - Y2h4[idx4]);
        double num3 = fabs(Y3h[i] - Y3h2[idx2]);
        double den3 = fabs(Y3h2[idx2] - Y3h4[idx4]);

        Q1[i] = (den1 > 1e-12) ? log(num1 / den1) / log(2.0) : 0.0;
        Q2[i] = (den2 > 1e-12) ? log(num2 / den2) / log(2.0) : 0.0;
        Q3[i] = (den3 > 1e-12) ? log(num3 / den3) / log(2.0) : 0.0;

        printf("%-10d %-15lf %-15lf %-15lf %-15lf\n", i, Xh[i], Q1[i], Q2[i], Q3[i]);
    }

    // Guardar resultados
    save_in_txt("results_Q1.txt", Xh, Q1, n1);
    save_in_txt("results_Q2.txt", Xh, Q2, n1);
    save_in_txt("results_Q3.txt", Xh, Q3, n1);

    rename("results_Q1.txt", filename1);
    rename("results_Q2.txt", filename2);
    rename("results_Q3.txt", filename3);

    // Generar gráfico si existe script Python
    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    } else {
        printf("⚠️  Warning: 'graph_convergence.py' not found. Skipping graph generation.\n");
    }
}

void convergence_factor_euler_4(int n1, double h1, double X0, double Xf,
                                double Y10, double Y20, double Y30, double Y40,
                                double (*f1)(double, double, double, double, double),
                                double (*f2)(double, double, double, double, double),
                                double (*f3)(double, double, double, double, double),
                                double (*f4)(double, double, double, double, double),
                                const char *filename1, const char *filename2,
                                const char *filename3, const char *filename4) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    // Arrays
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE * 2 + 1], Xh4[MAX_SIZE * 4 + 1];
    double Y1h[MAX_SIZE + 1], Y2h[MAX_SIZE + 1], Y3h[MAX_SIZE + 1], Y4h[MAX_SIZE + 1];
    double Y1h2[MAX_SIZE * 2 + 1], Y2h2[MAX_SIZE * 2 + 1], Y3h2[MAX_SIZE * 2 + 1], Y4h2[MAX_SIZE * 2 + 1];
    double Y1h4[MAX_SIZE * 4 + 1], Y2h4[MAX_SIZE * 4 + 1], Y3h4[MAX_SIZE * 4 + 1], Y4h4[MAX_SIZE * 4 + 1];
    double Q1[MAX_SIZE + 1], Q2[MAX_SIZE + 1], Q3[MAX_SIZE + 1], Q4[MAX_SIZE + 1];

    // Inicialización
    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Y1h[0] = Y1h2[0] = Y1h4[0] = Y10;
    Y2h[0] = Y2h2[0] = Y2h4[0] = Y20;
    Y3h[0] = Y3h2[0] = Y3h4[0] = Y30;
    Y4h[0] = Y4h2[0] = Y4h4[0] = Y40;

    // ---------- Euler con paso h ----------
    for (int i = 0; i < n1; i++) {
        Xh[i+1] = Xh[i] + h1;
        Y1h[i+1] = Y1h[i] + h1 * f1(Xh[i], Y1h[i], Y2h[i], Y3h[i], Y4h[i]);
        Y2h[i+1] = Y2h[i] + h1 * f2(Xh[i], Y1h[i], Y2h[i], Y3h[i], Y4h[i]);
        Y3h[i+1] = Y3h[i] + h1 * f3(Xh[i], Y1h[i], Y2h[i], Y3h[i], Y4h[i]);
        Y4h[i+1] = Y4h[i] + h1 * f4(Xh[i], Y1h[i], Y2h[i], Y3h[i], Y4h[i]);
    }

    // ---------- Euler con paso h/2 ----------
    for (int i = 0; i < 2*n1; i++) {
        Xh2[i+1] = Xh2[i] + h2;
        Y1h2[i+1] = Y1h2[i] + h2 * f1(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i], Y4h2[i]);
        Y2h2[i+1] = Y2h2[i] + h2 * f2(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i], Y4h2[i]);
        Y3h2[i+1] = Y3h2[i] + h2 * f3(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i], Y4h2[i]);
        Y4h2[i+1] = Y4h2[i] + h2 * f4(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i], Y4h2[i]);
    }

    // ---------- Euler con paso h/4 ----------
    for (int i = 0; i < 4*n1; i++) {
        Xh4[i+1] = Xh4[i] + h3;
        Y1h4[i+1] = Y1h4[i] + h3 * f1(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i], Y4h4[i]);
        Y2h4[i+1] = Y2h4[i] + h3 * f2(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i], Y4h4[i]);
        Y3h4[i+1] = Y3h4[i] + h3 * f3(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i], Y4h4[i]);
        Y4h4[i+1] = Y4h4[i] + h3 * f4(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i], Y4h4[i]);
    }

    // ---------- Calcular factores Q ----------
    printf("\n%-10s %-15s %-15s %-15s %-15s %-15s\n", "i", "x_i", "Q1_i", "Q2_i", "Q3_i", "Q4_i");
    printf("--------------------------------------------------------------------------------------\n");

    Q1[0] = Q2[0] = Q3[0] = Q4[0] = 0.0;

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;
        int idx4 = 4*i;

        double num1 = fabs(Y1h[i] - Y1h2[idx2]);
        double den1 = fabs(Y1h2[idx2] - Y1h4[idx4]);
        double num2 = fabs(Y2h[i] - Y2h2[idx2]);
        double den2 = fabs(Y2h2[idx2] - Y2h4[idx4]);
        double num3 = fabs(Y3h[i] - Y3h2[idx2]);
        double den3 = fabs(Y3h2[idx2] - Y3h4[idx4]);
        double num4 = fabs(Y4h[i] - Y4h2[idx2]);
        double den4 = fabs(Y4h2[idx2] - Y4h4[idx4]);

        Q1[i] = (den1 > 1e-12) ? log(num1 / den1) / log(2.0) : 0.0;
        Q2[i] = (den2 > 1e-12) ? log(num2 / den2) / log(2.0) : 0.0;
        Q3[i] = (den3 > 1e-12) ? log(num3 / den3) / log(2.0) : 0.0;
        Q4[i] = (den4 > 1e-12) ? log(num4 / den4) / log(2.0) : 0.0;

        printf("%-10d %-15lf %-15lf %-15lf %-15lf %-15lf\n", i, Xh[i], Q1[i], Q2[i], Q3[i], Q4[i]);
    }

    // Guardar resultados
    save_in_txt("results_Q1.txt", Xh, Q1, n1);
    save_in_txt("results_Q2.txt", Xh, Q2, n1);
    save_in_txt("results_Q3.txt", Xh, Q3, n1);
    save_in_txt("results_Q4.txt", Xh, Q4, n1);

    rename("results_Q1.txt", filename1);
    rename("results_Q2.txt", filename2);
    rename("results_Q3.txt", filename3);
    rename("results_Q4.txt", filename4);

    // Generar gráfico si existe script Python
    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    } else {
        printf("⚠️  Warning: 'graph_convergence.py' not found. Skipping graph generation.\n");
    }
}

void convergence_factor_heun_2(int n1, double h1, double X0, double Xf,
                               double Y10, double Y20,
                               double (*f1)(double, double, double),
                               double (*f2)(double, double, double),
                               const char *filename1,
                               const char *filename2) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    // Array values
    double Y1h[MAX_SIZE + 1], Y2h[MAX_SIZE + 1];
    double Y1h2[MAX_SIZE * 2 + 1], Y2h2[MAX_SIZE * 2 + 1];
    double Y1h4[MAX_SIZE * 4 + 1], Y2h4[MAX_SIZE * 4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE * 2 + 1], Xh4[MAX_SIZE * 4 + 1];
    double Q1[MAX_SIZE + 1], Q2[MAX_SIZE + 1];

    // Initialization
    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Y1h[0] = Y1h2[0] = Y1h4[0] = Y10;
    Y2h[0] = Y2h2[0] = Y2h4[0] = Y20;

    // ---------- Heun with step h ----------
    for (int i = 0; i < n1; i++) {
        double k1_y1 = f1(Xh[i], Y1h[i], Y2h[i]);
        double k1_y2 = f2(Xh[i], Y1h[i], Y2h[i]);

        double predictor_y1 = Y1h[i] + h1 * k1_y1;
        double predictor_y2 = Y2h[i] + h1 * k1_y2;

        double k2_y1 = f1(Xh[i] + h1, predictor_y1, predictor_y2);
        double k2_y2 = f2(Xh[i] + h1, predictor_y1, predictor_y2);

        Y1h[i+1] = Y1h[i] + (h1 / 2.0) * (k1_y1 + k2_y1);
        Y2h[i+1] = Y2h[i] + (h1 / 2.0) * (k1_y2 + k2_y2);
        Xh[i+1] = Xh[i] + h1;
    }

    // ---------- Heun with step h/2 ----------
    for (int i = 0; i < 2 * n1; i++) {
        double k1_y1 = f1(Xh2[i], Y1h2[i], Y2h2[i]);
        double k1_y2 = f2(Xh2[i], Y1h2[i], Y2h2[i]);

        double predictor_y1 = Y1h2[i] + h2 * k1_y1;
        double predictor_y2 = Y2h2[i] + h2 * k1_y2;

        double k2_y1 = f1(Xh2[i] + h2, predictor_y1, predictor_y2);
        double k2_y2 = f2(Xh2[i] + h2, predictor_y1, predictor_y2);

        Y1h2[i+1] = Y1h2[i] + (h2 / 2.0) * (k1_y1 + k2_y1);
        Y2h2[i+1] = Y2h2[i] + (h2 / 2.0) * (k1_y2 + k2_y2);
        Xh2[i+1] = Xh2[i] + h2;
    }

    // ---------- Heun with step h/4 ----------
    for (int i = 0; i < 4 * n1; i++) {
        double k1_y1 = f1(Xh4[i], Y1h4[i], Y2h4[i]);
        double k1_y2 = f2(Xh4[i], Y1h4[i], Y2h4[i]);

        double predictor_y1 = Y1h4[i] + h3 * k1_y1;
        double predictor_y2 = Y2h4[i] + h3 * k1_y2;

        double k2_y1 = f1(Xh4[i] + h3, predictor_y1, predictor_y2);
        double k2_y2 = f2(Xh4[i] + h3, predictor_y1, predictor_y2);

        Y1h4[i+1] = Y1h4[i] + (h3 / 2.0) * (k1_y1 + k2_y1);
        Y2h4[i+1] = Y2h4[i] + (h3 / 2.0) * (k1_y2 + k2_y2);
        Xh4[i+1] = Xh4[i] + h3;
    }

    // ---------- Calculate convergence factors Q ----------
    printf("\n%-10s %-15s %-15s %-15s\n", "i", "x_i", "Q1_i", "Q2_i");
    printf("----------------------------------------------------------\n");

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2 * i;
        int idx4 = 4 * i;

        double num1 = fabs(Y1h[i] - Y1h2[idx2]);
        double den1 = fabs(Y1h2[idx2] - Y1h4[idx4]);
        double num2 = fabs(Y2h[i] - Y2h2[idx2]);
        double den2 = fabs(Y2h2[idx2] - Y2h4[idx4]);

        Q1[i] = (den1 > 1e-12) ? log(num1 / den1) / log(2.0) : 0.0;
        Q2[i] = (den2 > 1e-12) ? log(num2 / den2) / log(2.0) : 0.0;

        printf("%-10d %-15lf %-15lf %-15lf\n", i, Xh[i], Q1[i], Q2[i]);
    }

    // ---------- Save results ----------
    save_in_txt("results_Q1.txt", Xh, Q1, n1);
    save_in_txt("results_Q2.txt", Xh, Q2, n1);
    rename("results_Q1.txt", filename1);
    rename("results_Q2.txt", filename2);

    // ---------- Generate graph ----------
    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    } else {
        printf("⚠️  Warning: 'graph_convergence.py' not found. Skipping graph generation.\n");
    }
}

void convergence_factor_heun_3(int n1, double h1, double X0, double Xf,
                               double Y10, double Y20, double Y30,
                               double (*f1)(double, double, double, double),
                               double (*f2)(double, double, double, double),
                               double (*f3)(double, double, double, double),
                               const char *filename1,
                               const char *filename2,
                               const char *filename3) {
    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    double Y1h[MAX_SIZE + 1], Y2h[MAX_SIZE + 1], Y3h[MAX_SIZE + 1];
    double Y1h2[MAX_SIZE * 2 + 1], Y2h2[MAX_SIZE * 2 + 1], Y3h2[MAX_SIZE * 2 + 1];
    double Y1h4[MAX_SIZE * 4 + 1], Y2h4[MAX_SIZE * 4 + 1], Y3h4[MAX_SIZE * 4 + 1];
    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE * 2 + 1], Xh4[MAX_SIZE * 4 + 1];
    double Q1[MAX_SIZE + 1], Q2[MAX_SIZE + 1], Q3[MAX_SIZE + 1];

    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Y1h[0] = Y1h2[0] = Y1h4[0] = Y10;
    Y2h[0] = Y2h2[0] = Y2h4[0] = Y20;
    Y3h[0] = Y3h2[0] = Y3h4[0] = Y30;

    // --- Heun with step h ---
    for (int i = 0; i < n1; i++) {
        double k11 = f1(Xh[i], Y1h[i], Y2h[i], Y3h[i]);
        double k12 = f2(Xh[i], Y1h[i], Y2h[i], Y3h[i]);
        double k13 = f3(Xh[i], Y1h[i], Y2h[i], Y3h[i]);

        double y1_pred = Y1h[i] + h1 * k11;
        double y2_pred = Y2h[i] + h1 * k12;
        double y3_pred = Y3h[i] + h1 * k13;

        double k21 = f1(Xh[i] + h1, y1_pred, y2_pred, y3_pred);
        double k22 = f2(Xh[i] + h1, y1_pred, y2_pred, y3_pred);
        double k23 = f3(Xh[i] + h1, y1_pred, y2_pred, y3_pred);

        Y1h[i+1] = Y1h[i] + (h1 / 2.0) * (k11 + k21);
        Y2h[i+1] = Y2h[i] + (h1 / 2.0) * (k12 + k22);
        Y3h[i+1] = Y3h[i] + (h1 / 2.0) * (k13 + k23);
        Xh[i+1] = Xh[i] + h1;
    }

    // --- Heun with step h/2 ---
    for (int i = 0; i < 2 * n1; i++) {
        double k11 = f1(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i]);
        double k12 = f2(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i]);
        double k13 = f3(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i]);

        double y1_pred = Y1h2[i] + h2 * k11;
        double y2_pred = Y2h2[i] + h2 * k12;
        double y3_pred = Y3h2[i] + h2 * k13;

        double k21 = f1(Xh2[i] + h2, y1_pred, y2_pred, y3_pred);
        double k22 = f2(Xh2[i] + h2, y1_pred, y2_pred, y3_pred);
        double k23 = f3(Xh2[i] + h2, y1_pred, y2_pred, y3_pred);

        Y1h2[i+1] = Y1h2[i] + (h2 / 2.0) * (k11 + k21);
        Y2h2[i+1] = Y2h2[i] + (h2 / 2.0) * (k12 + k22);
        Y3h2[i+1] = Y3h2[i] + (h2 / 2.0) * (k13 + k23);
        Xh2[i+1] = Xh2[i] + h2;
    }

    // --- Heun with step h/4 ---
    for (int i = 0; i < 4 * n1; i++) {
        double k11 = f1(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i]);
        double k12 = f2(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i]);
        double k13 = f3(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i]);

        double y1_pred = Y1h4[i] + h3 * k11;
        double y2_pred = Y2h4[i] + h3 * k12;
        double y3_pred = Y3h4[i] + h3 * k13;

        double k21 = f1(Xh4[i] + h3, y1_pred, y2_pred, y3_pred);
        double k22 = f2(Xh4[i] + h3, y1_pred, y2_pred, y3_pred);
        double k23 = f3(Xh4[i] + h3, y1_pred, y2_pred, y3_pred);

        Y1h4[i+1] = Y1h4[i] + (h3 / 2.0) * (k11 + k21);
        Y2h4[i+1] = Y2h4[i] + (h3 / 2.0) * (k12 + k22);
        Y3h4[i+1] = Y3h4[i] + (h3 / 2.0) * (k13 + k23);
        Xh4[i+1] = Xh4[i] + h3;
    }

    printf("\n%-10s %-15s %-15s %-15s %-15s\n", "i", "x_i", "Q1_i", "Q2_i", "Q3_i");
    printf("--------------------------------------------------------------------------\n");

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2 * i;
        int idx4 = 4 * i;

        double num1 = fabs(Y1h[i] - Y1h2[idx2]);
        double den1 = fabs(Y1h2[idx2] - Y1h4[idx4]);
        double num2 = fabs(Y2h[i] - Y2h2[idx2]);
        double den2 = fabs(Y2h2[idx2] - Y2h4[idx4]);
        double num3 = fabs(Y3h[i] - Y3h2[idx2]);
        double den3 = fabs(Y3h2[idx2] - Y3h4[idx4]);

        Q1[i] = (den1 > 1e-12) ? log(num1 / den1) / log(2.0) : 0.0;
        Q2[i] = (den2 > 1e-12) ? log(num2 / den2) / log(2.0) : 0.0;
        Q3[i] = (den3 > 1e-12) ? log(num3 / den3) / log(2.0) : 0.0;

        printf("%-10d %-15lf %-15lf %-15lf %-15lf\n", i, Xh[i], Q1[i], Q2[i], Q3[i]);
    }

    save_in_txt("results_Q1.txt", Xh, Q1, n1);
    save_in_txt("results_Q2.txt", Xh, Q2, n1);
    save_in_txt("results_Q3.txt", Xh, Q3, n1);
    rename("results_Q1.txt", filename1);
    rename("results_Q2.txt", filename2);
    rename("results_Q3.txt", filename3);

    if (system("test -f graph_convergence.py") == 0)
        system("python3 graph_convergence.py");
    else
        printf("⚠️  Warning: 'graph_convergence.py' not found. Skipping graph generation.\n");
}


void convergence_factor_midpoint_2(int n1, double h1, double X0, double Xf,
                                   double Y10, double Y20,
                                   double (*f1)(double, double, double),
                                   double (*f2)(double, double, double),
                                   const char *fileQ1,
                                   const char *fileQ2) {

    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE*2 + 1], Xh4[MAX_SIZE*4 + 1];
    double Y1h[MAX_SIZE + 1], Y1h2[MAX_SIZE*2 + 1], Y1h4[MAX_SIZE*4 + 1];
    double Y2h[MAX_SIZE + 1], Y2h2[MAX_SIZE*2 + 1], Y2h4[MAX_SIZE*4 + 1];
    double Q1[MAX_SIZE + 1], Q2[MAX_SIZE + 1];

    // Inicialización
    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Y1h[0] = Y1h2[0] = Y1h4[0] = Y10;
    Y2h[0] = Y2h2[0] = Y2h4[0] = Y20;

    // ---------- Midpoint paso h ----------
    for (int i = 0; i < n1; i++) {
        double k11 = f1(Xh[i], Y1h[i], Y2h[i]);
        double k12 = f2(Xh[i], Y1h[i], Y2h[i]);
        double k21 = f1(Xh[i] + h1/2.0, Y1h[i] + (h1/2.0)*k11, Y2h[i] + (h1/2.0)*k12);
        double k22 = f2(Xh[i] + h1/2.0, Y1h[i] + (h1/2.0)*k11, Y2h[i] + (h1/2.0)*k12);
        Y1h[i+1] = Y1h[i] + h1 * k21;
        Y2h[i+1] = Y2h[i] + h1 * k22;
        Xh[i+1] = Xh[i] + h1;
    }

    // ---------- Midpoint paso h/2 ----------
    for (int i = 0; i < 2*n1; i++) {
        double k11 = f1(Xh2[i], Y1h2[i], Y2h2[i]);
        double k12 = f2(Xh2[i], Y1h2[i], Y2h2[i]);
        double k21 = f1(Xh2[i] + h2/2.0, Y1h2[i] + (h2/2.0)*k11, Y2h2[i] + (h2/2.0)*k12);
        double k22 = f2(Xh2[i] + h2/2.0, Y1h2[i] + (h2/2.0)*k11, Y2h2[i] + (h2/2.0)*k12);
        Y1h2[i+1] = Y1h2[i] + h2 * k21;
        Y2h2[i+1] = Y2h2[i] + h2 * k22;
        Xh2[i+1] = Xh2[i] + h2;
    }

    // ---------- Midpoint paso h/4 ----------
    for (int i = 0; i < 4*n1; i++) {
        double k11 = f1(Xh4[i], Y1h4[i], Y2h4[i]);
        double k12 = f2(Xh4[i], Y1h4[i], Y2h4[i]);
        double k21 = f1(Xh4[i] + h3/2.0, Y1h4[i] + (h3/2.0)*k11, Y2h4[i] + (h3/2.0)*k12);
        double k22 = f2(Xh4[i] + h3/2.0, Y1h4[i] + (h3/2.0)*k11, Y2h4[i] + (h3/2.0)*k12);
        Y1h4[i+1] = Y1h4[i] + h3 * k21;
        Y2h4[i+1] = Y2h4[i] + h3 * k22;
        Xh4[i+1] = Xh4[i] + h3;
    }

    // ---------- Calcular Q1 y Q2 ----------
    printf("\n%-10s %-15s %-15s %-15s\n", "i", "x_i", "Q1_i", "Q2_i");
    printf("------------------------------------------------------------\n");

    Q1[0] = Q2[0] = 0.0;

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;
        int idx4 = 4*i;

        double num1 = fabs(Y1h[i] - Y1h2[idx2]);
        double den1 = fabs(Y1h2[idx2] - Y1h4[idx4]);
        double num2 = fabs(Y2h[i] - Y2h2[idx2]);
        double den2 = fabs(Y2h2[idx2] - Y2h4[idx4]);

        Q1[i] = (den1 > 1e-12) ? log(num1 / den1) / log(2.0) : 0.0;
        Q2[i] = (den2 > 1e-12) ? log(num2 / den2) / log(2.0) : 0.0;

        printf("%-10d %-15lf %-15lf %-15lf\n", i, Xh[i], Q1[i], Q2[i]);
    }

    save_in_txt("results_Q1.txt", Xh, Q1, n1);
    save_in_txt("results_Q2.txt", Xh, Q2, n1);
    rename("results_Q1.txt", fileQ1);
    rename("results_Q2.txt", fileQ2);

    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    } else {
        printf("⚠️  Warning: 'graph_convergence.py' not found. Skipping graph generation.\n");
    }
}

void convergence_factor_midpoint_3(int n1, double h1, double X0, double Xf,
                                   double Y10, double Y20, double Y30,
                                   double (*f1)(double, double, double, double),
                                   double (*f2)(double, double, double, double),
                                   double (*f3)(double, double, double, double),
                                   const char *fileQ1,
                                   const char *fileQ2,
                                   const char *fileQ3) {

    double h2 = h1 / 2.0;
    double h3 = h1 / 4.0;

    double Xh[MAX_SIZE + 1], Xh2[MAX_SIZE*2 + 1], Xh4[MAX_SIZE*4 + 1];
    double Y1h[MAX_SIZE + 1], Y1h2[MAX_SIZE*2 + 1], Y1h4[MAX_SIZE*4 + 1];
    double Y2h[MAX_SIZE + 1], Y2h2[MAX_SIZE*2 + 1], Y2h4[MAX_SIZE*4 + 1];
    double Y3h[MAX_SIZE + 1], Y3h2[MAX_SIZE*2 + 1], Y3h4[MAX_SIZE*4 + 1];
    double Q1[MAX_SIZE + 1], Q2[MAX_SIZE + 1], Q3[MAX_SIZE + 1];

    // Inicialización
    Xh[0] = Xh2[0] = Xh4[0] = X0;
    Y1h[0] = Y1h2[0] = Y1h4[0] = Y10;
    Y2h[0] = Y2h2[0] = Y2h4[0] = Y20;
    Y3h[0] = Y3h2[0] = Y3h4[0] = Y30;

    // ---------- Midpoint paso h ----------
    for (int i = 0; i < n1; i++) {
        double k11 = f1(Xh[i], Y1h[i], Y2h[i], Y3h[i]);
        double k12 = f2(Xh[i], Y1h[i], Y2h[i], Y3h[i]);
        double k13 = f3(Xh[i], Y1h[i], Y2h[i], Y3h[i]);

        double k21 = f1(Xh[i] + h1/2.0, Y1h[i] + (h1/2.0)*k11, Y2h[i] + (h1/2.0)*k12, Y3h[i] + (h1/2.0)*k13);
        double k22 = f2(Xh[i] + h1/2.0, Y1h[i] + (h1/2.0)*k11, Y2h[i] + (h1/2.0)*k12, Y3h[i] + (h1/2.0)*k13);
        double k23 = f3(Xh[i] + h1/2.0, Y1h[i] + (h1/2.0)*k11, Y2h[i] + (h1/2.0)*k12, Y3h[i] + (h1/2.0)*k13);

        Y1h[i+1] = Y1h[i] + h1 * k21;
        Y2h[i+1] = Y2h[i] + h1 * k22;
        Y3h[i+1] = Y3h[i] + h1 * k23;
        Xh[i+1] = Xh[i] + h1;
    }

    // ---------- Midpoint paso h/2 ----------
    for (int i = 0; i < 2*n1; i++) {
        double k11 = f1(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i]);
        double k12 = f2(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i]);
        double k13 = f3(Xh2[i], Y1h2[i], Y2h2[i], Y3h2[i]);

        double k21 = f1(Xh2[i] + h2/2.0, Y1h2[i] + (h2/2.0)*k11, Y2h2[i] + (h2/2.0)*k12, Y3h2[i] + (h2/2.0)*k13);
        double k22 = f2(Xh2[i] + h2/2.0, Y1h2[i] + (h2/2.0)*k11, Y2h2[i] + (h2/2.0)*k12, Y3h2[i] + (h2/2.0)*k13);
        double k23 = f3(Xh2[i] + h2/2.0, Y1h2[i] + (h2/2.0)*k11, Y2h2[i] + (h2/2.0)*k12, Y3h2[i] + (h2/2.0)*k13);

        Y1h2[i+1] = Y1h2[i] + h2 * k21;
        Y2h2[i+1] = Y2h2[i] + h2 * k22;
        Y3h2[i+1] = Y3h2[i] + h2 * k23;
        Xh2[i+1] = Xh2[i] + h2;
    }

    // ---------- Midpoint paso h/4 ----------
    for (int i = 0; i < 4*n1; i++) {
        double k11 = f1(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i]);
        double k12 = f2(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i]);
        double k13 = f3(Xh4[i], Y1h4[i], Y2h4[i], Y3h4[i]);

        double k21 = f1(Xh4[i] + h3/2.0, Y1h4[i] + (h3/2.0)*k11, Y2h4[i] + (h3/2.0)*k12, Y3h4[i] + (h3/2.0)*k13);
        double k22 = f2(Xh4[i] + h3/2.0, Y1h4[i] + (h3/2.0)*k11, Y2h4[i] + (h3/2.0)*k12, Y3h4[i] + (h3/2.0)*k13);
        double k23 = f3(Xh4[i] + h3/2.0, Y1h4[i] + (h3/2.0)*k11, Y2h4[i] + (h3/2.0)*k12, Y3h4[i] + (h3/2.0)*k13);

        Y1h4[i+1] = Y1h4[i] + h3 * k21;
        Y2h4[i+1] = Y2h4[i] + h3 * k22;
        Y3h4[i+1] = Y3h4[i] + h3 * k23;
        Xh4[i+1] = Xh4[i] + h3;
    }

    // ---------- Calcular Q1, Q2, Q3 ----------
    printf("\n%-10s %-15s %-15s %-15s %-15s\n", "i", "x_i", "Q1_i", "Q2_i", "Q3_i");
    printf("--------------------------------------------------------------------------\n");

    Q1[0] = Q2[0] = Q3[0] = 0.0;

    for (int i = 1; i <= n1; i++) {
        int idx2 = 2*i;
        int idx4 = 4*i;

        double num1 = fabs(Y1h[i] - Y1h2[idx2]);
        double den1 = fabs(Y1h2[idx2] - Y1h4[idx4]);
        double num2 = fabs(Y2h[i] - Y2h2[idx2]);
        double den2 = fabs(Y2h2[idx2] - Y2h4[idx4]);
        double num3 = fabs(Y3h[i] - Y3h2[idx2]);
        double den3 = fabs(Y3h2[idx2] - Y3h4[idx4]);

        Q1[i] = (den1 > 1e-12) ? log(num1 / den1) / log(2.0) : 0.0;
        Q2[i] = (den2 > 1e-12) ? log(num2 / den2) / log(2.0) : 0.0;
        Q3[i] = (den3 > 1e-12) ? log(num3 / den3) / log(2.0) : 0.0;

        printf("%-10d %-15lf %-15lf %-15lf %-15lf\n", i, Xh[i], Q1[i], Q2[i], Q3[i]);
    }

    save_in_txt("results_Q1.txt", Xh, Q1, n1);
    save_in_txt("results_Q2.txt", Xh, Q2, n1);
    save_in_txt("results_Q3.txt", Xh, Q3, n1);
    rename("results_Q1.txt", fileQ1);
    rename("results_Q2.txt", fileQ2);
    rename("results_Q3.txt", fileQ3);

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



void convergence_factor_rk4_3(int n1, double h1, double X0,
                              double Y10, double Y20, double Y30,
                              double (*f1)(double,double,double,double),
                              double (*f2)(double,double,double,double),
                              double (*f3)(double,double,double,double),
                              const char *filename1,
                              const char *filename2,
                              const char *filename3) {
    double h2 = h1/2.0;
    double h3 = h1/4.0;

    double Y1h[MAX_SIZE+1], Y2h[MAX_SIZE+1], Y3h[MAX_SIZE+1];
    double Y1h2[MAX_SIZE*2+1], Y2h2[MAX_SIZE*2+1], Y3h2[MAX_SIZE*2+1];
    double Y1h4[MAX_SIZE*4+1], Y2h4[MAX_SIZE*4+1], Y3h4[MAX_SIZE*4+1];
    double Xh[MAX_SIZE+1], Xh2[MAX_SIZE*2+1], Xh4[MAX_SIZE*4+1];
    double Q1[MAX_SIZE+1], Q2[MAX_SIZE+1], Q3[MAX_SIZE+1];

    // Inicialización
    Xh[0]=Xh2[0]=Xh4[0]=X0;
    Y1h[0]=Y1h2[0]=Y1h4[0]=Y10;
    Y2h[0]=Y2h2[0]=Y2h4[0]=Y20;
    Y3h[0]=Y3h2[0]=Y3h4[0]=Y30;

    // RK4 paso h
    for(int i=0;i<n1;i++){
        rk4_step3(Xh[i], &Y1h[i], &Y2h[i], &Y3h[i], h1, f1,f2,f3);
        Xh[i+1] = Xh[i]+h1;
    }
    // RK4 paso h/2
    for(int i=0;i<2*n1;i++){
        rk4_step3(Xh2[i], &Y1h2[i], &Y2h2[i], &Y3h2[i], h2, f1,f2,f3);
        Xh2[i+1] = Xh2[i]+h2;
    }
    // RK4 paso h/4
    for(int i=0;i<4*n1;i++){
        rk4_step3(Xh4[i], &Y1h4[i], &Y2h4[i], &Y3h4[i], h3, f1,f2,f3);
        Xh4[i+1] = Xh4[i]+h3;
    }

    printf("\n%-10s %-15s %-15s %-15s %-15s\n","i","x_i","Q1_i","Q2_i","Q3_i");
    printf("-------------------------------------------------------------\n");

    for(int i=1;i<=n1;i++){
        double xi = X0 + i*h1;                // x exacto
        int idx2 = (int)((xi - X0)/h2 + 0.5); // índice en h/2
        int idx4 = (int)((xi - X0)/h3 + 0.5); // índice en h/4

        double num1 = fabs(Y1h[i]-Y1h2[idx2]);
        double den1 = fabs(Y1h2[idx2]-Y1h4[idx4]);
        double num2 = fabs(Y2h[i]-Y2h2[idx2]);
        double den2 = fabs(Y2h2[idx2]-Y2h4[idx4]);
        double num3 = fabs(Y3h[i]-Y3h2[idx2]);
        double den3 = fabs(Y3h2[idx2]-Y3h4[idx4]);

        Q1[i] = (den1>1e-12)? log(num1/den1)/log(2.0) : 0.0;
        Q2[i] = (den2>1e-12)? log(num2/den2)/log(2.0) : 0.0;
        Q3[i] = (den3>1e-12)? log(num3/den3)/log(2.0) : 0.0;

        char s1[64], s2[64], s3[64];
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
        if (Q3[i] != 0.0) {
            snprintf(s3, sizeof(s3), "%15.6f", Q3[i]);
        } else {
            snprintf(s3, sizeof(s3), "N/A");
        }

        printf("%-10d %-15lf %-15s %-15s %-15s\n", i, xi, s1, s2, s3);
    }

    // Guardar resultados en TXT
    save_in_txt(filename1, Xh, Q1, n1);
    save_in_txt(filename2, Xh, Q2, n1);
    save_in_txt(filename3, Xh, Q3, n1);

    // Generar gráfico si existe script Python
    if (system("test -f graph_convergence.py") == 0) {
        system("python3 graph_convergence.py");
    } else {
        printf("⚠️  Warning: 'graph_convergence.py' not found. Skipping graph generation.\n");
    }
}




