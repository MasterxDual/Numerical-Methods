#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 1000

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
    f(t, v) = 10 - (0.01 * pow(v, 1.5))
    X0 = 0, Xf= 6
    v(0) = 0
    h = 0.05
    We cant calculate v(t) because it's a difficult ODE to solve analytically. Anyway, this exercise doesn't ask for it.
    So we have the following results:

    t[i]            v[i]
    0.000000        0.000000
    0.050000        0.500000
    0.100000        0.999823
    0.150000        1.499323
    0.200000        1.998405
    0.250000        2.496993
    0.300000        2.995020
    0.350000        3.492428
    0.400000        3.989165
    0.450000        4.485181
    0.500000        4.980432
    0.550000        5.474875
    0.600000        5.968469
    0.650000        6.461179
    0.700000        6.952967
    0.750000        7.443800
    0.800000        7.933645
    0.850000        8.422472
    0.900000        8.910251
    0.950000        9.396952
    1.000000        9.882549
    1.050000        10.367015
    1.100000        10.850326
    1.150000        11.332455
    1.200000        11.813381
    1.250000        12.293079
    1.300000        12.771528
    1.350000        13.248707
    1.400000        13.724595
    1.450000        14.199173
    1.500000        14.672420
    1.550000        15.144319
    1.600000        15.614852
    1.650000        16.084000
    1.700000        16.551748
    1.750000        17.018078
    1.800000        17.482976
    1.850000        17.946426
    1.900000        18.408412
    1.950000        18.868922
    2.000000        19.327940
    2.050000        19.785454
    2.100000        20.241450
    2.150000        20.695916
    2.200000        21.148840
    2.250000        21.600211
    2.300000        22.050016
    2.350000        22.498246
    2.400000        22.944889
    2.450000        23.389935
    2.500000        23.833374
    2.550000        24.275197
    2.600000        24.715396
    2.650000        25.153960
    2.700000        25.590882
    2.750000        26.026153
    2.800000        26.459766
    2.850000        26.891712
    2.900000        27.321986
    2.950000        27.750579
    3.000000        28.177486
    3.050000        28.602699
    3.100000        29.026214
    3.150000        29.448023
    3.200000        29.868121
    3.250000        30.286504
    3.300000        30.703166
    3.350000        31.118102
    3.400000        31.531308
    3.450000        31.942780
    3.500000        32.352513
    3.550000        32.760503
    3.600000        33.166748
    3.650000        33.571243
    3.700000        33.973986
    3.750000        34.374974
    3.800000        34.774203
    3.850000        35.171672
    3.900000        35.567378
    3.950000        35.961319
    4.000000        36.353493
    4.050000        36.743899
    4.100000        37.132534
    4.150000        37.519398
    4.200000        37.904489
    4.250000        38.287806
    4.300000        38.669349
    4.350000        39.049117
    4.400000        39.427109
    4.450000        39.803326
    4.500000        40.177767
    4.550000        40.550431
    4.600000        40.921320
    4.650000        41.290434
    4.700000        41.657773
    4.750000        42.023337
    4.800000        42.387128
    4.850000        42.749147
    4.900000        43.109394
    4.950000        43.467870
    5.000000        43.824578
    5.050000        44.179518
    5.100000        44.532693
    5.150000        44.884103
    5.200000        45.233752
    5.250000        45.581639
    5.300000        45.927769
    5.350000        46.272143
    5.400000        46.614763
    5.450000        46.955632
    5.500000        47.294752
    5.550000        47.632126
    5.600000        47.967757
    5.650000        48.301648
    5.700000        48.633801
    5.750000        48.964220
    5.800000        49.292908
    5.850000        49.619868
    5.900000        49.945103
    5.950000        50.268617
    6.000000        50.590414

    The final result v(6) = 50.590414

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
        case 1: {
            for(int i = 1; i <= n; i++) {
                // Compute next X incrementally to avoid rounding surprises
                X[i] = X[i-1] + h; // X[i] = X0 + i*h; equivalently
                // Euler method
                Y[i] = Y[i-1] + h * f(X[i-1], Y[i-1]); 
            }
            break;
        }
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

double f(double t, double v) {
    return (10 - (0.01 * pow(v, 1.5)));
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

