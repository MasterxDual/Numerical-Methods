#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_SIZE 500

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

/* 

Non-trivial solution:
y(t) = pow((2.0/3.0) * t, 1.5); with y(0) = 1e-16

h = 0.1, X0 = 0, Xf = 1, Y0 = 1e-16

i          X[i]            Exact Y         Euler Y         Exact Error     Local Trunc. Err
-------------------------------------------------------------------------------------------
0          0.000000        0.000000        0.000000        0.000000        0.000000       
1          0.100000        0.017213        0.000000        0.017213        0.006455       
2          0.200000        0.048686        0.000775        0.047912        0.004564       
3          0.300000        0.089443        0.009959        0.079484        0.003727       
4          0.400000        0.137706        0.031474        0.106232        0.003227       
5          0.500000        0.192450        0.063047        0.129403        0.002887       
6          0.600000        0.252982        0.102848        0.150135        0.002635       
7          0.700000        0.318794        0.149700        0.169094        0.002440       
8          0.800000        0.389492        0.202798        0.186694        0.002282       
9          0.900000        0.464758        0.261549        0.203209        0.002152       
10         1.000000        0.544331        0.325501        0.218830        0.002041       
X[i]            Y[i]
0.000000        0.000000
0.100000        0.000000
0.200000        0.000775
0.300000        0.009959
0.400000        0.031474
0.500000        0.063047
0.600000        0.102848
0.700000        0.149700
0.800000        0.202798
0.900000        0.261549
1.000000        0.325501



h = 0.01, X0 = 0, Xf = 0.5, Y0 = 1e-16
i          X[i]            Exact Y         Euler Y         Exact Error     Local Trunc. Err
-------------------------------------------------------------------------------------------
0          0.000000        0.000000        0.000000        0.000000        0.000000       
1          0.010000        0.000544        0.000000        0.000544        0.000204       
2          0.020000        0.001540        0.000036        0.001504        0.000144       
3          0.030000        0.002828        0.000366        0.002462        0.000118       
4          0.040000        0.004355        0.001082        0.003273        0.000102       
5          0.050000        0.006086        0.002108        0.003978        0.000091       
6          0.060000        0.008000        0.003390        0.004610        0.000083       
7          0.070000        0.010081        0.004892        0.005189        0.000077       
8          0.080000        0.012317        0.006590        0.005727        0.000072       
9          0.090000        0.014697        0.008465        0.006232        0.000068       
10         0.100000        0.017213        0.010503        0.006710        0.000065       
11         0.110000        0.019859        0.012693        0.007166        0.000062       
12         0.120000        0.022627        0.015026        0.007602        0.000059       
13         0.130000        0.025514        0.017493        0.008021        0.000057       
14         0.140000        0.028514        0.020089        0.008425        0.000055       
15         0.150000        0.031623        0.022808        0.008815        0.000053       
16         0.160000        0.034837        0.025643        0.009194        0.000051       
17         0.170000        0.038154        0.028592        0.009561        0.000050       
18         0.180000        0.041569        0.031650        0.009919        0.000048       
19         0.190000        0.045081        0.034813        0.010268        0.000047       
20         0.200000        0.048686        0.038079        0.010608        0.000046       
21         0.210000        0.052383        0.041443        0.010940        0.000045       
22         0.220000        0.056169        0.044903        0.011266        0.000044       
23         0.230000        0.060042        0.048458        0.011584        0.000043       
24         0.240000        0.064000        0.052104        0.011896        0.000042       
25         0.250000        0.068041        0.055839        0.012203        0.000041       
26         0.260000        0.072164        0.059661        0.012504        0.000040       
27         0.270000        0.076368        0.063568        0.012799        0.000039       
28         0.280000        0.080649        0.067559        0.013090        0.000039       
29         0.290000        0.085008        0.071632        0.013376        0.000038       
30         0.300000        0.089443        0.075785        0.013658        0.000037       
31         0.310000        0.093952        0.080017        0.013935        0.000037       
32         0.320000        0.098534        0.084326        0.014208        0.000036       
33         0.330000        0.103189        0.088711        0.014478        0.000036       
34         0.340000        0.107915        0.093171        0.014744        0.000035       
35         0.350000        0.112711        0.097705        0.015006        0.000035       
36         0.360000        0.117576        0.102310        0.015265        0.000034       
37         0.370000        0.122508        0.106987        0.015521        0.000034       
38         0.380000        0.127508        0.111735        0.015774        0.000033       
39         0.390000        0.132575        0.116551        0.016023        0.000033       
40         0.400000        0.137706        0.121436        0.016270        0.000032       
41         0.410000        0.142902        0.126388        0.016514        0.000032       
42         0.420000        0.148162        0.131406        0.016756        0.000031       
43         0.430000        0.153485        0.136490        0.016995        0.000031       
44         0.440000        0.158870        0.141639        0.017231        0.000031       
45         0.450000        0.164317        0.146852        0.017465        0.000030       
46         0.460000        0.169824        0.152128        0.017697        0.000030       
47         0.470000        0.175392        0.157466        0.017926        0.000030       
48         0.480000        0.181019        0.162866        0.018153        0.000029       
49         0.490000        0.186706        0.168327        0.018378        0.000029       
50         0.500000        0.192450        0.173848        0.018602        0.000029       
X[i]            Y[i]
0.000000        0.000000
0.010000        0.000000
0.020000        0.000036
0.030000        0.000366
0.040000        0.001082
0.050000        0.002108
0.060000        0.003390
0.070000        0.004892
0.080000        0.006590
0.090000        0.008465
0.100000        0.010503
0.110000        0.012693
0.120000        0.015026
0.130000        0.017493
0.140000        0.020089
0.150000        0.022808
0.160000        0.025643
0.170000        0.028592
0.180000        0.031650
0.190000        0.034813
0.200000        0.038079
0.210000        0.041443
0.220000        0.044903
0.230000        0.048458
0.240000        0.052104
0.250000        0.055839
0.260000        0.059661
0.270000        0.063568
0.280000        0.067559
0.290000        0.071632
0.300000        0.075785
0.310000        0.080017
0.320000        0.084326
0.330000        0.088711
0.340000        0.093171
0.350000        0.097705
0.360000        0.102310
0.370000        0.106987
0.380000        0.111735
0.390000        0.116551
0.400000        0.121436
0.410000        0.126388
0.420000        0.131406
0.430000        0.136490
0.440000        0.141639
0.450000        0.146852
0.460000        0.152128
0.470000        0.157466
0.480000        0.162866
0.490000        0.168327
0.500000        0.173848

As we can see, if we reduce h, the error keeps growing, so this may indicate that the not trivial solution is not reproduced 
by Euler's method, at least for this range of h.

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

    printf("Insert the method to use: 1. Euler's method\n");
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
    return cbrt(y);
}

double y(double x) {
    return pow((2.0/3.0) * x, 1.5);
}


double fprima(double x, double y) {
    if (y <= 0.0) return 0.0; 
    return (1.0 / 3.0) * pow(y, -1.0 / 3.0);
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

