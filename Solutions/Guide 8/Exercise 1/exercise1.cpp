#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 20
#define MAX_SIZE 100

/**
 * Function to define function f(x)
 * @param x The point at which to evaluate the function (in our case, X̂)
 * @return The value of the function at x
 */
double f(double x);

/* 
Using sustitution in the integral we reached the following integral:
    ∫(sin(x+1)/x+1)dx from -1 to 1. And using aproximation of two points with Gauss we finally have the result:
    The integral is: 1.604454
Notice that we can use the following integral too:
    ∫((sin(x)/x))dx from 0 to 2 and we have the result:
    The integral is: 1.604454
Also, we did it writing in paper the sustitution and the aproximation and we reached the same result.
*/

int main(int argc, char const *argv[]) {
    // Weighting factor
    double c0, c1, c2, c3, c4, c5;
    // Function arguments
    double x0, x1, x2, x3, x4, x5;
    // Limits of integration
    double a, b;
    // Number of points
    int number_of_points;
    // Result of the integral
    double integral;

    printf("Insert the limits of integration:\n");
    scanf("%lf %lf", &a, &b);    

    printf("Insert the number of points (between 2 and 6)\n");
    scanf("%d", &number_of_points);

    switch(number_of_points) {
        case 2: 
            // Implement 2-point Gauss-Legendre quadrature
            c0 = 1.0;
            c1 = 1.0;
            x0 = -0.577350269;
            x1 = 0.577350269;
            integral = (((b-a)/2) * (c0 * f((((b-a)*x0) + (b+a))/2)) + (c1 * f((((b-a)*x1) + (b+a))/2)));
            break;
        case 3:
            // Implement 3-point Gauss-Legendre quadrature
            c0 = 0.5555556;
            c1 = 0.8888889;
            c2 = 0.5555556;
            x0 = -0.774596669;
            x1 = 0.0;
            x2 = 0.774596669;
            integral = (((b-a)/2) * (c0 * f((((b-a)*x0) + (b+a))/2)) + (c1 * f((((b-a)*x1) + (b+a))/2)) + (c2 * f((((b-a)*x2) + (b+a))/2)));
            break;
        case 4:
            // Implement 4-point Gauss-Legendre quadrature
            c0 = 0.3478548;
            c1 = 0.6521452;
            c2 = 0.6521452;
            c3 = 0.3478548;
            x0 = -0.861136312;
            x1 = -0.339981044;
            x2 = 0.339981044;
            x3 = 0.861136312;
            integral = (((b-a)/2) * (c0 * f((((b-a)*x0) + (b+a))/2)) + (c1 * f((((b-a)*x1) + (b+a))/2)) + (c2 * f((((b-a)*x2) + (b+a))/2)) + (c3 * f((((b-a)*x3) + (b+a))/2)));
            break;
        case 5:
            // Implement 5-point Gauss-Legendre quadrature
            c0 = 0.2369269;
            c1 = 0.4786287;
            c2 = 0.5688889;
            c3 = 0.4786287;
            c4 = 0.2369269;
            x0 = -0.906179846;
            x1 = -0.538469310;
            x2 = 0.0;
            x3 = 0.538469310;
            x4 = 0.906179846;
            integral = (((b-a)/2) * (c0 * f((((b-a)*x0) + (b+a))/2)) + (c1 * f((((b-a)*x1) + (b+a))/2)) + (c2 * f((((b-a)*x2) + (b+a))/2)) + (c3 * f((((b-a)*x3) + (b+a))/2)) + (c4 * f((((b-a)*x4) + (b+a))/2)));
            break;
        case 6:
            // Implement 6-point Gauss-Legendre quadrature
            c0 = 0.1713245;
            c1 = 0.3607616;
            c2 = 0.4679139;
            c3 = 0.4679139;
            c4 = 0.3607616;
            c5 = 0.1713245;
            x0 = -0.932469514;
            x1 = -0.661209386;
            x2 = -0.238619186;
            x3 = 0.238619186;
            x4 = 0.661209386;
            x5 = 0.932469514;
            integral = (((b-a)/2) * (c0 * f((((b-a)*x0) + (b+a))/2)) + (c1 * f((((b-a)*x1) + (b+a))/2)) + (c2 * f((((b-a)*x2) + (b+a))/2)) + (c3 * f((((b-a)*x3) + (b+a))/2)) + (c4 * f((((b-a)*x4) + (b+a))/2)) + (c5 * f((((b-a)*x5) + (b+a))/2)));
            break;
        default:
            printf("Error: Number of points must be between 2 and 6\n");
            exit(0);
    }

    printf("The integral is: %lf\n", integral);
    return 0;
}


double f(double x) {
    return (sin(x)/(x));
}
