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
Reading the results of Composition Trapeze and Simpson from exercise 3 guide 7 we have this results:
--> f(x) = sin(2*x) * exp(-x)
    a = 0, b = pi = 3.14159
        Composition Trapeze:
            subintervals = 10
            The integral is: 0.366951
        Composition Simpson:
            subintervals = 5 panels * 2 = 10
            The integral is: 0.382793
Using Gauss quadrature method with two and three points we have:
    2 points:
        The integral is: 0.703456
    3 points:
       The integral is: 0.376098 
Conclusion:
--> Undersampling: 2 points cannot capture a function that oscillates twice as much.
--> Systematic overestimation: Gaussian points fall into regions where the function is relatively high.
--> Theretical limitation: 2-point Gaussian is accurate only for polynomials of degree ≤ 3, but the function of our problem
is transcendental and oscillatory.

Summary:
    2 points:
        Exact for: Polynomials of degree ≤ 3
        Ideal functions:
        f(x) = ax³ + bx² + cx + d
        Very smooth and monotonic functions
        Simple exponentials in small intervals
        When it fails: Oscillating functions, abrupt changes
        Example: f(x) = x³ + 2x² → Exact result
    3 points:
        Exact for: Polynomials of degree ≤ 5
        Ideal functions:
        Moderate oscillations (like your sin(2x)*e^(-x))
        Functions with 1-2 local extrema
        Exponentials with complex behavior
        Optimal trade-off: Accuracy vs. efficiency
        Example: Your function → Very good result
    4 points:
        Exact for: Polynomials of degree ≤ 7
        Ideal functions:
        Greater curvature and complexity
        Trigonometric functions with medium frequency
        Combinations of exponential and trigonometric functions
        When to use: You need high precision without being excessive
    5 points:
        Exact for: Polynomials of degree ≤ 9
        Ideal functions:
        Highly oscillatory (e.g., sin(10x))
        Multiple local extrema
        Functions with smooth singularities
        Cost-benefit: Begins to be expensive
    6 points:
        Exact for: Polynomials of degree ≤ 11
        Ideal functions:
        Extremely complex
        Many oscillations in the interval
        When 5 points do not converge
        Warning: Diminishing returns
        
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
    return sin(2*x) * exp(-x);
}
