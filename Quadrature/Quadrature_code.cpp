// ---------------------------------------------- Code for Numerical Integration using Quadratures --------------------------------------------

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

const double PI = 2*acos(0);
const int DP = 6;                   // No of decimal places round off

// Function to return value for analytical function 
double analytical_interpolation_function(double argument)
{
    return sin(PI*argument/2);
}

// Function to round off values to desired number of decimal places
double round_off(double value, int decimal_places)
{
    double factor = pow(10,decimal_places);
    return round(value*factor)/factor;
}

// Function to bubble sort arrays
void bubble_sort(int N, double points[], double weights[])
{
    for(int i=0;i<N+1;i++)
    {
        for(int j=0;j<N-i;j++)
        {
            if(points[j]>points[j+1])
            {
                double temp = points[j+1];
                points[j+1] = points[j];
                points[j] = temp;

                // The following segment swaps the weights also 
                double weight_temp = weights[j+1];
                weights[j+1] = weights[j];
                weights[j] = weight_temp;
            }
        }
    }
}

// Function to generate legendre points and corresponding weights
void generate_legendre_points(int N, double points[], double weights[])
{
    double tolerance = pow(10,-6);                                  // Convergence tolerance
    for(int i=0;i<N+1;i++)
    {
        double initial_guess = cos((2.0*i+1.0)*PI/(2.0*N+2.0));     // Chebyshev point as initial guess 
        double residue = 1;
        double derivative_phi_n;
        
        while(residue>tolerance)
        {
            double phi_0 = 1;
            double phi_1 = initial_guess; 
            double phi_n;
            double old_phi_n_minus_1;
            
            for(int j=1;j<N+1;j++)
            {
                phi_n = ((2.0*(j+1)-1)/(j+1))*initial_guess*phi_1 - (((j+1)-1.0)/(j+1))*phi_0;
                old_phi_n_minus_1 = phi_1;
                phi_0 = phi_1;
                phi_1 = phi_n;
            }
            derivative_phi_n = (N+1)*(old_phi_n_minus_1 - initial_guess*phi_n)/(1.0 - pow(initial_guess,2));
            // I have N+1 th order legendre polynomial, whose roots are x_i

            double x_update = initial_guess - phi_n/derivative_phi_n;
            residue = fabs(x_update - initial_guess);

            initial_guess = x_update;
        }
        points[i] = round_off(initial_guess,DP);
        weights[i] = 2/((1-pow(points[i],2))*pow(derivative_phi_n,2));
    }
    bubble_sort(N,points,weights);
}

// Function to generate lobatto points and corresponding weights
void generate_lobatto_points(int N, double points[], double weights[])
{
    double tolerance = pow(10,-6);
    points[0] = -1;
    points[N] = 1;

    if (N>1)
    {
        for(int i=1;i<N;i++)
        {
            double initial_guess = cos((2.0*i+1.0)*PI/(2.0*N+2.0));
            double residue = 1;
            double store_phi_n;

            while(residue>tolerance)
            {
                double phi_0 = 1;
                double phi_1 = initial_guess;
                double phi_n;
                double old_phi_n_minus_1;

                for(int j=1;j<N;j++)
                {
                    phi_n = ((2.0*(j+1)-1)/(j+1))*initial_guess*phi_1 - (((j+1)-1.0)/(j+1))*phi_0;
                    old_phi_n_minus_1 = phi_1;
                    phi_0 = phi_1;
                    phi_1 = phi_n;
                }

                double derivative_phi_n = N*(old_phi_n_minus_1 - initial_guess*phi_n)/(1.0 - pow(initial_guess,2));
                double second_derivative_phi_n = (2.0*derivative_phi_n - N*(N+1.0)*phi_n)/(1 - pow(initial_guess,2));

                double lobatto_pol_phi_n = (1-pow(initial_guess,2))*derivative_phi_n;
                double lobatto_pol_derivative_phi_n = (1-pow(initial_guess,2))*second_derivative_phi_n - 2*initial_guess*derivative_phi_n;

                double x_update = initial_guess - lobatto_pol_phi_n/lobatto_pol_derivative_phi_n;
                residue = fabs(x_update - initial_guess);

                store_phi_n = phi_n;
                initial_guess = x_update;
            }
            points[i] = round_off(initial_guess,DP);
            weights[i] = 2/(N*(N+1)*pow(store_phi_n,2));
        }
    }
    weights[0] = 2.0/(N*(N+1.0));
    weights[N] = 2.0/(N*(N+1.0));
    bubble_sort(N,points,weights);
}

double interpolation_function_value(int N, double argument, double points[])
{
    double sum = 0;
    for(int i=0;i<N+1;i++)
    {
        double product = 1;
        for(int j=0;j<N+1;j++)
        {
            if(i!=j)
            {
                product = product*(argument - points[j])/(points[i] - points[j]);
            }
        }
        sum = sum + product*analytical_interpolation_function(points[i]);
    }
    return sum;
}

// Function to evaluate numerical integral value using quadrature
double quadrature(int N, double points[], double weights[])
{
    double integral = 0;
    for(int i=0;i<N+1;i++)
    {
        integral = integral + weights[i]*interpolation_function_value(N,points[i],points);
    }
    return integral;
}

// Function to store result in a CSV file
void put_in_csv(const std::string& filename, int count, double result[])
{
    std::ofstream file(filename);
    for (int i=0;i<count;i++)
    {
        file << result[i];
        if (i<count-1)
        {
            file << ",";
        }
    }
    file << std::endl;
    file.close();
}

int main()
{
    
    /*
    Testing:
    // double x_start = -1;
    // double x_end = 1;
    // int N = 4;                      // N = Interpolation Order, (N+1) should be number of points used for interpolation
        for(int i=0;i<N+1;i++)
        {
            cout << legendre_points[i] << " " << legendre_weights[i] << "\n";
        }

        for(int i=0;i<N+1;i++)
        {
            cout << lobatto_points[i] << " " << lobatto_weights[i] << "\n";
        }
    */

    //cout << interpolation_function_value(N,0,lobatto_points);
    //cout << quadrature(N,lobatto_points,lobatto_weights) << "\n";

    int order_of_polynomial = 10;
    int N[order_of_polynomial];
    for(int i=0;i<order_of_polynomial;i++)
    {
        N[i] = i+1;
    }
    double integral_result_array_lobatto[order_of_polynomial];
    double integral_result_array_legendre[order_of_polynomial];

    for(int i=0;i<order_of_polynomial;i++)
    {
        double legendre_points[N[i]+1];
        double legendre_weights[N[i]+1];
        double lobatto_points[N[i]+1];
        double lobatto_weights[N[i]+1];

        generate_legendre_points(N[i],legendre_points,legendre_weights);
        generate_lobatto_points(N[i],lobatto_points,lobatto_weights);
        integral_result_array_legendre[i] = quadrature(N[i],legendre_points,legendre_weights);
        integral_result_array_lobatto[i] = quadrature(N[i],lobatto_points,lobatto_weights);
    }

    put_in_csv("legendre_results.csv",order_of_polynomial,integral_result_array_legendre);
    put_in_csv("lobatto_results.csv",order_of_polynomial,integral_result_array_lobatto);
}    
