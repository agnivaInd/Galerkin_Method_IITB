#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

const double PI = 2*acos(0);
const int DP = 6;                   // No of decimal places round off

double analytical_interpolation_function(double argument)
{
    return sin(PI*argument/2);
}

double analytical_derivative_interpolation_function(double argument)
{
    return (PI/2)*cos(PI*argument/2);
}

double round_off(double value, int decimal_places)
{
    double factor = pow(10,decimal_places);
    return round(value*factor)/factor;
}

void bubble_sort(int N, double points[])
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
            }
        }
    }
}

void generate_equally_spaced_points(double x_start, double x_end, int N, double points[])
{
    double delta_x = (x_end - x_start)/N;
    double counter = x_start;
    for(int i=0;i<N+1;i++)
    {
        points[i] = round_off(counter,DP);
        counter = counter + delta_x;
    }
}

void generate_legendre_points(int N, double points[])
{
    double tolerance = pow(10,-6);
    for(int i=0;i<N+1;i++)
    {
        double initial_guess = cos((2.0*i+1.0)*PI/(2.0*N+2.0));
        double residue = 1;
        
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
            double derivative_phi_n = (N+1)*(old_phi_n_minus_1 - initial_guess*phi_n)/(1.0 - pow(initial_guess,2));
            
            double x_update = initial_guess - phi_n/derivative_phi_n;
            residue = fabs(x_update - initial_guess);

            initial_guess = x_update;
        }
        points[i] = round_off(initial_guess,DP);
    }
    bubble_sort(N,points);
}

void generate_lobatto_points(int N, double points[])
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

                initial_guess = x_update;
            }
            points[i] = round_off(initial_guess,DP);
        }
    }
    bubble_sort(N,points);
}

// Function to construct lagrange polynomial and evaluate interpolation value
void interpolation_result(int number_of_int_points, int N, double points[], double final_interpolation_result_points[], double final_interpolation_result[])
{
    for(int k=0;k<number_of_int_points;k++)
    {
        double sum = 0;
        for(int j=0;j<N+1;j++)
        {
            double product = 1;
            for(int i=0;i<N+1;i++)
            {
                if(i!=j)
                {
                    product = product*(final_interpolation_result_points[k] - points[i])/(points[j] - points[i]);
                }
            }
            sum = sum + product*analytical_interpolation_function(points[j]);
        }
        final_interpolation_result[k] = sum;
    }
}

void derivative_interpolation(int number_of_int_points, int N, double points[], double final_interpolation_result_points[], double derivative_interpolation_result[])
{
    for(int m=0;m<number_of_int_points;m++)
    {
        double main_sum = 0;
        for(int i=0;i<N+1;i++)
        {
            double inner_sum = 0;
            for(int k=0;k<N+1;k++)
            {
                if(k!=i)
                {
                    double product = 1;
                    for(int j=0;j<N+1;j++)
                    {
                        if((j!=i)&&(j!=k))
                        {
                            product = product*(final_interpolation_result_points[m] - points[j])/(points[i] - points[j]);
                        }
                    }
                    inner_sum = inner_sum + product/(points[i] - points[k]);        //  dL_i/dx value
                }
            }
            main_sum = main_sum + inner_sum*analytical_interpolation_function(points[i]);
        }
        derivative_interpolation_result[m] = main_sum;
    }
}

double L2_norm(int number_of_int_points, double numerical_result[], double final_interpolation_result_points[])
{
    double numerator_sum = 0;
    double denominator_sum = 0;
    for(int k=0;k<number_of_int_points;k++)
    {
        numerator_sum = numerator_sum + pow(numerical_result[k] - analytical_interpolation_function(final_interpolation_result_points[k]),2);
        denominator_sum = denominator_sum + pow(analytical_interpolation_function(final_interpolation_result_points[k]),2);
    }
    return sqrt(numerator_sum/denominator_sum);     
}

void put_in_csv(const std::string& filename, int number_of_int_points, double result[])
{
    std::ofstream file(filename);
    for (int i=0;i<number_of_int_points;i++)
    {
        file << result[i];
        if (i < number_of_int_points-1)
        {
            file << ",";
        }
    }
    file << std::endl;
    file.close();
}

int main()
{
    double x_start = -1;
    double x_end = 1;
    int N = 3;                      // N = Interpolation Order
    
    double equally_spaced_points[N+1];       // N+1 should be number of points used for interpolation
    double legendre_points[N+1];
    double lobatto_points[N+1];
    
    int number_of_int_points = 50;  // Final interpolation number of grid points
    double final_interpolation_result_points[number_of_int_points];
    double global_counter = x_start;
    double global_delta_x = (x_end - x_start)/(number_of_int_points - 1);

    for(int i=0;i<number_of_int_points;i++)
    {
        final_interpolation_result_points[i] = global_counter;
        global_counter = global_counter + global_delta_x;
    }

    double result_equally_spaced_points_interpolation[number_of_int_points];
    double result_legendre_points_interpolation[number_of_int_points];
    double result_lobatto_points_interpolation[number_of_int_points];

    generate_equally_spaced_points(x_start,x_end,N,equally_spaced_points);
    generate_legendre_points(N,legendre_points);
    generate_lobatto_points(N,lobatto_points);

    interpolation_result(number_of_int_points,N,equally_spaced_points,final_interpolation_result_points,result_equally_spaced_points_interpolation);
    interpolation_result(number_of_int_points,N,legendre_points,final_interpolation_result_points,result_legendre_points_interpolation);
    interpolation_result(number_of_int_points,N,lobatto_points,final_interpolation_result_points,result_lobatto_points_interpolation);
    
    double result_derivative_equally_spaced_interpolation[number_of_int_points];
    double result_derivative_legendre_interpolation[number_of_int_points];
    double result_derivative_lobatto_interpolation[number_of_int_points];

    derivative_interpolation(number_of_int_points,N,equally_spaced_points,final_interpolation_result_points,result_derivative_equally_spaced_interpolation);
    derivative_interpolation(number_of_int_points,N,legendre_points,final_interpolation_result_points,result_derivative_legendre_interpolation);
    derivative_interpolation(number_of_int_points,N,lobatto_points,final_interpolation_result_points,result_derivative_lobatto_interpolation);

    /*
        Testing:

        cout << "x-value" << " " << "Equal" << " " << "Legendre" << " " << "Lobatto" << "\n";

        for(int k=0;k<number_of_int_points;k++)
        {
            cout << final_interpolation_result_points[k] << " " << result_equally_spaced_points_interpolation[k] << " "<< result_legendre_points_interpolation[k] << " " << result_lobatto_points_interpolation[k] << "\n";
        }
        cout << "\n\n";
        for(int k=0;k<number_of_int_points;k++)
        {
            cout << analytical_derivative_interpolation_function(final_interpolation_result_points[k]) << " " << result_derivative_interpolation[k] << "\n";
        }
    */

    put_in_csv("equally_spaced_interpolation.csv",number_of_int_points,result_equally_spaced_points_interpolation);
    put_in_csv("legendre_interpolation.csv",number_of_int_points,result_legendre_points_interpolation);
    put_in_csv("lobatto_interpolation.csv",number_of_int_points,result_lobatto_points_interpolation);

    put_in_csv("derivative_equally_spaced.csv",number_of_int_points,result_derivative_equally_spaced_interpolation);
    put_in_csv("derivative_legendre.csv",number_of_int_points,result_derivative_legendre_interpolation);
    put_in_csv("derivative_lobatto.csv",number_of_int_points,result_derivative_lobatto_interpolation);

    cout << "L2 norm for Equally Spaced Points: " << L2_norm(number_of_int_points,result_equally_spaced_points_interpolation,final_interpolation_result_points) << "\n";
    cout << "L2 norm for Legendre Points: " << L2_norm(number_of_int_points,result_legendre_points_interpolation,final_interpolation_result_points) << "\n";
    cout << "L2 norm for Lobatto Points: " << L2_norm(number_of_int_points,result_lobatto_points_interpolation,final_interpolation_result_points) << "\n";

    cout << "Interpolation of Function and derivative evaluated successfully\n"; 
}