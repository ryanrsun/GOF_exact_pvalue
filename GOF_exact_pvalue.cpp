//
//  main.cpp
//  testScripts
//
//  Created by Ryan Sun on 11/26/15.
//  Copyright © 2015 Ryan Sun. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <vector>
#include <armadillo>
#include <ctime>
#include "expandBinary.h"
#include "mvtnorm.h"

using namespace arma;

double factorial(int x);

double exact_calc(const std::vector<double> &bounds,
                  const std::vector<double> &cor_vec);



// Simple factorial function (integers only).
double factorial(int x)
{
    if (x == 0)
        return 1;
    else
        return x*factorial(x-1);
}
 

double exact_calc(const std::vector<double> &bounds,
                  const std::vector<double> &cor_vec)
{
    int d = bounds.size();
    double num_outer_loops = factorial(d);
    double num_inner_loops = pow(2.0, d);
    
    // To hold results
    std::vector<double> loop_sums(num_outer_loops);
    std::vector<double> loop_errs(num_outer_loops);
    
    // Make the upper and lower bounds
    std::vector<double> lower_bound((2*d-1), 0);
    std::vector<double> upper_bound((2*d-1), 999);
    
    for (int temp_it=1; temp_it<d; ++temp_it)
    {
        lower_bound[temp_it] = -999;
    }
    for (int temp_it=0; temp_it<d; ++temp_it)
    {
        upper_bound[temp_it] = bounds[temp_it];
    }
    
    // Random seed before integration starts
    srand(time(NULL));
    
    // All permutations of order of observations
    std::vector<int> class_A(d);
    for (int temp_it=0; temp_it<d; ++temp_it)
    {
        class_A[temp_it] = temp_it + 1;
    }
    
    // Outer loop shuffling the order of observations
    int new_iii;            // for keeping track of shuffling order
    int new_jjj;
    for (double outer_it=0; outer_it<num_outer_loops; ++outer_it)
    {
        // class_A is a vector holding our 'permuted key'
        
        std::vector<double> new_cor_A(cor_vec.size());
        
        // Shuffle the variance matrix elements to account for new order
        for (int iii=2; iii<=d; ++iii) {
            for (int jjj=1; jjj<=(iii-1); ++jjj) {
                
                if (class_A[iii-1] > class_A[jjj-1])
                {
                    new_iii = class_A[iii-1];
                    new_jjj = class_A[jjj-1];
                }
                else
                {
                    new_iii = class_A[jjj-1];
                    new_jjj = class_A[iii-1];
                }
                
                new_cor_A[jjj + ((iii-2)*(iii-1))/2 - 1] =
                cor_vec[new_jjj + ((new_iii-2)*(new_iii-1))/2 - 1];
            }
        }		// Done permuting variance matrix
        
        // Make the class_S
        std::vector<int> build_vec;
        std::vector<std::vector<int> >class_S;
        int expanded = expand_binary(d, build_vec, class_S);
        
        if (expanded!=0) {
            std::cout << "Error creating class_S!";
            return 1;
        }
        
        // Inner loop over all permutations of S
        double temp_sum_prob = 0;
        double temp_sum_err = 0;
        
        for (double inner_it=0; inner_it<num_inner_loops; ++inner_it)
        {
            
            std::vector<double> final_cor_A_S(cor_vec.size());
            
            // Perform the +/- permutation of Sigma_A
            // At the same time, build the variance matrix back up
            mat v_mat(d,d);
            
            for (int iii=2; iii<=d; ++iii) {
                for (int jjj=1; jjj<=(iii-1); ++jjj) {
                    
                    if ( (class_S[inner_it][iii-1] + class_S[inner_it][jjj-1])==1 )
                    {
                        final_cor_A_S[jjj + ((iii-2)*(iii-1))/2 - 1] =
                        -1 * new_cor_A[jjj + ((iii-2)*(iii-1))/2 - 1];
                    } else
                    {
                        final_cor_A_S[jjj + ((iii-2)*(iii-1))/2 - 1] =
                        1 * new_cor_A[jjj + ((iii-2)*(iii-1))/2 - 1];
                    }
                    
                    // build
                    v_mat(iii-1,jjj-1) = final_cor_A_S[jjj + ((iii-2)*(iii-1))/2 - 1];
                    v_mat(jjj-1,iii-1) = final_cor_A_S[jjj + ((iii-2)*(iii-1))/2 - 1];
                }
            }
            
            //  build diagonals
            for (int temp_it=0; temp_it<d; ++temp_it)
            {
                v_mat(temp_it, temp_it) = 1;
            }
            
            // Account for delta_p
            mat Delta_p(d-1,d);
            Delta_p.zeros();
            
            for (int temp_it=0; temp_it<(d-1); ++temp_it)
            {
                Delta_p(temp_it, temp_it) = -1;
                Delta_p(temp_it, temp_it+1) = 1;
            }
        
            
            // bottom left addition
            mat bottom_left = Delta_p*v_mat;
            // bottom right
            mat bottom_right = bottom_left*Delta_p.t();
            
            // The variances of the (p*1)x(p*1) lower right corner
            // Divide to make these unity
            std::vector<double> lower_vars(d-1);
            for (int temp_it=0; temp_it<(d-1); ++temp_it)
                lower_vars[temp_it] = bottom_right(temp_it, temp_it);
            
            // Add to final_cor_s_a
            // Divide to get the variances of 1
            // It's clear that we need to divide!!! or else the variance matrix is not p.s.d
            for (int temp_it=0; temp_it<(d-1); ++temp_it)
            {
                for (int second_it=0; second_it<d; ++second_it)
                    final_cor_A_S.push_back(bottom_left(temp_it, second_it) / sqrt(lower_vars[temp_it]));
                    //final_cor_A_S.push_back(bottom_left(temp_it, second_it));
                
                for (int second_it=0; second_it<temp_it; ++second_it)
                    final_cor_A_S.push_back(bottom_right(temp_it, second_it)
                                            / (sqrt(lower_vars[temp_it])*sqrt(lower_vars[second_it])));
                    //final_cor_A_S.push_back(bottom_right(temp_it, second_it));
            }
            
            // Perform the integration
            double integration_error;
            double lower_array[lower_bound.size()];
            double upper_array[upper_bound.size()];
            double cor_array[final_cor_A_S.size()];
            std::copy(lower_bound.begin(), lower_bound.end(), lower_array);
            std::copy(upper_bound.begin(), upper_bound.end(), upper_array);
            std::copy(final_cor_A_S.begin(), final_cor_A_S.end(), cor_array);
            
            /*
            double sigma[] = {0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1,  0.1, -0.9,  0.9,  0.0,  0.0,  0.0,  0.0, -0.9,  0.9,  0.0,  0.0, -0.9, 0.0, 0.0, -0.9,  0.9,  0.0,  0.0, -0.9,  0.0,  0.0,  0.0, -0.9,  0.9,  0.0,  0.0, -0.9};
            //double lower[] = {0, -999, -999, -999, -999, 0,0,0,0};
            //double upper[] = {0.4325330, 0.7598041, 1.0777005, 1.4462875, 1.9759605, 999, 999, 999, 999};
            double sigma[] = {1.0/10, -9.0/(10*sqrt(1.8)),  9.0/(10*sqrt(1.8))};
            double lower[] = {0, -999, 0};
            double upper[] = {1.081541, 1.776954, 999 };
            double error;
            double ret = pmvnorm_B(3, lower, upper, sigma, &error);
            */
            
            double integral_result = pmvnorm_B((2*d-1), lower_array, upper_array, cor_array, &integration_error);
            //double integral_result = pmvnorm_B(3, lower_bound, upper_bound, cor_vec, &integration_error);
            
            temp_sum_prob += integral_result;
            temp_sum_err += integration_error;
            
        }		// inner loop, for each S	
        
        loop_sums[outer_it] = temp_sum_prob;
        loop_errs[outer_it] = temp_sum_err;
        
        // Permute Again
        std::next_permutation(class_A.begin(), class_A.end());
    } // outer loop, for each A
    
    
    // The p-value
    double p_value = 1;
    for(std::vector<double>::iterator pval_it = loop_sums.begin(); pval_it != loop_sums.end(); ++pval_it) {
        /* std::cout << *it; ... */
        p_value -= loop_sums[*pval_it];
    }
    
    return p_value;
}





int main(int argc, const char * argv[]) {
    
    std::vector<double> bounds {0.4325330, 0.7598041, 1.0777005, 1.4462875, 1.9759605};
    std::vector<double> cor_vec {0.01912024, 0.15498892, 0.12003188, 0.21301551, 0.14181763, 0.21795605, 0.12711216, 0.23954297, 0.17157866, 0.10482516};
    //std::vector<double> cor_vec {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    
    clock_t begin = std::clock();
    
    double result = exact_calc(bounds, cor_vec);
    
    clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout << elapsed_secs << std::endl;
    std::cout << result << std::endl;
    
    
    
    
    
    
    
    
    
    

return 0;
}

/*

int main(int argc, const char * argv[]) {
    
 
    std::vector<int> myvec;     // just pass something in
    int p =4;
    std::vector<std::vector<int>> class_S;
    class_S.reserve(pow(2,p));
    mutate(p, myvec, class_S);
    
    for(int jjj=0; jjj<pow(2,p); ++jjj)
    {
        for(int kkk=0; kkk<p; ++kkk)
        {
            std::cout << class_S[jjj][kkk];
        }
        std::cout << std::endl;
    }
    return 0;
 
    
    for (int i=0; i<0; ++i)
        std::cout << i << std::endl;
    
    for (int i=0; i<1; ++i)
        std::cout << i << std::endl;
    
    
    cout << "Armadillo version: " << arma_version::as_string() << endl;
    
    mat A(3,3);  // directly specify the matrix size (elements are uninitialised)
    
    // endr indicates "end of row"
    A << 0.99427048 << 0.04304764 << 0.09784227 << endr
    << 0.04304764 << 0.98774438 << 0.15002645 << endr
    << 0.09784227 << 0.15002645 << 0.98382872 << endr;
    
    A(1,1) = 1;
    A(2,2) = 2;
    A.print("A:");

    return 0;
    
    double num_sims = 10000000;
    
    srand(time(NULL));
    
    mat B(num_sims, 3, fill::randn);
    
    std::cout << "hi" << std::endl;
    
    mat E = B*A;
    
    std::cout << "hi" << std::endl;
    
    double yescounter = 0;
    double notcounter = 0;
    std::vector<double> lower(3);
    std::vector<double> upper(3);
    lower[0]= -1;
    lower[1] = -2;
    lower[2] = -3;
    upper[0]= 1;
    upper[1] = 2;
    upper[2] = 3;
    
    for(double i=0; i<num_sims; ++i)
    {
        bool inside = true;
        for(int j=0; j<3; ++j)
        {
            //std::cout << E(i,j) << " ";
            if (E(i,j) < lower[j] || E(i,j) > upper[j]) {
                inside = false;
            }
        }
        //std::cout << std::endl;
        if (inside==true)
        {
            yescounter = yescounter + 1;
        } else {notcounter = notcounter + 1;}
    }
    
    std::cout << yescounter/num_sims << std::endl;
    std::cout << notcounter/num_sims << std::endl;
    
    return 0;
}
 
*/

