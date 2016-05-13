//
//  GOF_exact_pvalue.cpp
//
//  Created by Ryan Sun on 11/26/15.
//  Copyright © 2015 Ryan Sun. All rights reserved.
//


// The syntax for this is './GOF_exact_pvalue [num_bounds] [bounds_file] [cor_file] [method]
// where method is an integer 0-4.
// Method 0 is the standard exact calculation, do not use with num_bounds>5.
// Method 1 approximates the correlation matrix as exchangeable, saves d! permutations of order.
// Method 2 is Method 1 AND we don't calculate all the terms in the |Z| PDF - instead
// just do 2^d times the first term.
// Method 3 doesn't permute the ordering of the terms (aka doesn't permute the correlation matrix)
// at all, saves d! permutations.
// Method 4 is Method 3 AND we don't calculate all the terms in the |Z| PDF.

// Update on 5/12/16 - We get rid of all the matrix algebra (no longer need armadillo) and calculate
// all the relevant correlation matrices without the help of preexisting functions). Before this change,
// a p-value for d=5 took approximately 50 seconds on local mac.
// However after the update it's even slower! Approx 1 min. :(
// Then we updated it again to split up the reordering and +/- loops, now down to ~48 seconds for d=5.


#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <ctime>
#include "expandBinary.h"
#include "mvtnorm.h"

double factorial(int x);
double exact_calc(const std::vector<double> &bounds,
                  const std::vector<double> &cor_vec);
double exact_calc_exch(const std::vector<double> &bounds,
		                  const std::vector<double> &cor_vec);
double exact_calc_exch_noS(const std::vector<double> &bounds,
		                      const std::vector<double> &cor_vec);
double exact_calc_noperm(const std::vector<double> &bounds,
		                      const std::vector<double> &cor_vec);
double exact_calc_noperm_noS(const std::vector<double> &bounds,
		                      const std::vector<double> &cor_vec);

// Simple factorial function (integers only).
double factorial(int x)
{
    if (x == 0)
        return 1;
    else
        return x*factorial(x-1);
}
 
// The exact p-value calculation, given bounds and correlation
// vector where element ( j + (i-2)(i-1)/2 ) is the (i,j) element
// of the covariance matrix (default from upper.tri()).
double exact_calc(const std::vector<double> &bounds,
                  const std::vector<double> &cor_vec)
{
    int d = bounds.size();
    double num_outer_loops = factorial(d);
    double num_inner_loops = pow(2.0, d);
    
    // Hold result (and err) of each ordered integration
    std::vector<double> loop_sums(num_outer_loops);
    std::vector<double> loop_errs(num_outer_loops);
    
    // Make the upper and lower integration bounds from
	// the GOF bounds.
	// We use 999 for infinity.
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
   
	// Make the class_S, all permutations of +/-
	std::vector<int> build_vec;
	std::vector<std::vector<int> >class_S;
	int expanded = expand_binary(d, build_vec, class_S);
						
	if (expanded!=0) {
		std::cout << "Error creating class_S!";
	    return 1;
	}
   

	// Set the first permutation order - 1,2,....,d
    std::vector<int> class_A(d);
    for (int temp_it=0; temp_it<d; ++temp_it)
    {
        class_A[temp_it] = temp_it + 1;
    }
    
    // Outer loop shuffling the order of observations.
	// We shuffle at the end of the loop.
    int new_iii;            // for keeping track of shuffling order
    int new_jjj;
    double plus_term;		// for multiplication by Delta_p
    double minus_term;
    for (double outer_it=0; outer_it<num_outer_loops; ++outer_it)
    {
      	double temp_sum_prob = 0;
      	double temp_sum_err = 0;
      	
      	// Correlation vector for new ordering 
      	// We will reference this always when doing the +/- permutation
		std::vector<double> new_cor_A(cor_vec.size());
      	
      	// Shuffle the variance matrix elements to account for new order
		for (int iii=2; iii<=d; ++iii) 
		{
			for (int jjj=1; jjj<=(iii-1); ++jjj) 
			{          
				// In our labeling system, i>j always.
				// So if choosing between (1,3) and (3,1), look for the (3,1) element.	
				new_iii = std::max(class_A[iii-1], class_A[jjj-1]);
				new_jjj = std::min(class_A[iii-1], class_A[jjj-1]);
					
				// Swap in correct element.	
                new_cor_A[jjj + ((iii-2)*(iii-1))/2 - 1] =
                cor_vec[new_jjj + ((new_iii-2)*(new_iii-1))/2 - 1];
        	}		// Done permuting variance matrix
        }
        
      	
		// Inner loop shuffles the permutation of +/-
		for (double inner_it=0; inner_it<num_inner_loops; ++inner_it)
		{
			// New correlation vector for +/-
			std::vector<double> final_cor_A_S(cor_vec.size());
        
			// Shuffle the variance matrix elements to account for new order
			for (int iii=2; iii<=d; ++iii) 
			{
				for (int jjj=1; jjj<=(iii-1); ++jjj) 
				{          

                	// Check if it's a + or - term
                	if ( (class_S[inner_it][iii-1] + class_S[inner_it][jjj-1]) == 1)
                	{
                		final_cor_A_S[jjj + ((iii-2)*(iii-1))/2 - 1] = -1 * new_cor_A[jjj + ((iii-2)*(iii-1))/2 - 1];
                	} else {
                		final_cor_A_S[jjj + ((iii-2)*(iii-1))/2 - 1] = new_cor_A[jjj + ((iii-2)*(iii-1))/2 - 1];
                	}
            	}
        	}		// Done permuting variance matrix
        
        
        	// Build bottom left and bottom right matrices
        	std::vector<std::vector<double>> bottom_left((d-1), std::vector<double>(d));
        	std::vector<std::vector<double>> bottom_right((d-1), std::vector<double>(d-1));
        	 
        	// Do the multiplication by Delta_p to get bottom_left
        	// Delta_p looks like [-1	1	0	0
        	//						0	-1	1	0
        	//						0	0	-1	1]
        	// And we want \Delta_p %*% \Sigma_A_S
        	
        	for (int iii=1; iii<=(d-1); ++iii) 
			{
				for (int jjj=1; jjj<=d; ++jjj) 
				{          
					// Get the minus term first
					new_iii = std::max(iii, jjj);
					new_jjj = std::min(iii, jjj);
					
					// Diagonal elements are 1
					if (new_iii == new_jjj) {
						minus_term = -1;
					} else {
						minus_term = -1 * final_cor_A_S[new_jjj + ((new_iii-2)*(new_iii-1))/2 - 1];
					}
					
					// Now get the plus term
					new_iii = std::max(iii+1, jjj);
					new_jjj = std::min(iii+1, jjj);
					if (new_iii == new_jjj) {
						plus_term = 1;
					} else {
						plus_term = final_cor_A_S[new_jjj + ((new_iii-2)*(new_iii-1))/2 - 1];
					}
					
					// Remember the matrix starts from index [0,0]
					bottom_left[iii-1][jjj-1] = plus_term + minus_term;
            	}
        	}		// Done with bottom_left
        	
        	// Multiply bottom_left by t(Delta_p) to get bottom_right
        	for (int iii=0; iii<(d-1); ++iii)
        	{
        		for (int jjj=0; jjj<=iii; ++jjj)
        		{
        			bottom_right[iii][jjj] = bottom_left[iii][jjj+1] - bottom_left[iii][jjj];
        		}
        	}
        	
        	 // Add to final_cor_s_a
            // Divide to get the variances of 1
            // It's clear that we need to divide!!! or else the variance matrix is not p.s.d
            for (int temp_it=0; temp_it<(d-1); ++temp_it)
            {
                for (int second_it=0; second_it<d; ++second_it)
                    final_cor_A_S.push_back(bottom_left[temp_it][second_it] / sqrt(bottom_right[temp_it][temp_it]));
                
                for (int second_it=0; second_it<temp_it; ++second_it)
                    final_cor_A_S.push_back(bottom_right[temp_it][second_it]
                                            / (sqrt(bottom_right[temp_it][temp_it])*sqrt(bottom_right[second_it][second_it])));
            }
        	
        	
            // Perform the integration
            double integration_error;
            double* lower_array = &lower_bound[0];
            double* upper_array = &upper_bound[0];
            double* cor_array = &final_cor_A_S[0];
            // std::copy(lower_bound.begin(), lower_bound.end(), lower_array);
            // std::copy(upper_bound.begin(), upper_bound.end(), upper_array);
            // std::copy(final_cor_A_S.begin(), final_cor_A_S.end(), cor_array);
           
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
    for(int temp_it=0; temp_it<num_outer_loops; ++temp_it) {
        p_value -= loop_sums[temp_it];
    }
    
    return p_value;
}





////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
// Main function, choose which calculation to use and 
// read in parameters.

int main(int argc, const char * argv[]) {

		// first trailing argument is the number of time points t_k
    // second is the file containing the bounds
    // third file is the file containing the pairwise correlations
		// fourth file is the calculation method (exact or one of the approximations)
    int d = atoi(argv[1]);
    long num_cor = d*(d-1)/2;

		int method = atoi(argv[4]);

    // Dynamically allocate arrays for time points and all pairwise correlations
    std::vector<double> boundaryPts;
    std::vector<double> cors_vec;

    // Reserve memory to prevent fragmentation, easy since we know the exact size
    boundaryPts.reserve(d);

    // Open the boundary file for reading, last two trailing arguments give the name
    std::ifstream bounds_f(argv[2]);
    if (!bounds_f)
    {
        std::cerr << "Can't open bounds file" << std::endl;
        exit(1);
    }

    // Put the boundary pts into array, line by line
    double temp_input;
    for (int iii=0; iii<d; ++iii)
    {
        bounds_f >> temp_input;
        boundaryPts.push_back(temp_input);
    }

	  // If it's a filename, then reserve space for the correlation vector
    cors_vec.reserve(num_cor);

    // Open the correlations file
    std::ifstream cor_f(argv[3]);
    if (!cor_f)
    {
        std::cerr << "Can't open correlations file" << std::endl;
        exit(1);
    }

    // Put the correlations into array, line by line
    for (long iii=0; iii<num_cor; ++iii)
    {
        cor_f >> temp_input;
        cors_vec.push_back(temp_input);
    }

    //std::vector<double> bounds {0.4325330, 0.7598041, 1.0777005, 1.4462875, 1.9759605};
    //std::vector<double> cor_vec {0.01912024, 0.15498892, 0.12003188, 0.21301551, 0.14181763, 0.21795605, 0.12711216, 0.23954297, 0.17157866, 0.10482516};
    //std::vector<double> cor_vec {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
    
		
		//clock_t begin = std::clock();
		double result;

		// Choose which calculation method to use
		switch (method) 
		{
		case 0:
			result = exact_calc(boundaryPts, cors_vec);
			break;

		default:
			std::cout << "Invalid choice of method!";
			return 1;
		}
  	
    //clock_t end = std::clock();
    //double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    //std::cout << elapsed_secs << std::endl;
    std::cout << result << std::endl;
    
		return 0;
}


