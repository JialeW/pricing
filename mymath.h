//
//  Header.h
//  triniomial_tree
//
//  Created by jiale on 3/13/16.
//  Copyright © 2016 jiale. All rights reserved.
//

#ifndef Header_h
#define Header_h
#include <cmath>
#include "matrix.h"
double erf(double x)
{
    double y = 1.0 / ( 1.0 + 0.3275911 * x);
    return 1 - (((((
                    + 1.061405429  * y
                    - 1.453152027) * y
                   + 1.421413741) * y
                  - 0.284496736) * y
                 + 0.254829592) * y)
    * exp (-x * x);      
}

double cdf(double x, double mu, double sigma)
{
    return 0.5 * (1 + erf((x - mu) / (sigma * sqrt(2.))));
}

double normalcdf(const double& x)
{
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));
    
    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - normalcdf(-x);
    }
}


#endif /* Header_h */
