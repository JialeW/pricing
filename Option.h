
//  Option.h
//  triniomial_tree
//
//  Created by jiale on 3/11/16.
//  Copyright © 2016 jiale. All rights reserved.
//

#ifndef Option_h
#define Option_h

#include <vector>
#include "matrix.h"
#include "mymath.h"
using  std::vector;
using  std::string;



double trinomial_price(double s0, double K, double Time,double r,double d,double sigma,int N,double dx,string call_put_flag)
{
    
    double dt = Time/N;
    double nu = r - d - sigma*sigma/2;
    dx =  sqrt(sigma*sigma*dt + nu*dt*nu*dt);
    double edx = exp(dx);
    double pu = 0.5 * ( (sigma*sigma*dt  + nu*nu*dt*dt) /(dx*dx) + nu*dt/dx );
    double pm = 1.0 - (sigma*sigma*dt  + nu*nu*dt*dt) / (dx*dx) + r*dt;
    double pd = 0.5 * ( (sigma*sigma*dt  + nu*nu*dt*dt) /(dx*dx) - nu*dt/dx );
    double disc = exp(-r*dt);
    // initialise asset prices at maturity
    vector <double> s_t;
    vector <double> c; // option values at maturity
    s_t.push_back(s0 * exp(0.0 - N*dx));
    for (int i = 0; i < 2*N ; ++i)
    {
        c.push_back( call_put_flag=="call"? (0 > s_t[i]-K ? 0:s_t[i]-K)  :  (0 > K-s_t[i] ? 0:K-s_t[i]) );
        s_t.push_back(s_t.back() * edx);
    }
    c.push_back( call_put_flag=="call"? (0 > s_t.back()-K ? 0:s_t.back()-K)  :  (0 > K-s_t.back() ? 0:K-s_t.back()) );
    // option values
    for (int j = N - 1; j >= 0; --j)
    {
        for (int i = 0; i <= 2*j; ++i)
        {
            c[i] = disc*(pd*c[i] + pm*c[i+1] + pu*c[i+2]);
        }
    }
    return  c[0];
}


double Explicit(double s0, double K,double Time,double r,double d,double sigma,int N,int Nj,double dx,string call_put_flag)
{
    double dt = Time/N;
    double nu = r - d - sigma*sigma/2;
    dx =  sqrt(sigma*sigma*dt + nu*dt*nu*dt);
    double edx = exp(dx);
    double pu = 0.5 * dt * ( (sigma/dx)*(sigma/dx) + nu/dx );
    double pm = 1.0 - dt*(sigma/dx)*(sigma/dx) - r*dt;
    double pd = 0.5 * dt * ( (sigma/dx)*(sigma/dx) - nu/dx );
    double disc = exp(-r*dt);
    // initialise asset prices at maturity
    vector <double> s_t;
    vector <double> c; // option values at maturity
    s_t.push_back(s0 * exp(0.0 - N*dx));
    for (int i = 0; i < 2*N ; ++i)
    {
        c.push_back( call_put_flag=="call"? (0 > s_t[i]-K ? 0:s_t[i]-K)  :  (0 > K-s_t[i] ? 0:K-s_t[i]) );
        s_t.push_back(s_t.back() * edx);
    }
    c.push_back( call_put_flag=="call"? (0 > s_t.back()-K ? 0:s_t.back()-K)  :  (0 > K-s_t.back() ? 0:K-s_t.back()) );
    vector<double> temp_c = c;
    Nj = N;
    // option values
    for (int i = N - 1; i >= 0; --i)
    {
        for (int j = 1; j < 2*N; ++j)
        {
            c[j] = disc*(pd*temp_c[j-1] + pm*temp_c[j] + pu*temp_c[j+1]);
        }
        if(call_put_flag=="call")
        {
        c[0] = c[1];
        c[2*N] = c[2*N-1] + s_t[2*N] - s_t[2*N-1];
        }
        else
        {
            c[2*N] = c[2*N-1];
            c[0] = c[1] + s_t[1] - s_t[0];
        }
        temp_c = c;
    }
    return  c[N];
}


matrix<float> inv_coefficient_matrix(int N,float pu,float pm,float pd)
{
    int Nj = 2*N + 1;
    vector<float> a = {1.0,-1.0};
    vector<float> b(Nj-2,0.0);
    vector<float> c = {pu,pm,pd};
    c.insert(c.end(),b.begin(),b.end());// pu,pm,pd,0*(Nj-2)
    vector<float> coefficient_vector;
    coefficient_vector.insert(coefficient_vector.end(), a.begin(), a.end());// add 1,-1
    coefficient_vector.insert(coefficient_vector.end(), b.begin(), b.end());// add 0*(Nj-2) , finished the first row
    for (int i = 0; i < Nj -2; ++ i)
    {
        coefficient_vector.insert(coefficient_vector.end(), c.begin(), c.end());
    }
    coefficient_vector.insert(coefficient_vector.end(), a.begin(), a.end());
    matrix<float> coefficient_m(coefficient_vector, Nj,Nj);
    return coefficient_m.inverse();
}

float solve_implicit_tridiagonal_system(int n,int K, vector<float> c,float pu,float pm,float pd,float lambda_L,float lambda_U,string call_put_flag)
{
    matrix<float> m = inv_coefficient_matrix(n, pu,pm, pd);
    matrix<float> M ((int)c.size(),0);
    std::reverse(c.begin(), c.end());
    M.add_col(0,c);
    vector<float> lambda_u ={lambda_U};
    vector<float> lambda_l ={lambda_L};
    for(int j = 0; j < n; ++j)
    {
// delete the first and last row of original option values matrix M
        M.del_row(0);
        M.del_row(M.get_row() - 1);
// add the lambda_U and lambda_L into the option values matrix M
        M.add_row(0, lambda_u);
        M.add_row_InTheEnd(lambda_l);
        M = m * M;
    }
        return M.m[n][0];
}

float Implicit(float s0, float K, float Time,float r,float d,float sigma,int N,int Nj,float dx,string call_put_flag)
{
    Nj = N;
    float dt = Time/N;
    float nu = r - d - sigma*sigma/2;
    dx =  sqrt(sigma*sigma*dt + nu*dt*nu*dt);
    float edx = exp(dx);
    float pu = 0.0 - 0.5 * dt * ( (sigma/dx)*(sigma/dx) + nu/dx );
    float pm = 1.0 + dt*(sigma/dx)*(sigma/dx) - r*dt;
    float pd = 0.0 - 0.5 * dt * ( (sigma/dx)*(sigma/dx) - nu/dx );
    // initialise asset prices at maturity
    vector <float> st;
    vector <float> c; // option values at maturity
    st.push_back(s0 * exp(0.0 - Nj*dx));
    // derivative boundary condition
    float lambda_U = 0.0;
    float lambda_L = 0.0;
    // option values at maturity
    
    if ( call_put_flag=="call")
    {
        for (int i = 0; i < 2*Nj ; ++i)
        {
            c.push_back(0 > st[i]-K ?  0:st[i]-K);
            st.push_back(st[i] * edx);
        }
        c.push_back(0 > st.back()-K ?  0:st.back()-K);
        lambda_U = st[2*Nj] -  st[2*Nj-1];
    }
    else
    {
        for (int i = 0; i < 2*Nj ; ++i)
        {
            c.push_back(0 > K-st[i] ? 0:K-st[i]);
            st.push_back(st.back() * edx);
        }
        c.push_back(0 > K-st.back() ? 0:K-st.back());
        lambda_L = st[0] - st[1];
    }
    // step back through lattice
    return  solve_implicit_tridiagonal_system( N, K, c, pu, pm, pd, lambda_L, lambda_U,call_put_flag);
}
// ============================================================================
matrix<float> inv_coefficient_matrix2(int N,float pu,float pm,float pd)
{
    int Nj = 2*N + 1;
    vector<float> b(Nj-1,0.0);
    vector<float> coefficient_vector;
    coefficient_vector.push_back(1);// add 1,-1
    coefficient_vector.insert(coefficient_vector.end(), b.begin(), b.end());// add 0*(Nj-2) , finished the first row
    b.pop_back();
    vector<float> c = {-pu,-pm+2,-pd};
    c.insert(c.end(),b.begin(),b.end());// pu,pm,pd,0*(Nj-2)
    for (int i = 0; i < Nj -2; ++ i)
    {
        coefficient_vector.insert(coefficient_vector.end(), c.begin(), c.end());
    }
    coefficient_vector.push_back(0);
    coefficient_vector.push_back(1);
    matrix<float> coefficient_m(coefficient_vector, Nj,Nj);
    return coefficient_m;

}
float solve_CrankNicolson_tridiagonal_system(int n,int K, vector<float> c,float pu,float pm,float pd,float lambda_L,float lambda_U,string call_put_flag)
{
    matrix<float> m1 = inv_coefficient_matrix(n, pu,pm, pd);
    matrix<float> m2 = inv_coefficient_matrix2(n, pu,pm, pd);
    matrix<float> M ((int)c.size(),0);
    std::reverse(c.begin(), c.end());
    M.add_col(0,c);
    vector<float> lambda_u ={lambda_U};
    vector<float> lambda_l ={lambda_L};
    for(int j = 0; j < n; ++j)
    {
        // delete the first and last row of original option values matrix M
        M.del_row(0);
        M.del_row(M.get_row() - 1);
        // add the lambda_U and lambda_L into the option values matrix M
        M.add_row(0, lambda_u);
        M.add_row_InTheEnd(lambda_l);
        M = m1 * m2 * M;
    }
    return M.m[n][0];
}

float CrankNicolson(float s0, float K, float Time,float r,float d,float sigma,int N,int Nj,float dx,string call_put_flag)
{
    float dt = Time/N;
    float nu = r - d - sigma*sigma/2;
    dx =  sqrt(sigma*sigma*dt + nu*dt*nu*dt);
    float edx = exp(dx);
    float pu = - 0.25 * dt * ( (sigma/dx)*(sigma/dx) + nu/dx );
    float pm = 1.0 + 0.5*dt*(sigma/dx)*(sigma/dx) + 0.5*r*dt;
    float pd = 0.25 * dt * (  nu/dx - (sigma/dx)*(sigma/dx) );
    // initialise asset prices at maturity
    vector <float> st;
    vector <float> c; // option values at maturity
    st.push_back(s0 * exp(0.0 - Nj*dx));
    // derivative boundary condition
    float lambda_U = 0.0;
    float lambda_L = 0.0;
    // option values at maturity
    Nj = N;
    if ( call_put_flag=="call")
    {
        for (int i = 0; i < 2*Nj ; ++i)
        {
            c.push_back(0 > st[i]-K ?  0:st[i]-K);
            st.push_back(st[i] * edx);
        }
        c.push_back(0 > st.back()-K ?  0:st.back()-K);
        lambda_U = st[2*Nj] -  st[2*Nj-1];
    }
    else
    {
        for (int i = 0; i < 2*Nj ; ++i)
        {
            c.push_back(0 > K-st[i] ? 0:K-st[i]);
            st.push_back(st.back() * edx);
        }
        c.push_back(0 > K-st.back() ? 0:K-st.back());
        lambda_L = st[0] - st[1];
    }
    // step back through lattice
    return  solve_CrankNicolson_tridiagonal_system( N, K, c, pu, pm, pd, lambda_L, lambda_U,call_put_flag);
}


// M paths,N timesteps
float MC(float s0, float K, float Time,float r,float d,float sigma,int M,int N,string call_put_flag)
{
    float dt = Time/N;
    float nu = Time*(r - d - sigma*sigma/2);
    float sigsdt = sigma*sqrt(dt);
    float lns = log(s0);
    
    float sum_ct = 0;
    float sum_ct2 = 0;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);
    for (int i  = 0; i < M; ++i) {
        float lnst = 0;
        for (int j = 0; j < N; ++j) {
            lnst += distribution(generator);
        }
        float st = exp(lnst*sigsdt + lns + nu);
        float ct = call_put_flag =="call"?  (0 > st-K ? 0:st-K) : (0 > K-st ? 0:K-st);
        sum_ct += ct;
        sum_ct2 += ct*ct;
    }
    float option_value = exp(0.0 - r*Time) * sum_ct/M;
    float SD = sqrt( (sum_ct2 - sum_ct*sum_ct/M)*exp(0.0-2*r*Time)/(M-1) );
    std::cout << " option value is " << option_value << "\n";
    std::cout << " Standard deviation is " << SD << "\n";
    std::cout << " Standard error is " << SD/sqrt(M) << "\n";
    return option_value;
}


template <class T>
T MC2(T s0, T K, T Time,T r,T d,T sigma,int M,int N,string call_put_flag)
{
    T dt = Time/N;
    T nu = Time*(r - d - sigma*sigma/2);
    T sigsdt = sigma*sqrt(dt);
    T lns = log(s0);
    
    T sum_ct = 0;
    T sum_ct2 = 0;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);
    for (int i  = 0; i < M; ++i) {
        T lnst = 0;
        for (int j = 0; j < N; ++j) {
            lnst += distribution(generator);
        }
        T st = exp(lnst*sigsdt + lns + nu);
        T st2 = exp( lns + nu - lnst*sigsdt );
        T ct = call_put_flag =="call"?  (0 > st-K ? 0:st-K) : (0 > K-st ? 0:K-st);
        ct += call_put_flag =="call"?  (0 > st2-K ? 0:st2-K) : (0 > K-st2 ? 0:K-st2);
        ct /= 2;
        sum_ct += ct;
        sum_ct2 += ct*ct;
    }
    T option_value = exp(0.0 - r*Time) * sum_ct/M;
    T SD = sqrt( (sum_ct2 - sum_ct*sum_ct/M)*exp(0.0-2*r*Time)/(M-1) );
    std::cout << " option value is " << option_value << "\n\n";
    std::cout << " Standard deviation is " << SD << "\n";
    std::cout << " Standard error is " << SD/sqrt(M) << "\n";
    return option_value;
}
//  Black_Scholes to calcalate option
float Black_Scholes(float S, float K, float t,float T,float v,float r,float d,string call_put_flag)
{
    T = T - t;
    float T_sqrt = sqrt(T);
    float d1 = (log(float(S)/K)+((r-d)+v*v/2.)*T)/(v*T_sqrt);
    float d2 = d1-v*sqrt(T);
    if(call_put_flag=="call")
    {
        return S* cdf(d1,0,1)-K*exp(-r*T)* cdf(d2,0,1);
    }
  
    return  K*exp(-r*T)*cdf(-d2,0,1)-S*cdf(-d1,0,1);
}

// Black Scholes delta
float Black_Scholes_delta(float S,float K,float t,float T,float v,float r,float d,string call_put_flag)
{
    T = T - t;
    float T_sqrt = sqrt(T);
    float d1 = (log(float(S)/K)+((r-d)+v*v/2.)*T)/(v*T_sqrt);
    return  call_put_flag=="call"? normalcdf(d1):normalcdf(d1)-1 ;
}

float MC3(float s0, float K, float Time,float r,float d,float sigma,int M,int N,string call_put_flag)
{
    float dt = Time/N;
    float nudt = dt*(r - d - sigma*sigma/2);
    float sigsdt = sigma*sqrt(dt);
    float erddt = exp((r-d)*dt);
    float t = 0;
    float beta = 0.0 - 1;
    float delta1, stn1 = 0;
    
    float sum_ct = 0;
    float sum_ct2 = 0;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);
    for (int i  = 0; i < M; ++i)
    {
        float st1 = s0;
        float cv1 = 0;
        for (int j = 0; j <= N -1; ++j)
        {
            t = j*dt;
            delta1 = Black_Scholes_delta(st1, K, t, Time, sigma, r, d, call_put_flag);
            stn1 = st1 * exp(nudt + sigsdt*distribution(generator));
            cv1 += delta1*(stn1 - st1*erddt);
            st1 = stn1;
        }
        float ct = call_put_flag =="call"? ((0 > st1-K ? 0:st1-K) + beta*cv1) : ( (0 > K -st1 ? 0:K -st1) + beta*cv1 );
        sum_ct += ct;
        sum_ct2 += ct*ct;
    }
    float option_value = exp(0.0 - r*Time) * sum_ct/M;
    float SD = sqrt( (sum_ct2 - sum_ct*sum_ct/M)*exp(0.0-2*r*Time)/(M-1) );
    std::cout << " option value is " << option_value << "\n\n";
    std::cout << " Standard deviation is " << SD << "\n";
    std::cout << " Standard error is " << SD/sqrt(M) << "\n";
    return option_value;
}


// M paths,N timesteps
float MC4(float s0, float K, float Time,float r,float d,float sigma,int M,int N,string call_put_flag)
{
    float dt = Time/N;
    float nudt = dt*(r - d - sigma*sigma/2);
    float sigsdt = sigma*sqrt(dt);
    float erddt = exp((r-d)*dt);
    float t = 0;
    float beta = -1;
    float delta1, delta2, stn1 = 0, stn2 = 0;
    
    float sum_ct = 0;
    float sum_ct2 = 0;
    
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);
    for (int i  = 0; i < M; ++i)
    {
        float st1 = s0;
        float cv1 = 0;
        float st2 = s0;
        float cv2 = 0;
        for (int j = 0; j < N; ++j)
        {
            t = j*dt;
            delta1 = Black_Scholes_delta(st1, K, t, Time, sigma, r, d, call_put_flag);
            delta2 = Black_Scholes_delta(st2, K, t, Time, sigma, r, d, call_put_flag);
            float e = distribution(generator);
            stn1 = st1 * exp(nudt + sigsdt*e);
            stn2 = st2 * exp(nudt + sigsdt*(e));
            cv1 += delta1*(stn1 - st1*erddt);
            cv2 += delta2*(stn2 - st2*erddt);
            st1 = stn1;
            st2 = stn2;
        }
        float ct =  0 ;
        ct = call_put_flag =="call"? 0.5*( (0 > st1-K ? 0:st1-K) + beta*cv1 + (0 > st2-K ? 0:st2-K) +beta*cv2 ) : 0.5*( (0 > K - st1 ? 0:K - st1) + beta*cv1 + (0 > K - st2 ? 0:K - st2) +beta*cv2 );
        sum_ct += ct;
        sum_ct2 += ct*ct;
    }
    float option_value = exp(0.0 - r*Time) * sum_ct/M;
    float SD = sqrt( (sum_ct2 - sum_ct*sum_ct/M)*exp(0.0-2*r*Time)/(M-1) );
    std::cout << " option value is " << option_value << "\n";
    std::cout << " Standard deviation is " << SD << "\n";
    std::cout << " Standard error is " << SD/sqrt(M) << "\n";
    return option_value;
}


template <class T>
T binomial_price(T s0, T K, T Time,T r,T d,T sigma,int N,T dx,string call_put_flag)
{
    T dt = Time/N;
    T nu = r  - sigma*sigma/2;
    T edx = exp(dx);
    T pu = 0.5 + 0.5 * (nu*dt/dx );
    T pd = 1 - pu;
    T disc = exp(-r*dt);
    T dxu = sqrt(sigma*sigma*dt + nu*dt*nu*dt);
    T dxd = -dxu;
    // initialise asset prices at maturity
    vector <T> s_t;
    vector <T> c; // option values at maturity
    s_t.push_back(s0 * exp(N*dxd));
    for (int i = 0; i <=N ; ++i)
    {
        c.push_back( call_put_flag=="call"? (0 > s_t[i]-K ? 0:s_t[i]-K)  :  (0 > K-s_t[i] ? 0:K-s_t[i]) );
        s_t.push_back(s_t.back() * exp(dxu - dxd ));
    }
    c.push_back( call_put_flag=="call"? (0 > s_t.back()-K ? 0:s_t.back()-K)  :  (0 > K-s_t.back() ? 0:K-s_t.back()) );
    // option values
    for (int j = N - 1; j >= 0; --j)
    {
        for (int i = 0; i <= j; ++i)
        {
            c[i] = disc*(pd*c[i] + pu*c[i+1]);
        }
    }
    return  c[0];
}



#endif /* Option_h */














