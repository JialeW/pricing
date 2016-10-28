//
//  main.cpp
//  triniomial_tree
//
//  Created by jiale on 3/9/16.
//  Copyright © 2016 jiale. All rights reserved.
//

#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <time.h>
#include "Option.h"
#include "spy_sigmas.h" // sigma calculated form data1
#include "spy_day2.hpp" // data2
//#include "matrix.h"

using namespace  std;

vector<spy_day2> getData(); // for problem 10, match data1's sigma with data2. save 10 strike values each of SPY.
void problem10();

int main(int argc, const char * argv[]) {
    /*
    // trinomial tree   //first calculate a vector Sn[2*n +1];
    float s0 = 100.0;
    int N = 3;
    float K = 100;
    float Time = 1;
    float r = 0.06;
    float d = 0.0;
    float sigma = 0.2;
    float dx = 0.02;
    

    
    cout << "binomial_price\n";
    cout << binomial_price(s0,  K,  Time,r,d, sigma,1000,dx,"call") << "\n";
   
    
    cout << "==============================================  \n";
    cout << "Black_Scholes\n";
    cout << Black_Scholes(s0,  K, 0, Time,sigma,r, d,"call") <<"\n";
    cout << "==============================================  \n";
    cout << "Explicit\n";
    
    cout << Explicit(s0, K,Time,r,d,sigma, 200,30,dx,"call") <<endl;
   
    cout << "==============================================  \n";
    cout << " Implicit\n";
    cout << Implicit( s0,  K,  Time, r,d,sigma, N, 3, dx, "call") << "\n";
    
    cout << "==============================================  \n";
    cout << "CrankNicolson\n";
    cout << CrankNicolson( s0,  K,  Time, r,d,sigma, N, 3, dx, "call") << "\n";
    cout << "==============================================  \n";
    cout << "trinomial_price\n";
    cout << trinomial_price(s0, K,Time,r,d,sigma, N,dx,"call") << endl;
    cout << "==============================================  \n";
     cout << "MC \n";
    MC(s0, K,Time,r,d,sigma, 10, 10,"call");
    cout << "==============================================  \n";
     cout << "MC2 \n";
    MC2(s0, K,Time,r,d,sigma, 10, 10,"call");
    cout << "==============================================  \n";
     cout << "MC3 \n";
    MC3(s0, K,Time,r,d,sigma, 10, 10,"call");
    cout << "==============================================  \n";
     cout << "MC4 \n";
    MC4(s0, K,Time,r,d,sigma, 100,10,"call");
    
  
    
    // problem 7  Montr Carlo
    clock_t mc_start = clock();
    MC(s0, K,Time,r,d,sigma, 1000000, 300,"call");
    cout << "Time taken: " << (double)(clock() - mc_start)/CLOCKS_PER_SEC << "\n";
     
  
    */
    
   problem10();
    
    return 0;
}


vector<spy_day2> getData()
{
    vector<spy_sigmas> data1;
    vector<spy_day2> data2;
    
    // read from table1 and save data to data1
    fstream myfile1("/Users/jiale/Desktop/mytable_621_1.csv");
    if (!myfile1.is_open())      // check if the file exists
    {
        cout << "Can not find data1 !" << endl;
        system("Pause");// return with error code
    }
    string oneLine;
    getline(myfile1, oneLine);//jump tittle,read data from the second line
    // save 10 data from data1.csv
    data1.resize(10);
    for (int i = 0; i < 10; ++ i) {
        getline(myfile1, oneLine);
        data1[i].Save_Line(oneLine);
    }
    
    // read from table2 and save data to data2
    fstream myfile2("/Users/jiale/Desktop/mytable_621_2.csv");
    if (!myfile2.is_open())      // check if the file exists
    {
        cout << "Can not find data1 !" << endl;
        system("Pause");// return with error code
    }
    getline(myfile2, oneLine);//jump tittle,read data from the second line
    // save 10 data from data1.csv
    for (int i = 0; i < 1000; ++i)
    {
        if (myfile2.eof())
        {
            break;
        }
        getline(myfile2, oneLine);
        spy_day2 tempSpyDay2;
        tempSpyDay2.Save_Line(oneLine);
        data2.push_back(tempSpyDay2);
    }
    
    vector<spy_day2> final_data; // use this data to calculate option prices
    
    // both data1 and data2 are sorted by K, so we can match by K
    vector<spy_day2>::iterator it2 = data2.begin();
    for (vector<spy_sigmas>::iterator it1 = data1.begin(); it1 != data1.end(); ++it1)
    {
        // when all  K T and call_put_type are the same, save data1.sigma to data2
        while ( !( it2->getK() == it1->getK() && it2->getT() == it1->getT() && it2->getCall_put_type() == it1->getCall_put_type()) )
        {
            it2 ++;
        }
        it2->setSigma(it1->getSigma());
        // save this data to final_data
        final_data.push_back(*it2);
        it1 ++;
        it2 ++;
    }
    return final_data;
}

void problem10()
{
    ofstream fout; // write to a file
    fout.open("/Users/jiale/Desktop/FE621_hw2_p10.csv");
    fout << "Bid, Ask, Binomial tree, Trinomial tree, Explicit FiniteDifference, Implicit Finite Difference, Crank-Nicolson, Monte-Carlo\n";
    vector<spy_day2> final_data = getData();
    //Binomial tree, Trinomial tree, Explicit Finite
    //Difference, Implicit Finite Difference, Crank-Nicolson, Monte-Carlo
    matrix<float> price(10,6); // 6 method and 10 data
    float s0 = 186.63;
    float r = 0.0038;
    float d = 0;
    float dx = 0.2;
    int N = 10;
    for (int i = 0; i < 10; ++ i)
    {
        float K =  final_data[i].getK();
        float Time = (final_data[i].getT() - final_data[i].getQuote_Time())/365;
        float sigma = final_data[i].getSigma();
        string call_put_flag = final_data[i].getCall_put_type();
        price.m[i][0] = binomial_price(s0,  K, Time,r,d, sigma,100,dx,call_put_flag);
        price.m[i][1] = trinomial_price(s0, K,Time,r,d,sigma, N,dx,call_put_flag);
        price.m[i][2] = Explicit(s0, K,Time,r,d,sigma, N,100,dx,call_put_flag);
        
        price.m[i][3] = Implicit(s0, K,Time,r,d,sigma, N,100,dx,call_put_flag);
        cout <<  "price.m[i][3]"  <<  price.m[i][3] << "\n";
        price.m[i][4] = CrankNicolson( s0,  K,  Time, r,d,sigma, 3, 3, dx, call_put_flag);
        price.m[i][5] = MC(s0, K,Time,r,d,sigma, N, 100,call_put_flag);
        fout << final_data[i].getBidPrice() << "," << final_data[i].getAskPrice() << "," << price.m[i][0] << "," << price.m[i][1] << ","
        << price.m[i][2] << ","
        << price.m[i][3] << ","
        << price.m[i][4] << ","
        << price.m[i][5] << "\n";

    }
    fout.close();
}
