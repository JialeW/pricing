//
//  spy_day2.hpp
//  triniomial_tree
//
//  Created by jiale on 3/14/16.
//  Copyright © 2016 jiale. All rights reserved.
//

#ifndef spy_day2_hpp
#define spy_day2_hpp
#include <string>
class spy_day2
{
private:
    float K;
    float T;
    std::string call_put_type;
    float LastPrice;
    float BidPrice;
    float AskPrice;
    float Quote_Time;
    float sigma;
public:
    spy_day2();
    ~spy_day2();
    float getK();
    float getT();
    std::string getCall_put_type();
    float getLastPrice();
    float getBidPrice();
    float getAskPrice();
    float getQuote_Time();
    void setSigma(float mySigma);
    float getSigma();
    void Save_Line(std::string aLine);
};



#endif /* spy_day2_hpp */
