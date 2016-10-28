//
//  spy_day2.cpp
//  triniomial_tree
//
//  Created by jiale on 3/14/16.
//  Copyright © 2016 jiale. All rights reserved.
//

#include "spy_day2.hpp"
spy_day2::spy_day2()
{}
spy_day2::~spy_day2()
{}
float spy_day2::getK()
{
    return  K;
}
float spy_day2::getT()
{
    return T;
}
std::string spy_day2::getCall_put_type()
{
    return call_put_type;
}
float spy_day2::getLastPrice()
{
    return LastPrice;
}
float spy_day2::getBidPrice()
{
    return BidPrice;
}
float spy_day2::getAskPrice()
{
    return AskPrice;
}
float spy_day2::getQuote_Time()
{
    return Quote_Time;
}
void spy_day2::setSigma(float mySigma)
{
    sigma = mySigma;
}
float spy_day2::getSigma()
{
    return  sigma;
}

void spy_day2::Save_Line(std::string aLine)
{
    char *pc_line = (char *)aLine.c_str();
    K = std::stof(strtok(pc_line, ","));
    T = std::stof(strtok(NULL, ","));
    call_put_type = strtok(NULL, ",");
    LastPrice = std::stof(strtok(NULL, ","));
    BidPrice = std::stof(strtok(NULL, ","));
    AskPrice = std::stof(strtok(NULL, ","));
    std::string iv = strtok(NULL, ",");
    Quote_Time = std::stof(strtok(NULL, ","));
}

