//
//  spy_sigmas.cpp
//  triniomial_tree
//
//  Created by jiale on 3/14/16.
//  Copyright © 2016 jiale. All rights reserved.
//

#include "spy_sigmas.h"

spy_sigmas::spy_sigmas()
{}
spy_sigmas::~spy_sigmas()
{}
float spy_sigmas::getK()
{
    return  K;
}
float spy_sigmas::getT()
{
    return T;
}
std::string spy_sigmas::getCall_put_type()
{
    return call_put_type;
}
std::string spy_sigmas::getUnderlyinhg()
{
    return underlying;
}
float spy_sigmas::getSigma()
{
    return sigma;
}

void spy_sigmas::Save_Line(std::string aLine)
{
    char *pc_line = (char *)aLine.c_str();
    K = std::stof(strtok(pc_line, ","));
    T = std::stof(strtok(NULL, ","));
    call_put_type = strtok(NULL, ",");
    underlying = strtok(NULL, ",");
    sigma = std::stof(strtok(NULL, ","));
}

