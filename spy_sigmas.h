//
//  spy_sigmas.h
//  triniomial_tree
//
//  Created by jiale on 3/14/16.
//  Copyright © 2016 jiale. All rights reserved.
//

#ifndef spy_sigmas_h
#define spy_sigmas_h
#include <string>
class spy_sigmas
{
private:
    float K;
    float T;
    std::string call_put_type;
    std::string underlying;
    float sigma;
public:
    spy_sigmas();
    ~spy_sigmas();
    float getK();
    float getT();
    std::string getCall_put_type();
    std::string getUnderlyinhg();
    float getSigma();
    void Save_Line(std::string aLine);
};


#endif /* spy_sigmas_h */
