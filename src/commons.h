//
//  commons.h
//  PH-BB
//
//  Created by Semih Atakan on 7/5/15.
//  Copyright (c) 2015 3D Lab. All rights reserved.
//

#ifndef commons_h
#define commons_h


#include "Header.h"



const std::string currentDateTime();

double get_cpu_time();
double get_wall_time();

std::istream& safeGetline(std::istream& is, std::string& t);

#endif
