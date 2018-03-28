/*
** mpr and mapr: SciVis-2018 Project 2
** Copyright (C)  2018 University of Chicago. All rights reserved.
**
** This file is distributed to students in the 2018 CMSC 23710/33710
** ("SciVis") class, and no one else, for use in the SciVis class. This
** file is not open-source licensed, nor licensed for any other distribution.
** You should not allow this file to be copied or downloaded by someone who
** is not a 2018 SciVis student, as that would violate the copyright
** held by University of Chicago.
*/

#include "mpr.h"
#include "mprPrivate.h"

/*
******** mprQuantize(): quantize a value by assigning it an index
** into an interval that has been subdivided into a number of segments
**
** INPUT: min, val, max, num: The closed interval [min,max] (you can assume
** min < max) is divided evenly into num subintervals (you can assuming that
** num > 0). The value val (which you can assume is finite) may or may not
** be in the interval.
**
** OUTPUT: Returns the index of the subinterval (from 0 to num-1) containing
** val if val is in [min,max], or 0 if val < min, or num-1 if val > max. In
** any case, the value ret returned must satisfy 0 <= ret <= num-1. The
** subintervals of [min,max] have particular requirements for which ends are
** open and closed.  Here is a diagram for num=5:
**
**    index:       0    1    2    3    4 = num-1
** subintervals: [   )[   )[   )[   )[    ]
**               |----|----|----|----|----|
**    value:    min                      max
*/
uint
mprQuantize(real min, real val, real max, uint num) {
  uint ret=42;
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
  if(val<min)
  	ret = 0;
  else if(val>= max)
  	ret = (signed)num - 1;
  else{
  real interval = (max - min)/(signed)num;
  ret = floor((val - min)/interval);
  }
  // ( 9 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return ret;
}
