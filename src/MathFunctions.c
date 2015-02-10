/*************************************************************************
    > File Name: MathFunctions.c
    > Author: Ge Tan
    > Mail: gtan@me.com 
    > Created Time: Mon Aug  4 19:54:31 2014
 ************************************************************************/

#include<stdio.h>
#include "MathFunctions.h"
#include "Rdefines.h"
/*********************************************************************
 * Compute log(1+x), (the natural logarithm) for very small x.
 * The equivalent darwin code by Gaston Gonnet:
 * ln1x := proc( x:numeric  )
 *  local i;
 *  if |x| > 0.125 then ln(1+x)
 *  elif |x| > 0.01 then
 *    x - x^2/2 - sum( (-x)^i/i, i=3..18  )
 *  else
 *    x - (0.5-(1/3+(-1/4+(1/5+(-1/6+(1/7-x/8)*x)*x)*x)*x)*x)*x^2
 *  fi
 *  end:
 * ******************************************************************/
double log1x(double x){
  double absX = fabs(x);
  if(absX > 0.125){
    return log(1+x);
  }else if(absX > 0.01){
    double sum = 0;
    double p = -1;
    for(int i=1; i <= 18; i++){
      p *= -x;
      sum += p/i;
    }
    return sum;
  }else{
    return x - (0.5-(1/3+(-1/4+(1/5+(-1/6+(1/7-x/8)*x)*x)*x)*x)*x)*x*x;
  }
}

/*********************************************************************
 * computes exp(x)-1, for very small x
 * expx1 := proc( x:numeric  )
 *   if |x| > 0.6931 then exp(x)-1
 * elif |x| > 0.0433217 then
 *   t := procname(x/16);
 *   to 4 do t := 2*t+t^2 od;
 *   t;
 * else (1+(1+(1+(1+(1+(1+(1+x/8)*x/7)*x/6)*x/5)*x/4)*x/3)*x/2)*x fi
 * end:
 * ******************************************************************/
double exp1x(double x){
  double absX = fabs(x);
  if(absX < 0.0433217){
    return (1+(1+(1+(1+(1+(1+(1+x/8)*x/7)*x/6)*x/5)*x/4)*x/3)*x/2)*x;
  }else{
    return exp(x)-1;
  }
}

