# fOptimize

## Description
Fortran modules that provides simple Continuous Optimization methods. The main purpose is solve generic continuous functions
given by a Fortran function environment.

The routines has only academic purpose. Further applications will need improvements, of course.

## Features

* One-dimensional optimization methods:
  + Bisection
  + Golden Ratio
  + Newton
  + Secant
  
* Multidimensional unconstrained optimization methods:
  + Gradient (Cauchy)
  + Newton
  + Self-scaling Davidon--Fletcher--Powell (DFP)
  + Self-scaling Broyden--Fletcher--Goldfarb--Shanno (BFGS)
  + Self-scaling Broyden family (DFP+BFGS)
  + Box Method (Direct coordinate search)
  + Hooke & Jeeves

It is important to say that, except for the derivative-free methods, the remainder methods do exact linear search with Bisection method.
  
One can find futher information about the above methods in Classical Optimization books like Luenberger, Nocedal, Bazaraa
and so on. I based my implementations in Luenberger and Bazaraa books. Also, I recommend Introduction To Derivative-Free Optimization (Andrew Conn et al.). 

## Contact
Diego T. Volpatto: volpatto@lncc.br
