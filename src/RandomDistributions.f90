
! Bayesian hierarchical models for complex trait analysis using a mixture of 
! normal distributions of SNP effects 
! Copyright (C) 2014 Gerhard Moser
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License!
!    along with this program.  If not, see <http://www.gnu.org/licenses

module RDistributions
! Code downloaded from Sukhbinder Singh blog 
! http://sukhbinder.wordpress.com/fortran-random-number-generation/
! Non uniform random Number Generators in Fortran


implicit none

double precision, parameter :: PI=3.141592653589793238462

contains

  function rand_uniform(a,b) result(c)
    double precision :: a,b,c,temp
    call random_number(temp)
    c= a+temp*(b-a)
  end function rand_uniform
!
! Random Sample from normal (Gaussian) distribution
!
  function rand_normal(mean,stdev) result(c)
    double precision :: mean,stdev,c,r,theta,temp(2)
    if(stdev <= 0.0d0) then
       !Write(*,*) "Standard Deviation must be +ve"
       c=mean
    else
       call random_number(temp)
       r=(-2.0d0*log(temp(1)))**0.5d0
       theta = 2.0d0*PI*temp(2)
       c= mean+stdev*r*sin(theta)
    end if
  end function rand_normal
!
!  Random smaple from an exponential distribution
!
  function rand_exponential(mean) result(c)
    double precision :: mean,c,temp
    if (mean <= 0.0d0) then
       write(*,*) "mean must be positive"
    else
       call random_number(temp)
       c=-mean*log(temp)
    end if
  end function rand_exponential
!
!  Return a random sample from a gamma distribution
!
  recursive function rand_gamma(shape, scale) result(ans)
    double precision shape,scale,ans,u,w,d,c,x,xsq,g,v
    if (shape <= 0.0d0) then
       write(*,*) "Shape parameter must be positive"
    end if
    if (scale <= 0.0d0) then
       write(*,*) "Scale parameter must be positive"
    end if
!
!    ## Implementation based on "A Simple Method for Generating Gamma Variables"
!    ## by George Marsaglia and Wai Wan Tsang.
!    ## ACM Transactions on Mathematical Software
!    ## Vol 26, No 3, September 2000, pages 363-372.
!
    if (shape >= 1.0d0) then
       d = shape - 1.0d0/3.0d0
       c = 1.0d0/(9.0d0*d)**0.5d0
       do while (.true.)
          x = rand_normal(0.0d0, 1.0d0)
          v = 1.0d0 + c*x
          do while (v <= 0.0d0)
             x = rand_normal(0.0d0, 1.0d0)
             v = 1.0d0 + c*x
          end do
          v = v*v*v
          call random_number(u)
          xsq = x*x
          if ((u < 1.0d0 -.0331d0*xsq*xsq) .or. (dlog(u) < 0.5d0*xsq + d*(1.0d0 - v + dlog(v))) )then
             ans=scale*d*v
             return
          end if
       end do
    else
        g = rand_gamma(shape+1.0d0, 1.0d0)
        call random_number(w)
        ans=scale*g*(w)**(1.0d0/shape)
        return
     end if
   end function rand_gamma
!
! ## return a random sample from a chi square distribution
! ## with the specified degrees of freedom
!
   function rand_chi_square(dof) result(ans)
     double precision ans,dof
     ans=rand_gamma(0.5d0*dof, 2.0d0)
   end function rand_chi_square
!
! ## return a random sample from a scaled 
!    inverse chi square distribution with
!    df and scale parameter
!
   function rand_scaled_inverse_chi_square(dof,scale) result(ans)
     double precision ans,dof,scale
     ans=rand_inverse_gamma(dble(0.5)*dof, dble(0.5)*dof*scale)
   end function rand_scaled_inverse_chi_square

! ## return a random sample from an inverse gamma random variable
!
   function rand_inverse_gamma(shape, scale) result(ans)
     double precision shape,scale,ans
!    ## If X is gamma(shape, scale) then
!    ## 1/Y is inverse gamma(shape, 1/scale)
     ans= 1.0d0 / rand_gamma(shape, 1.0d0 / scale)
   end function rand_inverse_gamma
!
!## return a sample from a Weibull distribution
!
   function rand_weibull(shape, scale) result(ans)
     double precision shape,scale,temp,ans
     if (shape <= 0.0d0) then
        write(*,*) "Shape parameter must be positive"
     end if
     if (scale <= 0.0d0) then
        write(*,*) "Scale parameter must be positive"
     end if
     call random_number(temp)
     ans= scale * (-log(temp))**(1.0d0 / shape)
   end function rand_weibull
!
!## return a random sample from a Cauchy distribution
!
   function rand_cauchy(median, scale) result(ans)
     double precision ans,median,scale,p
     if (scale <= 0.0d0) then
        write(*,*) "Scale parameter must be positive"
     end if
     call random_number(p)
     ans = median + scale*tan(PI*(p - 0.5d0))
   end function rand_cauchy
!
!## return a random sample from a Student t distribution
!
   function rand_student_t(dof) result(ans)
     double precision ans,dof,y1,y2
     if (dof <= 0.d0) then
        write(*,*) "Degrees of freedom must be positive"
     end if
!
! ## See Seminumerical Algorithms by Knuth
      y1 = rand_normal(0.0d0, 1.0d0)
      y2 = rand_chi_square(dof)
      ans= y1 / (y2 / dof)**0.50d0
!
    end function rand_student_t
!
!## return a random sample from a Laplace distribution
!## The Laplace distribution is also known as the double exponential distribution.
!
    function rand_laplace(mean, scale)  result(ans)
      double precision ans,mean,scale,u
      if (scale <= 0.0d0) then
        write(*,*) "Scale parameter must be positive"
     end if
     call random_number(u)
     if (u < 0.5d0) then
        ans = mean + scale*log(2.0d0*u)
     else
        ans = mean - scale*log(2.0d0*(1.0d0-u))
     end if
   end function rand_laplace
!
! ## return a random sample from a log-normal distribution
!
   function rand_log_normal(mu, sigma) result(ans)
     double precision ans,mu,sigma
     ans= exp(rand_normal(mu, sigma))
   end function rand_log_normal
!
! ## return a random sample from a beta distribution
!
   function rand_beta(a, b) result(ans)
     double precision a,b,ans,u,v
     if ((a <= 0.0d0) .or. (b <= 0.0d0)) then
        write(*,*) "Beta parameters must be positive"
     end if
!    ## There are more efficient methods for generating beta samples.
!    ## However such methods are a little more efficient and much more complicated.
!    ## For an explanation of why the following method works, see
!    ## http://www.johndcook.com/distribution_chart.html#gamma_beta
     u = rand_gamma(a, 1.0d0)
     v = rand_gamma(b, 1.0d0)
     ans = u / (u + v)
   end function rand_beta

! Based on rdirichlet function in R {gtools}     
   function rdirichlet(n,irx) result(x)
     integer :: n, i
     double precision :: sx
     double precision, dimension(n) :: irx, x
     do i=1,n
        x(i)=rand_gamma(irx(i),1.0d0)
     enddo
     sx=sum(x)
     x=x/sx
   end function rdirichlet
   
 end module RDistributions





