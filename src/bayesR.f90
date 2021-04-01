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
!    along with this program.  If not, see <http://www.gnu.org/licenses/>

program bayesR
  use parz
  use bays
  use cmd_parser
  use routinez
  use RDistributions
#ifdef block
  use blockutilz
  use omp_lib
#endif  


  implicit none
  integer :: i, k, jj, counter, nm
  double precision :: ii, wtime
  character (len=8)  :: cdate
  character (len=10) :: ctime

#ifdef block
  wtime = omp_get_wtime()
#endif
  
  call date_and_time(date=cdate,time=ctime)
  call parse
  call get_size
  call load_phenos
  call allocate_data
  call parse_covar
  
  if(mcmc) then
     call init_model
     call write_log('train',cdate,ctime)
     call init_random_seed
     gstore=0d0
     indiststore=0d0
     snptracker(:,3)=2
     !starting values
     mu=sum(why, mask=trains==0)/dble(nt)
     g=dsqrt(sum(varcomp%vara)/(0.5*dble(nloci)))
     if(snpout) then
        open(unit=14,file=locfil,status='unknown',action='write')
     endif
#ifdef block
     call omp_set_num_threads(nthreads)
     call omp_set_dynamic(.false.) 
     call init_block
     allocate(permvec(n_blocks))
     do i=1,n_blocks
        permvec(i)=i
     enddo
#else
     allocate(permvec(nloci))
     do i=1,nloci
        permvec(i)=i
     enddo
     call load_snp
     call xcenter
     call compute_residuals
#endif
     vare=dot_product(yadj,yadj)/dble(nt)*0.5d0
     counter=0
     each_cycle : do rep=1,numit
        if(.not. VCE) then
           vare=dot_product(yadj,yadj)/rand_chi_square(dble(nt)+3.0d0)
        endif
        yadj=yadj+mu
        mu=rand_normal(sum(yadj)/dble(nt), dsqrt(vare/dble(nt)))
        yadj=yadj-mu
        !Fixed effects
        if(covar) then
           call update_fixed
        endif
        !SNP
#ifdef block
        call update_bayesR_block
#else        
        call update_bayesR
#endif
        call update_varcomp
        call update_segments

        if(rep>burnin .and. mod(rep,thin)==0) then
           counter=counter+1
           gstore=gstore+g
           mu_s=mu_s+mu
           vare_s=vare_s+vare
           
           if(covar) then
              ii=dble(counter-1)/dble(counter)
              alphastore=alphastore*ii+alpha/dble(counter)
           end if
           do i=1,nseg
              nm=segments%nmix(i)
              segments_s%p(i,1:nm)=segments_s%p(i,1:nm)+segments%p(i,1:nm)
              segments_s%snpinseg(i,1:nm)=segments_s%snpinseg(i,1:nm)+segments%snpinseg(i,1:nm)
              segments_s%varinseg(i,1:nm)=segments_s%varinseg(i,1:nm)+segments%varinseg(i,1:nm)
           enddo
           do i=1,ncomp
              varcomp_s%vara(i)=varcomp_s%vara(i)+varcomp%vara(i)
           enddo
           do i=1,nloci
              jj=snptracker(i,3)
              indiststore(i,jj)=indiststore(i,jj)+1
           enddo
           call output_hyp
           if(snpout) call output_snploc
        end if
     enddo each_cycle
   
     !posterior means
     gstore=gstore/dble(counter)
     mu_s=mu_s/dble(counter)
     vare_s=vare_s/dble(counter)
     segments_s%p=segments_s%p/dble(counter)
     segments_s%snpinseg=segments_s%snpinseg/dble(counter)
     segments_s%varinseg=segments_s%varinseg/dble(counter)
     varcomp_s%vara= varcomp_s%vara/dble(counter)
     indiststore=indiststore/dble(counter) 

     call output_model
     mu=mu_s
     call compute_dgv
     call write_predval(mbvfil,pred)
     if(covar) then
        call compute_fitted
        call write_predval(fittfil,pred)
     endif
  else   ! prediction
#ifdef block
     call init_block
     allocate(permvec(n_blocks))
     do i=1,n_blocks
      permvec(i)=i                                                                            
     enddo                                                                                       
#else        
     allocate(permvec(nloci))
     do i=1,nloci                                                                                
        permvec(i)=i                                                                             
     enddo
     call load_snp                                                                               
     call xcenter                                                                                
     call compute_residuals                                                                      
#endif                                                                                           
     call write_log('test',cdate,ctime)
     call load_param
     call compute_dgv
     call write_predval(mbvfil,pred)
     if(covar) then
        call compute_fitted
        call write_predval(fittfil,pred)
     end if
  end if

call date_and_time(date=cdate,time=ctime)
#ifdef block
 wtime= omp_get_wtime() -wtime
 write(21,'(a,t32,f20.6)')  'user time = ', wtime/60.0d0
#endif
call write_log('end',cdate,ctime)
end program bayesR
