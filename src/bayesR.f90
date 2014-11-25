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
use cmd_parser
use routinez
use RDistributions

implicit none
integer :: i, j, k, kk, jj,snploc
character (len=8)  :: cdate
character (len=10) :: ctime, ci, ca
logical :: overflow

call date_and_time(date=cdate,time=ctime)
call parse

call get_size
call load_phenos_plink
call allocate_data
call load_snp_binary
call xcenter 
call parse_priors
call init_random_seed
if(mcmc) then
   open(unit=21,file=logfil,status='unknown',form='formatted')
   write(21,901) 'Program BayesR'
   write(21,908) 'Run started at',cdate(1:4),cdate(5:6),cdate(7:8),ctime(1:2),ctime(3:4),ctime(5:6)
   write(21,902) 'Prefix for input files',trim(inprefix)
   write(21,902) 'Prefix for output files',trim(outprefix)
   write(21,903) 'Phenotype column',trait_pos
   write(21,903) 'No. of loci',nloci
   write(21,903) 'No. of individuals',nind
   write(21,903) 'No. of training individuals',nt
   write(21,906) 'Prior Vara', vara, dfvara
   write(21,906) 'Prior Vare', vare, dfvare
   write(21,903) 'Model size',msize
   write(21,903) 'No. of cycles',numit
   write(21,903) 'Burnin ',burnin
   write(21,903) 'Thinning rate',thin
   write(21,903) 'No. of mixtures',ndist
   write(21,905) 'Variance of dist ', gpin
   write(21,905) 'Dirichlet prior', delta
   write(21,903) 'Seed ', seed1
   write(21,903) 'SNP output ', snpout
   call flush(21)

   nnind=dble(nt)
   if(snpout) then
      open(unit=14,file=locfil,status='unknown',action='write')
   endif

   open(unit=25,file=hypfil,status='unknown',form='formatted')
   write(25,'(2(A10,1x),2(A12,1x),A7)',advance='no') 'Replicate','Nsnp','Va','Ve',' '
   do i=1,ndist
      write(ci,'(I8)') i
      ci=adjustl(ci)
      ca="Nk"//trim(ci)
      write(25,'(A10,1x)',advance="no") ca
   end do
   !write(25,'(A7)',advance='no') ' '
   do i=1,ndist
      write(ci,'(I8)') i
      ci=adjustl(ci)
      ca="Vk"//trim(ci)
      write(25,'(A12)',advance="no") ca
   end do
   write(25,*)
   !Calculate vara from h2 or use apriori estimates
   ! df=-2 will produce 'flat' (improper priors)
   ! df < -2 sets vara = h2*vary
   if(dfvara < -2) then
      VCE=.false.
      yhat=sum(why, mask=trains==0)/nnind
      vary= sum((why-yhat)*(why-yhat),mask=trains==0)/(nnind-1.0d0)
      vara=vara*vary
   else
      VCE=.true.
      vara_ap=vara
      vare_ap=vare
      if(dfvara == -2) then
         vara_ap=0d0
      endif
      if(dfvare == -2) then
         vare_ap=0d0
      endif
   endif

   !initialize
   pstore=0d0
   gstore=0d0
   mu_vare_store=0
   snpstore=0d0
   indiststore=0d0
   snptracker=2
   xpx=0d0
   do i=1,nloci
      xpx(i)=dot_product(X(:,i),X(:,i))
   enddo

   !starting values
   mu=1.0d0
   yadj=0.0d0
   yhat=sum(why, mask=trains==0)/nnind
   vary= sum((why-yhat)*(why-yhat),mask=trains==0)/(nnind-1.0d0)
   gp=gpin*vara
   scale=0.0d0
   p(1)=0.5d0
   p(2:ndist)=1.0d0/gpin(2:ndist)
   p(2:ndist)=0.5*p(2:ndist)/sum(p(2:ndist))
   g=dsqrt(vara/(0.5*dble(nloci)))
   do k=1,nloci
      permvec(k)=k
   enddo
   call compute_residuals
   each_cycle : do rep=1,numit
      included=0
      if(.not. VCE) then
         vare=dot_product(yadj,yadj)/rand_chi_square(nnind+3.0d0)
      endif
      yadj=yadj+mu
      mu=rand_normal(sum(yadj)/nnind, dsqrt(vare/nnind))
      yadj=yadj-mu
      log_p(1)=dlog(p(1))
      do i=2,ndist
         log_p(i)=dlog(p(i))
         log_gp(i)=dlog(gp(i))
         vare_gp(i)=vare/gp(i)
      enddo
      if(permute) then
         call permutate(permvec,nloci)
      endif
      do k=1,nloci
         snploc=permvec(k)
         z => X(:,snploc)
         zz=xpx(snploc)
         zz_vare=zz/vare
         gk=g(snploc)
         if(snptracker(snploc) > 1) then
            yadj=yadj+z*gk
         endif
         rhs= dot_product(yadj,z)
         lhs=zz/vare
         s(1)=log_p(1)
         do kk=2,ndist
            logdetV=dlog(gp(kk)*zz_vare+1.0d0)
            uhat=rhs/(zz+vare_gp(kk))
            s(kk)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(kk)
         enddo
         stemp=0.0d0
         do kk=1,ndist
            skk=s(kk)
            sk=0.0d0
            overflow=.false.
            do j=1,ndist
               if(j==kk) cycle
               clike=s(j)-skk
               if(clike .lt. -700) then !undeflow
                  cycle
               else if (clike .gt. 700) then 
                  overflow=.true.
                  exit
               endif
               sk=sk+dexp(clike)
            enddo
            if (overflow .eqv. .true.) then
               stemp(kk) = 0.0
            else
               stemp(kk)=1.0d0/(1.0d0+sk)
            endif
         enddo
         ssculm=0.0d0
         call random_number(r)
         indistflag=1
         do kk=1,ndist
            ssculm=ssculm+stemp(kk)
            if (r<ssculm) then
               indistflag=kk
               exit
            endif
         enddo
         snptracker(snploc)=indistflag
         if(indistflag==1) then
            gk=0.0d0
         else
            v1=zz+vare/gp(indistflag)
            gk=rand_normal(rhs/v1, dsqrt(vare/v1))
            yadj=yadj-z*gk  
            included=included+1
         endif
         g(snploc)=gk
         if(msize>0 .and. rep>mrep) then
            if(included>=msize) exit
         endif
      enddo ! each loci

      do i=1,ndist
         snpindist(i)=count(snptracker==i)
         varindist(i)=sum(g*g, mask= snptracker==i)
      enddo
      included=nloci-snpindist(1)

      if(VCE) then
         scale=(dble(included)*sum(g**2) + vara_ap*dfvara)/(dfvara+dble(included))
         vara=rand_scaled_inverse_chi_square(dble(included)+dfvara,scale)
         gp=gpin*vara
         vare=(dot_product(yadj,yadj)+vare_ap*dfvare)/ (nnind+dfvare)
         vare=rand_scaled_inverse_chi_square(nnind+dfvare,vare)
      endif

      dirx=dble(snpindist)+delta
      p=rdirichlet(ndist,dirx)

      if(rep>burnin .and. mod(rep,thin)==0) then
         counter=counter+1
         gstore=gstore+g
         pstore=pstore+p
         mu_vare_store(1)=mu_vare_store(1)+mu
         mu_vare_store(2)=mu_vare_store(2)+included
         mu_vare_store(3)=mu_vare_store(3)+vara
         mu_vare_store(4)=mu_vare_store(4)+vare
         varstore=varstore+varindist
         snpstore=snpstore+snpindist
         do i=1,nloci
            jj=snptracker(i)
            indiststore(i,jj)=indiststore(i,jj)+1
         enddo
         write(25,'(i10,1x,i10,1x,2(E15.7,1x),20(i10,1x))',advance='no')  rep, included , & 
              vara, vare, snpindist
         write(25,'(20E15.7,1x)') varindist
         call flush(25)
         if(snpout) call output_snploc
      end if
      ! re-calibrate residuals
      if(mod(rep,1000)==0) then
         !call compute_residuals
      endif
   enddo each_cycle

   !posterior means
   gstore=gstore/counter
   pstore=pstore/counter
   mu_vare_store=mu_vare_store/counter
   varstore=varstore/counter
   snpstore=snpstore/counter
   do i=1,nloci
      indiststore(i,:)=indiststore(i,:)/counter 
   enddo
   call output_model
   mu=mu_vare_store(1)
   call compute_dgv
   call write_dgv
else   ! end mcmc 
   open(unit=21,file=logfil,status='unknown',form='formatted')
   write(21,901) 'Program BayesR'
   write(21,908) 'Run started at',cdate(1:4),cdate(5:6),cdate(7:8),ctime(1:2),ctime(3:4),ctime(5:6)
   write(21,902) 'Prefix for input files',trim(inprefix)
   write(21,902) 'Prefix for output files',trim(outprefix)
   write(21,903) 'Phenotype column',trait_pos
   write(21,903) 'No. of loci',nloci
   write(21,903) 'No. of individuals',nind
   write(21,903) 'No. of individuals to predict',nt
   call load_param
   call compute_dgv
   call write_dgv
end if

call date_and_time(date=cdate,time=ctime)
write(21,908) 'Run ended at',cdate(1:4),cdate(5:6),cdate(7:8),ctime(1:2),ctime(3:4),ctime(5:6)
close(21)

901 format(a)
902 format(a,t30,': ',a)
903 format(a,t30,'= ',i8)
904 format(a,t30,'= ',f20.6)
905 format(a,t30,'= ',10f10.5)
906 format(a,t30,'= ',2f10.6)
907 format(a,t30,'= ',f10.2,a)
908 format(a20,1x,a4,'-',a2,'-',a2,' ',a2,':',a2':',a2)
end program bayesR
