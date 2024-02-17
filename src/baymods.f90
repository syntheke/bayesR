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

 module parz
  implicit none
  !Global variables
  integer :: nloci, nind, nt, ndist, numit, seed1, ios, rep, &
             burnin, thin, trait_pos, nfixed
  logical :: mcmc, snpout, permute, covar, plink
  !File names
  character(len=200) :: genfil, phenfil, bimfil, inprefix, outprefix, logfil, freqfil, &
       mbvfil, hypfil, locfil, modfil, paramfil, alphafil, covarfil, fittfil, &
       varcompfil,segmentfil, snpmodfil
  ! Data
  double precision, dimension(:), allocatable :: why, pred, freqstore
  double precision, target, dimension(:,:), allocatable :: W
  ! Genetic variables
  integer, dimension(:,:), allocatable ::snpindist, snptracker
  double precision :: vare, vare_ap, varb, dfvare, mu, scale
  double precision, dimension(:), allocatable :: g, yadj, dirx, alpha
  integer,  dimension(:), allocatable :: alphanull
  ! BayesR specific
  logical :: VCE
  integer:: included, msize, mrep
  double precision, dimension(:), allocatable :: xpx, wpw
  integer, allocatable, dimension(:)::permvec
  !Auxiliary variables
  integer, dimension(:), allocatable ::trains
  double precision :: vary
  ! Storage
  double precision :: mu_s, vare_s
  double precision, dimension(:), allocatable :: gstore, alphastore
  double precision, dimension(:,:), allocatable :: indiststore
  ! V and segment dimensions
  integer, parameter :: mxdist = 9
  integer :: ncomp, nseg
  ! block related
  integer :: block_size, nthreads 
  double precision, dimension(:,:), allocatable :: geno_mean
#ifdef block
  integer(kind=1), dimension(:,:), allocatable :: X
#else
  double precision, target, dimension(:,:), allocatable :: X
#endif

end module parz

module multiprior
  use parz
  implicit none

  type variance
     double precision, allocatable, dimension(:) :: vara
     double precision, allocatable, dimension(:) :: dfvara
     double precision, allocatable, dimension(:) :: vara_ap
  end type variance
  
  type variance_s
     double precision, allocatable, dimension(:) :: vara
  end type variance_s
  
  type mixture
     integer, allocatable, dimension(:) :: nmix 
     integer, allocatable, dimension (:) :: varcomp
     double precision, allocatable, dimension(:,:) :: gpin
     double precision, allocatable, dimension(:,:) :: delta
     double precision, allocatable, dimension (:,:) :: gp
     double precision, allocatable, dimension(:,:) :: p
     integer, allocatable, dimension (:,:) :: snpinseg
     double precision, allocatable, dimension (:,:) :: varinseg
  end type mixture
  
  type mixture_s
     double precision, allocatable, dimension(:,:) :: p
     double precision, allocatable, dimension (:,:) :: snpinseg
     double precision, allocatable, dimension (:,:) :: varinseg
  end type mixture_s

  type(variance) :: varcomp
  type(variance_s) :: varcomp_s
  type(mixture):: segments
  type(mixture_s):: segments_s

contains
  
 subroutine alloc_varcomp
    allocate(varcomp%vara(ncomp))
    allocate(varcomp%dfvara(ncomp))
    allocate(varcomp%vara_ap(ncomp))
    varcomp%vara=0.0d0
    varcomp%dfvara=0.0d0
    varcomp%vara_ap=0.0d0
  end subroutine alloc_varcomp

  subroutine alloc_varcomp_s
    allocate(varcomp_s%vara(ncomp))
    varcomp_s%vara=0.0d0
  end subroutine alloc_varcomp_s

  subroutine alloc_segments
    allocate(segments%nmix(nseg))
    allocate(segments%varcomp(nseg))
    allocate(segments%gpin(nseg,mxdist))
    allocate(segments%delta(nseg,mxdist))
    allocate(segments%gp(nseg,mxdist))
    allocate(segments%p(nseg,mxdist))
    allocate(segments%snpinseg(nseg,mxdist))
    allocate(segments%varinseg(nseg,mxdist))
    segments%nmix=0
    segments%varcomp=0.0d0
    segments%gpin=0.0d0
    segments%delta=0.0d0
    segments%gp=0.0d0
    segments%p=0.0d0
    segments%snpinseg=0
    segments%varinseg=0.0d0
  end subroutine alloc_segments
  
  subroutine alloc_segments_s
    allocate(segments_s%p(nseg,mxdist))
    allocate(segments_s%snpinseg(nseg,mxdist))
    allocate(segments_s%varinseg(nseg,mxdist))
    segments_s%p=0.0d0
    segments_s%snpinseg=0.0d0
    segments_s%varinseg=0.0d0
  end subroutine alloc_segments_s

end module multiprior

module cmd_parser
  use parz
  use multiprior
implicit None

integer :: narg, nopt
character(len=1024) :: arg
character(len=1024), allocatable, dimension(:) :: cmd_line

integer, parameter :: &
     a_int     = 1, &
     a_float   = 2, &
     a_char    = 3, &
     a_flag    = 4 

type param
   character(len=20), allocatable, dimension(:) :: key
   character(len=20), allocatable, dimension(:) :: argtype
   character(len=100), allocatable,dimension(:) :: desc
   integer, allocatable, dimension(:) :: kind
   character(len=200), allocatable,dimension(:) :: default
end type param

type(param) :: register

contains

  subroutine alloc_register
    allocate(register%key(nopt))
    allocate(register%argtype(nopt))
    allocate(register%desc(nopt))
    allocate(register%kind(nopt))
    allocate(register%default(nopt))
  end subroutine alloc_register

  subroutine init_register
#ifdef block
    nopt=30
#else    
    nopt=28
#endif    
    call alloc_register
    call include_option(1, '-bfile'  ,'[prefix]',  'prefix PLINK binary files'               ,a_char,     '')
    call include_option(2, '-file'   ,'[prefix]',  'prefix flat input files'                 ,a_char,     '')
    call include_option(3, '-out'    ,'[prefix]',  'prefix for output'                       ,a_char,     '')
    call include_option(4, '-n'      ,'[num]',     'phenotype column'                        ,a_int,      '1')
    call include_option(5, '-vara'   ,'[num]',     'SNP variance prior'                      ,a_float, '0.01')
    call include_option(6, '-vare'   ,'[num]',     'error variance prior'                    ,a_float, '0.01')
    call include_option(7, '-dfvara' ,'[num]',     'degrees of freedom Va'                   ,a_float, '-2.0')
    call include_option(8, '-dfvare' ,'[num]',     'degrees of freedom Ve'                   ,a_float, '-2.0')
    call include_option(9,'-delta'  ,'[num]',     'prior for Dirichlet'                     ,a_float,   '1,1,1,1')
    call include_option(10, '-msize'  ,'[num]',     'number of SNPs in reduced update'        ,a_int,      '0')
    call include_option(11, '-mrep'   ,'[num]',     'number of full cycles in reduced update' ,a_int,   '5000')
    call include_option(12,'-numit' ,'[num]',     'length of MCMC chain'                    ,a_int,   '50000')
    call include_option(13,'-burnin' ,'[num]',     'burnin steps'                            ,a_int,   '20000')
    call include_option(14,'-thin'   ,'[num]',     'thinning rate'                           ,a_int,   '10')
    call include_option(15,'-ndist'  ,'[num]',     'number of mixture distributions'         ,a_int,     '4')
    call include_option(16,'-gpin'   ,'[num]',     'effect sizes of mixtures (% x Va)'       ,a_float,  '0.0,0.0001,0.001,0.01')
    call include_option(17,'-seed'   ,'[num]',     'initial value for random number'         ,a_int,      '0')
    call include_option(18,'-predict','[flag]',    'perform prediction'                      ,a_flag,     'f')
    call include_option(19,'-snpout' ,'[flag]',    'output detailed SNP info'                ,a_flag,     'f')
    call include_option(20,'-permute','[flag]',     'permute order of SNP'                    ,a_flag,   'f')
    call include_option(21,'-model'  ,'[filename]','model summary file (for prediction) '    ,a_char,   '')
    call include_option(22,'-freq'   ,'[filename]','SNP frequency file (for prediction)'     ,a_char,   '')
    call include_option(23,'-param'  ,'[filename]','SNP effect file (for prediction)'        ,a_char,   '')
    call include_option(24,'-covar'  ,'[filename]','design matrix for fixed effects'         ,a_char,   '')
    call include_option(25,'-alpha'   ,'[filename]','Fixed effects estimates file (prediction)',a_char,   '')
    call include_option(26,'-snpmodel','[filename]','grouped effects SNP model'              ,a_char,   '')
    call include_option(27,'-varcomp','[filename]','vara priors file (when vara >1)'          ,a_char,   '')
    call include_option(28,'-segments','[filename]','segment priors file (when nseg >1)'      ,a_char,   '')
#ifdef block
    call include_option(29,'-blocksize' ,'[num]'  ,'Number of SNP in rhs block'              ,a_int,    '4')
    call include_option(30,'-nthreads' ,'[num]'   ,'Number of threads'                       ,a_int,    '4')
#endif
  end subroutine init_register

 subroutine include_option(pos,key,argtype,desc,kind,default)
   implicit none
   integer :: pos
   character(len=*) :: key
   character(len=*) :: desc
   character(len=*) :: argtype
   integer :: kind
   character(len=*) :: default
   
   register%key(pos)=trim(key)
   register%argtype(pos)=trim(argtype)
   register%desc(pos)=trim(desc)
   register%kind(pos)=kind
   register%default(pos)=trim(default)
 end subroutine include_option

 logical function str_match(str1,str2)
   implicit none
   character(*) :: str1, str2
   integer :: comp
   !exact match
   comp = index(trim(str1),trim(str2))*index(trim(str2),trim(str1))
   str_match = .false.
   if (comp /= 0) str_match = .true.
 end function str_match

  integer function ntokens(line,separator)
    !http://www.tek-tips.com/viewthread.cfm?qid=1688013
    character,intent(in):: line*(*)
    integer i, n, toks
    character(len=1), dimension(:) :: separator
    i = 1;
    n = len_trim(line)
    toks = 0
    ntokens = 0
    do while(i <= n) 
        do while(any(separator==line(i:i)))
          i = i + 1
          if (n < i) return
       enddo
       toks = toks + 1
       ntokens = toks
       do
          i = i + 1
          if (n < i) return
          if (any(separator==line(i:i))) exit
       enddo
    enddo
  end function ntokens

  subroutine tokenize(str,sep,dim,val)
    character(len=*), intent(in) :: str
    character(len=1), intent(in) :: sep
    integer, intent(in) :: dim
    character(len=128),dimension(dim), intent(out) :: val
    integer :: pos1, pos2, n

    pos1=1
    n=0
    do
       pos2 = index(str(pos1:), sep)
       if (pos2 == 0) then
          n=n+1
          val(n) = str(pos1:)
          exit
       endif
       n=n+1
       val(n) = str(pos1:pos1+pos2-2)
       pos1 = pos2+pos1
    enddo
  end subroutine tokenize

  subroutine get_cmdLine
    integer ::i 

    narg = command_argument_count()  
    allocate(cmd_line(narg))
    do i=1,narg
       call get_command_argument(i,arg)
       cmd_line(i)=arg
    end do
  end subroutine get_cmdLine

  integer function cast_int(value)
    character(len=*) :: value
    read(value,*)cast_int
  end function cast_int
  
  real function cast_float(value)
    character(len=*) :: value
    read(value,*)cast_float
  end function cast_float

  logical function cast_logical(value)
    character(len=*) :: value
    cast_logical=.false.
    if(str_match(trim(value),'t')) cast_logical=.true.
  end function cast_logical

  logical function is_key(key)
    character(len=*) :: key
    character(len=2),parameter :: str='-'
    is_key=.false.
    !if(str_match(trim(key(1:1)),str)) is_key=.true.
    if(str_match(trim(key(1:1)),str) .and..not. is_numeric(trim(key(2:2)))) is_key=.true.
  end function is_key

  logical function is_numeric(key)
    character(len=*) :: key
    integer :: value,ios
    is_numeric=.true.
    read(key,*,iostat=ios) value
    if(ios /=0)is_numeric=.false.
  end function is_numeric

  subroutine update_register
    integer ::  i, k, kind
    character(len=1024) :: key
    character(len=100) :: err1, err2
    logical :: valid

    err1='Unknown command line option : '
    err2='Missing argument for :'

    i=1
    do while (i.le.narg)
       key=trim(cmd_line(i))
       k=1
       valid=.false.
       do while(k.le.nopt) 
          if(str_match(key,trim(register%key(k)))) then
             kind=register%kind(k)
             valid=.true.
             if(kind==4) then
                register%default(k)='t'
                i=i+1
             else if(i==narg) then
                print *, trim(err2),trim(key),i
                stop 'ERROR: Problem parsing the command line arguments'
                valid=.false.
             else if(is_key(trim(cmd_line(i+1)))) then
                print *, trim(err2),trim(key),i
                stop 'ERROR: Problem parsing the command line arguments'
             else 
                register%default(k)=trim(cmd_line(i+1))
                valid=.true.
                i=i+2
             endif
             exit
          end if
          k=k+1
          valid=.false.
       enddo
       if(.not.valid) then
          print *, trim(err1),trim(key),i
          stop 'ERROR: Problem parsing the command line arguments'
       end if
    enddo
  end subroutine update_register

  subroutine parse_help
   integer :: i,k
   character(len=1024) :: key
   if(narg==0) then
      write(*,'(a)')adjustl('argument  type        description &
           &                                default')
      do k=1,nopt
         write(*,111) adjustl(register%key(k)),adjustl(register%argtype(k)), &
              adjustl(register%desc(k)), adjustl(register%default(k))
      end do
      stop
   else
      do i=1,narg
         key=trim(cmd_line(i))
         if(str_match(trim(key(1:2)),'-h')) then
            write(*,'(a)')adjustl('argument  type        description &
                 &                                default')
            do k=1,nopt
               write(*,111) adjustl(register%key(k)),adjustl(register%argtype(k)), &
                    adjustl(register%desc(k)), adjustl(register%default(k))
            end do
            stop
         end if
111      format(a10,a12,a45,a25)
      end do
   end if
 end subroutine parse_help

subroutine parse_out
  integer::i
  do i=1,nopt
    if(str_match(trim(register%key(i)),'-out')) then
        outprefix = trim(register%default(i))
        logfil = trim(outprefix)//'.log'
        freqfil= trim(outprefix)//'.frq'
        mbvfil= trim(outprefix)//'.gv'
        hypfil= trim(outprefix)//'.hyp'
        locfil= trim(outprefix)//'.snp'
        modfil=trim(outprefix)//'.model'
        paramfil=trim(outprefix)//'.param'
        fittfil=trim(outprefix)//'.fitted'
        alphafil=trim(outprefix)//'.alpha'
        exit
     end if
  end do
end subroutine parse_out

subroutine parse_plink
  integer::i
  logical:: fileExist

  plink=.false.
  do i=1,nopt
     if(str_match(trim(register%key(i)),'-bfile')) then
        inprefix=trim(register%default(i))
        if(.not.inprefix=="") then
           plink=.true.
           genfil=trim(inprefix)//'.bed'
           inquire(file=genfil,exist=fileExist)
           if(.not.fileExist)then
              write(*,*)'file ',trim(genfil),' not found'
              stop
           end if
           phenfil=trim(inprefix)//'.fam'
           inquire(file=phenfil,exist=fileExist)
           if(.not.fileExist)then
              write(*,*)'file ',trim(phenfil),' not found'
              stop
           end if
           bimfil=trim(inprefix)//'.bim'
           inquire(file=bimfil,exist=fileExist)
           if(.not.fileExist)then
              write(*,*)'file ',trim(bimfil),' not found'
              stop
           end if
           exit
        end if
     end if
  enddo
end subroutine parse_plink

subroutine init_model
  call parse_snpmodel
  call parse_varcomp
  call parse_segments
  call set_varcomp
  call set_segments
  call init_mixture
  call hyp_header
end subroutine init_model
            
 subroutine parse_snpmodel
   integer :: i, j
   character(len=200) :: fname
   logical:: fileExist

   do i=1,nopt
     if(str_match(trim(register%key(i)),'-snpmodel')) then
        fname=trim(register%default(i))
        if(.not.fname=="") then
           snpmodfil=trim(fname)
           inquire(file=snpmodfil,exist=fileExist)
           if(.not.fileExist)then
               write(*,*)'file ',trim(snpmodfil),' not found'
               stop
            end if
            open(37,file=trim(snpmodfil),status='old',form='formatted')
            do j=1,nloci
               read(37,*,iostat=ios) snptracker(j,1:2)
            enddo
            close(37,status='keep')
            nseg=maxval(snptracker(:,1))
            ncomp=maxval(snptracker(:,2))
            ! some basic input checking
            do j=1,nseg
               if(count(snptracker(:,1)==j)==0) then
                  write(*,*) 'Number of segments in snpmodel < ', nseg
                  stop
               endif
            enddo
            do j=1,ncomp
               if(count(snptracker(:,2)==j)==0) then
                  write(*,*) 'Number of varcomps in snpmodel < ', ncomp
                  stop
               endif
            enddo
         end if
      endif
   enddo
 end subroutine parse_snpmodel

subroutine parse_varcomp
  integer::i,j,nitem
  logical:: fileExist
  character(len=200) :: fname

  do i=1,nopt
     if(str_match(trim(register%key(i)),'-varcomp')) then
        fname=trim(register%default(i))
        if(.not.fname=="") then
           varcompfil=trim(fname)
           inquire(file=varcompfil,exist=fileExist)
           if(.not.fileExist)then
               write(*,*)'file ',trim(varcompfil),' not found'
               stop
           end if
           open(33,file=trim(varcompfil),status='old',form='formatted')
           ! variance component
           !ncomp=  no of vara components in snpmodel (parse_snpmodel)
           call alloc_varcomp
           call alloc_varcomp_s
           nitem=0
           open(33,file=trim(varcompfil),status='old',form='formatted')
           do
              read(33,*,iostat=ios)
              if (ios.ne.0) exit
              nitem=nitem+1
           enddo
           rewind(33)
           if(nitem /= ncomp) then
              write(*,'(a)') 'Mismatch in number of variance components'
              write(*,'(2(a,i4,2x))')trim(snpmodfil),ncomp,trim(varcompfil),nitem
              stop
           endif
           do j=1,ncomp
              read(33,*,iostat=ios) varcomp%vara(j), varcomp%dfvara(j)
              nitem=nitem+1
           enddo
           close(33)
        else
           ncomp=1
           call alloc_varcomp
           call alloc_varcomp_s
           do j=1,nopt
              if(str_match(trim(register%key(j)),'-vara')) then
                 varcomp%vara(1) = cast_float(register%default(j))
              else if(str_match(trim(register%key(j)),'-dfvara')) then
                 varcomp%dfvara(1) = cast_float(register%default(j))
              endif
           end do
        end if
     else if(str_match(trim(register%key(i)),'-vare')) then
        vare = cast_float(register%default(i))
     else if(str_match(trim(register%key(i)),'-dfvare')) then
        dfvare = cast_float(register%default(i))
     end if
  end do
end subroutine parse_varcomp

subroutine parse_segments
  integer::i,j,k, nitem
  logical:: fileExist
  character (len=128) :: lstring
  character(len=200) :: fname
  character(len=128), allocatable, dimension(:) :: c_string

  do i=1,nopt
     if(str_match(trim(register%key(i)),'-segments')) then
        fname=trim(register%default(i))
        if(.not.fname=="") then
           segmentfil=trim(fname)
           inquire(file=segmentfil,exist=fileExist)
           if(.not.fileExist)then
              write(*,*)'file ',trim(segmentfil),' not found'
              stop
           end if
           open(33,file=trim(segmentfil),status='old',form='formatted')
           ! segement information
           !no of segments from snpmodel (parse_snpmodel)
           call alloc_segments
           call alloc_segments_s
           nitem=0
           open(33,file=trim(segmentfil),status='old',form='formatted')
           do
              read(33,*,iostat=ios) lstring
              if (ios.ne.0) exit
              nitem=nitem+1
           enddo
           rewind(33)
           if(nitem /= nseg) then
              write(*,'(a)') 'Mismatch in number of segments'
              write(*,'(2(a,i4,2x))')trim(snpmodfil),ncomp,trim(segmentfil),nitem
              stop
           endif
           do j=1,nseg
              read(33,*) lstring
              read(lstring,*,iostat=ios) k
              backspace(33)
              read(33,*) segments%nmix(j),segments%gpin(j,1:k), segments%delta(j,1:k)
           enddo
           close(33)
        else
           nseg=1
           call parse_ndist
           allocate(c_string(ndist))
           call alloc_segments
           call alloc_segments_s
           segments%nmix(1)=ndist
           do j=1,nopt
              if(str_match(trim(register%key(j)),'-gpin')) then
                 nitem=ntokens(register%default(j),(/','/))
                 if(nitem /= ndist) then
                    print *, 'Error: Number of mixture classes ',ndist
                    print *,'        but ',nitem,'effect sizes specified (gpin)'
                    stop
                 else
                    call tokenize(trim(register%default(j)),',',ndist,c_string)
                    do k=1,ndist
                       segments%gpin(1,k)=cast_float(trim(c_string(k)))
                    enddo
                 endif
              else if(str_match(trim(register%key(j)),'-delta')) then
                 nitem=ntokens(register%default(j),(/','/))
                 if(nitem /= ndist) then
                    print *, 'Error: Number of mixture classes',ndist
                    print *,'        but',nitem,'prior values for Dirichlet specified (delta)'
                    stop
                 else
                    call tokenize(trim(register%default(j)),',',ndist,c_string)
                    do k=1,ndist
                       segments%delta(1,k)=cast_float(trim(c_string(k)))
                    enddo
                 endif
              end if
           end do
        end if
     end if
  enddo
end subroutine parse_segments

subroutine set_varcomp
  integer :: i
  double precision :: yhat

  yhat=sum(why, mask=trains==0)/dble(nt)
  vary= sum((why-yhat)*(why-yhat),mask=trains==0)/(dble(nt)-1.0d0)
  !Calculate vara from h2 or use apriori estimates
  do i=1,ncomp
     if(varcomp%dfvara(i) < -2) then
        varcomp%vara(i)=varcomp%vara(i)*vary
     else if(varcomp%dfvara(i)== -2.0d0) then
        varcomp%vara_ap(i)=0.0d0
        VCE=.true.
     else
        varcomp%vara_ap(i)=varcomp%vara(i)
     endif
  end do
  if(dfvare == -2) then
     vare_ap=0.0d0
  else
     vare_ap=vare
  endif
end subroutine set_varcomp

subroutine set_segments()
  integer :: i, s
  
  if(nseg==1 .and. ncomp==1) then
     snptracker(:,1)=1
     snptracker(:,2)=1
     segments%varcomp(1)=1
  else
     do i=1,nloci
        s=snptracker(i,1)
        segments%varcomp(s)=snptracker(i,2)
     end do
  end if
end subroutine set_segments

subroutine init_mixture
  !Not optimal yet
  integer :: i,nm

  do i=1,nseg
     nm=segments%nmix(i)
     segments%gp(i,1:nm)=segments%gpin(i,1:nm)*(sum(varcomp%vara)/dble(nseg))
     segments%p(i,1)=0.5d0
     segments%p(i,2:nm)=1.0d0/segments%gpin(i,2:nm)
     segments%p(i,2:nm)=0.5*segments%p(i,2:nm)/sum(segments%p(i,2:nm))
  end do
end subroutine init_mixture

subroutine parse_flatinput
  integer::i
  logical:: fileExist

  if(plink.eqv..false.) then
     do i=1,nopt
        if(str_match(trim(register%key(i)),'-file')) then
           inprefix=trim(register%default(i))
           if(.not.inprefix=="") then
              genfil=trim(inprefix)//'.txt'
              inquire(file=genfil,exist=fileExist)
              if(.not.fileExist)then
                 write(*,*)'file ',trim(genfil),' not found'
                 stop
              end if
              phenfil=trim(inprefix)//'.phe'
              inquire(file=phenfil,exist=fileExist)
              if(.not.fileExist)then
                 write(*,*)'file ',trim(phenfil),' not found'
                 stop
              end if
              bimfil=trim(inprefix)//'.map'
              inquire(file=bimfil,exist=fileExist)
              if(.not.fileExist)then
                 write(*,*)'file ',trim(bimfil),' not found'
                 stop
              end if
              exit
           end if
        endif
     enddo
  endif
end subroutine parse_flatinput

subroutine parse_trait_pos
  integer::i

  do i=1,nopt
     if(str_match(trim(register%key(i)),'-n')) then
        trait_pos = cast_int(register%default(i))
     end if
  end do
end subroutine parse_trait_pos

subroutine parse_predict
  integer::i
  character(len=200) :: fname
  logical :: flag, fileExist
  mcmc=.true.
  do i=1,nopt
     if(str_match(trim(register%key(i)),'-predict')) then
        flag = cast_logical(register%default(i))
        if(flag) mcmc=.false. 
     end if
  end do
  if(.not.mcmc) then
     do i=1,nopt
        if(str_match(trim(register%key(i)),'-model')) then
           fname=trim(register%default(i))
           if(str_match(trim(fname),' ')) stop 'Error: No model file specified'
           modfil=trim(fname)
           inquire(file=modfil,exist=fileExist)
           if(.not.fileExist)then
              write(*,*)'file ',trim(modfil),' not found'
              stop
           end if
        else if(str_match(trim(register%key(i)),'-freq')) then
           fname=trim(register%default(i))
           if(str_match(trim(fname),' ')) stop 'Error: No SNP frequency file specified'
           freqfil=trim(fname)
           inquire(file=freqfil,exist=fileExist)
           if(.not.fileExist)then
              write(*,*)'file ',trim(freqfil),' not found'
              stop
           end if
        else if(str_match(trim(register%key(i)),'-param')) then
           fname=trim(register%default(i))
           if(str_match(trim(fname),' ')) stop 'Error: No param file specified'
           paramfil=trim(fname)
           inquire(file=paramfil,exist=fileExist)
           if(.not.fileExist)then
              write(*,*)'file ',trim(paramfil),' not found'
              stop
           end if
        endif
     end do
  end if
end subroutine parse_predict


subroutine parse_covar
  integer:: i,ii, k, good, count
  character(len=128), allocatable, dimension(:) ::tmp 
  logical :: fileExist
  character(len=1024):: str,fname 

  covar=.false.
  do i=1,nopt
     if(str_match(trim(register%key(i)),'-covar')) then
        fname=trim(register%default(i))
        if(.not.fname=="") then
           covar=.true.
           covarfil=trim(fname)
           inquire(file=covarfil,exist=fileExist)
           if(.not.fileExist)then
              write(*,*)'file ',trim(covarfil),' not found'
              stop
           end if
        end if
        exit
     end if
  enddo
  
  if(covar) then
     varb=dble(1e10)
     count=0
     open(55,file=trim(covarfil),status='old',form='formatted')
     read(55,'(a)') str
     backspace 55
     nfixed = ntokens(str,(/' '/))
     allocate(W(nt,nfixed),tmp(nfixed),wpw(nfixed),alpha(nfixed), &
              alphanull(nfixed),alphastore(nfixed))
     do i=1,nind
        read(55,'(a)') str
        if(trains(i)==0) then
           count=count+1
           good=index(str, "NA")
           if(good /= 0) then
              write(*,*)'missing covariate for individual #:',i
              stop
           else
              call tokenize(trim(str),' ',nfixed,tmp)
              do k=1,nfixed
                 W(count,k)=cast_float(trim(tmp(k)))
              enddo
           endif
        endif
     end do
     deallocate(tmp)
     if(mcmc) then
        alphanull=1
        wpw=0.d0
        do i=1,nfixed
           if(all(W(:,i)==W(1,i))) then
              alphanull(i)=0
              wpw(i)=0.0
           else 
              wpw(i)=dot_product(W(:,i),W(:,i))
           endif
        enddo
        alpha=0.0d0
        alphastore=0.0d0
     else
        do i=1,nopt
           if(str_match(trim(register%key(i)),'-alpha')) then
              fname=trim(register%default(i))
              if(fname=="") stop 'Error: No alpha file specified'
              !read alpha from training into alphastore
              open(14,file=trim(fname),status='old',form='formatted',iostat=ios)
              if(ios/=0) then
                 write(21,*) 'Error opening alpha file'
                 stop
              end if
              do ii=1,nfixed
                 read(14,*) alphastore(ii)
              enddo
              close(14,status='keep')
              exit
           end if
        end do
     end if
  end if
 end subroutine parse_covar

subroutine parse_ndist
  integer::i
  do i=1,nopt
     if(str_match(trim(register%key(i)),'-ndist')) then
        ndist = cast_int(register%default(i))
        exit
     end if
  end do
end subroutine parse_ndist

!subroutine parse_priors
!  integer::i,k, nitem
!  character(len=128), dimension(ndist) :: c_string
!  
!  do i=1,nopt
!     if(str_match(trim(register%key(i)),'-vara')) then
!        vara = cast_float(register%default(i))
!     else if(str_match(trim(register%key(i)),'-vare')) then
!        vare = cast_float(register%default(i))
!     else if(str_match(trim(register%key(i)),'-varb')) then
!        varb = cast_float(register%default(i))
!     else if(str_match(trim(register%key(i)),'-dfvara')) then
!        dfvara = cast_float(register%default(i))
!     else if(str_match(trim(register%key(i)),'-dfvare')) then
!        dfvare = cast_float(register%default(i))
!     else if(str_match(trim(register%key(i)),'-gpin')) then
!        nitem=ntokens(register%default(i),(/','/))
!        if(nitem /= ndist) then
!           if(ndist>2) then
!              print *, 'Error: Number of mixtue classes ',ndist
!              print *,'        but ',nitem,'effect sizes specified (gpin)'
!              stop
!           end if
!        else
!           call tokenize(trim(register%default(i)),',',ndist,c_string)
!           do k=1,ndist
!              gpin(k)=cast_float(trim(c_string(k)))
!           enddo
!        endif
!     else if(str_match(trim(register%key(i)),'-delta')) then
!        nitem=ntokens(register%default(i),(/','/))
!        if(nitem==1) then
!           call tokenize(trim(register%default(i)),',',1,c_string)
!           delta=cast_float(trim(c_string(1)))
!        else if(nitem /= ndist) then
!           print *, 'Error: Number of mixtue classes',ndist
!           print *,'        but',nitem,'prior values for Dirichlet pecified (delta)'
!           stop
!        else
!           call tokenize(trim(register%default(i)),',',ndist,c_string)
!           do k=1,ndist
!              delta(k)=cast_float(trim(c_string(k)))
!           enddo
!        endif
!     end if
!  end do
!end subroutine parse_priors

subroutine parse_initialise
  integer::i
  do i=1,nopt
     if(str_match(trim(register%key(i)),'-numit')) then
        numit = cast_int(register%default(i))
     else if(str_match(trim(register%key(i)),'-burnin')) then
        burnin= cast_int(register%default(i))
     else if(str_match(trim(register%key(i)),'-thin')) then
        thin= cast_int(register%default(i))
     else if(str_match(trim(register%key(i)),'-seed')) then
        seed1= cast_int(register%default(i))
     else if(str_match(trim(register%key(i)),'-snpout')) then
        snpout= cast_logical(register%default(i))
     else if(str_match(trim(register%key(i)),'-permute')) then
        permute = cast_logical(register%default(i))
     else if(str_match(trim(register%key(i)),'-msize')) then
        msize = cast_int(register%default(i))
        permute=.true.
     else if(str_match(trim(register%key(i)),'-mrep')) then
        mrep = cast_int(register%default(i))
     endif
#ifdef block
     if(str_match(trim(register%key(i)),'-blocksize')) then
        block_size = cast_int(register%default(i))
        nthreads=block_size
     else if(str_match(trim(register%key(i)),'-nthreads')) then
        nthreads = cast_int(register%default(i))
     endif
#endif     
  end do
end subroutine parse_initialise

subroutine hyp_header
  integer :: i,j,sx
  character (len=10) :: ci, ca

  open(unit=25,file=hypfil,status='unknown',form='formatted')
  !write(25,'(2(A10,1x),A7)',advance='no') 'Replicate','Nsnp',' '
  write(25,'(2(A10,1x))',advance='no') 'Replicate','Nsnp'
  do i=1,ncomp
     write(ci,'(I8)') i
     ci=adjustl(ci)
     ca="Va"//trim(ci)
     if(ncomp==1) ca="Va"
     write(25,'(A15,1x)',advance="no") trim(ca)
  end do
  write(25,'(A15,1x)',advance="no") 'Ve'
  do i=1,nseg
     sx=segments%nmix(i)
     do j=1,sx
        if(nseg==1) then
           write(ci,'(I8)') j
        else
           write(ci,'(I8)') i*10+j
        endif
        ci=adjustl(ci)
        ca="Nk"//trim(ci)
        write(25,'(A10,1x)',advance="no") trim(ca)
     end do
  end do
  do i=1,nseg
     sx=segments%nmix(i)
     do j=1,sx
        if(nseg==1) then
           write(ci,'(I8)') j
        else
           write(ci,'(I8)') i*10+j
        endif
        ci=adjustl(ci)
        ca="Vk"//trim(ci)
        write(25,'(A15,1x)',advance="no") trim(ca)
     end do
  end do
  write(25,*)
end subroutine hyp_header

subroutine parse
  call init_register
  call get_cmdLine
  call parse_help
  call update_register
  call parse_out
  call parse_plink
  call parse_flatinput
  call parse_trait_pos
  call parse_predict
  call parse_ndist
  call parse_initialise
end subroutine parse

end module cmd_parser

module perm
  implicit none
contains
subroutine permutate(p,n)
  integer :: n,p(n), k,j,i,ipj,itemp,m
  double precision :: u(100)
  do i=1,n
     p(i)=i
  enddo
  ! generate up to 100 u(0,1) numbers at a time.
  do i=1,n,100
     m=min(n-i+1,100)
     call random_number(u)
     do j=1,m
        ipj=i+j-1
        k=int(u(j)*(n-ipj+1))+ipj
        itemp=p(ipj)
        p(ipj)=p(k)
        p(k)=itemp
     enddo
  enddo
end subroutine permutate
end module perm

module bays
  use RDistributions
  use parz
  use perm
  use multiprior
  implicit none

contains

  subroutine update_fixed
    integer :: k
    double precision, dimension(:), pointer :: z
    double precision :: rhs, gk,v1,zz

    do k=1,nfixed
       if(alphanull(k)==1) then
          z=>W(:,k)
          zz=wpw(k)
          gk=alpha(k)
          rhs=dot_product(yadj,z)/vare
          rhs=rhs + zz*gk/vare
          v1=zz/vare+1.0d0/varb
          alpha(k)=rand_normal(rhs/v1, dsqrt(1.0d0/v1))
          gk=alpha(k)-gk
          yadj=yadj-z*gk
       end if
    end do
  end subroutine update_fixed

#ifndef block
  subroutine update_bayesR
    integer :: snploc, indistflag, sidx, nx
    double precision, dimension(nseg,mxdist) :: log_p, vare_gp
    double precision, dimension(mxdist) :: s, stemp
    double precision, dimension(:), pointer :: z
    double precision :: zz, zz_vare, gk, rhs, logdetV, uhat, sk, skk, &
         clike, ssculm, v1, rd
    logical :: overflow
    integer :: i, j, k, kk

    included=0
    log_p(:,1)=dlog(segments%p(:,1))
    do i=1,nseg
       nx=segments%nmix(i)
       do j=2,nx
          log_p(i,j)=dlog(segments%p(i,j))
          vare_gp(i,j)=vare/segments%gp(i,j)
       enddo
    enddo
    if(permute) then
       call permutate(permvec,nloci)
    endif
    do k=1,nloci
       snploc=permvec(k)
       sidx=snptracker(snploc,1)
       indistflag=snptracker(snploc,3)
       nx=segments%nmix(sidx)
       z => X(:,snploc)
       zz=xpx(snploc)
       gk=g(snploc)
       if(indistflag > 1) then
          yadj=yadj+z*gk
       endif
       rhs= dot_product(yadj,z)
       
       !s(1)=log_p(1)
       s(1)=log_p(sidx,1)
       do kk=2,nx
          !logdetV=dlog(gp(kk)*zz_vare+1.0d0)
          logdetV=dlog(segments%gp(sidx,kk)*zz/vare+1.0d0)
          !uhat=rhs/(zz+vare_gp(kk))
          !uhat=rhs/(zz+vare/segments%gp(sidx,kk))
          uhat=rhs/(zz+vare_gp(sidx,kk))
          !s(kk)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(kk)
          !s(kk)=-0.5d0*(logdetV-(rhs*uhat/vare))+dlog(segments%p(sidx,kk))
          s(kk)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(sidx,kk)
       enddo
       stemp=0.0d0
       do kk=1,nx
          skk=s(kk)
          sk=0.0d0
          overflow=.false.
          do j=1,nx
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
       call random_number(rd)
       indistflag=1
       do kk=1,nx
          ssculm=ssculm+stemp(kk)
          if (rd<ssculm) then
             indistflag=kk
             exit
          endif
       enddo
       snptracker(snploc,3)=indistflag
       snpindist(snploc,indistflag)= snpindist(snploc,indistflag)+1
       if(indistflag==1) then
          gk=0.0d0
       else
          !v1=zz+vare/gp(indistflag)
          v1=zz+vare/segments%gp(sidx,indistflag)
          gk=rand_normal(rhs/v1, dsqrt(vare/v1))
          yadj=yadj-z*gk  
          included=included+1
       endif
       g(snploc)=gk
       if(msize>0 .and. rep>mrep) then
          if(included>=msize) exit
       endif
    enddo ! each loci
  end subroutine update_bayesR
#endif
  
  subroutine update_varcomp
    integer :: i, inc

    varcomp%vara=0.0d0
    do i=1,ncomp
       inc=count(snptracker(:,3) > 1 .and. snptracker(:,2)==i)
       scale=(inc*sum(g*g, mask= snptracker(:,2)==i) +varcomp%vara_ap(i)*varcomp%dfvara(i)) / &
            (varcomp%dfvara(i)+inc)
       varcomp%vara(i)=rand_scaled_inverse_chi_square(inc+varcomp%dfvara(i),scale)
    enddo
    vare=(dot_product(yadj,yadj)+vare_ap*dfvare)/ (dble(nt)+dfvare)
    vare=rand_scaled_inverse_chi_square(dble(nt)+dfvare,vare)
  end subroutine update_varcomp
    
 subroutine update_segments
    integer :: i, sidx, nm, ns, indistflag

    segments%snpinseg(:,:)=0
    segments%varinseg(:,:)=0.0d0
    do i=1,nloci
       sidx=snptracker(i,1)
       indistflag=snptracker(i,3)
       segments%snpinseg(sidx,indistflag)=segments%snpinseg(sidx,indistflag)+1
       segments%varinseg(sidx,indistflag)=segments%varinseg(sidx,indistflag)+g(i)*g(i)
    enddo

    do i=1,nseg
       nm=segments%nmix(i)
       ns=segments%varcomp(i)
       segments%gp(i,1:nm)=segments%gpin(i,1:nm)*varcomp%vara(ns)
       dirx(1:nm)=dble(segments%snpinseg(i,1:nm))+segments%delta(i,1:nm)
       segments%p(i,1:nm)=rdirichlet(nm,dirx(1:nm))
    enddo
  end subroutine update_segments
  
end module bays

module blockutilz
  use RDistributions
  use parz
  use perm
  use multiprior
  !use cmd_parser
  !use routinez
  !use cmd_parser
  !varibales for block update
  integer :: rhs_ncol, n_blocks, rest_size
  integer, parameter :: ngeno=4 
  integer, dimension(:,:),allocatable :: rhs_block, rhs_ind, block_info
  double precision, dimension(:,:), allocatable :: rhs_count
  integer, dimension(:), allocatable :: lookup, map

contains

  subroutine block_dimensions
    allocate(geno_mean(0:3,nloci))
    rhs_ncol=ngeno**block_size
    n_blocks=int(nloci/block_size)
    rest_size=nloci-(n_blocks*block_size)
    if(rest_size>0) then 
       n_blocks=n_blocks+1
    end if
    allocate(rhs_block(rhs_ncol,block_size), rhs_ind(nt,n_blocks), &
         rhs_count(rhs_ncol,n_blocks),block_info(n_blocks,block_size), &
         map(block_size), stat=ios)
    if( ios /= 0 )then
       stop 'Unable to allocate required storage for rhs block'
    endif
  end subroutine block_dimensions
 
  subroutine init_block
    call block_dimensions
    call build_blocks
    call rhs_group
    call mapfun
    call process_snp
  end subroutine init_block
  
  subroutine process_snp
    integer, dimension(nt,block_size) :: Xblock
    integer :: i, j,k, tr, bl, nl, b_size, m
    integer*1 :: b1
    integer(kind=1), dimension(0:3) :: igen 
    integer(kind=1) :: val

    igen(0)=0; igen(1)=3;igen(2)=1;igen(3)=2
    open (unit=41,file=trim(genfil),status='old',access='stream',form='unformatted')
    read(41)b1           !plink magic number 1
    if (.not.btest(b1,0).and..not.btest(b1,1).and.btest(b1,2).and.btest(b1,3).and.&
         & .not.btest(b1,4).and.btest(b1,5).and.btest(b1,6).and..not.btest(b1,7)) then
       write(21,'(a)') 'Plink magic number 1 - ok'
    else
       write(*,*)  'Binary genotype file may not be a plink file'
       write(21,'(a)') 'Binary genotype file may not be a plink file'
       call flush(21)
    end if

    read(41)b1           !plink magic number 2
    if (btest(b1,0).and.btest(b1,1).and..not.btest(b1,2).and.btest(b1,3).and.&
         & btest(b1,4).and..not.btest(b1,5).and..not.btest(b1,6).and..not.btest(b1,7)) then
       write(21,'(a)') 'Plink magic number 2 - ok'
    else
       write(*,*)  'Binary genotype file may not be a plink file'
       write(21,'(a)') 'Binary genotype file may not be a plink file'
       call flush(21)
    end if

    read(41)b1           !mode 
    if (btest(b1,0)) then
       write(21,'(a)') 'SNP-major mode for PLINK .bed file'
    else
       write(*,'(a)') 'should not be individual mode - check >>>'
       write(21,'(a)') 'SNP file should not be in individual mode'
       call flush(21)
    end if

    bl=1
    b_size=block_size-count(block_info(bl,:)==0)
    nl=1 !#loci in block
    do j=1,nloci
       k=0
       tr=1
       do i=1,nind
          if(k==0) read(41)b1
          val=igen(ibits(b1,k,2))
          if(trains(i)==0) then
             Xblock(tr,nl)=val
             tr=tr+1
          endif
          k=mod(k+2,8)
       enddo
       nl=nl+1
       if((nl-1)==b_size) then
          call genotype_mean
          call rhs_ct
          call compute_residuals
          if(j==nloci) exit
          bl=bl+1
          b_size=block_size-count(block_info(bl,:)==0)
          nl=1
       endif
    enddo
    close(41,status='keep')

  contains
    
    subroutine genotype_mean
      integer :: i,jj,pos,m, ios
      double precision :: q,qtest,stdev,mean
      integer, dimension(ngeno) :: gc

      do jj=1,b_size
         pos=(j-b_size)+jj

         if(mcmc) then
            gc=0
            nomiss=count(Xblock(:,jj) < 3)
            q=dble(sum(Xblock(:,jj),mask= Xblock(:,jj) < 3)) / (2.0d0*nomiss)
            mean=2.0d0*q
            stdev=dsqrt(mean*(1.0d0-q))
            if(q==1.0d0 .or. q==0.0d0) then
               geno_mean(:,pos)=0.0d0
               xpx(pos)=0.0d0
            else 
               do m=1,ngeno
                  gc(m)=count(Xblock(:,jj)==(m-1))
               enddo
               geno_mean(0,pos)=(0.0d0-mean)/stdev
               geno_mean(1,pos)=(1.0d0-mean)/stdev
               geno_mean(2,pos)=(2.d0-mean)/stdev
               geno_mean(3,pos)=0.0d0
               xpx(pos)=dot_product(dble(gc),geno_mean(:,pos)**2)
            endif
            freqstore(pos)=q
            if(pos==nloci) then
               open(45,file=trim(freqfil),status='unknown')
               do m=1,nloci
                  write(45,'(F10.6)') freqstore(m)
               enddo
               close(45,status='keep')
            endif
         else
          if(pos==1) then
             open(45,file=trim(freqfil),status='unknown')
             do m=1,nloci
                read(45,'(E15.7)') freqstore(m)
             enddo
             close(45,status='keep')
          endif
          q=freqstore(pos)
          mean=2.0d0*q
          stdev=dsqrt(mean*(1.0d0-q))
          if(q==1.0d0 .or. q==0.0d0) then
             geno_mean(:,pos)=0.0d0
          else
             nomiss=count(Xblock(:,jj) < 3)
             if(nomiss==0) then
                geno_mean(0:3,pos)=0.0d0
             else
                qtest=dble(sum(Xblock(:,jj),mask= Xblock(:,jj) < 3)) / (2.0d0*nomiss)
                geno_mean(0,pos)=(0.0d0-mean)/stdev
                geno_mean(1,pos)=(1.0d0-mean)/stdev
                geno_mean(2,pos)=(2.d0-mean)/stdev
                geno_mean(3,pos)=(2.0d0*qtest-mean)/stdev
             end if
          end if
       endif
    end do
  end subroutine genotype_mean

  subroutine rhs_ct
    integer:: i, ii,idx
    integer, dimension(rhs_ncol) :: countvec, countvec_temp
    
    countvec=0
    !$OMP parallel private(i,idx,countvec_temp) if(nt >5000)
    countvec_temp=0
    !$OMP DO
    do i=1,nt
       idx=lookup(sum(Xblock(i,1:b_size)*map(1:b_size))+1)
       !rhs_count(pos,i)=rhs_count(pos,i)+1
       countvec_temp(idx)=countvec_temp(idx)+1
       rhs_ind(i,bl)=idx
    enddo
    !$OMP end do
    !$OMP critical
    countvec=countvec+countvec_temp
    !$OMP end critical
    !$OMP barrier
    rhs_count(:,bl)=countvec
    !$OMP end parallel
  end subroutine rhs_ct

  subroutine compute_residuals()
    !used only for stanard implementation
    integer :: tr, i, ii, pos

    do jj=1,b_size
       pos=(j-b_size)+jj
       tr=0
       do i=1,nind
          if(trains(i)==0) then
             tr=tr+1
             if(pos==1) then
                yadj(tr)=why(i)-mu
             endif
             yadj(tr)=yadj(tr)-geno_mean(Xblock(tr,jj),pos)*g(pos)
          endif
       enddo
    enddo
  end subroutine compute_residuals
end subroutine process_snp

subroutine mapfun
    integer :: i,n,n1, tp, id
    integer, dimension(block_size) :: temp
    
      temp(:)=3
      map(1)=1
      do i=2,block_size
         map(i)=ngeno**i
      enddo
      n1=sum(temp*map)+1
      allocate(lookup(n1))
      lookup(:)=0
      do i=1,rhs_ncol
         n=sum(rhs_block(i,:)*map)+1
         lookup(n) = i
      enddo
    end subroutine mapfun
  
  subroutine build_blocks
    integer :: i,bs,start,end

    block_info=0
    !do i=1,nloci
    !   permvec(i)=i
    !enddo
    start=0
    end=0
    do i=1,n_blocks
       if(i==n_blocks .and. rest_size>0) then
          bs=rest_size
       else 
          bs=block_size
       end if
       start=end+1
       end=end+bs
       block_info(i,1:bs)=[(id,id=start,end)]
    enddo
  end subroutine build_blocks

  subroutine rhs_group
    integer :: i,j
    integer, dimension(block_size) :: ks
    integer, dimension((block_size-1)) :: rs
    
    ks=0
    rs=0
    do i=1,block_size-1
       rs(i)=ngeno**i
    enddo
    do i=1,rhs_ncol
       do j=1,block_size
          rhs_block(i,j)=ks(j)
       enddo
       ks(1)=ks(1)+1
       do j=1,block_size-1
          if(mod(i,rs(j))==0) then
             ks(j)=0
             ks(j+1)=ks(j+1)+1
          endif
       end do
    end do
  end subroutine rhs_group

  subroutine update_bayesR_block
   integer :: i, j,k, kk, bl, b_size, x1, pos, snploc, b, &
       indistflag, sidx, nx
  double precision :: gk, go, zz, rhs, logdetV, uhat, sk, skk, &
         clike, ssculm, v1, rd
  integer, dimension(block_size) :: s_loc
  double precision, dimension(block_size) :: xe, gbl
  double precision, dimension(rhs_ncol) :: resvec, resvec_temp, delta_e
  double precision, dimension(ngeno,block_size) :: bres
  double precision, dimension(mxdist) :: s, stemp
  double precision, dimension(nseg,mxdist) :: log_p, vare_gp
  logical :: changed, overflow

  included=0
  log_p(:,1)=dlog(segments%p(:,1))
  do i=1,nseg
     nx=segments%nmix(i)
     do j=2,nx
        log_p(i,j)=dlog(segments%p(i,j))
        vare_gp(i,j)=vare/segments%gp(i,j)
     enddo
  enddo
  
  if(permute) then
     call permutate(permvec,n_blocks)
  endif

  eachblock : do b=1,n_blocks
     bl=permvec(b)
     if(bl==n_blocks .and. rest_size>0) then
        b_size=rest_size
     else
        b_size=block_size
     end if
     s_loc(1:b_size)=block_info(bl,1:b_size)
     resvec=0.0d0
     !$OMP parallel private(k,x1,resvec_temp)
     resvec_temp=0.0d0
     !$OMP DO
     do k=1,nt
        x1=rhs_ind(k,bl)
        resvec_temp(x1)=resvec_temp(x1)+yadj(k)
     enddo
     !$OMP END DO
     !$omp critical
     resvec=resvec+resvec_temp
     !$OMP END critical
     !$OMP barrier
     !$OMP DO
     do k=1,b_size
        bres(1,k)=sum(resvec,mask = rhs_block(:,k)== 0)
        bres(2,k)=sum(resvec,mask = rhs_block(:,k)== 1)
        bres(3,k)=sum(resvec,mask = rhs_block(:,k)== 2)
        bres(4,k)=sum(resvec,mask = rhs_block(:,k)== 3)
     enddo
     !$OMP end do
     !$OMP end parallel
     do k=1,b_size
        xe(k)=dot_product(bres(:,k),geno_mean(:,s_loc(k)))
     enddo
     delta_e=0.0d0 
     gbl=0.0d0
     changed=.false.
     within_block : do pos=1,b_size
        rhs=xe(pos)
        if(pos>1) then
           if(changed) call xe_update
        end if
        snploc=s_loc(pos)
        sidx=snptracker(snploc,1)
        indistflag=snptracker(snploc,3)
        nx=segments%nmix(sidx)
        zz=xpx(snploc)
        go=g(snploc)
        if(indistflag > 1) rhs=rhs+zz*go
        s(1)=log_p(sidx,1)
        do kk=2,nx
           logdetV=dlog(segments%gp(sidx,kk)*zz/vare+1.0d0)
           uhat=rhs/(zz+vare_gp(sidx,kk))
           s(kk)=-0.5d0*(logdetV-(rhs*uhat/vare))+log_p(sidx,kk)
        enddo
        stemp=0.0d0
        do kk=1,nx
           skk=s(kk)
           sk=0.0d0
           overflow=.false.
           do j=1,nx
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
        call random_number(rd)
        indistflag=1
        do kk=1,nx
           ssculm=ssculm+stemp(kk)
           if (rd<ssculm) then
              indistflag=kk
              exit
           endif
        enddo
        snpindist(snploc,indistflag)= snpindist(snploc,indistflag)+1
        if(indistflag==1) then
           gk=0.0d0
        else
           v1=zz+vare/segments%gp(sidx,indistflag)
           gk=rand_normal(rhs/v1, dsqrt(vare/v1))
           included=included+1
        end if
        if(indistflag>1 .or. snptracker(snploc,3)>1) then
           changed=.true.
        end if
        gbl(pos)=gk-go
        g(snploc)=gk
        snptracker(snploc,3)=indistflag
     end do within_block

     call res_update_rhs

     if(msize>0 .and. rep>mrep) then
        if(included>=msize) exit
     endif
  end do eachblock

contains
     
  subroutine xe_update
    integer ::i, loc, s1, s2
    double precision, dimension(0:(ngeno-1)) :: ny2
    
    loc=s_loc(pos-1)
    ny2=0.0d0
    do i=1,rhs_ncol
       if(rhs_count(i,bl)==0.0d0) cycle
       s1=rhs_block(i,pos-1)
       s2=rhs_block(i,pos)
       delta_e(i)=delta_e(i)+rhs_count(i,bl)*geno_mean(s1,loc)*gbl(pos-1)
       ny2(s2)=ny2(s2)+delta_e(i)
    enddo
    rhs = xe(pos) -dot_product(geno_mean(:,s_loc(pos)),ny2)
  end subroutine xe_update


  subroutine res_update_rhs
    integer :: i,k,s
    double precision, dimension(rhs_ncol) :: delta_res

    delta_res=0.0d0
    do k=1,b_size
       if(gbl(k) /= 0.0d0) then
          do i=1,rhs_ncol
             s=rhs_block(i,k)
             delta_res(i)=delta_res(i)+geno_mean(s,s_loc(k))*gbl(k)
          end do
       end if
    end do
    if(any(gbl /= 0.0d0)) then
       !$OMP PARALLEL DO private(i)
       do i=1,nt
          yadj(i)=yadj(i)-delta_res(rhs_ind(i,bl))
       end do
       !$OMP end parallel do 
    end if
  end subroutine res_update_rhs
  
end subroutine update_bayesR_block

end module blockutilz

module routinez
  use parz
  use perm
  use multiprior
#ifdef block
  use blockutilz
#endif
  implicit none

contains

subroutine get_size()
  nind=0
  open(35,file=trim(phenfil),status='old',form='formatted')
  do 
     read(35,*,iostat=ios)
     if (ios.ne.0) exit
     nind=nind+1
  enddo
  close(35,status='keep')
  nloci=0
  open(36,file=trim(bimfil),status='old',form='formatted')
  do 
     read(36,*,iostat=ios)
     if (ios.ne.0) exit
     nloci=nloci+1
  enddo
  close(36,status='keep')
end subroutine get_size

subroutine load_phenos()
  character(len=1024):: str 
  character(len=20) :: why_str
  integer ::pos1, pos2,n ,i, incr

  if(plink) then
     incr=5
  else
     incr=0
  endif
  allocate(trains(nind), why(nind))
  open(31,file=trim(phenfil),status='old',form='formatted')
  do i=1,nind
    read(31,'(a)') str
    pos1=1
    n=0
    do
       pos2 = index(str(pos1:), " ")
       if (pos2 == 0) then
          n=n+1
          why_str = str(pos1:)
          exit
       endif
       n=n+1
       why_str = str(pos1:pos1+pos2-2)
       if (n== trait_pos+incr) then
          exit
       endif
       pos1 = pos2+pos1
    enddo
    why_str=trim(why_str)
    if(why_str /=  'NA') then
       trains(i)=0
       read(why_str,*) why(i)
    else
       trains(i)=1
       why(i)=-9999
    endif
  enddo
  close(unit=31,status='keep')
end subroutine load_phenos

subroutine load_snp
 if(plink) then
     call load_snp_binary
  else
     call load_snp_flat
  end if
end subroutine load_snp

subroutine load_snp_flat()
  integer :: i
  open(31,file=trim(genfil),status='old',form='formatted')
  do i=1,nind
     read(31,'(10000000f1.0)') X(i,1:nloci)
  enddo
  close(unit=31,status='keep')
end subroutine load_snp_flat

subroutine load_snp_binary()
integer :: i, j,k, tr,ios
integer*1 :: b1
double precision, dimension(0:3) :: igen 
double precision :: val
igen(0)=0.0d0; igen(1)=3.0d0;igen(2)=1.0d0;igen(3)=2.0d0

 open (unit=41,file=trim(genfil),status='old',access='stream',form='unformatted')
 read(41)b1           !plink magic number 1
 if (.not.btest(b1,0).and..not.btest(b1,1).and.btest(b1,2).and.btest(b1,3).and.&
      & .not.btest(b1,4).and.btest(b1,5).and.btest(b1,6).and..not.btest(b1,7)) then
    write(21,'(a)') 'Plink magic number 1 - ok'
 else
    write(*,*)  'Binary genotype file may not be a plink file'
    write(21,'(a)') 'Binary genotype file may not be a plink file'
    call flush(21)
 end if

 read(41)b1           !plink magic number 2
 if (btest(b1,0).and.btest(b1,1).and..not.btest(b1,2).and.btest(b1,3).and.&
      & btest(b1,4).and..not.btest(b1,5).and..not.btest(b1,6).and..not.btest(b1,7)) then
    write(21,'(a)') 'Plink magic number 2 - ok'
 else
    write(*,*)  'Binary genotype file may not be a plink file'
    write(21,'(a)') 'Binary genotype file may not be a plink file'
    call flush(21)
 end if

 read(41)b1           !mode 
 if (btest(b1,0)) then
    write(21,'(a)') 'SNP-major mode for PLINK .bed file'
 else
    write(*,'(a)') 'should not be individual mode - check >>>'
    write(21,'(a)') 'SNP file should not be in individual mode'
    call flush(21)
 end if

 do j=1,nloci
    k=0
    tr=1
    do i=1,nind
       if(k==0) read(41)b1
       val=igen(ibits(b1,k,2))
       if(trains(i)==0) then
          X(tr,j)=val
          tr=tr+1
       endif
       k=mod(k+2,8)
    enddo
 enddo
 close(41,status='keep')
end subroutine load_snp_binary

subroutine init_random_seed()
! one integer for seed in parameter file
! seed1 >0: initialized with seed1
  integer :: i, n, clock
  integer,dimension(:), allocatable :: seed

  call random_seed(size = n)
  allocate(seed(n))
  
  if (seed1 /= 0) then
     do i=1,n
        seed(i)=abs(seed1)+(i-1)
     enddo
     write(21,'(a,i8)') 'New seeds generated from seed1 ',seed1
  else
     write(21,*)  'new seeds generated from clock'
     call system_clock(count=clock)
     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  endif
  call random_seed(put = seed)
  deallocate(seed)
end subroutine init_random_seed


subroutine xcenter()
  !RESCALE DESIGN MATRIX SO THAT GENOTYPES ARE CHANGED FROM BEING CODED
  !AS 0,1,2 TO BEING CODED AS DEVIATIONS FROM ZERO
  integer :: j, nomiss
  double precision :: q, qtest
  double precision, dimension(nt) :: xtemp
  if(mcmc) then
     xpx=0.0d0
     do j=1,nloci
        xtemp=X(:,j)
        nomiss=count(xtemp < 3.0d0)
        q=sum(xtemp,mask= xtemp < 3.0d0) / (2.0d0*nomiss)
        if(q==1.0d0 .or. q==0.0d0) then
           X(:,j)=0.0d0
        else 
           where (xtemp >2.0d0) xtemp=2.0d0*q
           X(:,j)=(xtemp-2.0d0*q)/sqrt(2.0d0*q*(1.0d0-q))
        endif
        freqstore(j)=q
        xpx(j)=dot_product(X(:,j),X(:,j))
     enddo
     open(45,file=trim(freqfil),status='unknown')
     do j=1,nloci
        write(45,'(F10.6)') freqstore(j)
     enddo
     close(45,status='keep')
  else
     open(45,file=trim(freqfil),status='unknown')
      do j=1,nloci
        read(45,'(E15.7)') freqstore(j)
     enddo
     close(45,status='keep')
     do j=1,nloci
        q=freqstore(j)
        if(q==1.0d0 .or. q==0.0d0) then
           X(:,j)=0.0d0
        else
           xtemp=X(:,j)
           nomiss=count(xtemp < 3)
           qtest=dble(sum(xtemp,mask= xtemp < 3)) / (2.0d0*nomiss)
           where (xtemp >2.0d0) xtemp=2.0d0*qtest
           X(:,j)=(xtemp-2.0d0*q)/sqrt(2.0d0*q*(1.0d0-q))
        endif
     enddo
  endif
end subroutine xcenter

subroutine xcenter2()
  !RESCALE DESIGN MATRIX SO THAT GENOTYPES ARE CHANGED FROM BEING CODED                            
  !AS 0,1,2 TO BEING CODED AS DEVIATIONS FROM ZERO                                                 
  integer :: j
  double precision :: xbar, xvar
  double precision, dimension(nt) :: xtemp
  if(mcmc) then
     do j=1,nloci
        xtemp=X(:,j)
        xbar = sum(xtemp,mask= xtemp < 3.0d0)/dble(nt)
        where (xtemp >2.0d0) xtemp=xbar
        xvar=sum(xtemp**2)/dble(nt)
        X(:,j)=(xtemp-xbar)/dsqrt(xvar)
     enddo
  endif
end subroutine xcenter2

subroutine allocate_data
  !trains keeps track of training (0) and test(1) individuals
  !swap around for predicition 
  if(.not.mcmc) then
     where(trains==0) trains=3
     where(trains==1) trains=0
     where(trains==3) trains=1
  end if
  nt=count(trains==0)
  allocate(X(nt,nloci),pred(nind), dirx(ndist), g(nloci), yadj(nt), snpindist(nloci,mxdist), &
       gstore(nloci), indiststore(nloci,ndist), freqstore(nloci), xpx(nloci), &
       snptracker(nloci,3), stat=ios)
  if( ios /= 0 )then
     write(21,'(a50)') 'ERROR :: Unable to allocate required storage'
     stop 'Unable to allocate required storage for data'
     call flush(21)
  endif
end subroutine allocate_data

subroutine compute_fitted
  integer :: i, tr
  tr=0
  do i=1,nind
     if(trains(i)==0) then
        tr=tr+1
        pred(i)=pred(i)+dot_product(W(tr,1:nfixed),alphastore(1:nfixed))
     endif
  end do
end subroutine compute_fitted

subroutine write_dgv()
 integer :: i
 character(len=2) :: missvalue

 missvalue='NA'
 open(unit=61,file=mbvfil,status='unknown',form='formatted')
 do i=1,nind
    if (trains(i)==0) then
       !write (rchar, '(E15.7)') pred(i)
       !write(61,*)adjustl(rchar)
       write(61,'(E15.7e3)') pred(i)
    else
       write(61,'(a2)') missvalue
    end if
 enddo
 close(unit=61,status='keep')
end subroutine write_dgv

subroutine write_predval(infile,fitted)
  character(len=200), intent(in) :: infile
  double precision , dimension(:), intent(in) :: fitted
  integer :: i
  character(len=2) :: missvalue

 missvalue='NA'
 open(unit=61,file=trim(infile),status='unknown',form='formatted')
 do i=1,nind
    if (trains(i)==0) then
       !write (rchar, '(E15.7)') pred(i)
       !write(61,*)adjustl(rchar)
       write(61,'(E15.7e3)') fitted(i)
    else
       write(61,'(a2)') missvalue
    end if
 enddo
 close(unit=61,status='keep')
end subroutine write_predval

subroutine load_param()
  character(len=10):: dummy
  integer :: i, nc, ios
  double precision, dimension(20):: gtemp
  
  open(31,file=paramfil,status='old',form='formatted')
  read(31,*)
  nc=ndist+1
  do i=1,nloci
     read(31,*,iostat=ios) gtemp(1:nc)
     gstore(i)=gtemp(nc)
     if(ios/=0) then
        write(21,'(a)') 'Error opening file ',adjustl(paramfil)
     call flush(21)
     end if
  enddo
  close(31)
  open(32,file=modfil,status='old',form='formatted')
  read(32,*) dummy, mu
  close(32)
end subroutine load_param

subroutine output_hyp
  integer :: i, mx

  write(25,'(2(i10,1x))', advance='no') rep, included                                             
  write(25,'(200(E15.7e3,1x))',advance='no')  (varcomp%vara(i),i=1,ncomp), vare                         
  do i=1,nseg
     mx=segments%nmix(i)
     write(25,'(200(i10,1x))',advance='no')  segments%snpinseg(i,1:mx)
  enddo
  do i=1,nseg
     mx=segments%nmix(i)
     write(25,'(200(E15.7e3,1x))',advance='no')  segments%varinseg(i,1:mx)
  enddo
  write(25,*)
  call flush(25)
end subroutine output_hyp


subroutine write_log(sel,cdat,ctim)
  integer :: i,nm
  character(len=*), intent(in) :: sel
  character(len=10), intent(in) :: ctim
  character(len=8), intent(in) :: cdat

  if(trim(sel) /= 'end') then
     open(unit=21,file=logfil,status='unknown',form='formatted')
     write(21,901) 'Program BayesR'
     write(21,908) 'Run started at',cdat(1:4),cdat(5:6),cdat(7:8), &
          ctim(1:2),ctim(3:4),ctim(5:6)
     write(21,902) 'Prefix for input files',trim(inprefix)
     write(21,902) 'Prefix for output files',trim(outprefix)
     
     if(trim(sel)=='train') then
        write(21,903) 'Phenotype column',trait_pos
        write(21,903) 'No. of loci',nloci
        write(21,903) 'No. of individuals',nind
        write(21,903) 'No. of training individuals',nt
        if(ncomp==1) then
           write(21,906) 'Prior Vara', varcomp%vara(1), varcomp%dfvara(1)
        else
           write(21,902) 'SNP model file',trim(snpmodfil)
           write(21,902) 'Variance model file',trim(varcompfil)
           write(21,901) 'Prior Vara '
           do i=1,ncomp
              write(21,'(i4,2E15.7)') i,varcomp%vara(i),varcomp%dfvara(i)
           enddo
        endif
        if(nseg==1) then
           write(21,903) 'No. of mixtures',ndist
           write(21,905) 'Variance of dist ', segments%gpin(1,1:ndist)
           write(21,905) 'Dirichlet prior', segments%delta(1,1:ndist)
        else
           write(21,902) 'Segment model file',trim(segmentfil)
           write(21,901) 'Prior segements (No. of mixtures, Variance of dist, Dirichlet prior)'
           do i=1,nseg
              nm=segments%nmix(i)
              write(21,910) i,nm,segments%gpin(i,1:nm),segments%delta(i,1:nm)
           enddo
        end if
        write(21,906) 'Prior Vare', vare, dfvare
        if(msize > 0) then
           write(21,903) 'Model size (reduced update)',msize
           write(21,903) 'Full cycles (reduced update)', mrep
        endif
        write(21,903) 'No. of cycles',numit
        write(21,903) 'Burnin ',burnin
        write(21,903) 'Thinning rate',thin
#ifdef block
        write(21,903) 'Block size',block_size
        write(21,903) 'Number of threads', nthreads
#endif
        write(21,903) 'Seed ', seed1
        write(21,909) 'SNP output ', snpout
     else 
        write(21,902) 'Model file (training)',trim(modfil)
        write(21,902) 'Parameter file (training)',trim(paramfil)
        write(21,902) 'SNP frequency file (training)',trim(freqfil)
        write(21,903) 'Phenotype column',trait_pos
        write(21,903) 'No. of loci',nloci
        write(21,903) 'No. of individuals',nind
        write(21,903) 'No. of individuals to predict',nt
     endif
     flush(21)
  else
     write(21,908) 'Run ended at',cdat(1:4),cdat(5:6),cdat(7:8), &
          ctim(1:2),ctim(3:4),ctim(5:6)
     close(21)
  endif
901 format(a)
902 format(a,t34,': ',a)
903 format(a,t34,'= ',i8)
904 format(a,t34,'= ',f20.6)
905 format(a,t34,'= ',10f10.5)
906 format(a,t34,'= ',2f10.6)
907 format(a,t34,'= ',f10.2,a)
908 format(a,1x,a4,'-',a2,'-',a2,' ',a2,':',a2':',a2)
909 format(a,t34,'= ',l)
910 format(2i8,25f10.5)
end subroutine write_log

subroutine output_model
  integer :: i,ns,j
  character(LEN=10) :: ci,ca

  ns=maxval(segments%nmix)
  open(unit=10,file=paramfil,status='unknown',iostat=ios)
  if(ios/=0) then
     write(21,*) 'Error opening ',trim(paramfil)
     stop
  end if
  do i=1,ns
     write(ci,'(I8)') i
     ci=adjustl(ci)
     ca="PIP"//trim(ci)
     write(10,'(A15,1x)',advance="no") trim(ca)
  end do
  write(10,'(A15)') 'alpha'
  do i=1,nloci
     write(10,'(50(E15.7,1X))') indiststore(i,:), gstore(i)
  enddo
  close(10,status='keep')

  open(unit=12,file=modfil,status='unknown',iostat=ios)
  if(ios/=0) then
     write(21,*) 'Error opening ',trim(modfil)
     stop
  end if
  write(12,800) 'Mean',mu_s
  write(12,800) 'Nsnp',sum(segments_s%snpinseg(:,2:mxdist))
  write(12,800) 'Va',sum(varcomp_s%vara)
  write(12,800) 'Ve',vare_s
  if(ncomp>1) then
     do i=1,ncomp
        write(ci,'(I8)') i
        ci=adjustl(ci)
        ca="Va"//trim(ci)
        write(12,800) ca, varcomp_s%vara(i)
     enddo
  endif
  do i=1,nseg
     ns=segments%nmix(i)
     do j=1,ns
        if(nseg==1) then
           write(ci,'(I8)') j
        else
           write(ci,'(I8)') i*10+j
        endif
        ci=adjustl(ci)
        ca="Nk"//trim(ci)
        write(12,800) ca,segments_s%snpinseg(i,j)
     end do
  end do
  do i=1,nseg
     ns=segments%nmix(i)
     do j=1,ns
        if(nseg==1) then
           write(ci,'(I8)') j
        else
           write(ci,'(I8)') i*10+j
        endif
        ci=adjustl(ci)
        ca="Pk"//trim(ci)
        write(12,800) ca,segments_s%p(i,j)
     end do
  end do
  do i=1,nseg
     ns=segments%nmix(i)
     do j=1,ns
        if(nseg==1) then
           write(ci,'(I8)') j
        else
           write(ci,'(I8)') i*10+j
        endif
        ci=adjustl(ci)
        ca="Vk"//trim(ci)
        write(12,800) ca, segments_s%varinseg(i,j)
     end do
  end do
  close(12,status='keep')
  if(covar) then
     open(unit=12,file=alphafil,status='unknown',iostat=ios)
     if(ios/=0) then
        write(21,*) 'Error opening ',trim(alphafil)
        stop
     end if
     do i=1,nfixed
        write(12,'(E15.7e3)') alphastore(i)
     enddo
     close(12,status='keep')
  endif
800 format(a,t10,E15.7)
end subroutine output_model

subroutine output_snploc
  integer :: i,j
  character(LEN=20) :: ci, cg, ce
  character(len=1) :: double = ":"

  ! convert int to char
  !write(ci,'(I8)') rep
  !ci=adjustl(ci)

!  ! write SNP effects in sparse format
  do i=1,nloci
     j= snptracker(i,3)
     if(j>1) then
        write(ci,'(I8)') i
        ci=adjustl(ci)
        write(cg,'(I8)') j
        cg=adjustl(cg)
        write(ce,'(e15.6)') g(i)**2
        ce=adjustl(ce)
        write(14,'(1X,a,a,a,a,a)',advance='no') trim(cg),double,trim(ci),double,trim(ce)
     endif
  enddo
  write(14,*)
end subroutine output_snploc

#ifndef block
  subroutine compute_residuals()
    !used only for stanard implementation
    !blocks use compute_residuals in blockutilz
    integer :: tr, i, j, b_size,gidx,jj
    tr=0
    do i=1,nind
       if(trains(i)==0) then
          tr=tr+1
          yadj(tr)=why(i)-mu-dot_product(X(tr,1:nloci),g)
       endif
    enddo
  end subroutine compute_residuals
#endif


  subroutine compute_dgv
    integer :: i, j, jj, tr, bl ,tidx,b_size
    integer, dimension(block_size) :: s_loc
    
    tr=0
    pred=-9999.0d0
    do i=1,nind
       if(trains(i)==0) then
          tr=tr+1
#ifndef block
          pred(i)=mu+dot_product(X(tr,1:nloci),gstore(1:nloci))
#else
          pred(i)=mu
          do bl=1,n_blocks
             if(bl==n_blocks .and. rest_size>0) then
                b_size=rest_size
             else
                b_size=block_size
             end if
             s_loc(1:b_size)=block_info(bl,1:b_size)
             tidx=rhs_ind(tr,bl)
             do j=1,b_size
                jj=s_loc(j)
                pred(i)=pred(i)+gstore(jj)*geno_mean(rhs_block(tidx,j),jj)
             enddo
          enddo
#endif        
       endif
    end do
  end subroutine compute_dgv

end module routinez
