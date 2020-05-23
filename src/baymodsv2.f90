! Bayesian hierarchical models for complex trait analysis using a mixture of 
! normal distributions of SNP effects 
! Copyright (C) 2014 Gerhard Moser
!
!    This program is free software: you can redistribute itg and/or modify
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
  integer :: nloci, nind, nt, ndist, numit, seed1, clock, ios, rep, &
             burnin, thin, ntrain, ntest,counter, trait_pos
  logical :: mcmc, snpout, permute
  !File names
  character(len=200) :: genfil, phenfil, bimfil, inprefix, outprefix, logfil, freqfil, &
                        mbvfil, hypfil, locfil, modfil, paramfil
  ! Data
  double precision, dimension(:), allocatable :: why, pred, freqstore
  integer(kind=1), dimension(:,:), allocatable :: X
  ! Genetic variables
  integer, dimension(:), allocatable ::snpindist, snptracker
  double precision :: vara, vare, mu, scale
  double precision, dimension(:), allocatable :: p, gp, gpin, delta, g, g_old, res, &
       dirx, varindist, res_old, res_null
  ! BayesR specific
  logical :: VCE
  integer:: included, msize, mrep
  double precision :: logdetV, uhat, zz, lhs, v1, total_ssq, rhs, gk, go, vara_ap, vare_ap, dfvara, dfvare
  double precision, dimension(:), allocatable :: xpx, s, stemp
  integer, allocatable, dimension(:)::permvec
! rhs update
  integer :: block_size, rhs_ncol, n_blocks, t_blocks, rest_size, pos, bl, &
             t_size, b_size, nthreads
  integer, parameter :: ngeno=4 
  integer, dimension(:,:),allocatable :: rhs_block, rhs_ind, block_info
  double precision, dimension(:,:), allocatable :: geno_mean, rhs_count,bres
  double precision, dimension(:), allocatable :: class_sum, ny2, delta_e,delta_res, gbl, xe, &
                       resvec, resvec_temp
  integer, dimension(:), allocatable :: s_loc, shufflevec, lookup,map,countvec, countvec_temp
  logical :: shuffle
  !Auxiliary variables
  integer:: indistflag, burn
  integer, dimension(:), allocatable ::trains
  double precision :: xhat, yhat, vary, sk, skk, r, ssculm, nnind, clike
  !Summary statistics
  integer:: nref
  double precision :: msep,bhat, ahat,corr
  ! Storage
  double precision, dimension(:), allocatable :: gstore, pstore, mu_vare_store, &
       snpstore, varstore
  double precision, dimension(:,:), allocatable :: indiststore

  double precision, dimension(:), allocatable :: log_p, log_gp, vare_gp
  double precision :: zz_vare
end module parz

module cmd_parser
 use parz
implicit None

integer :: narg, nopt
character(len=1024) :: arg
character(len=1024), allocatable, dimension(:) :: cmd_line, tmp_array

integer, parameter :: &
     a_int     = 1, &
     a_float   = 2, &
     a_char    = 3, &
     a_flag    = 4 

type param
   integer, allocatable ,dimension(:) :: pos 
   character(len=20), allocatable, dimension(:) :: key
   character(len=20), allocatable, dimension(:) :: argtype
   character(len=100), allocatable,dimension(:) :: desc
   integer, allocatable, dimension(:) :: kind
   character(len=200), allocatable,dimension(:) :: default
end type param

type(param) :: register

contains

  subroutine alloc_register
    allocate(register%pos(nopt))
    allocate(register%key(nopt))
    allocate(register%argtype(nopt))
    allocate(register%desc(nopt))
    allocate(register%kind(nopt))
    allocate(register%default(nopt))
  end subroutine alloc_register

  subroutine init_register
    nopt =25
    call alloc_register
    call include_option(1, '-bfile'  ,'[prefix]',  'prefix PLINK binary files'               ,a_char,     '')
    call include_option(2, '-out'    ,'[prefix]',  'prefix for output'                       ,a_char,     '')
    call include_option(3, '-n'      ,'[num]',     'phenotype column'                        ,a_int,      '1')
    call include_option(4, '-vara'   ,'[num]',     'SNP variance prior'                      ,a_float, '0.01')
    call include_option(5, '-vare'   ,'[num]',     'error variance prior'                    ,a_float, '0.01')
    call include_option(6, '-dfvara' ,'[num]',     'degrees of freedom Va'                   ,a_float, '-2.0')
    call include_option(7, '-dfvare' ,'[num]',     'degrees of freedom Ve'                   ,a_float, '-2.0')
    call include_option(8,'-delta'  ,'[num]',      'prior for Dirichlet'                     ,a_float,   '1.0')
    call include_option(9, '-msize'  ,'[num]',     'number of SNPs in reduced update'        ,a_int,      '0')
    call include_option(10, '-mrep'   ,'[num]',    'number of full cycles in reduced update' ,a_int,   '5000')
    call include_option(11,'-numit' ,'[num]',      'length of MCMC chain'                    ,a_int,   '50000')
    call include_option(12,'-burnin' ,'[num]',     'burnin steps'                            ,a_int,   '20000')
    call include_option(13,'-thin'   ,'[num]',     'thinning rate'                           ,a_int,   '10')
    call include_option(14,'-ndist'  ,'[num]',     'number of mixture distributions'         ,a_int,     '4')
    call include_option(15,'-gpin'   ,'[num]',     'effect sizes of mixtures (% x Va)'       ,a_float,  '0.0,0.0001,0.001,0.01')
    call include_option(16,'-seed'   ,'[num]',     'initial value for random number'         ,a_int,      '0')
    call include_option(17,'-predict','[flag]',    'perform prediction'                      ,a_flag,     'f')
    call include_option(18,'-snpout' ,'[flag]',    'output detailed SNP info'                ,a_flag,     'f')
    call include_option(19,'-permute','[flag]',    'permute order of SNP'                    ,a_flag,   'f')
    call include_option(20,'-model'  ,'[filename]','model summary file (for prediction) '    ,a_char,   '')
    call include_option(21,'-freq'   ,'[filename]','SNP frequency file (for prediction)'     ,a_char,   '')
    call include_option(22,'-param'  ,'[filename]','SNP effect file (for prediction)'        ,a_char,   '')
    call include_option(23,'-blocksize' ,'[num]'  ,'Number of SNP in rhs block'              ,a_int,    '4')
    call include_option(24,'-nthreads' ,'[num]'   ,'Number of threads'                       ,a_int,    '4')
    call include_option(25,'-shuffle' , '[flag]'  ,'permute order of blocks'                 ,a_flag,   'f') 
  end subroutine init_register

  subroutine include_option(pos,key,argtype,desc,kind,default)
    implicit none
    integer :: pos
    character(len=*) :: key
    character(len=*) :: desc
    character(len=*) :: argtype
    integer :: kind
    character(len=*) :: default
   
    register%pos(pos)=pos
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

  integer function ntokens(line)
    !http://www.tek-tips.com/viewthread.cfm?qid=1688013
    character,intent(in):: line*(*)
    integer i, n, toks
    character(len=1) :: separator
    separator=','
    i = 1;
    n = len_trim(line)
    toks = 0
    ntokens = 0
    do while(i <= n) 
       do while(line(i:i) == separator) 
          i = i + 1
          if (n < i) return
       enddo
       toks = toks + 1
       ntokens = toks
       do
          i = i + 1
          if (n < i) return
          if (line(i:i) == separator) exit
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
  
  double precision function cast_float(value)
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
    integer :: isnum
    character(len=2),parameter :: str='-'
    is_key=.false.
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
    err2='Missing argumnet for :'

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
       111 format(a10,a12,a45,a25)
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
          exit
       end if
    end do
  end subroutine parse_out

  subroutine parse_plink
    integer::i
    logical:: fileExist

    do i=1,nopt
       if(str_match(trim(register%key(i)),'-bfile')) then
          inprefix=trim(register%default(i))
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
    enddo
  end subroutine parse_plink

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
             inprefix=trim(register%default(i))
             if(str_match(trim(inprefix),' ')) stop 'Error: No model file specified'
             modfil=trim(inprefix)
             inquire(file=modfil,exist=fileExist)
             if(.not.fileExist)then
                write(*,*)'file ',trim(modfil),' not found'
                stop
             end if
          else if(str_match(trim(register%key(i)),'-freq')) then
             inprefix=trim(register%default(i))
             if(str_match(trim(inprefix),' ')) stop 'Error: No SNP frequency file specified'
             freqfil=trim(inprefix)
             inquire(file=freqfil,exist=fileExist)
             if(.not.fileExist)then
                write(*,*)'file ',trim(freqfil),' not found'
                stop
             end if
          else if(str_match(trim(register%key(i)),'-param')) then
             inprefix=trim(register%default(i))
             if(str_match(trim(inprefix),' ')) stop 'Error: No param file specified'
             paramfil=trim(inprefix)
             inquire(file=freqfil,exist=fileExist)
             if(.not.fileExist)then
                write(*,*)'file ',trim(paramfil),' not found'
                stop
             end if
          end if
       end do
    end if
  end subroutine parse_predict

  subroutine parse_ndist
    integer::i
    do i=1,nopt
       if(str_match(trim(register%key(i)),'-ndist')) then
          ndist = cast_int(register%default(i))
          exit
       end if
    end do
  end subroutine parse_ndist

  subroutine parse_priors
    integer::i,k, nitem
    character(len=128), dimension(ndist) :: c_string

    do i=1,nopt
       if(str_match(trim(register%key(i)),'-vara')) then
          vara = cast_float(register%default(i))
       else if(str_match(trim(register%key(i)),'-vare')) then
          vare = cast_float(register%default(i))
       else if(str_match(trim(register%key(i)),'-dfvara')) then
          dfvara = cast_float(register%default(i))
       else if(str_match(trim(register%key(i)),'-dfvare')) then
          dfvare = cast_float(register%default(i))
       else if(str_match(trim(register%key(i)),'-gpin')) then
          nitem=ntokens(register%default(i))
          if(nitem /= ndist) then
             print *, 'Error: Number of mixtue classes ',ndist
             print *,'        but ',nitem,'effect sizes specified (gpin)'
             stop
          else
             call tokenize(trim(register%default(i)),',',ndist,c_string)
             do k=1,ndist
                gpin(k)=cast_float(trim(c_string(k)))
             enddo
          endif
       else if(str_match(trim(register%key(i)),'-delta')) then
          nitem=ntokens(register%default(i))
          if(nitem==1) then
             call tokenize(trim(register%default(i)),',',1,c_string)
             delta=cast_float(trim(c_string(1)))
          else if(nitem /= ndist) then
             print *, 'Error: Number of mixtue classes',ndist
             print *,'        but',nitem,'prior values for Dirichlet pecified (delta)'
             stop
          else
             call tokenize(trim(register%default(i)),',',ndist,c_string)
             do k=1,ndist
                delta(k)=cast_float(trim(c_string(k)))
             enddo
          endif
       end if
    end do
  end subroutine parse_priors

  subroutine parse_initialise
    integer::i
    do i=1,nopt
       if(str_match(trim(register%key(i)),'-msize')) then
          msize = cast_int(register%default(i))
          permute=.true.
       else if(str_match(trim(register%key(i)),'-mrep')) then
          mrep = cast_int(register%default(i))
       else if(str_match(trim(register%key(i)),'-numit')) then
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
       else if(str_match(trim(register%key(i)),'-blocksize')) then
          block_size = cast_int(register%default(i))
          nthreads=block_size
       else if(str_match(trim(register%key(i)),'-nthreads')) then
          nthreads = cast_int(register%default(i))
       else if(str_match(trim(register%key(i)),'-shuffle')) then
          shuffle = cast_logical(register%default(i)) 
       endif
    end do
  end subroutine parse_initialise

  subroutine parse
    call init_register
    call get_cmdLine
    call parse_help
    call update_register
    call parse_out
    call parse_plink
    call parse_trait_pos
    call parse_predict
    call parse_ndist
    call parse_initialise
  end subroutine parse

end module cmd_parser

module routinez
use parz
implicit none

contains

  subroutine get_size()
    character(len=3) :: cdum
    nind=0
    open(35,file=trim(phenfil),status='old',form='formatted')
    do 
       read(35,*,iostat=ios) cdum
       if (ios.ne.0) exit
       nind=nind+1
    enddo
    close(35,status='keep')
    nloci=0
    open(36,file=trim(bimfil),status='old',form='formatted')
    do 
       nloci=nloci+1
       read(36,*,iostat=ios) cdum
       if (ios.ne.0) exit
    enddo
    close(36,status='keep')
    nloci=nloci-1
  end subroutine get_size

  subroutine load_phenos_plink()
    character(len=1024):: str 
    character(len=20) :: why_str
    integer ::pos1, pos2,n ,i 

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
          if (n== trait_pos+5) then
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
  end subroutine load_phenos_plink

  subroutine load_snp_binary()
    integer :: i, j,k, tr
    integer(kind=1) :: b1
    integer(kind=1), dimension(0:3) :: igen 
    integer(kind=1) :: val
    igen(0)=0
    igen(1)=3
    igen(2)=1
    igen(3)=2

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

  subroutine process_snp
    integer :: j,k, nomiss
    double precision :: q,qtest,stdev,mean
    integer, dimension(ngeno) :: gc
    integer, dimension(nt) :: xtemp
    if(mcmc) then
       xpx=0.0d0
       do j=1,nloci
          xtemp=X(:,j)
          gc=0
          nomiss=count(xtemp < 3)
          q=dble(sum(xtemp,mask= xtemp < 3)) / (2.0d0*nomiss)
          mean=2.0d0*q
          stdev=dsqrt(mean*(1.0d0-q))
          if(q==1.0d0 .or. q==0.0d0) then
             geno_mean(:,j)=0.0d0
             xpx(j)=0.0d0
          else 
             do k=1,ngeno
                gc(k)=count(xtemp==(k-1))
             enddo
             geno_mean(0,j)=(0.0d0-mean)/stdev
             geno_mean(1,j)=(1.0d0-mean)/stdev
             geno_mean(2,j)=(2.d0-mean)/stdev
             geno_mean(3,j)=0.0d0
             xpx(j)=dot_product(dble(gc),geno_mean(:,j)**2)
          end if
          freqstore(j)=q
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
          mean=2.0d0*q
          stdev=dsqrt(mean*(1.0d0-q))
          if(q==1.0d0 .or. q==0.0d0) then
             geno_mean(:,j)=0.0d0
          else
             xtemp=X(:,j)
             nomiss=count(xtemp < 3)
             qtest=dble(sum(xtemp,mask= xtemp < 3)) / (2.0d0*nomiss)
             geno_mean(0,j)=(0.0d0-mean)/stdev
             geno_mean(1,j)=(1.0d0-mean)/stdev
             geno_mean(2,j)=(2.d0-mean)/stdev
             geno_mean(3,j)=(2.0d0*qtest-mean)/stdev
          end if
       end do
    endif
  end subroutine process_snp

  subroutine allocate_data
    !trains keeps track of training (0) and test(1) individuals
    !swap around for predicition 
    if(.not.mcmc) then
       where(trains==0) trains=3
       where(trains==1) trains=0
       where(trains==3) trains=1
    end if
    nt=count(trains==0)
    allocate(pred(nind), gpin(ndist), gp(ndist), p(ndist), geno_mean(0:3,nloci), &
         X(nt,nloci), delta(ndist),dirx(ndist), g(nloci), g_old(nloci), res(nt), & 
         snpindist(ndist), varindist(ndist), s(ndist), stemp(ndist), &
         xpx(nloci), gstore(nloci), snpstore(ndist), varstore(ndist), pstore(ndist), &
         indiststore(nloci,ndist), mu_vare_store(4), freqstore(nloci), permvec(nloci),&
         snptracker(nloci), log_p(ndist), log_gp(ndist),vare_gp(ndist), stat=ios)  
    if( ios /= 0 )then
       write(21,'(a50)') 'ERROR :: Unable to allocate required storage'
       stop 'Unable to allocate required storage for data'
       call flush(21)
    endif
  end subroutine allocate_data

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

  subroutine compute_dgv_rhs
    integer :: i,j, tr
    double precision :: pr

    tr=0
    pred=-9999.0d0
    do i=1,nind
       if(trains(i)==0) then
          tr=tr+1
          pr=mu
          do j=1,nloci
             pr=pr+geno_mean(X(tr,j),j)*gstore(j)
          enddo
          pred(i)=pr
       endif
    end do
  end subroutine compute_dgv_rhs

  subroutine write_dgv()
    integer :: i
    character(len=2) :: missvalue

    missvalue='NA'
    open(unit=61,file=mbvfil,status='unknown',form='formatted')
    do i=1,nind
       if (trains(i)==0) then
          !write (rchar, '(E15.7)') pred(i)
          !write(61,*)adjustl(rchar)
          write(61,'(E15.7)') pred(i)
       else
          write(61,'(a2)') missvalue
       end if
    enddo
    close(unit=61,status='keep')
  end subroutine write_dgv

  subroutine load_param()
    character(len=100):: str 
    character(len=10):: dum
    integer :: i, nc, ios
    double precision, dimension(20):: gtemp

    open(31,file=paramfil,status='old',form='formatted')
    read(31,'(a)') str
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
    read(32,*) dum, mu
    close(32)
  end subroutine load_param

  subroutine output_model
    integer :: i
    character(LEN=10) :: ci,ca

    open(unit=10,file=paramfil,status='unknown',iostat=ios)
    if(ios/=0) then
       write(21,*) 'Error opening ',trim(paramfil)
       stop
    end if
    do i=1,ndist
       write(ci,'(I8)') i
       ci=adjustl(ci)
       ca="PIP"//trim(ci)
       write(10,'(2(A5))',advance="no") ' ',ca
    end do
    write(10,*) '    beta'
    do i=1,nloci
       write(10,'(50(E15.7,1X))') indiststore(i,:), gstore(i)
    enddo
    close(10,status='keep')

    open(unit=12,file=modfil,status='unknown',iostat=ios)
    if(ios/=0) then
       write(21,*) 'Error opening ',trim(modfil)
       stop
    end if
    write(12,800) 'Mean',mu_vare_store(1)
    write(12,800) 'Nsnp',mu_vare_store(2)
    write(12,800) 'Va',mu_vare_store(3)
    write(12,800) 'Ve', mu_vare_store(4)
    do i=1,ndist
       write(ci,'(I8)') i
       ci=adjustl(ci)
       ca="Nk"//trim(adjustl(ci))
       write(12,800) ca,snpstore(i)
    end do
    do i=1,ndist
       write(ci,'(I8)') i
       ci=adjustl(ci)
       ca="Pk"//trim(adjustl(ci))
       write(12,800) ca,pstore(i)
    end do
    do i=1,ndist
       write(ci,'(I8)') i
       ci=adjustl(ci)
       ca="Vk"//trim(adjustl(ci))
       write(12,800) ca,varstore(i)
    end do
    close(12,status='keep')

    800 format(a,t10,E15.7)
  end subroutine output_model

  subroutine compute_residuals_rhs()
    integer :: tr, i,j
    tr=0
    do i=1,nind
       if(trains(i)==0) then
          tr=tr+1
          res(tr)=why(i)-mu
          do j=1,nloci
             res(tr)=res(tr)-geno_mean(X(tr,j),j)*g(j)
          end do
       endif
    enddo
  end subroutine compute_residuals_rhs

  subroutine output_snploc
    integer :: i,j
    character(LEN=20) :: ci, cg, ce
    character(len=1) :: double = ":"

    do i=1,nloci
       j= snptracker(i)
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

  subroutine block_dimensions
    rhs_ncol=ngeno**block_size
    n_blocks=int(nloci/block_size)
    rest_size=nloci-(n_blocks*block_size)
    if(rest_size>0) then 
       n_blocks=n_blocks+1
    end if
    allocate(rhs_block(rhs_ncol,block_size),xe(block_size),s_loc(block_size),ny2(0:(ngeno-1)), &
         gbl(block_size), rhs_ind(nt,n_blocks),res_old(nt), res_null(nt), delta_res(rhs_ncol), &
         rhs_count(rhs_ncol,n_blocks),delta_e(rhs_ncol),block_info(n_blocks,block_size),&
         shufflevec(n_blocks),resvec(rhs_ncol),resvec_temp(rhs_ncol),bres(ngeno,block_size), &
         countvec(rhs_ncol),countvec_temp(rhs_ncol),stat=ios)
    if( ios /= 0 )then
       stop 'Unable to allocate required storage for rhs block'
    endif
  end subroutine block_dimensions

  subroutine rhs_group
    integer :: i,j
    integer, dimension(:), allocatable :: ks,rs 

    allocate(ks(block_size),rs(block_size-1))
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

  subroutine build_blocks
    integer :: i,bs,start,end

    do i=1,nloci 
       permvec(i)=i
    enddo
    if(permute) then
       call permutate(permvec,nloci)
    endif
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
       block_info(i,1:bs)=permvec(start:end)
    enddo
  end subroutine build_blocks

  subroutine mapfun(bsize)
    integer :: i,n
    integer, intent(in) :: bsize 
    integer, dimension(:), allocatable :: temp

    if(allocated(lookup)) deallocate(lookup)
    if(allocated(map)) deallocate(map)
    allocate(map(bsize),temp(bsize))

    temp(:)=3
    map(1)=1
    do i=2,bsize
       map(i)=ngeno**i
    enddo
    n=sum(temp*map)+1
    allocate(lookup(n))
    do i=1,rhs_ncol
       n=sum(rhs_block(i,:)*map)+1
       lookup(n) = i
    enddo
  end subroutine mapfun

  subroutine rhs_ct
    integer:: i,j,bs, pos
    integer, dimension(9) :: sindex

    call mapfun(block_size)
    do i=1,n_blocks
       if(i==n_blocks .and. rest_size>0) then
          bs=rest_size
       else
          bs=block_size
       end if
       sindex(1:bs)=block_info(i,1:bs)
       countvec=0
       !$OMP parallel private(j,pos,countvec_temp) if(nt >5000)
       countvec_temp=0
       !$OMP DO
       do j=1,nt
          pos=lookup(sum(X(j,sindex(1:bs))*map(1:bs))+1)
          !rhs_count(pos,i)=rhs_count(pos,i)+1
          countvec_temp(pos)=countvec_temp(pos)+1
          rhs_ind(j,i)=pos
       enddo
       !$OMP end do
       !$OMP critical
       countvec=countvec+countvec_temp
       !$OMP end critical
       !$OMP barrier
       rhs_count(:,i)=countvec
       !$OMP end parallel
    end do
  end subroutine rhs_ct

  subroutine xe_update
     integer ::i, loc, s1, s2

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
         res(i)=res(i)-delta_res(rhs_ind(i,bl))
      end do
      !$OMP end parallel do 
   end if
 end subroutine res_update_rhs

end module routinez
