! 	print*, '>>>>>>>>>>>>>>>>          CRYOFF           <<<<<<<<<<<<<<'
! 	print*, '>>>>>>>>>>>>>> CReate Your Own Force Field <<<<<<<<<<<<<<'
! 	print*, '>>>>>>>>>>>>>>>>    Group of Feng Wang     <<<<<<<<<<<<<<'
! 	print*, '>>>>>>>>>>>>>>>>        Jicun    Li        <<<<<<<<<<<<<<'
! 	print*, '>>>>>>>>>>>>>>>>        Zhonghua Ma        <<<<<<<<<<<<<<'
! 	print*, '>>>>>>>>>>>>>>>>        Ying     Yuan      <<<<<<<<<<<<<<'
! 	print*, '>>>>>>>>>>>>>>>>        Feng     Wang      <<<<<<<<<<<<<<'
! 	print*, '>>>>>>>>>>>>>>>>    2022S-01-21 13:14:00   <<<<<<<<<<<<<<'
! 	print*, '>> Usage:  CRYOFF File.ff'
! 	print*
!
!    Notes:
!    1. For Ewald, only cubic box is supported (Ewald is not working. In the process of fixing it)
!    2. DE part (should be working but not fully tested )
!
! 2022-11-29: Support for hybrid fit where atomic forces and net forces and torques are fitted together
! 2021-12-24: New printing of error bar of parameters in case of KCV (Ying)
! 2021-11-23: Implementation of the EQV (equivalence keyword.) (Ying)
! 2021-11-18: Improved support to cross validation (Ying)
! 2021-11-16: add charge matrix decomposition through spectral analysis. (CMD) (Ying)
! 2021-09-17: added coupled dihedral support (Ying) 
! 2021-09-17: reactivated BND and BD3 (bonded 3 body)  (Ying/Wang) 
! 2021-08-12: reduced memory footprint for MPI jobs. 
! 2021-08-07: add support for harmonic restraining non-linear parameters.
! 2021-08-04: add support for equivalent molecules. 
! 2019-09-04: implment support for cross validation value calculation 
! 2019-08-19: a bug was fixed when there are more than 1000 frames in the ref file. 
! 2019-05-09: merged TT damping and cross terms
! 2019-05-05: added Thole damping for Coulomb
! 2018-11-23: merged new function type
! 2018-10-30  built-in NetF and Torq calculation
! 2018-07-27: add atom type for performance
! 2017-11-30: New Ewald code
! 2017-11-19: dErr return Chisq
! 2017-02-27: PBC for VDW
! 2014-09-05: Output for DrvCRYOFF should use weighted RMSE
! 2014-02-12: Add piecewise spline function, order 0 or 1
! 2013-11-12: MPI
! 2013-06-28: Serial
! 2013-06-28: Output including original forces
! 2013-06-12: Testing Pass
! 2013-06-03: Add exclusion part
! 2013-06-01: Add error for parameters, R^2 for parameters
! 2013-05-29: Test All potential type
! 2013-05-27: Add Dihedral term
! 2013-05-24: Force and Matrix in ONE subroutine for Bond and Angle interaction
! 2013-05-21: Force and Matrix in ONE subroutine for Pair interaction
! 2013-05-21: A bug for damping Buckingham force fixed

! steps for the calc
! 1. readin the top file so we can know the FF parameters
! 2. readin the ref file so we can know the Reference forces, torques
! 3. build the SVD matrix
! 4. solve the SVD
! 5. print the results

#define noMPI


module extrafuncs
contains

!integer function findloc(array,value,dim)
! integer, dimension(:) :: array
! integer :: value, dim
!   
! findloc=minloc((array-value)**2,dim=1)
!End function findloc

integer function findatm(size, TagIatm, array)
implicit none
integer :: size
integer :: i
logical :: lfound
character*20 :: TagIatm
character*20, dimension (:) :: array
    lfound=.false.
    do i=1,size
      if (TagIatm .eq. array(i)) then
          lfound=.true.
          exit
      endif
    enddo 
    if (lfound==.false.) then 
        call ErrStop('Atom Name Not Found '//trim(TagIatm), 1)
    endif
    findatm=i
     
End function findatm 

end module extrafuncs
    

Module CRYOFFmod
#ifdef MPI
        include "mpif.h"
#endif

    integer, parameter:: &
        & Lout=2,      & ! output file unit
        & Lbug=5,      & ! output file unit
        & Ltop=3,      & ! input .ff file unit
        & Lref=4,      & ! input .ref file unit
        & MaxPara=10,  & ! maximum number of parameters for a potential type
        & MaxChg=32,   & ! max terms  for charge constraints
        & MaxAdj=5,    & ! max number of adjacient atoms for intramolecular potential
        & MaxLev=10,    & ! max level for .ff fi
        & Maxpptp=10,    &! Maximum number of Pair Interaction Per Type-Pair
        & MaxLab=25      ! Maximum Length for label Lab
    
    real*8, parameter:: &
        & Pi=4.d0*atan(1.d0), TwoPi=2.d0*Pi, &
        & Deg2Rad=Pi/180.d0, Rad2Deg=180.d0/Pi, &
        & WgtTol=1.d-10,      &
        & Reps=epsilon(1.d0), &
        & FcovCou = 332.063713741257D0, &
        & bdwarn = 2.5 !warning if bond lenght exceeds this. 
    
    character*1, parameter:: CR=achar(13), LF=achar(10)

    type:: Mol
    integer  Natm, Ntyp, Nvir, NmolTot ! No. of atom, exclusion, pot, vsite, molecules_of_this_type
    integer  EquivID
    integer, pointer:: Ityp(:), &  ! pot type
        & Npar(:), &              ! no. params
        & Nlnr(:), Ilnr(:), &     ! no. and idx of linear params
        & Ngrp(:), Iatm(:,:,:), & ! no. groups and atom idx
        & IatmVir(:), NatmVir(:), Ivir(:,:) ! vsite

    	real*8   FudgeQQ, FudgeVDW
    	real*8,  pointer:: Para(:,:), Rvir(:,:),ParaRMS(:,:)  
    	logical, pointer:: YesFit(:), YesExc(:,:), YesPair(:,:)
    character*20  Name
    character*20, pointer:: Ttyp(:), Satm(:), VDWatm(:), COUatm(:)
    end type Mol
    type(Mol), allocatable, target:: Mole(:)

    type Par
	real*8, pointer:: P
    end type Par
    type(Par), allocatable:: Popt(:)

    integer IsnpTot, NsnpTot, NmolTyp, Npair, Ncst, MaxAtm, MaxVir, &
        & NlnrTot, NrefTot, NatmRefTot, Mdim, NrefAtmSnp, &
        & Npar, Nopt, Iwgt, NatmTyp, NcouTyp, NvdwTyp, NvirTot

    ! NvirTot Number of virutal sites in the reffile.

    integer, allocatable:: &
        & NatmSnpTot(:), IbgnAtmSnp(:), &
        & NrefSnpTot(:), NmolRefTot(:), &
        & ImolRefTot(:), IbgnMolRef(:), IbgnAtmRef(:), &
        & IatmAtmTot(:), IatmRefTot(:), &
        & NatmRefSnp(:), IatmRefSnp(:), &
        & NparPair(:), NlnrPair(:), IlnrPair(:), IatmRef(:), &
        & IatmAtm(:), NchgGrp(:), Nchg(:,:), Ichg(:,:), &
        & IcouPair(:,:), IvdwPair(:,:), Ivdw(:), Icou(:), &
        & IcouTot(:), IvdwTot(:), Iuniq(:), IkcvFrm(:), &
        & npptp(:,:,:), cnpptp(:,:,:)

    real*8, allocatable:: XatmTot(:), YatmTot(:), ZatmTot(:), &
        & FatmTot(:), FwgtTot(:), &
        & Xatm(:), Yatm(:), Zatm(:), &
        & Fatm(:), Fwgt(:), Fchg(:), Wchg(:), Ffix(:), FfixTot(:), &
        & Ai(:,:), Aj(:,:), Ak(:,:), Al(:,:), Ar(:,:), &
        & Amat(:,:), AmatTot(:,:), Fmat(:,:), FmatTot(:,:), Asav(:,:),  Fsav(:,:), Anet(:,:), Ator(:,:), &
        & Pmat(:), ParaPot(:), Asig(:), Chisq(:,:), &
        & Pinf(:), Psup(:), RminPair(:), RmaxPair(:), &
        & PRC(:), PRP(:)     !Center and Penalty for the Non-linear Parameter Harmonic Restraint. 
    
    logical, allocatable :: PLR(:)       !Logical for Harmonic Restraint
    
    integer Nff ! Number of lines in txtFF. (noncomment lines in ff file).
    real*8 :: wgtfac, deps, optstp
    integer :: IniSnp, LstSnp, MaxSnp  ! Initial and Last snapshot to read. MaxSnp=LstSnp-IniSnp
    integer :: Iopt, Pgen, MaxIter, Mopt
    integer :: NatmTot
    integer :: ItypKCV !CV type 0 means sequential 1 means random
    integer :: ipnt, Mfit, Irank, Mcst, MnetF, Mtorq, Matm
    character*800 text
    real*8, allocatable,target:: ParaPair(:,:)

    character*20, allocatable:: NameTot(:), SatmTot(:), &
        & Name(:), Satm(:),  &
        & TtypPair(:), COUatmTyp(:), VDWatmTyp(:)
    character*MaxLab, allocatable:: LabAtm(:), LabAtmTot(:), TatmPair(:)
    character*80 Fout
    character*1024, allocatable:: txtFF(:)   !initial guess array for nonlinear optimzation is large

    character*256 :: AllPairType
    data AllPairType /'COU THC GLJ BUC FDB STR EXP PEX GEX POW TTP SRD FDP'/
    logical:: YesBug=.false., YesOut(3)=.false., YesInt=.false., YesPBC=.false., YesHyb=.false.
    logical, allocatable:: YesNetF(:), YesTorq(:), YesFitPair(:)
    logical, allocatable :: YesExc(:, :)  !the working (preframe) YesExc only useful for ewald
    
    real*8 :: WssqFatm, RmsFatm,  WssqNetF,  RmsNetF, WssqTorq, RmsTorq
    real*8 :: WRmsFatm, WRmsNetF, WRmsTorq
!    for printing parameter index    EQV
    logical :: YesPotIdx=.false., YesEQV=.false.
    
!   for Charge product Matrix Decomposition (CMD)
    logical :: YesCMD=.false.
    integer CMDchgNum  !numer of CMD charge atom types
    character*800 CMDatmlist
    character*MaxLab, allocatable:: TatmIpairQQ(:), CMDatm(:) !CMD charge pair and CMD charge atom types
    
!   for Ewald
    logical :: YesCub, YesQQ, YesEwa=.false.

    integer Imax, Jmax, Kmax, NpairQQ, Nqq
    real*8  Etol, Rcou, Rvdw, Vol, beta , Rcon, &
        & Avec(3), Bvec(3), Cvec(3), Abox, Bbox, Cbox, &
        & Uvec(3), Vvec(3), Wvec(3), Ubox, Vbox, Wbox

    integer, allocatable:: ImaxTot(:), JmaxTot(:), KmaxTot(:), IpairQQ(:)
    ! RcouTot removed. It does not really make sense to use different rcou for each frame even considering volume changes
    real*8,  allocatable::  VolTot(:),  &
        & AboxTot(:),   BboxTot(:),   CboxTot(:),   &
        & AvecTot(:,:), BvecTot(:,:), CvecTot(:,:), &
        & UvecTot(:,:), VvecTot(:,:), WvecTot(:,:) , &
        & QQsin(:,:,:,:), QQcos(:,:,:,:)
    character*20, allocatable:: tagQQ(:)
    
!   for Simplex
    real*8, allocatable:: P0(:), dP(:), P(:,:), Xi(:,:),  Y(:)
    
!   for CV
    integer Nkcv, Iseed
    integer, allocatable:: IndxKCV(:), IkcvTot(:)
    real*8, allocatable:: PmatKcv(:,:)
    real*8, allocatable :: ParaRMS(:,:)

!   for MPI
    integer comm_world, comm_pop, comm_force, &
        & np_world, np_pop, np_force, &
        & MyID, MyID_pop, MyID_force
    integer, allocatable:: IdspRef(:), NcntRef(:)

! 	for DE
    integer dei(6)
    real*8 der(3)
!	integer Nrest
!	logical YesRest
End Module CRYOFFmod
!
!
Program CRYOFF
        use CRYOFFmod

        use extrafuncs
        
        implicit none

        integer i, j, k, ii, Ierr, &
		& Iref, Narg, Natm, Iatm, &
		& NatmSnp, Itmp, Ntyp, &
		& Iter, Iu(6), &
		& Level
        
        integer*8 nmem
        
	integer, allocatable:: Ifrm(:)
	real*8, allocatable:: X(:)
    
	real*8 Rtmp, Wgt, Yret, dErr, &
		&  XYZi(3), Fi(3), CPUini, CPUend

	real*8, external:: getValue, getBetaKmax, getEdir, getErec

	character*80  tag, tag1, tag2, TagIatm
	character*80  Ftop, Fref 
	character*1024 txt


	external Level

    
#ifdef MPI
	integer Irow, Icol
#endif
	integer  method(3)
	real*8   bestval

! ------------------ set up MPI --------------------
	MyID=0
	MyID_pop=0
	MyID_force=0
	np_world=1
#ifdef MPI
	call MPI_Init(Ierr)                            ! starts MPI
	comm_world= MPI_COMM_WORLD
	call MPI_Comm_size(comm_world, np_world, Ierr) ! get number of procs
	call MPI_Comm_rank(comm_world, MyID, Ierr)     ! get current proc ID
#endif

! ------------------ obtain input file name --------------------

	if(MyID==0) then
		Narg = Iargc()
		if(Narg>0) then
			call getarg(1, Ftop)
			i = index(Ftop, '.ff', .true.)
			if(i/=0) Ftop = Ftop(1:i-1)
		else
! 		get input file name
		print*, '>> Type Input File Name (*.ff/Exit):'
		read(*, *) Ftop
		Ftop=trim(adjustl(Ftop))
		txt=Ftop; call upcase(txt)
		if(index(txt, 'EXIT')==1) stop '>> Normal Exit of CRYOFF.'

! 		try to open input files
		open(unit=Ltop, file=trim(adjustl(Ftop))//'.ff', status='old', IOstat=Ierr)
		do while(Ierr/=0)
			print*, '>> Input File  '//trim(adjustl(Ftop))//'.ff  does NOT  exist!'
            close(Ltop)
            print *, '>> Please enter input file name or type exit:'
			read(*, *) Ftop
!			Ftop=trim(adjustl(Ftop))
			txt=trim(adjustl(Ftop)); call upcase(txt)
			if(index(txt, 'EXIT')==1) stop '>> Normal Exit of CRYOFF.'
			open(unit=Ltop, file=trim(adjustl(Ftop))//'.ff', status='old', IOstat=Ierr)
		end do
		close(Ltop)
        endif
    endif
        
! ------------------ Initial Output --------------------
    if(MyID==0) then
 	print*, '>>>>>>>>>>>>>>>>          CRYOFF           <<<<<<<<<<<<<<'
 	print*, '>>>>>>>>>>>>>> CReate Your Own Force Field <<<<<<<<<<<<<<'
 	print*, '>>>>>>>>>>>>>>>>    Group of Feng Wang     <<<<<<<<<<<<<<'
 	print*, '>>>>>>>>>>>>>>>>        Jicun    Li        <<<<<<<<<<<<<<'
 	print*, '>>>>>>>>>>>>>>>>        Zhonghua Ma        <<<<<<<<<<<<<<'
 	print*, '>>>>>>>>>>>>>>>>        Ying     Yuan      <<<<<<<<<<<<<<'
 	print*, '>>>>>>>>>>>>>>>>        Feng     Wang      <<<<<<<<<<<<<<'
 	print*, '>>>>>>>>>>>>>>>>       Version 3.0.0         <<<<<<<<<<<<<<'
    print*, LF,LF
		if(np_world>1) print*, '>> Running in MPI Mode over ',np_world,' CPUs'
		if(np_world==1) print*, '>> Running in Serial Mode'
		print*, '   Name of Input File:            '//trim(adjustl(Ftop))//'.ff'
!       This is the maximum number of adjacient atoms for intramolecular potential. 
!		print*, '   Maximum number of Atoms:        ', MaxAdj
		print*, '   Maximum number of parameters for each potential:    ', MaxPara
		print*, '   Maximum number of terms for Charge Constraints: ', MaxChg, LF

! ------------------ Read Input File to txtFF and determine ref file and output file name--------------------
! determine Nff (number of noncomment lines) Nopt (number of non-linear parameters) 
! Read input file into buffer txtFF
        
        call readff(Ftop)
        
        Nopt=0 ! number of nonlinear parameters to optimize
        NmolTyp=0

! Read input file buffer first time to get NmolTyp, Nopt, and filename for .ref and .off
        
        do i=1, Nff
		read(txtFF(i), '(A)', IOstat=Ierr) txt

		if(index(txt, ' FIX')==0) then
        k=Len_trim(txt)
				do j=1, k
					if(txt(j:j+1)=='_[') Nopt = Nopt+1
				end do
        end if
        
        tag=txt; call upcase(tag)
		if(tag(1:3)=='MOL') NmolTyp=NmolTyp+1
		if(tag(1:3)=='FIL') then
				Fref=''; Fout=''
				read(txt, *) tag, tag1, tag2
                tag=tag1; call upcase(tag)
				if(index(tag, '.REF', .true.)/=0) Fref=tag1
				if(index(tag, '.OFF', .true.)/=0) Fout=tag1
				tag=tag2; call upcase(tag)
				if(index(tag, '.REF', .true.)/=0) Fref=tag2
				if(index(tag, '.OFF', .true.)/=0) Fout=tag2
				if(len_trim(Fref)==0) read(txt, *) tag, Fref
				if(len_trim(Fout)==0) read(txt, *) tag, tag, Fout
                
                tag=Fref; call upcase(tag)
				j = index(tag, '.REF', .true.)
				if(j/=0) Fref = Fref(1:j-1)
                
                tag=Fout; call upcase(tag)
				j = index(tag, '.OFF', .true.)
				if(j/=0) Fout = Fout(1:j-1)
        end if 
        enddo

        
		print*, '   Name of Ref   File:           '//trim(Fref)//'.ref'
		print*, '   Name of Output File:           '//trim(Fout)//'.off'
		print*, '   Number of Molecular Types:         ', NmolTyp, LF
    end if   !MyID=0
    
#ifdef MPI
! Transmit input file buffer to other nodes and allocate arrays in other nodes.
	call MPI_Bcast(Nff,     1, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(Nopt,    1, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(NmolTyp, 1, MPI_INTEGER4, 0, comm_world, Ierr)

	if(MyID/=0) then
		allocate(txtFF(Nff), stat=Ierr)		!For MYID==0, txtFF allocation occured in readff subroutine
		call ErrStop('Memory Allocation Failure Popt', Ierr)
	end if

	call MPI_Bcast(txtFF, 1024*Nff, MPI_CHARACTER, 0, comm_world, Ierr)
#endif
        		
    allocate(PLR(Nopt), PRC(Nopt), PRP(Nopt), stat=Ierr)
    allocate(Popt(Nopt), P0(Nopt), dP(Nopt), Pinf(Nopt), Psup(Nopt), &
		&	Mole(NmolTyp), YesNetF(NmolTyp), YesTorq(NmolTyp), stat=Ierr)
    call ErrStop('Memory Allocation Failure Popt', Ierr)

! Read input file buffer again to set system options. 

	if(MyID==0)	print*, '>> Parsing System Options'
    call SetSystemOptions
    allocate (Chisq(0:4,Nkcv),stat=Ierr)
    call ErrStop('chisq allocation error', Ierr)

!#ifdef MPI

	np_force=np_world/np_pop
    
    allocate(IdspRef(np_force), NcntRef(np_force), stat=Ierr)
    call ErrStop("MPI NcntRef allocation failure.", Ierr)
    IdspRef(1)=0

#ifdef MPI
    
	irow=MyID/np_force
	icol=mod(MyID,np_force)
    
    if (np_pop > 1) then
	call MPI_Comm_split(comm_world, icol, irow, comm_pop,Ierr)
	call MPI_Comm_split(comm_world, irow, icol, comm_force,Ierr)
	call MPI_Comm_rank(comm_pop, MyID_pop, Ierr)
    call MPI_Comm_rank(comm_force, MyID_force, Ierr)
    else
        comm_force=comm_world
        MyID_force=MyID
    endif 
    
	
#endif
! Should print out the final interpretation of each input read. 


! Third pass of reading input file. 
    call ReadMolecularDef()
    
    ! Start reading Intermolecular interactions.     

    ! I am against the philosophy have to indpendent vdw atom type and coulombic atom type. 
    ! The model would then be confusing. 
	allocate(COUatmTyp(NatmTyp), VDWatmTyp(NatmTyp))

	! Probably a bad habit to recycle variable like this. However, the code did check out to be correct. 
    NatmTyp=0
	do i=1, NmolTyp
		do j=1, Mole(i)%Natm
			NatmTyp=NatmTyp+1
			COUatmTyp(NatmTyp)=Mole(i)%COUatm(j)
			VDWatmTyp(NatmTyp)=Mole(i)%VDWatm(j)
		end do
	end do

	call uniq(NatmTyp, NcouTyp, COUatmTyp)
	call uniq(NatmTyp, NvdwTyp, VDWatmTyp)
    ! NcouTyp and NvdwTyp could be less than NatmTyp
    if (MyID==0) then
    print*, "A total of ,",NvdwTyp," atom types read in for short range interactions."
    write(*,'(10A6)')  VDWatmTyp(1:NvdwTyp)
    print*, "A total of ,",NcouTyp," atom types read in for partial charges."
    write(*,'(10A6)') COUatmTyp(1:NcouTyp)
    print *
    endif
    
    allocate (npptp(NvdwTyp,NvdwTyp,Maxpptp)) !npair per type-pair for speeding up vdw loop in InteFrcMat
    allocate (cnpptp(NcouTyp,NcouTyp,Maxpptp))
    npptp(:,:,:)=-1
    cnpptp(:,:,:)=-1
    
    call ReadIntermolecularTerms()

    call ReadChargeConstraint()
    
    ! Done reading Input files. 
    
    ! Set up Simplex
	allocate(X(Nopt), stat=Ierr)
	call ErrStop('Memory Allocation Error X(Nopt)', Ierr)

	do i=1, Nopt
		P0(i)=Popt(i)%P
		X(i)=P0(i)
	end do

	if(Iopt==1) then
		select case(Pgen)
        case(0) ! initial simplex generated using read data
            if (MyID==0) print *, "Initial simplex set up using read data"
!            P(Mopt, 1:Nopt) = P0(1:Nopt)
        case(1:2) ! initial simplex generated
			P(1, 1:Nopt) = P0(1:Nopt)            
			do i=1, Nopt
				P(i+1, 1:Nopt) = P0(1:Nopt)
				P(i+1, i) = P(i+1, i)+dP(i)
			end do
        end select
	else if(Iopt==2) then ! initial directions generated using unit vectors
		Xi(:, :) = 0.d0
        select case (Pgen)
        case(0)
                if (MyID==0) print *, "Initial Powell set up using read data"       
        case(1:2)
		do i=1, Nopt
			Xi(i, i) = dP(i)
        end do
        end select
	else if(Iopt==3) then
		method=(/dei(4),dei(5),dei(6)/)
	end if

#ifdef MPI
	if(Iopt==1) call MPI_Bcast(P,  Nopt*Mopt, MPI_REAL8, 0, comm_world, Ierr)
	if(Iopt==2) call MPI_Bcast(Xi, Nopt*Nopt, MPI_REAL8, 0, comm_world, Ierr)
#endif
    ! Finish setting up simplex

!---------------write .off  file------------
	if(MyID==0) then
		open(unit=Lout, file=trim(adjustl(Fout))//'.off', status='new', IOstat=Ierr)
        call ErrStop("will not overwrite existing .off file",Ierr) 
        
        write(Lout, '(A)') 'Force Field/Parameter File: '//trim(Ftop)//'.ff'
        
        if (size(txtFF)/=Nff) call ErrStop("Internal Consistency Error with txtFF", Nff)
        
		do i=1, Nff
			read(txtFF(i), *) txt
			j=Level('['//trim(txt)//']')
			if(j<3) then
				ii=index(txtFF(i), ' ')
                ! put back []. no longer modifying txtFF at this stage.
				txt='[ '//trim(adjustl(txtFF(i)(:ii)))//' ]'//trim(txtFF(i)(ii+1:))
            else
                txt=trim(txtFF(i))
            end if
            ! format alignment
			write(tag, '(I0, A)') 2*min(j, 4), 'X'
			write(Lout, '('//trim(tag)//', A)') trim(txt)
		end do

		write(Lout, '(/,A)') 'Atom Types:'
		write(Lout, '(A)') '  idx  COU                 VDW'
		do i=1, max(NcouTyp, NvdwTyp)
			write(txt, '(2X, I0)') i
			write(Lout, '(A8, A, A)') txt, COUatmTyp(i), trim(VDWatmTyp(i))
        end do

        if(YesPotIdx) then				!EQV
            write(Lout, '(/,A)') 'Job Type: List fitting potential index.'//LF
		elseif(Nopt==0) then
			write(Lout, '(/,A)') 'Job Type: Linear Parameter Optimization Only'//LF
		else
			write(Lout, '(/,A)') 'Job Type: Non-Linear Optimization'//LF
			if (Iopt==1) write(Lout, '(A)') '  Optimization Method:    simplex'
			if (Iopt==2) write(Lout, '(A)') '  Optimization Method:    Powell'
			if (Iopt==3) write(Lout, '(A)') '  Optimization Method:    Differential Evolution'

			write(Lout, '(A, ES14.6)')      '  |Ymax-Ymin| Criterion: ', deps
			write(Lout, '(A, I11, A)')      '  Max No. of Iterations: ', MaxIter, LF
			write(Lout, '(A)') '  Initial Parameter Sets:'
            write(Lout, *) '   Parameter    InitialValue  LowerBnd/Rest.  UpperBnd/Rest.Penalty   Stepsize'
         	do i=1, Nopt
                if (PLR(i)) then
                    write(Lout, '(4X, A2, I2, A, 4F14.6)')  'P(', i, ')    ', P0(i), PRC(i), PRP(i), dP(i)
                else
                   write(Lout, '(4X, A2, I2, A, 4F14.6)')  'P(', i, ')    ', P0(i), Pinf(i), Psup(i), dP(i)
                endif
            enddo
            
! Text variable recycled 
            text = '    V\P'
			write(txt, '(I10)') 1
			text = trim(text)//trim(txt)
			do i=2, Nopt
				write(txt, '(I20)') i
				text = trim(text)//trim(txt)
			end do

			write(Lout, *) LF
			if(Iopt==1) then
				write(Lout, '(A, \)') '  The vertices of initial simplex are generated using: '
				if(Pgen==0) write(Lout, *) 'read data.'
				if(Pgen==1) write(Lout, *) 'uniform vectors. (will be overwritten by ff file specifications if present)'
				if(Pgen==2) write(Lout, *) 'factional scaling. (will be overwritten by ff file specifications if present)'
                if (YesOut(3)) then
				write(Lout, '(T32, A)') 'Vertices of Initial Simplex'
				write(Lout, '(A)') trim(text)
				do i=1, Mopt
					write(Lout, '(2X, I3, 600F20.9)') i, P(i, 1:Nopt)
                end do
                endif !YesOut(3)
			else if(Iopt==2) then
				write(Lout, *) ' The vectors of initial directions are generated using '
				if(Pgen==0) write(Lout, *) 'read data.'
				if(Pgen==1) write(Lout, *) 'uniform vectors. (will be overwritten by ff file specifications if present)'
				if(Pgen==2) write(Lout, *) 'fractional scaling. (will be overwritten by ff file specifications if present)'
                if (YesOut(3)) then
				write(Lout, '(T32, A)') 'Vertices of the initial simplex'
				write(Lout, '(A)') trim(text)
				do i=1, Mopt
					write(Lout, '(2X, I3, 600F20.9)') i, P(i, 1:Nopt)
                end do
                endif !YesOut(3)
			else if(Iopt==3) then
				write(Lout, *) ' The initial parameters:'
				write(Lout, *) ' NP=', dei(1), ' strategy=',dei(2), &
					&        ' refresh=',dei(3),' method=(/',dei(4),',',&
					&        dei(5),',',dei(6),'/)'
				write(Lout,*) ' CR_XC=',der(1),' F_XC=',der(2), 'F_CR=', der(3)
			end if
        end if
                write(Lout, *) LF        
    end if

    ! determine number of linear paramters with the following two loops over intra and inter
	NlnrTot = 0
	do i=1, NmolTyp
		Ntyp = Mole(i)%Ntyp
		do j=1, Ntyp
			if(Mole(i)%YesFit(j)) then
                if(Mole(i)%EquivID == i) then
				Mole(i)%Ilnr(j) = NlnrTot+1       
				NlnrTot = NlnrTot+Mole(i)%Nlnr(j)
                else
                Mole(i)%Ilnr(j) =Mole(Mole(i)%EquivID)%Ilnr(j)
                endif
                    
			end if
		end do
	end do

	do i=1, Npair
		if(YesFitPair(i)) then
			IlnrPair(i) = NlnrTot+1
			NlnrTot = NlnrTot+NlnrPair(i)
		end if
    end do   

! Take care of EQV:
    call ReadSetEQVPotential() 
    
    if(YesPotIdx .or. YesEQV) then
        call printPotIndex
        if(YesPotIdx) stop
    endif
    
	if(MyID==0) then
        write(Lout, *), 'Condition number for SVD set to ', Rcon
		print*, '>> Optimization Setting'
		print*, '   No. of Linear Parameters:     ', NlnrTot
        write(Lout, *), '   No. of Linear Parameters:     ', NlnrTot
		print*, '   No. of Non-Linear Parameters: ', Nopt, LF
        write(Lout, *), '   No. of Non-Linear Parameters: ', Nopt, LF
	end if

	allocate( Pmat(NlnrTot), PmatKcv(NlnrTot,NKcv), ParaPot(NlnrTot), Asig(NlnrTot),  &
		& Ai(3, MaxPara), Aj(3, MaxPara), Ak(3, MaxPara), Al(3, MaxPara), Ar(3, MaxPara), &
		& Anet(3, NlnrTot), Ator(3, NlnrTot), stat=Ierr)
	call ErrStop('Memory Allocation Failure Pmat block', Ierr)

    call cpu_time(CPUini)
	if(MyID==0) then ! call ReadRefFile first time to get dimenions needed for allocating memory. 
        call ReadRefFile(Fref,1)
    endif 

#ifdef MPI
	call MPI_Bcast(NsnpTot,    1, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(NatmTot,    1, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(MaxAtm,     1, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(NrefTot,    1, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(NatmRefTot, 1, MPI_INTEGER4, 0, comm_world, Ierr)
#endif

	MaxSnp=NsnpTot
	Mdim = 3*NatmRefTot+Ncst
    i = max(Mdim,Npair*np_force)

	allocate( IatmRef(MaxAtm), IatmAtm(MaxAtm),  &
			& Name(MaxAtm), Satm(MaxAtm), Icou(MaxAtm), Ivdw(MaxAtm), &
			& Xatm(MaxAtm), Yatm(MaxAtm), Zatm(MaxAtm), &
			& Fatm(Mdim), Fwgt(Mdim), LabAtm(Mdim), Ffix(Mdim),  &
            & Fmat(Mdim, 1),  Fsav(i, 1),  &
			& NatmSnpTot(NsnpTot), & ! Number of atoms in the snapshot
			& IbgnAtmSnp(NsnpTot), & ! Index for the first atom in the snapshot
			& NmolRefTot(NsnpTot), & ! Number of reference molecule in the snapshot
			& IbgnMolRef(NsnpTot), & ! Index for the first reference molecule in the snapshot
			& NatmRefSnp(NsnpTot), & ! Number of reference atoms in the snapshot
			& IatmRefSnp(NsnpTot), &
			& ImolRefTot(NrefTot), & ! Molecular type for each molecule
			& IbgnAtmRef(NrefTot), & ! The begin atom index (within each frame) of each reference molecule
			& IatmRefTot(NatmTot), & ! The reference atom index(global) of each atom (global) (0 means not reference atom) 
			& IatmAtmTot(NatmTot), & ! For each atom, the index of this atom in the molecule 
			& NameTot(NatmTot),   SatmTot(NatmTot), IcouTot(NatmTot), IvdwTot(NatmTot), &
			& XatmTot(NatmTot),   YatmTot(NatmTot),   ZatmTot(NatmTot), &
			& FatmTot(Mdim), IndxKCV(Mdim), IkcvTot(Mdim),    FwgtTot(Mdim),      LabAtmTot(Mdim),  &
			& RminPair(Npair), RmaxPair(Npair), stat=Ierr)
	call ErrStop('Memory allocation error when allocating: Fmat!', Ierr)
    
         
            if(MyID==0) then
                allocate (Amat(Mdim, NlnrTot), stat=Ierr)
                allocate (AmatTot(Mdim, NlnrTot), FmatTot(Mdim, 1), FfixTot(Mdim), Asav(Mdim, NlnrTot), stat=Ierr)
                call ErrStop('Memory allocation error when allocating: AmatTot! at node 0', Ierr)
                
                nmem=MaxAtm*10+NatmTot*10+Mdim*10+Npair*10
                nmem=nmem+Mdim*NlnrTot*3
                nmem=nmem*8
                nmem=nmem+Nff*1024+Mdim*3*MaxLab
			                
                print *, LF
                print *, "Estimated memory requirement for the head process in bytes: ", nmem
                print *, "For MPI jobs, all the slaves combined needs about the same amount of the memory as the head."
                print *, LF
            else
                ii=3*(NsnpTot/np_force+1)*MaxAtm+ncst
!               This may not be optimal but should be good enough. 
                
                allocate (Amat(ii, NlnrTot), stat=Ierr)
                call ErrStop('Memory allocation error when allocating: Amat! in slave', Ierr)
            endif

        if (YesPBC) then 
        allocate( AboxTot(NsnpTot),   BboxTot(NsnpTot),   CboxTot(NsnpTot), &
            & YesExc(MaxAtm, MaxAtm), VolTot(NsnpTot),    &
			& ImaxTot(NsnpTot),   JmaxTot(NsnpTot),   KmaxTot(NsnpTot), &
			& AvecTot(NsnpTot,3), BvecTot(NsnpTot,3), CvecTot(NsnpTot,3), &
			& UvecTot(NsnpTot,3), VvecTot(NsnpTot,3), WvecTot(NsnpTot,3), stat=Ierr) 
        call ErrStop('Memory allocation error when allocating Ewald arrays', Ierr)
        endif
                    
	if(MyID==0) then
		allocate(IkcvFrm(NsnpTot), stat=Ierr)
		do i=1, NsnpTot
			IkcvFrm(i)=mod(i-1, Nkcv)+1
		end do
		if(ItypKCV==1) call Shuffle(NsnpTot, IkcvFrm, Iseed)
! 		do i=1, NsnpTot
! 			write(999, *) i, IkcvFrm(i)
! 		end do
        
        ! Read Lref file a second time. 
        call ReadRefFile(Fref,2)
        call cpu_time(CPUend)
        print *, "Walltime for reading the ref file:", CPUend-CPUini, "seconds."
	end if ! MyID 0

#ifdef MPI

	call MPI_Bcast(NatmSnpTot, NsnpTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(IbgnAtmSnp, NsnpTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(NmolRefTot, NsnpTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(IbgnMolRef, NsnpTot, MPI_INTEGER4, 0, comm_world, Ierr)

	call MPI_Bcast(NatmRefSnp, NsnpTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(IatmRefSnp, NsnpTot, MPI_INTEGER4, 0, comm_world, Ierr)

	call MPI_Bcast(ImolRefTot, NrefTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(IbgnAtmRef, NrefTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(IatmRefTot, NatmTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(IatmAtmTot, NatmTot, MPI_INTEGER4, 0, comm_world, Ierr)

	call MPI_Bcast(NameTot, 20*NatmTot, MPI_CHAR, 0, comm_world, Ierr)
	call MPI_Bcast(Satmtot, 20*NatmTot, MPI_CHAR, 0, comm_world, Ierr)
	call MPI_Bcast(IvdwTot, NatmTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(IcouTot, NatmTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(XatmTot, NatmTot, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(YatmTot, NatmTot, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(ZatmTot, NatmTot, MPI_REAL8, 0, comm_world, Ierr)

	call MPI_Bcast(IkcvTot, Mdim, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(FatmTot, Mdim, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(FwgtTot, Mdim, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(LabAtmTot, MaxLab*Mdim, MPI_CHAR, 0, comm_world, Ierr)

	call MPI_Bcast(RminPair, Npair, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(RmaxPair, Npair, MPI_REAL8, 0, comm_world, Ierr)
    
    if (YesPBC) then 
    call MPI_Bcast(VolTot , NsnpTot, MPI_REAL8, 0, comm_world, Ierr)

	call MPI_Bcast(ImaxTot, NsnpTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(JmaxTot, NsnpTot, MPI_INTEGER4, 0, comm_world, Ierr)
	call MPI_Bcast(KmaxTot, NsnpTot, MPI_INTEGER4, 0, comm_world, Ierr)

	call MPI_Bcast(AboxTot, NsnpTot, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(BboxTot, NsnpTot, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(CboxTot, NsnpTot, MPI_REAL8, 0, comm_world, Ierr)

	call MPI_Bcast(AvecTot, NsnpTot*3, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(BvecTot, NsnpTot*3, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(CvecTot, NsnpTot*3, MPI_REAL8, 0, comm_world, Ierr)

	call MPI_Bcast(UvecTot, NsnpTot*3, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(VvecTot, NsnpTot*3, MPI_REAL8, 0, comm_world, Ierr)
	call MPI_Bcast(WvecTot, NsnpTot*3, MPI_REAL8, 0, comm_world, Ierr)
    endif
    

#endif

!------------------write on screen -----------------------
	if(MyID==0) then
		print*
		print*, '>> Beginning CRYOFF fitting ......'
    endif
        
		if(Iopt==0) then
			Iter=0
			if(MyID==0) print*, '   Linear Parameter Fitting'
		else
			if(MyID==0) then
            print*, '   Non-Linear Parameter Optimization'

			print*, '   Initial Parameter'
			write(*, '(1X, 600F12.6)') P(1,1:Nopt)

			write(*,*) '  Evaluating Initial Object Function Value and Simplex'
            endif

            ! Getting initial objective function for all simplex vertices   

    	if (Iopt==1) then
	    	do i=1, Mopt
		    	X(1:Nopt) = P(i, 1:Nopt)
			    Y(i) = dErr(X)
		    end do
	    else if (Iopt==2 .or. Iopt==3) then
		    Yret = dErr(X)
	    end if

         if(MyID==0) then   
			if(Iopt==1) then
				do i=1, Mopt
					write(*, '(1X, 600F12.6)') Y(i), P(i, 1:Nopt)
				end do
			else if (Iopt==2 .or. Iopt==3) then
				print*, Yret
			end if
			write(*, '(1X, 75("-"))')

			print *, "Optimization Starting"
            write(Lout, '(A, /)') 'Optimization Iteration:'
		end if
	end if ! MyID 0

#ifdef MPI
	if( (Iopt==1 .or. Iopt==2) .and. np_pop/=1) &
		& call ErrStop('Simplex/Powell NOT parallized yet', 1)
#endif

	YesBug=.true.
 	if(Iopt==0) then
        if(MyID==0) then 
        write(*, '(A, I5, A)') 'Non-iterative Linear Optimization'
		write(Lout, '(A, I5, A, /)') 'Non-iterative Linear Optimization'  
        write(Lout, '(/, A)') 'Fitting Results:'
        endif
        ! For linear optimization this is the actual fit. 
        ipnt=2  ! dErr print results with ipnt >= 1
       	Yret = dErr(X)
    end if
    
    if(Iopt==1) then 
        call Simplex(Lout, Mopt, Nopt, P, Y, dErr, &
		&					 deps, Iter, MaxIter, MyID)
        Yret=Y(1)
        X(1:Nopt) = P(1, 1:Nopt)
    endif 
    
	if(Iopt==2) then 
        call Powell(Lout, Nopt, X, Pinf, Psup, Xi, deps, Iter, MaxIter, Yret, MyID)
    endif 
	if(Iopt==3) call DE_Fortran90(derr, Nopt, Pinf, Psup, deps, dei(1),&
		&            MaxIter, der(2), der(1), dei(2), dei(3), Lout, X, Yret,&
		&            iter, der(3), method, .False.,1,MyID, np_world, comm_world, MyID_pop,&
		&            np_pop,comm_pop,MyID_force, np_force,comm_force)
    

	if(MyID==0 .and. Iopt > 0) then
		print*

		write(Lout, '(A, I8)') '  No. of Iterations: ', Iter
        if (MaxIter==0) then
           
            else
		if (Iter>=MaxIter) then
			write(*, '(A, I5, A)') '>> Convergence criterion NOT met after', MaxIter, '  iterations'
			write(Lout, '(A, I5, A, /)') '  Convergence criterion NOT met after', MaxIter, '  iterations'
		else
			write(*, '(A, I5, A)') '>> Convergence criterion met after', Iter, '  iterations'
			write(Lout, '(A, I5, A, /)') '  Convergence criterion met after', Iter, '  iterations'
        end if
        endif 

		if(Iopt==1) then
			write(*, '(/, 1X, A)') '>>Vertices of the final simplex and function values at the vertices'
			write(*, '(A)') trim(text)//'                   Y'
			do i=1, Mopt
				write(*, '(2X, I3, 999F20.9)') i, P(i, 1:Nopt), Y(i)
			end do

			write(Lout, '(2X, A)') 'Vertices of the final simplex and function values at the vertices:'
			write(Lout, '(A)') trim(text)//'                   Y'
			do i=1, Mopt
				write(Lout, '(2X, I3, 600F20.9)') i, P(i, :), Y(i)
			end do
		else if(Iopt==3) then
			write(*, '(/, 1X, A)') '>>The final DE results:'
			do  i = 1, Nopt
				write(*, *) i, X(i)
			end do
		end if

		write(Lout, '(/, A)') 'Fitting Results:'
    end if !MyID 0
    
    if (Iopt > 0) then
        ipnt=1  ! dErr print results with ipnt >= 1
   		Yret = dErr(X)
    endif
      

	if(MyID==0) then
		close(Lout)

		print*
		print*, '>>>>>>>> CRYOFF Fitting Finished Successfully <<<<<<<<'
		print*
	end if

#ifdef MPI
	call MPI_Finalize(Ierr) ! MPI finish up ...
#endif
End Program CRYOFF
!
!
Real*8 Function dErr(X)
    use CRYOFFmod
    implicit none
    integer  i, j, k, ii, Ierr, &
        & Ntmp, Itmp, Jtmp, Ityp, Isnp, Icv, Iref, &
        & Iatm, IbgnAtm, IendAtm, &
        & NatmSnp, Ntyp, Ilnr, Nref, Lwork,LIwork, &
        & IbgnSnp, IendSnp, ItmpA
    integer jto,kcv,mtot, ompnum
    character*8 ompvalue
    real*8, allocatable :: Work(:)
    integer, allocatable :: IWork(:)

    real CPUini, CPUend
    real*8 wgt, Fqq, Fvdw, XYZi(3), XYZj(3), &
        &  Fi(3), Fj(3), Rtmp, ABC(3,3), Atrs(3,3), &
        &  X(Nopt)
    character*MaxLab  TagIatm, Lab
    character*256 Ttyp
    logical YesCst
    logical :: firstrun=.true.
    integer :: ncalls=0
    save :: firstrun, ncalls
    
! build the SVD matrix and solve the SVD problem

!initialize nonlinear parameters
    do i=1, Nopt
        Popt(i)%P=X(i)
    end do
    Amat = 0.d0
    Asav = 0.d0
    Pmat = 0.d0
    PmatKcv = 0.d0
    Ffix = 0.d0

    Fatm=FatmTot
    Fwgt=FwgtTot
    LabAtm=LabAtmTot
    IndxKCV=IkcvTot
    
    IbgnSnp=1
    IendSnp=NsnpTot
    
    IbgnAtm=IatmRefSnp(IbgnSnp)
    IendAtm=IatmRefSnp(IendSnp)+NatmRefSnp(IendSnp)-1

#ifdef MPI
    Itmp=int(NsnpTot/np_force)
    Jtmp=mod(NsnpTot, np_force)
    do i=1, np_force
        if(i<=Jtmp) then
            Ntmp=Itmp+1
            Isnp=Ntmp*(i-1)+1
        else
            Ntmp=Itmp
            Isnp=(Ntmp+1)*Jtmp+(i-Jtmp-1)*Ntmp+1
        end if
    
        Iatm=IatmRefSnp(Isnp)
    
        NcntRef(i)=3*(IatmRefSnp(Isnp+Ntmp-1)+NatmRefSnp(Isnp+Ntmp-1)-Iatm)
        IdspRef(i)=3*(Iatm-1)
    
        if(MyID_force+1==i .and. Ntmp>=1) then
            IbgnSnp=Isnp
            IbgnAtm=Iatm
            IendSnp=Isnp+Ntmp-1
            IendAtm=IatmRefSnp(IendSnp)+NatmRefSnp(IendSnp)-1
        end if
    end do
#endif

    if(YesOut(3) .and. (MyID<6 .or. np_force-MyID<6)) &
        & write(*, '(6X, A, I4, 4(A, I6))') 'MPI-ID: ', MyID, &
        & '    Frame: ', IbgnSnp, ' to ', IendSnp, &
        & '    Ref. Atom: ',   IbgnAtm, ' to ', IendAtm

    call cpu_time(CPUini)
    do Isnp=IbgnSnp, IendSnp
        IsnpTot=Isnp   !used in IntrFrcMat to check bond length only for the first frame. 
        Icv=IndxKCV(Isnp)
        Itmp=IbgnAtmSnp(Isnp)
        NatmSnp=NatmSnpTot(Isnp)
        Ntmp=Itmp+NatmSnp-1
        Name(1:NatmSnp)=NameTot(Itmp:Ntmp)
        Satm(1:NatmSnp)=SatmTot(Itmp:Ntmp)
        Icou(1:NatmSnp)=IcouTot(Itmp:Ntmp)
        Ivdw(1:NatmSnp)=IvdwTot(Itmp:Ntmp)
        Xatm(1:NatmSnp)=XatmTot(Itmp:Ntmp)
        Yatm(1:NatmSnp)=YatmTot(Itmp:Ntmp)
        Zatm(1:NatmSnp)=ZatmTot(Itmp:Ntmp)
        IatmRef(1:NatmSnp)=IatmRefTot(Itmp:Ntmp)
        IatmAtm(1:NatmSnp)=IatmAtmTot(Itmp:Ntmp)

        if(YesPBC) then
            Vol  = VolTot( Isnp)
            Abox = AboxTot(Isnp); Uvec     = UvecTot(Isnp, :)
            Bbox = BboxTot(Isnp); Vvec     = VvecTot(Isnp, :)
            Cbox = CboxTot(Isnp); Wvec     = WvecTot(Isnp, :)
            Imax = ImaxTot(Isnp); ABC(:,1) = AvecTot(Isnp, :)
            Jmax = JmaxTot(Isnp); ABC(:,2) = BvecTot(Isnp, :)
            Kmax = KmaxTot(Isnp); ABC(:,3) = CvecTot(Isnp, :)
            call invA(ABC, Atrs)
        end if


!      get the forces of FIX potential
        Nref = NmolRefTot(Isnp)
        Iref = IbgnMolRef(Isnp)

        !Force-Matrix elements Evaluations. 
         ! maybe it is a good idea to create a PBC module? 

        !Intermolecular
        call InteFrcMat(NatmSnp, ABC, Atrs)

        ! Loop over each molecule
        do ii=Iref, Iref+Nref-1   
            !Intramolecular
            if (.not. YesInt) call IntrFrcMat(ii, ABC, Atrs)

            call VirtFrcMat(ii)
        end do
        
        ! Compute Molecular force (Fi) and torque (Fj) to contribution from fixed forces 

        if (YesInt .or. YesHyb) then ! Inter or Hybrid only
            Fi = 0.d0; Fj = 0.d0
            jto=-1
            do i=1, NatmSnp
                Iatm = IatmRef(i)
                if(Iatm/=0) then
                    Itmp = 3*Iatm
                    if (jto<i) then !torque center has not been found. 
                        j = i
!                    it1=IatmAtmTot(j)
                    ! scan until find the center (j) to be used for computation of torq 
                        do while(VDWatmTyp(Ivdw(j))/='TORQ' .and. j<NatmSnp)
                            j = j+1
                        end do
                        jto=j
                    endif
!                    if (IatmAtmTot(j)==1) then 
!                           j=j-1 !Make sure to stay in the same molecule
!                    endif

                    XYZi = [ Xatm(i)-Xatm(jto), Yatm(i)-Yatm(jto), Zatm(i)-Zatm(jto) ]

! 		    write(101, *) Ffix(Itmp-2:Itmp)

                    if(VDWatmTyp(Ivdw(i))/='NETF' .and. VDWatmTyp(Ivdw(i))/='TORQ') then
                        Fi = Fi + Ffix(Itmp-2:Itmp)
                        Fj(1) = Fj(1) + XYZi(2)*Ffix(Itmp  ) - XYZi(3)*Ffix(Itmp-1)
                        Fj(2) = Fj(2) + XYZi(3)*Ffix(Itmp-2) - XYZi(1)*Ffix(Itmp  )
                        Fj(3) = Fj(3) + XYZi(1)*Ffix(Itmp-1) - XYZi(2)*Ffix(Itmp-2)
                    else 
                        if(VDWatmTyp(Ivdw(i))=='NETF') then
                            Ffix(Itmp-2:Itmp) = Fi
                            Fi = 0.d0
                        else if(VDWatmTyp(Ivdw(i))=='TORQ') then
                            Ffix(Itmp-2:Itmp) = Fj
                            Fj = 0.d0
                        end if
                    end if
                end if
            end do
        end if !total force and total torque. 

        ! removed fixed force before SVD
        do i=1, NatmSnp
            Itmp = IatmRef(i)*3
! 	    write(101, *) Fatm(Jtmp-2:Jtmp), Ffix(Jtmp-2:Jtmp)
            if(Itmp/=0) Fatm(Itmp-2:Itmp) = Fatm(Itmp-2:Itmp)-Ffix(Itmp-2:Itmp)
        end do

! 		the SVD matrix was built using ALL frames with the following format
! 		[dU/dx(1)|P1  dU/dx(1)|P2  dU/dx(1)|P3  ...  dU/dx(1)|Pn] [P1]   [Fx(1)]
! 		[dU/dy(1)|P1  dU/dy(1)|P2  dU/dy(1)|P3  ...  dU/dy(1)|Pn] [P2]   [Fy(1)]
! 		[dU/dz(1)|P1  dU/dz(1)|P2  dU/dz(1)|P3  ...  dU/dz(1)|Pn] [P3]   [Fz(1)]
! 		[                                                       ] [. ] = [ ... ]
! 		[dU/dx(m)|P1  dU/dx(m)|P2  dU/dx(m)|P3  ...  dU/dx(m)|Pn] [. ]   [Fx(m)]
! 		[dU/dy(m)|P1  dU/dy(m)|P2  dU/dy(m)|P3  ...  dU/dy(m)|Pn] [. ]   [Fy(m)]
! 		[dU/dz(m)|P1  dU/dz(m)|P2  dU/dz(m)|P3  ...  dU/dz(m)|Pn] [Pn]   [Fz(m)]

        if(Ipnt>1.and.MyID==0) print*, '   MPI-ID=0 Calculating SVD Matrix for Frame', Isnp

        if(YesInt .or. YesHyb) then
            Anet(:,:) = 0.d0
            Ator(:,:) = 0.d0
            jto=-1
            do i=1, NatmSnp
                Itmp = 3*IatmRef(i)
                ItmpA = Itmp-IdspRef(MyID+1)
                TagIatm = VDWatmTyp(Ivdw(i))

                if(Itmp/=0) then
                
                if (jto<i) then ! center of torque has not been found
                    j = i
                ! scan until find the center (j) to be used for computation of torq 
                ! The molecule needs to be on the same size of PBC. Otherwise need to call image on XYZi
                    do while(VDWatmTyp(Ivdw(j))/='TORQ' .and. j<NatmSnp)
                        j = j+1
                    end do
!                if (IatmAtmTot(j)==1) then 
!                        j=j-1                 !Make sure to stay in the same molecule
!                endif
                    jto=j
                endif 
                XYZi = [ Xatm(i)-Xatm(jto), Yatm(i)-Yatm(jto), Zatm(i)-Zatm(jto) ]

                if(TagIatm/='NETF' .and. TagIatm/='TORQ') then
                    Anet(1:3, :) = Anet(1:3, :) + Amat(ItmpA-2:ItmpA, :)

                    Ator(1, :) = Ator(1, :) + XYZi(2)*Amat(ItmpA,   :) - XYZi(3)*Amat(ItmpA-1, :)
                    Ator(2, :) = Ator(2, :) + XYZi(3)*Amat(ItmpA-2, :) - XYZi(1)*Amat(ItmpA,   :)
                    Ator(3, :) = Ator(3, :) + XYZi(1)*Amat(ItmpA-1, :) - XYZi(2)*Amat(ItmpA-2, :)

! zero this for inter fit.  
                    if (.not. YesHyb) then    
                    Amat(ItmpA-2:ItmpA, :) = 0.d0
                    Fmat(Itmp-2:Itmp, 1) = 0.d0
                    Fwgt(Itmp-2:Itmp)    = 0.d0
                    endif 
                else
                   ! applying weight for total force and torque.                    
                    Wgt = Fwgt(Itmp)
                    if(TagIatm=='NETF') then
                        Amat(ItmpA-2:ItmpA, :) = Anet(1:3, :)*Wgt
                        Anet = 0.d0
                    else if(TagIatm=='TORQ') then
                        Amat(ItmpA-2:ItmpA, :) = Ator(1:3, :)*Wgt
                        Ator = 0.d0
                    end if
                    
!                   changed Nov. 21, 2022. Needs some testing
                    if(sum(abs(Amat(ItmpA-2:ItmpA,:)))<Reps) then 
                        Fwgt(Itmp-2:Itmp) = 0.d0
                        Wgt=0.d0
                    endif 
                    ! right side of the equations. 
                    Fmat(Itmp-2:Itmp, 1) = Fatm(Itmp-2:Itmp)*Wgt

                end if
                end if
            end do
        else
            ! This code path for atomic fit. The netf and torq will be zerored. Thus does not allow inter-intra to be fit together.
            do i=1, NatmSnp
                Itmp = 3*IatmRef(i)
                ItmpA = Itmp-IdspRef(MyID+1)
                if(Itmp/=0) then
                !Remove atoms if all Amat elements are small. 
                !The Amat elements can be small for certain interactions such as exponential repulsion where the distances are large. 
                !This is probably safe as such rows are numerically unstable will probably need to tiny singular values anyway.
                if(VDWatmTyp(Ivdw(i))=='NETF' .or. VDWatmTyp(Ivdw(i))=='TORQ' .or. sum(abs(Amat(ItmpA-2:ItmpA,:)))<Reps) then
!                    if(VDWatmTyp(Ivdw(i))=='NETF' .or. VDWatmTyp(Ivdw(i))=='TORQ') then
                            Fwgt(Itmp-2:Itmp)    = 0.d0
                            Fmat(Itmp-2:Itmp, 1) = 0.d0
                    else
                            Amat(ItmpA-2:ItmpA, :) = Amat(ItmpA-2:ItmpA, :)*Fwgt(Itmp)
                            Fmat(Itmp-2:Itmp, 1) = Fatm(Itmp-2:Itmp)*Fwgt(Itmp)
                    end if
                end if
            end do
        end if
        ! second loop for hybrid fits to apply the weight for atomic foces. 
        if (YesHyb) then
                do i=1, NatmSnp
                Itmp = 3*IatmRef(i)
                ItmpA = Itmp-IdspRef(MyID+1)
                if(Itmp/=0) then
                !Remove atoms if all Amat elements are small. 
                if(sum(abs(Amat(ItmpA-2:ItmpA,:)))<Reps) then     
                            Fwgt(Itmp-2:Itmp)    = 0.d0
                            Fmat(Itmp-2:Itmp, 1) = 0.d0
                endif 
                if(VDWatmTyp(Ivdw(i))/='NETF' .and. VDWatmTyp(Ivdw(i))/='TORQ') then  ! Do not touch the weight for NETF and TORQ
                           Amat(ItmpA-2:ItmpA, :) = Amat(ItmpA-2:ItmpA, :)*Fwgt(Itmp)
                           Fmat(Itmp-2:Itmp, 1) = Fatm(Itmp-2:Itmp)*Fwgt(Itmp)
                end if
                end if
                end do
        endif 
        
        end do ! loop Isnp

#ifdef MPI
!--------------collect amat, fmat, fwgt, ffix from other force processes-----
! int MPI_Gatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
!                 void* recvbuf, int *recvcounts, int *displs,
!               MPI_Datatype recvtype, int root, MPI_Comm comm)
!     IN     sendbuf      发送消息缓冲区的起始地址(可变)
!     IN     sendcount    发送消息缓冲区中的数据个数(整型)
!     IN     sendtype     发送消息缓冲区中的数据类型(句柄)
!     OUT    recvbuf      接收消息缓冲区的起始地址(可变,仅对于根进程)
!     IN     recvcounts   整型数组(长度为组的大小), 其值为从每个进程接收的数据个数(仅对于根进程)
!     IN     displs       整数数组,每个入口i表示相对于recvbuf的位移, 此位移处存放着从进程i中接收的输入数据(仅对于根进程)
!     IN     recvtype     接收消息缓冲区中数据类型(仅对于根进程)(句柄)
!     IN     root         接收进程的序列号(句柄)
!     IN     comm         通信子(句柄)

        if(IbgnSnp<=IendSnp) then
            IbgnAtm=3*IbgnAtm-2
            IendAtm=3*IendAtm
            Ntmp=IendAtm-IbgnAtm+1
            do i=1, NlnrTot
                Fsav(1:Ntmp, 1)=Amat(1:Ntmp, i)
                call MPI_GATHERV(Fsav(:,1), Ntmp, MPI_REAL8, Amat(:,i), NcntRef, IdspRef, MPI_REAL8, 0, comm_force, Ierr)
            end do
        
            Fsav(1:Ntmp, 1)=Fmat(IbgnAtm:IendAtm, 1)
            call MPI_GATHERV(Fsav(:,1), Ntmp, MPI_REAL8, Fmat(:,1), NcntRef, IdspRef, MPI_REAL8, 0, comm_force, Ierr)
            
            fsav(1:Ntmp, 1)=Ffix(IbgnAtm:IendAtm)
            call MPI_GATHERV(fsav,Ntmp,MPI_REAL8,ffix,NcntRef,IdspRef,MPI_REAL8, 0, comm_force,Ierr)
            
            fsav(1:Ntmp, 1)=Fwgt(IbgnAtm:IendAtm)
            call MPI_GATHERV(fsav,Ntmp,MPI_REAL8, Fwgt,NcntRef,IdspRef,MPI_REAL8, 0, comm_force,Ierr)
            
            if(ipnt>0) then !Do not need this if we are not printing RminPair and RmaxPair
                call MPI_GATHER(RminPair,Npair,MPI_REAL8, fsav, Npair, MPI_REAL8, 0, comm_force,Ierr)
                if(MyID==0) then
                do ii=1,Npair
                    do Ntmp=1,np_force-1
                        RminPair(ii)=min(RminPair(ii),fsav(Ntmp*Npair+ii,1))
                    enddo
                enddo
                endif

                call MPI_GATHER(RmaxPair,Npair,MPI_REAL8, fsav, Npair, MPI_REAL8, 0, comm_force,Ierr)        
                if(MyID==0) then
                do ii=1,Npair
                    do Ntmp=1,np_force-1
                    RmaxPair(ii)=max(RmaxPair(ii),fsav(Ntmp*Npair+ii,1))
                    enddo
                enddo
                endif
        
            endif !ipnt
        
! 		write(999, *) IbgnSnp, IendSnp, NdatRpair, IdspRpair, fsav(1:Ntmp, 1)
! 		write(999, *) IbgnSnp, IendSnp, fsav(1:Ntmp, 1)
        end if
#endif

        if(MyID_force==0) then
       !Take care of charge constraint (only needed for head node)      
		Itmp = 3*NatmrefTot
		Mcst=0
		do i=1, Ncst
			Wgt  = Wchg(i)
			Fwgt(Itmp+i)   = 0.d0
			Fmat(Itmp+i,1) = 0.d0
			LabAtmTot(Itmp+i) = 'ChargeConstraint'
			if(Wgt>WgtTol) then
				YesCst=.false.
				Ntmp = NchgGrp(i)
				do j=1, Ntmp
					Jtmp = IpairQQ(Ichg(i, j))
					if(YesFitPair(Jtmp)) then
						Ilnr = IlnrPair(Jtmp)
						Amat(Itmp+i, Ilnr) = dble(Nchg(i, j))*Wgt
						YesCst=.true.
                    else
                        ! Why we need this when we have charge constraint but not fitting the pair? 
						Fchg(i)=Fchg(i)-dble(Nchg(i, j))*ParaPair(Jtmp,1)
					end if
				end do

				if(YesCst) then
					Mcst=Mcst+1
					Fwgt(Itmp+i) = Wgt
					Fmat(Itmp+i, 1) = Fchg(i)*Wgt
				end if
			end if
        end do

	call cpu_time(CPUend)
	if(MyID==0 .and. Ipnt>1) write(*, '(6X, A, F10.1, /)') 'MPI-ID=0 Total Wall Time(s) for Matrix Construction: ', CPUend-CPUini
	CPUini=CPUend
    
!----------------------------------------------------------------------------------
! Now start fit 
    ! In CV computation, Tot is for the entire set including all the folds. 
        
FwgtTot=Fwgt
FfixTot=Ffix
AmatTot=Amat
FmatTot=Fmat

do kcv=1, Nkcv
	Fwgt=FwgtTot
	Ffix=FfixTot
	Amat=AmatTot
	Fmat=FmatTot
	LabAtm=LabAtmTot

		Mfit = 0
		do i=1, Mdim
		if(Fwgt(i)>WgtTol) then  
		if( (Nkcv==1 .and. IndxKCV(i)==kcv) .or. (Nkcv>1 .and. IndxKCV(i)/=kcv) .or. LabAtm(i)=='ChargeConstraint' ) then
        ! Need more checking, for CV computation, the matrix should only be the size of the N-1 fold of data.  
				Mfit = Mfit+1
 				if(Mfit/=i) then
					Fwgt(Mfit) = Fwgt(i)
					Fmat(Mfit, 1) = Fmat(i, 1)
					Amat(Mfit, 1:NlnrTot) = Amat(i, 1:NlnrTot)
					LabAtm(Mfit) = LabAtm(i)    
					Fwgt(i)    = 0.d0
					Fmat(i, 1) = 0.d0
					Amat(i, 1:NlnrTot) = 0.d0
					LabAtm(i) = ''

					Ffix(Mfit) = Ffix(i)
					Ffix(i) = 0.d0
 				end if
		end if  !index KCV
		end if
		end do

!   		print*, 'MDIM' , Mdim, Mfit !, Ffix
        ! only need the inital none-zero part to be saved? 
		Asav(1:Mfit,1:NlnrTot) = Amat(1:Mfit,1:NlnrTot)
		Fsav(1:Mfit,1) = Fmat(1:Mfit,1)

        ! Does YesOut control level of printing? 
		if(YesOut(1)) then
			open(unit=Lbug, file=trim(adjustl(Fout))//'.mat')
			do i=1, Mfit
                write(Lbug,*) i,"out of ", Mfit, LabAtm(i), Fmat(i,1)
				write(Lbug, '(100E26.16)') Amat(i,:)
			end do
			close(Lbug)
        end if

     ! Now we finished the building the matrix
    ! Now starts the actual fitting. 
!!$http://www.cs.berkeley.edu/~demmel/DOE2000/Report0100.html
!!$Note that DGELSD is only 3 to 5 times slower than the fastest algorithm DGELS
!!$, whereas the old algorithm DGELSS was 7 to 34 times slower.
!!$ A divide-and-conquer approach is used in DGELSD.
        if (firstrun) then 
        call GET_ENVIRONMENT_VARIABLE('OMP_NUM_THREADS',ompvalue)
        ompnum=0
        read (ompvalue,*,iostat=ierr) ompnum
        if (ompnum < 2) then 
            print *, "Warning: set OMP_NUM_THREADS to larger value to take advantage of parallel SVD library routine DGELSD"
        else 
            print *, "DGELSD will use ",ompnum, "CPUs in paralel" 
        endif 
            firstrun=.false.
        endif
        ! work space query
		allocate(Work(1), IWork(1))
		call dgelsd(Mfit, NlnrTot, 1, Amat, Mdim, Fmat, Mdim, &
			&		Asig, Rcon, Irank, Work, -1, IWork, Ierr)
		call ErrStop('SVD failed', Ierr)

		Lwork=int(Work(1))
		LIwork=max(1, int(IWork(1)))
		deallocate(Work, IWork)

        ! Actual SVD with DGELSD
		allocate(Work(Lwork), IWork(LIwork), stat=Ierr)
        call ErrStop('SVD Memory Allocation Failure', Ierr)
		call dgelsd(Mfit, NlnrTot, 1, Amat, Mdim, Fmat, Mdim, Asig, &
		&			Rcon, Irank, Work, Lwork, IWork, Ierr)
		call ErrStop('SVD failed', Ierr)
        ! On exit Amat has been destroyed. Fmat is the solution 

		deallocate(Work, IWork)
        
        call cpu_time(CPUend)
        if(MyID==0 .and. Ipnt>1) write(*, '(6X, A, F10.1, /)') 'MPI-ID=0 Wall Time(s) for SVD: ', CPUend-CPUini
        

! 		backup potential parameters and fitted forces
		Pmat(1:NlnrTot) = Fmat(1:NlnrTot, 1)
        PmatKcv(1:NlnrTot,kcv) = Pmat(1:NlnrTot)

!       Compute objective function. 
        
        Chisq(0:4,kcv) = 0.d0

		do i=1, Mfit
			Rtmp = dot_product(Asav(i, 1:NlnrTot), Pmat(1:NlnrTot))
			Fmat(i, 1) = Rtmp
			Rtmp=(Fsav(i,1)-Rtmp)*(Fsav(i,1)-Rtmp)
			Chisq(0,kcv) = Chisq(0,kcv)+Rtmp
        enddo 
        
        if (Nkcv>1) then 
            MnetF=0; Mtorq=0; Matm=0
			do i=1, Mdim
			if(FwgtTot(i)>WgtTol) then
			if( IndxKCV(i)==kcv) then
				Rtmp = dot_product(AmatTot(i, 1:NlnrTot), Pmat(1:NlnrTot))
				Rtmp=((FmatTot(i, 1)-Rtmp)/FwgtTot(i))**2

				Lab = trim(LabAtmTot(i))
				if(index(Lab, '-NETF')/=0) then
					MnetF = MnetF+1
					Chisq(1,kcv) = Chisq(1,kcv)+Rtmp
				else if(index(Lab, '-TORQ')/=0) then
					Mtorq = Mtorq+1
					Chisq(2,kcv) = Chisq(2,kcv)+Rtmp
				else if(index(Lab, 'ChargeConstraint')==0) then
					Matm = Matm+1
					Chisq(3,kcv) = Chisq(3,kcv)+Rtmp
				end if

			end if
			end if
            end do
            
        endif !Nkcv>1
        
end do !KCV



!Cross Vadiation computation is not supposed to be used with non-linear optimization. 
        !Return value of this routine is used to drive non-linear optimization. 
        !Note that the Chisq is already weighted by Wgt^2. For Iwgt=3, the weight is RMS force component. 
        !Need to divide by Mfit, so the convergence does not become harder as the number of forces increase. 
        dErr=Chisq(0,1)/dble(Mfit)
        
        ncalls=ncalls+1
        if (MyID==0) then
        if (mod(ncalls,100) == 0) print *, "Number of calls to dErr ", ncalls, "obj func: ", dErr
        endif
                
		do i=1, Nopt
            if (PLR(i)) then
                dErr = dErr+PRP(i)/dble(Mfit)*(X(i)-PRC(i))*(X(i)-PRC(i))
                else
				if(X(i)<Pinf(i)) then
					dErr = dErr+100.d0*(X(i)-Pinf(i))*(X(i)-Pinf(i))
				else if(X(i)>Psup(i)) then
					dErr = dErr+100.d0*(X(i)-Psup(i))*(X(i)-Psup(i))
				end if
            endif
        end do
        
        if(ipnt>0) call printfittingoutput

	end if !MyID_froce==0

! the following boardcast is required for Simplex method
#ifdef MPI
	call MPI_Bcast(dErr, 1, MPI_REAL8, 0, comm_force, Ierr)
#endif
    End Function dErr
    
Subroutine printfittingoutput
use CRYOFFmod
implicit none
real*8      RMSE, Para(MaxPara), Ppot(MaxPara),Prms(MaxPara),Rtmp
real*8      SumSqPmatKcv(NlnrTot), PmatKcvRMS(NlnrTot)
integer     i,j, k, Ntyp, Ityp, Ntmp, Ilnr, Nlnr,Nfixq
character*MaxLab  tag, TagIatm, Lab
character*256 txt, Ttyp
real*8 Wgt, Fl(3), Fref(3), Ffit(3), Fdif(3)
real*8,  allocatable:: RMS(:,:,:)
integer, allocatable:: Nrms(:,:)
integer    kcv, Itmp, ii, nfreq
real*8  sumkcv1, sumkcv2, sumkcv3, sumkcv4
logical   YesPrn

! Start printing here. 
    write(LOUT,*) "Effective Rank of SVD Matrix,", Irank, "for ",NlnrTot, "linear parameters."
    if (Irank<NlnrTot) write(LOUT,*) "Warning: Rank deficiency"
    
    if(YesOut(2)) open(unit=Lbug, file=trim(adjustl(Fout))//'.fit')
    
    if (Nkcv> 1) then 
        !write(Lout,*) "For KCV calculation, the fitted gradients are only based on parameters fitted to the last fold." 
        write(Lout,'(2X,A,I3,A)') "KCV calculation: The fitted gradients are based on parameters from the averaged SVD parameters with ", Nkcv," fold."  
        sumkcv1=0.d0 ;sumkcv2=0.d0; sumkcv3=0.d0;Pmat(1:NlnrTot)=0.d0;SumSqPmatKcv(1:NlnrTot)=0
        do ii=1, Nkcv
            sumkcv1=sumkcv1+chisq(1,ii)
            sumkcv2=sumkcv2+chisq(2,ii)
            sumkcv3=sumkcv3+chisq(3,ii)
            Pmat(1:NlnrTot)=Pmat(1:NlnrTot)+PmatKcv(1:NlnrTot,ii)
            SumSqPmatKcv(1:NlnrTot)=SumSqPmatKcv(1:NlnrTot)+PmatKcv(1:NlnrTot,ii)**2
        enddo 
        if (MnetF>0) print*, "Cross Validation Value: (NetF)",sqrt(sumkcv1*3.d0/(dble(MnetF*Nkcv)))  !Assume equal folds 
        if (Mtorq>0) print*, "Cross Validation Value: (Torq)",sqrt(sumkcv2*3.d0/(dble(Mtorq*Nkcv)))  !Assume equal folds
        if (Matm>0) print*, "Cross Validation Value: (Atomic)",sqrt(sumkcv3*3.d0/(dble(Matm*Nkcv)))  !Assume equal folds
        Pmat(1:NlnrTot)=Pmat(1:NlnrTot)/Nkcv
        PmatKcvRMS(1:NlnrTot)=sqrt(SumSqPmatKcv(1:NlnrTot)/Nkcv-Pmat(1:NlnrTot)**2)
    endif
    
!------------------------------Print out summary of forces ---------------------------------
! Compute statistics

    allocate(RMS(0:NmolTyp, 1:3, 3), Nrms(0:NmolTyp, 1:3))
    RMS=0.d0; Nrms=0
!   		Ntmp=(Mfit-Mcst)/3
    
    LabAtm=LabAtmTot
    Fwgt=FwgtTot
    Ffix=FfixTot
    Amat=AmatTot  !Amat has original index, Asav is condensed by weight. 
    Fmat=FmatTot 
    Ntmp=NatmRefTot
        
    ii = 0
    do i=1, Ntmp
        Itmp=i*3
        Wgt  = Fwgt(Itmp)
        Fref = FatmTot(Itmp-2:Itmp)
        if (Wgt > WgtTol) then
            ii = ii + 1
            Ffit(1) = dot_product(Amat(Itmp-2, 1:NlnrTot), Pmat(1:NlnrTot))
            Ffit(2) = dot_product(Amat(Itmp-1, 1:NlnrTot), Pmat(1:NlnrTot))
            Ffit(3) = dot_product(Amat(Itmp, 1:NlnrTot), Pmat(1:NlnrTot))
!       		Fref = Fmat(Itmp-2:Itmp, 1)/Wgt+Ffix(Itmp-2:Itmp)
            Ffit = Ffit/Wgt+Ffix(Itmp-2:Itmp)
        else
!                     Fref = Ffix(Itmp-2:Itmp)
            Ffit = Ffix(Itmp-2:Itmp)
        endif 

        Fdif = Fref-Ffit
        Fl = Ffix(Itmp-2:Itmp)
        
        Lab = trim(LabAtm(Itmp))

        do j=1, NmolTyp
            if(index(Lab, trim(Mole(j)%Name))/=0) exit
        end do
        k=1
        if(index(Lab, '-NETF')/=0) k=2
        if(index(Lab, '-TORQ')/=0) k=3
        Nrms(0, k) = Nrms(0,k)+1
        Nrms(j, k) = Nrms(j,k)+1
        Rtmp=dot_product(Fref, Fref)
        RMS(0,k,1) = RMS(0,k,1)+Rtmp
        RMS(j,k,1) = RMS(j,k,1)+Rtmp
        Rtmp=dot_product(Ffit, Ffit)
        RMS(0,k,2) = RMS(0,k,2)+Rtmp
        RMS(j,k,2) = RMS(j,k,2)+Rtmp
        Rtmp=dot_product(Fdif, Fdif)
        RMS(0,k,3) = RMS(0,k,3)+Rtmp
        RMS(j,k,3) = RMS(j,k,3)+Rtmp

        YesPrn=.false.
        if (YesOut(3)) YesPrn=.true.
        nfreq=max(1,Ntmp/180)
        if ( ii > 0 .and. ii<21) YesPrn=.true.  ! Print the first 20 record that does not have a weight of zero. 
        if (mod(i,nfreq)==0) YesPrn=.true.
        if (i==Ntmp) YesPrn=.true.
                 
        if(i==1) then
            write(Lout, '(/, A, I0, A)') '  Fitting, Weight, and Error:'
            write(Lout, '(A)') '  RefIdx                Frm.Grp.Atm         Weight' &
                & //'        Ref           Fit           Fit-Ref           Fix'
        endif 
        
        if (.not. YesHyb) then 
        if (.not. YesInt .and. index(Lab, '-NETF')/=0) YesPrn=.false.
        if (.not. YesInt .and. index(Lab, '-TORQ')/=0) YesPrn=.false.
        if (YesInt .and.  index(Lab, '-NETF')==0 .and. index(Lab, '-TORQ')==0) YesPrn=.false.
        endif
                 
        if(YesPrn) then
            write(txt, '(2X, I0)') i
            write(Lout, '(A9, A32, F14.6, F14.6, F14.6, F14.6, F14.6)') &
                & txt, trim(Lab)//'.x', Wgt, Fref(1), Ffit(1), Fdif(1), Fl(1)
        
            write(Lout, '(A41, 14X, F14.6, F14.6, F14.6, F14.6)') &
                & '.y',  Fref(2), Ffit(2), Fdif(2), Fl(2)
        
            write(Lout, '(A41, 14X, F14.6, F14.6, F14.6, F14.6)') &
                & '.z',  Fref(3), Ffit(3), Fdif(3), Fl(3)
        endif
        
        
    end do
    if(YesOut(2)) close(Lbug)

    do i=1, Mcst
        Itmp=3*Ntmp+i
        Wgt= Fwgt(Itmp)
        Fref = Fmat(Itmp, 1)/Wgt
        Rtmp = dot_product(Amat(Itmp, 1:NlnrTot), Pmat(1:NlnrTot))
        Ffit = Rtmp/Wgt
        Fdif = Fref-Ffit
        Lab = trim(LabAtm(Itmp))
        write(txt, *) Ntmp+i; txt='    '//trim(adjustl(txt))
        write(Lout, '(A9, A22, F14.6, F14.6, F14.6, F14.6)') &
            & txt, trim(Lab), Wgt, Fref(1), Ffit(1), Fdif(1)
    end do
       
        ! Note when doing cross valdiation Chisq(:,1) overwritten at this stage. 
    Chisq(1:4,1) = 0.d0
    MnetF=0; Mtorq=0; Matm=0
    do i=1, Mdim
        Wgt = Fwgt(i)
        if (Wgt > WgtTol) then
            Rtmp = dot_product(Amat(i, 1:NlnrTot), Pmat(1:NlnrTot))
            Rtmp=(Fmat(i,1)-Rtmp)*(Fmat(i,1)-Rtmp)
!	Chisq(0,1) = Chisq(0,1)+Rtmp

            Lab = trim(LabAtm(i))
            if(index(Lab, '-NETF')/=0) then
                MnetF = MnetF+1
                Chisq(1,1) = Chisq(1,1)+Rtmp
            else if(index(Lab, '-TORQ')/=0) then
                Mtorq = Mtorq+1
                Chisq(2,1) = Chisq(2,1)+Rtmp
            else if(index(Lab, 'ChargeConstraint')/=0) then
                Chisq(3,1) = Chisq(3,1)+Rtmp
            else
                Matm = Matm+1
                Chisq(4,1) = Chisq(4,1)+Rtmp
            end if
        endif
    end do

    write(Lout, '(/, A)') '  Objective Function:'
    write(Lout, '(A)') '    Type       No.          chi^2'
    if(Mfit>0)  write(Lout, '(A, I6, ES24.15)') '    Total  ', Mfit,  Chisq(0,1)/dble(Mfit) 
    if(Matm>0)  write(Lout, '(A, I6, ES24.15)') '    Fatm   ', Matm,  Chisq(4,1)/dble(Matm)
    if(MnetF>0) write(Lout, '(A, I6, ES24.15)') '    NetF   ', MnetF, Chisq(1,1)/dble(MnetF)
    if(Mtorq>0) write(Lout, '(A, I6, ES24.15)') '    Torq   ', Mtorq, Chisq(2,1)/dble(Mtorq)
    if(Mcst>0)  write(Lout, '(A, I6, ES24.15)') '    Charge ', Mcst,  Chisq(3,1)/dble(Mcst)


    write(Lout, '(/, A)') '  RMS Error: (Per particle)'
    write(Lout, '(A)') '    Mol.Type            No.       RMS|Ref|        RMS|Fit|         RMSD|Fit-Ref|'
    NvirTot = 0
    do ii=1, NmolTyp+1
        i=ii; if(ii>NmolTyp) i=0
        txt='ALL'; if(i>0) txt=trim(Mole(i)%Name)
        j=len_trim(txt); k=10
        if (i > 0) then
            NvirTot = NvirTot + Mole(i)%NmolTot*Mole(i)%Nvir
            Nrms(i,1)=Nrms(i,1) - Mole(i)%NmolTot*Mole(i)%Nvir
        else
            Nrms(i,1)= Nrms(i,1) - NvirTot
        endif
        if(Nrms(i,1)>0) RMS(i,1,:)=sqrt(RMS(i,1,:)/dble(Nrms(i,1)))
        if(Nrms(i,2)>0) RMS(i,2,:)=sqrt(RMS(i,2,:)/dble(Nrms(i,2)))
        if(Nrms(i,3)>0) RMS(i,3,:)=sqrt(RMS(i,3,:)/dble(Nrms(i,3)))
        !if(Nrms(i,1)>0 .and. Matm>0) then
        if(Nrms(i,1)>0 .and. (.not. YesInt)) then             
            write(Lout, '(4X, A16, I6, 3F16.6)')   trim(txt)//'.Fatm '//txt(j+1:j+k), Nrms(i,1), RMS(i,1,1), RMS(i,1,2), RMS(i,1,3)
        endif
        !if(Nrms(i,2)>0 .and. MnetF>0) then 
        if(Nrms(i,2)>0 .and. (YesInt .or. YesHyb) ) then 
            write(Lout, '(4X, A16, I6, 3F16.6)')   trim(txt)//'.NetF '//txt(j+1:j+k), Nrms(i,2), RMS(i,2,1), RMS(i,2,2), RMS(i,2,3)
        endif
        !if(Nrms(i,3)>0 .and. Mtorq>0) then 
        if(Nrms(i,3)>0 .and. (YesInt .or. YesHyb) ) then 
            write(Lout, '(4X, A16, I6, 3F16.6)')   trim(txt)//'.Torq '//txt(j+1:j+k), Nrms(i,3), RMS(i,3,1), RMS(i,3,2), RMS(i,3,3)
        endif 
    end do

        if (Nrms(i,1)>0 .and. (.not. YesInt) .and. Matm == 0 )then
            write (Lout, '(4X,A)') 'Warning: atomic force is not being fitted.'
        endif
        if(Nrms(i,2)>0 .and. (YesInt .or. YesHyb) .and.  MnetF == 0 ) then
            write (Lout,'(4X,A)') 'Warning: net force is not being fitted.'
        endif
        if(Nrms(i,3)>0 .and. (YesInt .or. YesHyb) .and.  Mtorq == 0 ) then
            write (Lout,'(4X,A)') 'Warning: net torque is not being fitted.'
        endif

    write(Lout, *)
        
    deallocate(RMS, Nrms)           
        
    write(Lout, '(/, A, I4)') '  Rank of SVD Matrix: ', Irank
    write(Lout, '(2X, A)') 'Singular Value Spectrum:' 

    write(Lout, '(5G14.6)') Asig(1:Irank)
        
    write(Lout, *)

    Matm =Matm /3
    MnetF=MnetF/3
    Mtorq=Mtorq/3
    
    ! collect force field parameters
    do i=1, Npair
        if(YesFitPair(i)) then
            Ttyp = TtypPair(i)
            Ilnr = IlnrPair(i)
            Nlnr = NlnrPair(i)
            Ntmp = Ilnr+Nlnr-1
            Para(1:Nlnr) = Pmat(Ilnr:Ntmp)
            call getPotPara(Ttyp(1:3), Para, Ppot)
            ParaPair(i,1:Nlnr) = Ppot(1:Nlnr)
            call getParaRMS(Ttyp(1:3),Para, PmatKcvRMS(Ilnr:Ntmp) ,Prms)
            ParaRMS(i,1:Nlnr) = Prms(1:Nlnr)
        end if
    end do
    
    do i=1, NmolTyp
        Ntyp = Mole(i)%Ntyp
        do j=1, Ntyp
            Ttyp = Mole(i)%Ttyp(j)
            Ityp = Mole(i)%Ityp(j)
            Ilnr = Mole(i)%Ilnr(j)
            Nlnr = Mole(i)%Nlnr(j)
            Ntmp = Ilnr+Nlnr-1
            if(Mole(i)%YesFit(j)) then
                Para(1:Nlnr) = Pmat(Ilnr:Ntmp)
                call getPotPara(Ttyp, Para, Ppot)
                Mole(i)%Para(j,1:Nlnr)=PPot(1:Nlnr)
                call getParaRMS(Ttyp(1:3),Para, PmatKcvRMS(Ilnr:Ntmp) ,Prms)
                Mole(i)%ParaRMS(j,1:Nlnr) = Prms(1:Nlnr)
            endif
            
            if((Ityp==3 .and. Ttyp(1:3)=='HAR') &
              & .or. (Ityp==4 .and. Ttyp(1:3)=='COS') &
              & .or. (Ityp==4 .and. Ttyp(1:3)=='HAR') ) &
              & Mole(i)%Para(j,1)=Mole(i)%Para(j,1)*Rad2Deg
            if(Ityp==4 .and. Ttyp(1:3)=="NCO") Mole(i)%Para(j,3)=Mole(i)%Para(j,3)*Rad2Deg
           ! if(Ityp==4 .and. Ttyp(1:3)=="COS") Mole(i)%Para(j,1)=Mole(i)%Para(j,1)*Rad2Deg
            if(Ityp==5 .and. Ttyp(1:3)=="CNC") Mole(i)%Para(j,4)=Mole(i)%Para(j,4)*Rad2Deg
            if(Ityp==5 .and. Ttyp(1:3)=="CCO") Mole(i)%Para(j,1)=Mole(i)%Para(j,1)*Rad2Deg
        end do
    end do
    
    !print CV
    if(Nkcv>1) then
        write(Lout, '(A)') '  KCV information:'
        call printKCV(PmatKcvRMS)
    endif
    
    !print force field
    write(Lout, '(/, A)') '  Force Field Potential:'
    if(Nkcv>1) write(Lout,'(A,I3,A)') "   KCV calculation: The parameters are based on the averaged SVD parameters with ", Nkcv," fold."                
    write(Lout, '(A)') '    Intra-Potential:'
    do i=1, NmolTyp
        do j=1, Mole(i)%Ntyp
            write(Lout, '(8X, A, I4, A, 9G18.8)') '[ '//trim(Mole(i)%Ttyp(j))//' ]', Mole(i)%Ngrp(j), '  FIxT  ', Mole(i)%Para(j,1:Mole(i)%Npar(j))
            do k=1, Mole(i)%Ngrp(j)
                if(Mole(i)%Ityp(j)==6) then   !bond-angle cross
                    write(Lout, *) '        ',Mole(i)%Iatm(j, k, 1:3)
                else
                    write(Lout, *) '        ', Mole(i)%Iatm(j, k, 1:Mole(i)%Ityp(j))
                endif
            end do
        end do
    end do
            
    write(Lout, '(A)') '    Inter-Potential:'
    do i=1, Npair
        if (RmaxPair(i)<Reps) RminPair(i) = 0.d0
            write(Lout, '(\, 8X, A, 99G18.8)') trim(TatmPair(i))//': '//trim(TtypPair(i)), ParaPair(i,1:NparPair(i))
        if(YesFitPair(i) .and. Nkcv > 1) then
            write(Lout, '(\A, 99G18.8)') ' RMSD:', ParaRMS(i,1:NlnrPair(i))
        end if
        write(Lout, '(2(A, F12.6))') '    Min:', RminPair(i), '    Max:', RmaxPair(i)
    end do

    if(YesCMD) then
        write(Lout, '(/, A)') '    Charge product Matrix Decomposition (CMD)' 
        call CMDcharge(ParaPair)
    endif
    
    write(Lout, '(/, A)') '    Molecular-Definition:'
    do i=1, NmolTyp
        write(Lout, '(4X, A)') trim(Mole(i)%Name)
        do j=1, Mole(i)%Ntyp
            call getBondedTag(i,j,txt)
            txt='#define    '//trim(txt)
            ii=Mole(i)%Ityp(j)
            if(ii==2) then
                if(trim(Mole(i)%Ttyp(j)) == "QUA") then !quartic bond
                    write(Lout, '(8X, A, 4G18.8)') trim(txt), Mole(i)%Para(j,1)/10., Mole(i)%Para(j,2)*4.184/0.1**2, &
                    & Mole(i)%Para(j,3)*4.184/0.1**3, Mole(i)%Para(j,4)*4.184/0.1**4
                else !harmonic bond
                    if(Mole(i)%YesFit(j) .and. Nkcv > 1 ) then
                        write(Lout, '(8X, A, G18.8, A,2G18.8,A,G18.8 )') trim(txt), Mole(i)%Para(j,1)/10.,'  RMSD:', Mole(i)%ParaRMS(j,1)/10, Mole(i)%Para(j,2)*4.184/0.1**2,' RMSD:', Mole(i)%ParaRMS(j,2)*4.184/0.1**2
                    else
                        write(Lout, '(8X, A, 2G18.8)') trim(txt), Mole(i)%Para(j,1)/10., Mole(i)%Para(j,2)*4.184/0.1**2
                    endif
                endif
            elseif(ii==3) then !harmonic angle
                if(Mole(i)%YesFit(j) .and. Nkcv > 1 ) then
                    write(Lout, '(8X, A, G18.8, A,2G18.8,A,G18.8)') trim(txt), Mole(i)%Para(j,1),' RMSD:', Mole(i)%ParaRMS(j,1)*Rad2Deg, Mole(i)%Para(j,2)*4.184,' RMSD:', Mole(i)%ParaRMS(j,2)*4.184
                else
                    write(Lout, '(8X, A, 2G18.8)') trim(txt), Mole(i)%Para(j,1), Mole(i)%Para(j,2)*4.184
                endif
            elseif(ii==6) then !bond-angle, bond-bond coupling
                if(trim(Mole(i)%Ttyp(j)) == "MUB") then
                    write(Lout, '(8X, A, 4G18.8)') trim(txt), Mole(i)%Para(j,1)*4.184/0.01, Mole(i)%Para(j,2)/10., Mole(i)%Para(j,3)/10., Mole(i)%Para(j,4)/10.    
                elseif(trim(Mole(i)%Ttyp(j)) == "QBB") then
                    write(Lout, '(8X, A, 5G18.8)') trim(txt), Mole(i)%Para(j,1)/10., Mole(i)%Para(j,2)*4.184/0.01, Mole(i)%Para(j,3)*4.184/0.01, &
                    & Mole(i)%Para(j,4)*4.184/0.1**3,Mole(i)%Para(j,5)*4.184/0.1**4 
                endif
            elseif(ii==4) then !dihedral
                if(trim(Mole(i)%Ttyp(j)) == "NCO") then
                    if(Mole(i)%YesFit(j) .and. Nkcv > 1 ) then
                        write(Lout, '(8X, A, 2G18.8,A, G18.8,I3)') trim(txt), Mole(i)%Para(j,3), Mole(i)%Para(j,1)*4.184, '  RMSD:', Mole(i)%ParaRMS(j,1)*4.184, int(Mole(i)%Para(j,2))
                    else
                        write(Lout, '(8X, A, 2G18.8,I3)') trim(txt), Mole(i)%Para(j,3), Mole(i)%Para(j,1)*4.184, int(Mole(i)%Para(j,2))
                    endif
                elseif (trim(Mole(i)%Ttyp(j)) == "COS") then
                    write(Lout, '(8X, A,2G18.8,I3, A)') trim(txt), Mole(i)%Para(j,1),Mole(i)%Para(j,2)*4.184,int(Mole(i)%Para(j,3))
                else
                    write(Lout, '(8X, A,2G18.8,I3, A)') trim(txt), Mole(i)%Para(j,1),Mole(i)%Para(j,2)*4.184,int(Mole(i)%Para(j,3))
                endif
            elseif(ii==5) then !phi-psi cross term
                if(trim(Mole(i)%Ttyp(j))=="CNC") then
                    write(Lout, '(8X, A, 2G18.8,2I3, A)') trim(txt), Mole(i)%Para(j,4), Mole(i)%Para(j,1)*4.184, &
                    & int(Mole(i)%Para(j,2)), int(Mole(i)%Para(j,3))
                elseif(trim(Mole(i)%Ttyp(j))=="CCO") then
                    write(Lout, '(8X, A, 2G18.8,2I3, A)') trim(txt), Mole(i)%Para(j,1), Mole(i)%Para(j,2)*4.184, &
                    & int(Mole(i)%Para(j,3)), int(Mole(i)%Para(j,4))
                else
                    write(Lout, '(8X, A, 2G18.8,2I3, A)') trim(txt), Mole(i)%Para(j,1), Mole(i)%Para(j,2)*4.184, &
                    & int(Mole(i)%Para(j,3)), int(Mole(i)%Para(j,4))
                endif
            endif
        end do
    end do
    
    write(Lout, '(/, A)') '    Table-Potential:'
    Iuniq=0
    do i=1, Npair
        if(Iuniq(i)==0) then
           Iuniq(i)=i
           tag=TatmPair(i)
           do j=i+1, Npair
               if(Iuniq(j)==0) then
                  txt  = TatmPair(j)
                  k=index(txt, '~')
                  if(trim(tag)==trim(txt) .or. trim(tag) &
                  & ==trim(txt(k+1:))//'~'//txt(1:k-1)) Iuniq(j)=i
               end if
           end do
       end if
    end do
    
    txt=''
    do i=1, Npair
        k=Iuniq(i)
        if(k==i) then
            write(Lout, '(\, 8X, A, 99G18.8)') trim(TatmPair(i))//': '//trim(TtypPair(i)), ParaPair(i,1:NparPair(i))
            do j=i+1, Npair
                if(Iuniq(j)==k) write(Lout, '(\, 4X, A, 99G18.8)') trim(TtypPair(j)), ParaPair(j,1:NparPair(j))
            end do
            write(Lout, *)
        end if
    end do
    
End subroutine printfittingoutput

    
Subroutine printPotIndex()
use CRYOFFmod
implicit none
integer   i,j, Ntyp
character*256 txt
    
    !print 
    write(Lout, '(A)') '    Inter-Potential Index:' ! nonlinear sometime later
    do i=1, Npair
        if(YesFitPair(i)) then
            write(Lout, '(8X, A, A, I5)') trim(TatmPair(i))//': '//trim(TtypPair(i)), '    Index:', IlnrPair(i)
        endif
    end do

    write(Lout, '(/, A)') '    Intra-Potential Index:'
    do i=1, NmolTyp
        write(Lout, '(4X, A)') trim(Mole(i)%Name)
        Ntyp = Mole(i)%Ntyp
        do j=1, Ntyp
            if(Mole(i)%YesFit(j)) then
                call getBondedTag(i,j,txt)
                write(Lout, '(8X, A, A, I5)') trim(txt),'    Index: ', (Mole(i)%Ilnr(j)) 
            endif
        end do
    end do
End subroutine printPotIndex     

    
Subroutine printKCV(kcvPmatRMS)
use CRYOFFmod
implicit none
real*8    KcvPmatRMS(NlnrTot), KcvPmat(MaxPara)
integer   i,j, k, Ntyp, Ityp, Ntmp, Ilnr, Nlnr
character*MaxLab  tag, TagIatm, Lab
character*256 txt, Ttyp
    
    !print 
    write(Lout, '(A)') '    Inter-Potential Pmat:'
    do i=1, Npair
        if(YesFitPair(i)) then
            Ttyp = TtypPair(i)
            Ilnr = IlnrPair(i)
            Nlnr = NlnrPair(i)
            Ntmp = Ilnr+Nlnr-1
            KcvPmat(1:Nlnr) = Pmat(Ilnr:Ntmp)
            do j=1, Nlnr
                write(Lout, '(8X, A\)') trim(TatmPair(i))//': '//trim(TtypPair(i))
                if(abs(KcvPmat(j))>Reps) then
                    write(Lout, '(A,G18.8, A, G18.8, A ,F12.2)') ' Pmat: ', KcvPmat(j),' rmsd:', KcvPmatRMS(Ilnr+j-1), &
                        & ' %rmsd:', abs(KcvPmatRMS(Ilnr+j-1)/KcvPmat(j))*100
                else
                    write(Lout, '(A,G18.8, A, G18.8, A)') ' Pmat: ', KcvPmat(j),' rmsd:', KcvPmatRMS(Ilnr+j-1), &
                        & ' %rmsd: NA'
                endif
            enddo
        endif
    end do

    write(Lout, '(/, A)') '    Intra-Potential Pmat:'
    do i=1, NmolTyp
        write(Lout, '(4X, A)') trim(Mole(i)%Name)
        Ntyp = Mole(i)%Ntyp
        do j=1, Ntyp
            if(Mole(i)%YesFit(j)) then
                Ilnr = Mole(i)%Ilnr(j)
                Nlnr = Mole(i)%Nlnr(j)
                Ntmp = Ilnr+Nlnr-1
                KcvPmat(1:Nlnr)= Pmat(Ilnr:Ntmp)
                call getBondedTag(i,j,txt)  
                write(Lout, '(8X, A\)') trim(txt)
                do k=1, Nlnr-1
                    if(abs(KcvPmat(k))>Reps) then
                        write(Lout, '(A,G18.8, A, G18.8, A ,F12.2\)') ' Pmat: ',KcvPmat(k),' rmsd:', KcvPmatRMS(Ilnr+k-1), &
                            & ' %rmsd:', abs(KcvPmatRMS(Ilnr+k-1)/KcvPmat(k))*100
                    else
                        write(Lout, '(A,G18.8, A, G18.8, A\)') ' Pmat: ', KcvPmat(k),' rmsd:', KcvPmatRMS(Ilnr+k-1),' %rmsd: NA'
                    endif
                enddo
                k=Nlnr
                if(abs(KcvPmat(k))>Reps) then
                    write(Lout, '(A,G18.8, A, G18.8, A ,F12.2)') ' Pmat: ',KcvPmat(k),' rmsd:', KcvPmatRMS(Ilnr+k-1), &
                        & ' %rmsd:', abs(KcvPmatRMS(Ilnr+k-1)/KcvPmat(k))*100
                else
                    write(Lout, '(A,G18.8, A, G18.8, A)') ' Pmat: ', KcvPmat(k),' rmsd:', KcvPmatRMS(Ilnr+k-1),' %rmsd: NA'
                endif
            endif
        end do
    end do

End subroutine printKCV


Subroutine getBondedTag(MoleID,NtypID,txt)
use CRYOFFmod
implicit none
integer   i,j, MoleID, NtypID, Ityp
character*MaxLab  tag, TagIatm
character*256 txt,Ttyp
    i=MoleID; j=NtypID
    Ttyp = Mole(i)%Ttyp(j)
    Ityp = Mole(i)%Ityp(j)
    if(Ityp==2) then
        TagIatm=Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 1))
        tag=Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 2))
        txt=trim(TagIatm)//'_'//trim(tag)
        if(TagIatm>tag) txt=trim(tag)//'_'//trim(TagIatm)
        txt=trim(txt)//'    '//trim(Ttyp)
    elseif(Ityp==3) then
        TagIatm=Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 1))
        tag=Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 3))
        txt=trim(TagIatm)//'_'//trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 2)))//'_'//trim(tag)
        if(TagIatm>tag) txt=trim(tag)//'_'//trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 2)))//'_'//trim(TagIatm)
        txt=trim(txt)//'    '//trim(Ttyp)                
    elseif(Ityp==6) then 
        TagIatm=Mole(i)%COUatm(Mole(i)%Iatm(j,1, 1))   
        tag=Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 3))
        txt=trim(TagIatm)//'_'//trim(Mole(i)%COUatm(Mole(i)%Iatm(j,1, 2)))//'_'//trim(tag)
        if(TagIatm>tag) txt=trim(tag)//'_'//trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 2)))//'_'//trim(TagIatm)
        txt=trim(txt)//'    '//trim(Ttyp)                
    elseif(Ityp==4) then
        TagIatm=Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 2))
        tag=Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 3))
        txt=trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 1)))//'_'&
        & //trim(TagIatm)//'_'//trim(tag)//'_'//trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 4)))
        if(TagIatm>tag) txt=trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 4)))//'_'&
        & //trim(tag)//'_'//trim(TagIatm)//'_'//trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1, 1)))
        txt=trim(txt)//'    '//trim(Ttyp)
    elseif(Ityp==5) then
        txt=trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1,1)))//'_'//trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1,2)))//'_' &
            & //trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1,3)))//'_'//trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1,4)))//'_' &
            & //trim(Mole(i)%COUatm(Mole(i)%Iatm(j, 1,5)))
        txt=trim(txt)//'    '//trim(Ttyp)                 
    endif   
End Subroutine getbondedtag    

    
Subroutine CMDcharge(CMDParaPair)
use CRYOFFmod
implicit none
integer :: ii,jj,kk,chargsign
real*8 :: Lambda, CMDParaPair(Npair, MaxPara)
character*80 :: tag1, tag2, tagIatm, tagJatm
integer, allocatable:: CMDpairindex(:,:)   !global or local?
real*8, allocatable:: CMDAmat(:,:), CMDchg(:)      !global or local?
integer :: INFO, LWORK
integer, parameter :: LWMAX=1000
real*8 :: WORK(LWMAX)
real*8, allocatable :: U(:,:), VT(:,:), S(:)
logical:: YesChgProd=.false., YesChg=.false.
    
    call getNumOftag(CMDatmlist,CMDchgNum)                              
    allocate(CMDatm(CMDchgNum), CMDpairindex(CMDchgNum, CMDchgNum), CMDAmat(CMDchgNum,CMDchgNum),CMDchg(CMDchgNum))
    allocate(U(CMDchgNum,CMDchgNum), VT(CMDchgNum,CMDchgNum), S(CMDchgNum))
    
    read(CMDatmlist,*) CMDatm(1:CMDchgNum)
    do ii=1,CMDchgNum
        do jj=ii,CMDchgNum
            YesChgProd=.false.
            tagIatm=trim(CMDatm(ii))
            tagJatm=trim(CMDatm(jj))
            tag1=trim(TagIatm)//'~'//trim(TagJatm)
            tag2=trim(TagJatm)//'~'//trim(TagIatm)
            do kk=1,NpairQQ
                if(tag1==trim(TatmIpairQQ(kk)) .or. tag2==trim(TatmIpairQQ(kk)) ) then 
                    YesChgProd=.true.       !charge product has been found. 
                    exit
                endif
            enddo
            if (YesChgProd) then
                CMDpairindex(ii,jj)=IpairQQ(kk)
                CMDAmat(ii,jj)=CMDParapair(CMDpairindex(ii,jj),1)
                CMDAmat(jj,ii)=CMDAmat(ii,jj)
            else
                write(Lout, '(4X, A, A, A)') '    CMD matrix can not be constructed because the charge product ',trim(tag1), ' is missing.'
                return
            endif
        enddo
    enddo

!   Query the optimal workspace.
    LWORK = -1
    call DGESVD( 'All', 'All', CMDchgNum, CMDchgNum, CMDAmat, CMDchgNum, S, U, CMDchgNum, VT, CMDchgNum, WORK, LWORK, INFO )
    LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!   Compute SVD.
    call DGESVD( 'All', 'All', CMDchgNum, CMDchgNum, CMDAmat, CMDchgNum, S, U, CMDchgNum, VT, CMDchgNum, WORK, LWORK, INFO )
!   Check for convergence.
    if( INFO.GT.0 ) then    
        write(*,*) 'CMD SVD failed to converge.'
        return
    endif
      
!Compute charges
    write(Lout, '(4X, A)') '    First three Singular Values for CMD'
    ii=min(3,CMDchgNum)
    write(Lout, '(4X, 3F14.6)') S(1:ii)

    do ii=1,CMDchgNum
        Lambda=U(1,ii)/VT(ii,1)
        if ((-1.d0-100*Reps) <= Lambda .and. Lambda <= (-1.d0+100*Reps)) then
            write(Lout, '(4X, A, F14.6 , A )') '    Singular Value ', S(ii) ,': skipped since it corresponds to a  product of imaginary charges '
            if (ii==CMDchgNum) write(Lout,'(4X, A )') '   No product of real charges can be found.'
        else if ((1.d0-100*Reps) <= Lambda .and. Lambda <= (1.d0+100*Reps)) then
            YesChg=.true.
            exit  
        else
            write(Lout, '(4X, A)') '    CMD internal consistency failure.'
            return
        endif
    enddo
    
    if(YesChg) then
        write(Lout, '(4X, A)') '    CMD charges'
        Lambda=U(1,ii)/abs(U(1,ii))
        if ((-1.d0-100*Reps) < Lambda .and. Lambda < (-1.d0+100*Reps)) then
            chargsign=-1
        else
            chargsign=1
        endif
        do jj=1,CMDchgNum
            CMDchg(jj)=U(jj,ii)*S(ii)*VT(ii,jj)
            if(chargsign==1) then
                write(Lout, '(8X, A, F12.5)') trim(CMDatm(jj)), sign(sqrt(CMDchg(jj)), U(jj,ii))
            else
                write(Lout, '(8X, A, F12.5)') trim(CMDatm(jj)), sign(sqrt(CMDchg(jj)), -U(jj,ii))
            endif
        enddo
    endif

End subroutine CMDcharge

    
Subroutine getNumOfCMDatm(tag,Qnum)
    use CRYOFFmod
    implicit none
	integer i, j, Qnum
	character(*) tag
    
    Qnum=0
    j=1
    if(len_trim(tag)>0) then
        do i=1,len_trim(tag)
            if (tag(i:i)==' ') then
                Qnum=Qnum+1
            endif
            if (i==len_trim(tag)) then
                Qnum=Qnum+1
            endif
        enddo
    endif
End Subroutine getNumOfCMDatm
    
        
Subroutine setWgt(Iwgt, Wgt, F, Frms, yestorq)
    use CRYOFFmod, only: wgtfac, Reps
    implicit none

    integer :: Iwgt
    real*8 :: Wgt, F(3), Frms   !Frms is weighted
    logical, intent (in) :: yestorq
	real*8:: Fabs,wgtf

    wgtf=wgtfac
    if (yestorq) wgtf=1.d0
    
    select case(Iwgt) 
        case (3)
	    if( Frms>Reps) Wgt=Wgt*wgtf/Frms
        case (4)
	    Fabs=norm2(F)            
	    if( Fabs>Reps) then 
        if (Fabs>wgtfac*Frms) then 
            Wgt=Wgt/Fabs
        else
            Wgt=Wgt/(wgtfac*Frms)
        endif 
        endif
        case default
        call ErrStop('Invalid Weight',max(Iwgt,1))
        end select
End Subroutine
!
!
Subroutine VirtFrcMat(Iref)
	use CRYOFFmod
	integer i, j, k, Jtmp, Ktmp, Iref, JtmpA
	real*8  Ftmp(3)
!   Redistribute forces on Virtual Site to atoms. 
	i=ImolRefTot(Iref)
	do j=1, Mole(i)%Nvir
		Jtmp = 3*IatmRef(IbgnAtmRef(Iref)+Mole(i)%IatmVir(j))
		if(Jtmp==0) cycle
        JtmpA = Jtmp-IdspRef(MyID+1)
		Ftmp = Ffix(Jtmp-2:Jtmp)
		Ffix(Jtmp-2:Jtmp) = 0.d0

		Anet(1:3, :) = Amat(JtmpA-2:JtmpA, :)
		Amat(JtmpA-2:JtmpA, :) = 0.d0

        ! The Mole%Ivir array is only used for distributing forces not used to compute position
		do k=1, Mole(i)%NatmVir(j)
			Ktmp = 3*IatmRef(IbgnAtmRef(Iref)+Mole(i)%Ivir(j, k))
			if(Ktmp==0) cycle
            ! It is indeed better to simply update both. 
            KtmpA = Ktmp-IdspRef(MyID+1)
			Ffix(Ktmp-2:Ktmp) = Ffix(Ktmp-2:Ktmp) + Mole(i)%Rvir(j, k)*Ftmp
			Amat(KtmpA-2:KtmpA, :) = Amat(KtmpA-2:KtmpA, :) + Mole(i)%Rvir(j, k)*Anet(1:3, :)
		end do
	end do
End Subroutine VirtFrcMat
!
!
Subroutine IntrFrcMat(Iref, ABC, Atrs)
	use CRYOFFmod
    implicit none
    
	integer i, j, k, Iatm, Jatm, Katm, Latm, Ratm, &
		&	Iref, Ityp, Ngrp, Nlnr, Ilnr, Ntmp, &
        &   IatmA, JatmA, KatmA, LatmA, RatmA

	real*8  Rtmp, Para(MaxPara), &
		&	Fi(3), Fj(3), Fk(3), Fl(3), Fr(3), &
		&	XYZi(3), XYZj(3), XYZk(3), XYZl(3), XYZr(3), &
		&	ABC(3,3), Atrs(3,3)
    
   
    real*8 Rtmp1,Rtmp2,Rtmp3,Rtmp4
!  Update Matrix Elements or compute the force contribution from fixed parameters (YesFrc) for Intramolecular terms.
    
	logical YesFrc
	character*256 Ttyp
    
	i=ImolRefTot(Iref)
	do j=1, Mole(i)%Ntyp
		Ityp = Mole(i)%Ityp(j)
		Ttyp = Mole(i)%Ttyp(j)
		Ngrp = Mole(i)%Ngrp(j)
		Para = Mole(i)%Para(j, :)
		YesFrc=.NOT. Mole(i)%YesFit(j)
		if(.NOT. YesFrc) then
			Nlnr = Mole(i)%Nlnr(j)
			Ilnr = Mole(i)%Ilnr(j)
			Ntmp = Ilnr+Nlnr-1
		end if

		do k=1, Ngrp
			Iatm = IbgnAtmRef(Iref)+Mole(i)%Iatm(j, k, 1)
			Jatm = IbgnAtmRef(Iref)+Mole(i)%Iatm(j, k, 2)
			Katm = 0; Latm = 0; Ratm = 0

			XYZi = [ Xatm(Iatm), Yatm(Iatm), Zatm(Iatm) ]
			XYZj = [ Xatm(Jatm), Yatm(Jatm), Zatm(Jatm) ]

			if(Ityp==2) then
				XYZj  = XYZj-XYZi
                if(YesPBC) call setPBC(Rtmp, XYZj, Abox, Bbox, Cbox, ABC, Atrs)
                Rtmp1=norm2(XYZj)
                if (IsnpTot == 1) then
                    if (Rtmp1 > bdwarn)    print*, "Warning, bond too long check .ff file", k, Iatm, Jatm
                endif
                        
                ! The distance has been precomputed. However, you still need the vector. 
				call BondFrcMat(YesFrc, Ttyp, XYZj, Rtmp1, Para, Ai, Aj, Fi, Fj)
			else if(Ityp==3) then
				Katm = IbgnAtmRef(Iref)+Mole(i)%Iatm(j, k, 3)
				XYZi = XYZi-XYZj
				XYZk = [ Xatm(Katm), Yatm(Katm), Zatm(Katm) ] - XYZj
				if(YesPBC) then
					call setPBC(Rtmp, XYZi, Abox, Bbox, Cbox, ABC, Atrs)
					call setPBC(Rtmp, XYZk, Abox, Bbox, Cbox, ABC, Atrs)
                end if
                Rtmp1=norm2(XYZi)
                Rtmp2=norm2(XYZk)
                if (IsnpTot == 1) then
                    if (Rtmp1 > bdwarn)    print*, "Warning, bond in angle too long check .ff file", k, Iatm, Jatm
                    if (Rtmp2 > bdwarn)    print*, "Warning, bond in angle too long check .ff file", k, Jatm, Katm
                endif
                                
				call AngleFrcMat(YesFrc, Ttyp, XYZi, XYZk, &
					&	Rtmp1, Rtmp2, Para, Ai, Aj, Ak, Fi, Fj, Fk)
			else if(Ityp==4) then
				Katm = IbgnAtmRef(Iref)+Mole(i)%Iatm(j, k, 3)
				Latm = IbgnAtmRef(Iref)+Mole(i)%Iatm(j, k, 4)

                XYZk = [ Xatm(Katm), Yatm(Katm), Zatm(Katm) ]
				XYZl = [ Xatm(Latm), Yatm(Latm), Zatm(Latm) ] - XYZk
				XYZk = XYZk - XYZj;                 Rtmp1=norm2(XYZk)
                XYZj = XYZj - XYZi

				if(YesPBC) then
					call setPBC(Rtmp, XYZj, Abox, Bbox, Cbox, ABC, Atrs)
					call setPBC(Rtmp, XYZk, Abox, Bbox, Cbox, ABC, Atrs)
					call setPBC(Rtmp, XYZl, Abox, Bbox, Cbox, ABC, Atrs)
                end if
                
                if (IsnpTot == 1) then
                    Rtmp=norm2(XYZj); if (Rtmp > bdwarn)    print*, "Warning, bond in dihedral too long check .ff file"
                    Rtmp=Rtmp1; if (Rtmp > bdwarn)    print*, "Warning, bond in dihedral too long check .ff file"
                    Rtmp=norm2(XYZl); if (Rtmp > bdwarn)    print*, "Warning, bond in dihedral too long check .ff file"
                endif

				call DihdrFrcMat(YesFrc, Ttyp, XYZj, XYZk, XYZl, &
					 & Rtmp1, Para, &
					 & Ai, Aj, Ak, Al, Fi, Fj, Fk, Fl)
                
            else if(Ityp==5) then
                Katm = IbgnAtmRef(Iref)+Mole(i)%Iatm(j, k, 3)
                Latm = IbgnAtmRef(Iref)+Mole(i)%Iatm(j, k, 4)
                Ratm = IbgnAtmRef(Iref)+Mole(i)%Iatm(j, k, 5)
                
                XYZk = [ Xatm(Katm), Yatm(Katm), Zatm(Katm) ]
                XYZl = [ Xatm(Latm), Yatm(Latm), Zatm(Latm) ]
                
                XYZr = [ Xatm(Ratm), Yatm(Ratm), Zatm(Ratm) ] - XYZl
                XYZl = XYZl - XYZk ;Rtmp2=norm2(XYZl)
                XYZk = XYZk - XYZj ;Rtmp1=norm2(XYZk)
                XYZj = XYZj - XYZi                 

                if(YesPBC) then
                    call setPBC(Rtmp, XYZj, Abox, Bbox, Cbox, ABC, Atrs)
                    call setPBC(Rtmp, XYZk, Abox, Bbox, Cbox, ABC, Atrs)
                    call setPBC(Rtmp, XYZl, Abox, Bbox, Cbox, ABC, Atrs)
                    call setPBC(Rtmp, XYZr, Abox, Bbox, Cbox, ABC, Atrs)
                end if
                
                if (IsnpTot == 1) then
                    Rtmp=norm2(XYZj); if (Rtmp > bdwarn)    print*, "Warning, bond in dihedral too long check .ff file"
                    Rtmp=Rtmp1; if (Rtmp > bdwarn)    print*, "Warning, bond in dihedral too long check .ff file"
                    Rtmp=Rtmp2; if (Rtmp > bdwarn)    print*, "Warning, bond in dihedral too long check .ff file"
                    Rtmp=norm2(XYZr); if (Rtmp > bdwarn)    print*, "Warning, bond in dihedral too long check .ff file"
                endif
                
                call CDihdrFrcMat(YesFrc, Ttyp, XYZj, XYZk, XYZl, XYZr, &
                    & Rtmp1, Rtmp2, Para, &
                    & Ai, Aj, Ak, Al, Ar, Fi, Fj, Fk, Fl, Fr)                
                
			else if(Ityp==6) then
				Katm = IbgnAtmRef(Iref)+Mole(i)%Iatm(j, k, 3)
				XYZk = [ Xatm(Katm), Yatm(Katm), Zatm(Katm) ]
				XYZi = XYZi-XYZj; Rtmp1=norm2(XYZi)
				XYZk = XYZk-XYZj; Rtmp2=norm2(XYZk)
				XYZj = XYZi-XYZk; Rtmp3=norm2(XYZj)
				if(YesPBC) then
					call setPBC(Rtmp, XYZi, Abox, Bbox, Cbox, ABC, Atrs)
					call setPBC(Rtmp, XYZk, Abox, Bbox, Cbox, ABC, Atrs)
					call setPBC(Rtmp, XYZj, Abox, Bbox, Cbox, ABC, Atrs)
                end if
                
                if (IsnpTot == 1) then
                    if (Rtmp1 > bdwarn)    print*, "Warning, bond in bb coupling too long check .ff file"
                    if (Rtmp2 > bdwarn)    print*, "Warning, bond in bb coupling too long check .ff file"
                    if (Rtmp3 > bdwarn)    print*, "Warning, bond in bb coupling too long check .ff file"
                endif
                
				call Bond3bFrcMat(YesFrc,Ttyp, XYZi, XYZk, XYZj, &
					& Rtmp1, Rtmp2, Rtmp3, Para, Ai, Aj, Ak, Fi, Fj, Fk)
			end if

! 				print*, Rtmp, ABC, Atrs, Fi, Fj

			Iatm = 3*IatmRef(Iatm)
			Jatm = 3*IatmRef(Jatm)
            IatmA = Iatm-IdspRef(MyID+1)
            JatmA = Jatm-IdspRef(MyID+1)
			if(YesFrc) then
				Ffix(Iatm-2:Iatm) = Ffix(Iatm-2:Iatm)+Fi
				Ffix(Jatm-2:Jatm) = Ffix(Jatm-2:Jatm)+Fj
			else
				if(Iatm/=0) Amat(IatmA-2:IatmA, Ilnr:Ntmp) = Amat(IatmA-2:IatmA, Ilnr:Ntmp)+Ai(1:3,1:Nlnr)
				if(Jatm/=0) Amat(JatmA-2:JatmA, Ilnr:Ntmp) = Amat(JatmA-2:JatmA, Ilnr:Ntmp)+Aj(1:3,1:Nlnr)
			end if

			if(Ityp==3 .or. Ityp==4 .or. Ityp==5 .or. Ityp==6 ) then
				Katm = 3*IatmRef(Katm)
                KatmA = Katm -IdspRef(MyID+1)
				if(YesFrc) then
					Ffix(Katm-2:Katm) = Ffix(Katm-2:Katm)+Fk
				else
					if(Katm/=0) Amat(KatmA-2:KatmA, Ilnr:Ntmp) = Amat(KatmA-2:KatmA, Ilnr:Ntmp)+Ak(1:3, 1:Nlnr)
				end if
			end if
			if(Ityp==4 .or. Ityp==5 ) then
				Latm = 3*IatmRef(Latm)
                LatmA = Latm-IdspRef(MyID+1)
				if(YesFrc) then
					Ffix(Latm-2:Latm) = Ffix(Latm-2:Latm)+Fl
				else
					if(Latm/=0) Amat(LatmA-2:LatmA, Ilnr:Ntmp) = Amat(LatmA-2:LatmA, Ilnr:Ntmp)+Al(1:3, 1:Nlnr)
				end if
            end if
            if(Ityp==5) then       
                Ratm = 3*IatmRef(Ratm)
                RatmA = Ratm-IdspRef(MyID+1)
                if(YesFrc) then
                    Ffix(Ratm-2:Ratm) = Ffix(Ratm-2:Ratm)+Fr
                else
                    if(Ratm/=0) Amat(RatmA-2:RatmA, Ilnr:Ntmp) = Amat(RatmA-2:RatmA, Ilnr:Ntmp)+Ar(1:3, 1:Nlnr)
                endif
            endif
		end do
	end do
End Subroutine IntrFrcMat
!
!
Subroutine EwaldFrcMat(NatmSnp, ABC, Atrs)
	use CRYOFFmod
	integer i, j, k, ii, jj, kk, Iatm, Jatm, NatmSnp, Ilnr, n, Jmin, Kmin
	real*8  Rtmp, Rfac, Para(MaxPara), &
		&	XYZi(3), XYZj(3), ABC(3,3), Atrs(3,3), &
		&	Ftmp(3), Kx, Ky, Kz, QQtmp, Rbet
	real*8, allocatable:: Frec(:,:)

	logical YesFit
	character*256 Ttyp
	character*MaxLab  TagIatm, TagJatm, Pair, Grp

	integer, parameter:: Nkmax=30
! 	real*8	Rx, Ry, Rz, Kx(0:Nkmax), Ky(0:Nkmax), Kz(0:Nkmax), &
! &		Ksq(0:Nkmax, 0:Nkmax, 0:Nkmax), &
! &		sinKx(0:Nkmax), sinKy(0:Nkmax), sinKz(0:Nkmax), &
! &		cosKx(0:Nkmax), cosKy(0:Nkmax), cosKz(0:Nkmax)

! 	direct part of Ewald
!   We already decided if we want to do ewald earlier. We probably do not need to keep determining. 
!   The original design is allowing COU and EWA to coexist, which does not make sense. 
    
    ! Note this code is not working the new code now only has part of the Amat in slave nodes. 
    
	YesEwa=.false.
	do k=1, Npair
		Ttyp = TtypPair(k)
		if(Ttyp(1:3)/='EWA' .and. Ttyp(1:3)/='PEW') cycle
		if(Ttyp(1:3)=='EWA') YesEwa=.true.

		Pair = TatmPair(k)
		Para = ParaPair(k, :)
		Ilnr = IlnrPair(k)
		YesFit=YesFitPair(k)
		do i=1, NatmSnp
			Grp  = Name(i)
			Iatm = 3*IatmRef(i)
			TagIatm = COUatmTyp(Icou(i))
			XYZi = [ Xatm(i), Yatm(i), Zatm(i) ]
			do j=1, i-1
				Jatm = 3*IatmRef(j)
				TagJatm = COUatmTyp(Icou(j))
				if(Iatm+Jatm==0) cycle
				if(trim(TagIatm)//'~'//TagJatm/=Pair &
				& .and. trim(TagJatm)//'~'//TagIatm/=Pair ) cycle

				XYZj = [ Xatm(j), Yatm(j), Zatm(j) ]-XYZi
				call setPBC(Kx, XYZj, Abox, Bbox, Cbox, ABC, Atrs)

				Ftmp=0.d0
				if(Kx<=Rcou) then
					Rbet=Kx*beta
					Rtmp=2.d0*Rbet*exp(-Rbet*Rbet)/sqrt(Pi)
					Ftmp = -(erfc(Rbet)+Rtmp)*XYZj/(Kx*Kx*Kx)
					if(YesExc(i, j)) Ftmp = Ftmp+XYZj/(Kx*Kx*Kx)
				end if

				if(Ttyp(1:3)=='PEW') then
					Rbet=-0.25d0/(beta*beta)
					Rfac=-4.d0*TwoPi/Vol
					if(.NOT. YesCub) then
						do ii=0, Imax
							Jmin=-Jmax; if(ii==0) Jmin=0
							do jj=Jmin, Jmax
								Kmin=-Kmax; if(ii==0 .and. jj==0) Kmin=1
								do kk=Kmin, Kmax
									Kx=dot_product( dble([ii,jj,kk]), [Uvec(1),Vvec(1),Wvec(1)] )
									Ky=dot_product( dble([ii,jj,kk]), [Uvec(2),Vvec(2),Wvec(2)] )
									Kz=dot_product( dble([ii,jj,kk]), [Uvec(3),Vvec(3),Wvec(3)] )
									Rtmp = Kx*Kx+Ky*Ky+Kz*Kz
									Rtmp = exp(Rbet*Rtmp)/Rtmp*sin(Kx*XYZj(1)+Ky*XYZj(2)+Kz*XYZj(3))
									Ftmp = Ftmp + Rfac*Rtmp*[Kx, Ky, Kz]
								end do
							end do
						end do
					else
! 						Ksq=0.d0; Kmax=max(Imax, Jmax, Kmax)
! 						do ii=0, Kmax
! 							Kx(ii)=dble(ii)*u(1)
! 							do jj=0, Kmax
! 								Ky(jj)=dble(jj)*v(2)
! 								do kk=0, Kmax
! 									if(ii+jj+kk/=0) then
! 										Kz(kk)=dble(kk)*w(3)
! 										Rtmp=Kx(ii)*Kx(ii)+Ky(jj)*Ky(jj)+Kz(kk)*Kz(kk)
! 										Ksq(kk, jj, ii) = exp(Rbet*Rtmp)/Rtmp
! 									end if
! 								end do
! 							end do
! 						end do

! 						do ii=1, Kmax
! 							sinKx(ii)=sin(Kx(ii)*dXYZ(1)); cosKx(ii)=cos(Kx(ii)*dXYZ(1))
! 							sinKy(ii)=sin(Ky(ii)*dXYZ(2)); cosKy(ii)=cos(Ky(ii)*dXYZ(2))
! 							sinKz(ii)=sin(Kz(ii)*dXYZ(3)); cosKz(ii)=cos(Kz(ii)*dXYZ(3))
! 						end do
! 						do ii=1, Kmax
! 							Rx=Kx(ii)*sinKx(ii)
! 							Ry=Ky(ii)*sinKy(ii)
! 							Rz=Kz(ii)*sinKz(ii)
! 							do jj=1, Kmax
! 								do kk=1, Kmax
! 									Ftmp(1) = Ftmp(1) + 4.d0*Rx*Ksq(kk,jj,ii)*cosKy(jj)*cosKz(kk)
! 									Ftmp(2) = Ftmp(2) + 4.d0*Ry*Ksq(kk,ii,jj)*cosKx(jj)*cosKz(kk)
! 									Ftmp(3) = Ftmp(3) + 4.d0*Rz*Ksq(ii,jj,kk)*cosKx(kk)*cosKy(jj)
! 								end do
! 								Ftmp(1) = Ftmp(1) + 2.d0*Rx*( Ksq(jj,0,ii)*cosKz(jj) + Ksq(0,jj,ii)*cosKy(jj) )
! 								Ftmp(2) = Ftmp(2) + 2.d0*Ry*( Ksq(jj,ii,0)*cosKz(jj) + Ksq(0,ii,jj)*cosKx(jj) )
! 								Ftmp(3) = Ftmp(3) + 2.d0*Rz*( Ksq(ii,jj,0)*cosKy(jj) + Ksq(ii,0,jj)*cosKx(jj) )
! 							end do
! 							Ftmp(1) = Ftmp(1) + Ksq(0, 0, ii)*Rx
! 							Ftmp(2) = Ftmp(2) + Ksq(0, ii, 0)*Ry
! 							Ftmp(3) = Ftmp(3) + Ksq(ii, 0, 0)*Rz
! 						end do
					end if
				end if

				if(YesFit) then
					if(Iatm/=0) Amat(Iatm-2:Iatm, Ilnr) = Amat(Iatm-2:Iatm, Ilnr)+Ftmp
					if(Jatm/=0) Amat(Jatm-2:Jatm, Ilnr) = Amat(Jatm-2:Jatm, Ilnr)-Ftmp
				else
					Rtmp=Para(1) !*FcovCou
					if(Iatm/=0) Ffix(Iatm-2:Iatm) = Ffix(Iatm-2:Iatm)+Ftmp*Rtmp
					if(Jatm/=0) Ffix(Jatm-2:Jatm) = Ffix(Jatm-2:Jatm)-Ftmp*Rtmp
				end if

			end do
		end do
	end do

	if(.NOT. YesEwa) return

! 	Optimized Ewald o(N^2), may also for other pair interactions
	allocate(QQsin(-Imax:Imax, -Jmax:Jmax, -Kmax:Kmax, Nqq), &
		&	 QQcos(-Imax:Imax, -Jmax:Jmax, -Kmax:Kmax, Nqq), Frec(3, Nqq) )

	QQsin = 0.d0; QQcos = 0.d0
	do ii=1, NatmSnp
		TagIatm=COUatmTyp(Icou(ii))
		do n=1, Nqq    ! find the charge idx
			if(trim(TagQQ(n))==trim(TagIatm)) then
				XYZi = [ Xatm(ii), Yatm(ii), Zatm(ii) ]
				do i=0, Imax
					Jmin=-Jmax; if(i==0) Jmin=0
					do j=Jmin, Jmax
						Kmin=-Kmax; if(i==0 .and. j==0) Kmin=1
						do k=Kmin, Kmax
							Ftmp=dble([i,j,k])
							Rtmp=XYZi(1)*( Ftmp(1)*Uvec(1)+Ftmp(2)*Vvec(1)+Ftmp(3)*Wvec(1) ) &
							&	+XYZi(2)*( Ftmp(1)*Uvec(2)+Ftmp(2)*Vvec(2)+Ftmp(3)*Wvec(2) ) &
							&	+XYZi(3)*( Ftmp(1)*Uvec(3)+Ftmp(2)*Vvec(3)+Ftmp(3)*Wvec(3) )
							QQsin(i,j,k, n) = QQsin(i,j,k, n)+sin(Rtmp)
							QQcos(i,j,k, n) = QQcos(i,j,k, n)+cos(Rtmp)
						end do
					end do
				end do
				exit
			end if
		end do
	end do

	Rbet=-0.25d0/(beta*beta)
	Rfac=-4.d0*TwoPi/Vol
	do ii=1, NatmSnp
		Iatm = 3*IatmRef(ii)
		if(Iatm==0) cycle

		TagIatm=COUatmTyp(Icou(ii))
		YesEwa=.false.
		do n=1, Nqq
			if(trim(TagQQ(n))==trim(TagIatm)) then
				YesEwa=.true.
				exit
			end if
		end do
		if(.not. YesEwa) cycle

		Frec=0.d0; Ftmp=0.d0
		XYZi = [ Xatm(ii), Yatm(ii), Zatm(ii) ]
		do n=1, Nqq
			do kk=1, Npair
				Ttyp = TtypPair(kk)
				Pair = TatmPair(kk)
				if(Ttyp(1:3)=='EWA' &
				& .and. (trim(TagIatm)//'~'//TagQQ(n)==Pair &
				&	.or. trim(TagQQ(n))//'~'//TagIatm==Pair) ) then
					Ilnr = IlnrPair(kk)
					Para = ParaPair(kk, :)
					YesFit=YesFitPair(kk)
					exit
				end if
			end do

			QQtmp=0.d0
			do i=0, Imax
				Jmin=-Jmax; if(i==0) Jmin=0
				do j=Jmin, Jmax
					Kmin=-Kmax; if(i==0 .and. j==0) Kmin=1
					do k=Kmin, Kmax
						XYZj=dble([i,j,k])
						Kx=dot_product( XYZj, [Uvec(1),Vvec(1),Wvec(1)] )
						Ky=dot_product( XYZj, [Uvec(2),Vvec(2),Wvec(2)] )
						Kz=dot_product( XYZj, [Uvec(3),Vvec(3),Wvec(3)] )

						Rtmp=Kx*XYZi(1)+Ky*XYZi(2)+Kz*XYZi(3)
						QQtmp = QQsin(i,j,k, n)*cos(Rtmp) &
						&	  - QQcos(i,j,k, n)*sin(Rtmp)

						Rtmp = Kx*Kx+Ky*Ky+Kz*Kz
						Rtmp = QQtmp*exp(Rbet*Rtmp)/Rtmp

						Frec(:, n) = Frec(:, n) + Rtmp*[Kx, Ky, Kz]
					end do
				end do
			end do

			if(.not. YesFit) Ftmp=Ftmp+Para(1)*Frec(:, n)
			if(YesFit) Amat(Iatm-2:Iatm, Ilnr) = Amat(Iatm-2:Iatm, Ilnr)+Frec(:, n)*Rfac
		end do
		if(.not. YesFit) Ffix(Iatm-2:Iatm) = Ffix(Iatm-2:Iatm)+Ftmp*Rfac !*FcovCou
	end do

	deallocate(QQsin, QQcos, Frec)

End Subroutine EwaldFrcMat
!
!
Subroutine InteFrcMat(NatmSnp, ABC, Atrs)
	use CRYOFFmod
    implicit none
	integer i, j, k, Iatm, Jatm, NatmSnp, Ilnr, Nlnr, Ntmp, K1, K2
    integer IatmA, JatmA
    integer Kic, Kjc, Kiw, Kjw
	real*8  Rtmp, Fudge, Para(MaxPara), &
		&	XYZi(3), XYZj(3), Fi(3), Fj(3), ABC(3,3), Atrs(3,3)
    real*8 Fvdw, Fqq
    integer ii, ij
	logical YesFrc
    integer cpptp,ccpptp   !current number of pairs per type pair for speed up compuation
	character*3 Ttyp
    character*MaxLab Lab
    
    ! It seems the code will do ewald without PBC!
    ! Note that Jicun recycle YesEwa locally in subroutines.
    ! One has to be careful. 
    if (YesEwa) call EwaldFrcMat(NatmSnp, ABC, Atrs)
    Para(:) = 0.d0

		do i=2, NatmSnp
            Lab  = Name(i)
			Iatm = 3*IatmRef(i) !Iatm is actually Icomponent in the gradient routine.
            IatmA = Iatm-IdspRef(MyID+1)
                      
			XYZi = [ Xatm(i), Yatm(i), Zatm(i) ]
			do j=1, i-1
                
             if (Iatm+IatmRef(j)==0) cycle
! Get fudge factor right
             
				Fvdw=1.d0;	Fqq=1.d0 
                
               ! Same molecule. determine fudge factor
				if(Name(j)==Lab) then
					do ii=1, NmolTyp
						if(index(Lab, trim(Mole(ii)%Name))/=0) exit
					end do

					if(Mole(ii)%YesExc( IatmAtm(i), IatmAtm(j))) then 
                        Fvdw=0.d0
                        Fqq=0.d0
                    endif 
! Looks like YesPair overwrites YesExc. Could be confusing if a pair can be in both states simultaneously. 
					if(Mole(ii)%YesPair(IatmAtm(i), IatmAtm(j))) then
						Fvdw=Mole(ii)%FudgeVDW
						Fqq =Mole(ii)%FudgeQQ
					end if

                end if
                
                XYZj = [ Xatm(j), Yatm(j), Zatm(j) ]-XYZi
                if(YesPBC) call setPBC(Rtmp, XYZj, Abox, Bbox, Cbox, ABC, Atrs)
                Rtmp = norm2(XYZj)
                XYZj=XYZj/Rtmp
                
                if(YesPBC) then 
                    ! Rvdw and Rcou only affect PBC calculations. 
                    ! Assume YesPBC will have coulombic thus yesEwald
 					YesExc(i, j)=Mole(ii)%YesExc(IatmAtm(i), IatmAtm(j))
					YesExc(j, i)=YesExc(i, j)                   
                    if(Rtmp>Rvdw) Fvdw=0.d0
                    if(Rtmp>Rcou) Fqq=0.d0                    
                endif
      
                Kic=Icou(i)
                Kjc=Icou(j)
                Kiw=Ivdw(i)
                Kjw=Ivdw(j)
                
            
            if((npptp(Kiw,Kjw,1)==-1).or.(cnpptp(Kic,Kjc,1)==-1)) then
!                if(.true.) then
            
            cpptp=0
            ccpptp=0
            
            do k=1, Npair
               
		    Ttyp  = TtypPair(k)(1:3)
		    Para = ParaPair(k, :)
		    YesFrc= .NOT. YesFitPair(k)
            
            Fudge=0.0
		    if(Ttyp=='COU' .or. Ttyp=='THC' ) then
!                if (Fqq<Reps) cycle
			    K1=IcouPair(k, 1)
			    K2=IcouPair(k, 2)
                if((Kic==K1 .and. Kjc==K2) .or. (Kic==K2 .and. Kjc==K1)) then 
                    Fudge=Fqq
                    ccpptp=ccpptp+1 
!                    if (ccpptp>Maxpptp) call ErrStop ("Increase maxpptp and recompile",ccpptp)
                    if (ccpptp>1) then
                    call ErrStop ("Having two different Coulombic Interactions between a pair of atoms do not make sense", ccpptp) 
                    endif
                    cnpptp(Kic,Kjc,ccpptp)=k
                    cnpptp(Kjc,Kic,ccpptp)=k
                endif 
            else  !vdw
!                if (Fvdw<Reps) cycle
		        K1=IvdwPair(k, 1)
			    K2=IvdwPair(k, 2)
				if ((Kiw==K1 .and. Kjw==K2) .or. (Kiw==K2 .and. Kjw==K1)) then 
                    Fudge = Fvdw
                    cpptp=cpptp+1
                    if (cpptp>Maxpptp) call ErrStop ("Increase maxpptp and recompile",cpptp)
                    npptp(Kiw,Kjw,cpptp)=k
                    npptp(Kjw,Kiw,cpptp)=k
                endif
            end if
            if (Fudge<Reps) cycle
  
              
                Jatm = 3*IatmRef(j)
                JatmA = Jatm-IdspRef(MyID+1)
               
                    ! Compute if there is a match
!					YesCub=.false.   !Maybe eventually needed for Ewald? Removed for now
                    ! Store minimum distance

					RminPair(k)=min(RminPair(k), Rtmp)
					RmaxPair(k)=max(RmaxPair(k), Rtmp)

                   
					call PairFrcMat(YesFrc, Ttyp, XYZj, Rtmp, Para, Ai, Aj, Fi, Fj)

					if(YesFrc) then
						if(Iatm/=0) Ffix(Iatm-2:Iatm) = Ffix(Iatm-2:Iatm)+Fi*Fudge
						if(Jatm/=0) Ffix(Jatm-2:Jatm) = Ffix(Jatm-2:Jatm)+Fj*Fudge
                    else
                        Ilnr = IlnrPair(k)
			            Nlnr = NlnrPair(k)
			            Ntmp = Ilnr+Nlnr-1
						if(Iatm/=0) Amat(IatmA-2:IatmA, Ilnr:Ntmp) = Amat(IatmA-2:IatmA, Ilnr:Ntmp)+Ai(1:3,1:Nlnr)*Fudge
						if(Jatm/=0) Amat(JatmA-2:JatmA, Ilnr:Ntmp) = Amat(JatmA-2:JatmA, Ilnr:Ntmp)+Aj(1:3,1:Nlnr)*Fudge
                    end if
              
            end do  !Enddo Npair 
            if (cpptp==0) then 
                npptp(Kiw,Kjw,1)=0
                npptp(Kjw,Kiw,1)=0
            endif
            if (ccpptp==0) then 
                cnpptp(Kic,Kjc,1)=0
                cnpptp(Kjc,Kic,1)=0
            endif
            
            else 
                
!            do ii=1, Maxpptp
                do ij=1, 1  ! Maximum one coulombic 
                k=cnpptp(Kic,Kjc,ij)
                if (k<=0) exit
                
		    Ttyp  = TtypPair(k)(1:3)
		    Para = ParaPair(k, :)
		    YesFrc= .NOT. YesFitPair(k)
            
            Fudge=Fqq
  
            if (Fudge<Reps) cycle

                Jatm = 3*IatmRef(j)
                JatmA = Jatm-IdspRef(MyID+1)

                      ! Compute if there is a match
!					YesCub=.false.   !Maybe eventually needed for Ewald? Removed for now
                    ! Store minimum distance

					RminPair(k)=min(RminPair(k), Rtmp)
					RmaxPair(k)=max(RmaxPair(k), Rtmp)

					call PairFrcMat(YesFrc, Ttyp, XYZj, Rtmp, Para, Ai, Aj, Fi, Fj)

					if(YesFrc) then
						if(Iatm/=0) Ffix(Iatm-2:Iatm) = Ffix(Iatm-2:Iatm)+Fi*Fudge
						if(Jatm/=0) Ffix(Jatm-2:Jatm) = Ffix(Jatm-2:Jatm)+Fj*Fudge
                    else
                        Ilnr = IlnrPair(k)
			            Nlnr = NlnrPair(k)
			            Ntmp = Ilnr+Nlnr-1
						if(Iatm/=0) Amat(IatmA-2:IatmA, Ilnr:Ntmp) = Amat(IatmA-2:IatmA, Ilnr:Ntmp)+Ai(1:3,1:Nlnr)*Fudge
						if(Jatm/=0) Amat(JatmA-2:JatmA, Ilnr:Ntmp) = Amat(JatmA-2:JatmA, Ilnr:Ntmp)+Aj(1:3,1:Nlnr)*Fudge
                    end if
              
            end do  !Enddo cnpptp
            
            do ij=1, Maxpptp
                k=npptp(Kiw,Kjw,ij)
                if (k<=0) exit
                
		    Ttyp  = TtypPair(k)(1:3)
		    Para = ParaPair(k, :)
		    YesFrc= .NOT. YesFitPair(k)
            
            Fudge=Fvdw
            if (Fudge<Reps) cycle
                Jatm = 3*IatmRef(j)
                JatmA = Jatm-IdspRef(MyID+1)
                    ! Compute if there is a match
!					YesCub=.false.   !Maybe eventually needed for Ewald? Removed for now
                    ! Store minimum distance

					RminPair(k)=min(RminPair(k), Rtmp)
					RmaxPair(k)=max(RmaxPair(k), Rtmp)

                   
					call PairFrcMat(YesFrc, Ttyp, XYZj, Rtmp, Para, Ai, Aj, Fi, Fj)

					if(YesFrc) then
						if(Iatm/=0) Ffix(Iatm-2:Iatm) = Ffix(Iatm-2:Iatm)+Fi*Fudge
						if(Jatm/=0) Ffix(Jatm-2:Jatm) = Ffix(Jatm-2:Jatm)+Fj*Fudge
                    else
                        Ilnr = IlnrPair(k)
			            Nlnr = NlnrPair(k)
			            Ntmp = Ilnr+Nlnr-1
						if(Iatm/=0) Amat(IatmA-2:IatmA, Ilnr:Ntmp) = Amat(IatmA-2:IatmA, Ilnr:Ntmp)+Ai(1:3,1:Nlnr)*Fudge
						if(Jatm/=0) Amat(JatmA-2:JatmA, Ilnr:Ntmp) = Amat(JatmA-2:JatmA, Ilnr:Ntmp)+Aj(1:3,1:Nlnr)*Fudge
                    end if
              
                        end do  !Enddo npptp
                        
            endif
                
            end do
        end do

End Subroutine InteFrcMat
!
!
Subroutine setPBC(Rij, XYZ, Abox, Bbox, Cbox, ABC, ABCinv)
	integer	i, j, k, Na, Nb, Nc
	real*8	Vol, Rtmp, Rij, Rmin, Abox, Bbox, Cbox, &
		&	XYZ(3), uvw(3), XYZmin(3), ABC(3,3), ABCinv(3,3)
! setPBC is being called at various places. 
! Not sure if it is simply computing distance RIJ
! I note XYZ vector is not intent in it is modified.
! Rij is the length of XYZ vector. 
! It looks to me the code makes two attemps to compute the XYZ vector to take care of image effect and get its length
! Performance problem since in many cases, we only need the distance. 

	uvw = matmul(ABCinv, XYZ)
	uvw = uvw - nint( uvw)
	XYZ = matmul(ABC, uvw)
	Rij = norm2(XYZ)

    ! Assuming cubic box? I believe volume has been computed. 
 	Vol=Abox*ABC(2,2)*ABC(3,3)

	Rtmp=Vol/sqrt( Bbox*Bbox*Cbox*Cbox - dot_product(ABC(:,2), ABC(:,3))**2 )
	Na=0; Rtmp=2.d0*Rij/Rtmp; if(Rtmp>1.d0) Na=ceiling(Rtmp)

	Rtmp=ABC(2,2)/sqrt(1.d0+(ABC(2,3)/ABC(3,3))**2)
	Nb=0; Rtmp=2.d0*Rij/Rtmp; if(Rtmp>1.d0) Nb=ceiling(Rtmp)

	Nc=0; Rtmp=2.d0*Rij/ABC(3,3); if(Rtmp>1.d0) Nc=ceiling(Rtmp)

	Rmin = Rij*Rij
	XYZmin = XYZ
	do i=-Na, Na
		do j=-Nb, Nb
			do k=-Nc, Nc
				XYZ = matmul(ABC, uvw+dble([i, j, k]))
				Rij = dot_product(XYZ, XYZ)
				if(Rij<Rmin) then
					Rmin   = Rij
					XYZmin = XYZ
				end if
			end do
		end do
	end do

	XYZ = XYZmin
	Rij = sqrt(Rmin)
End Subroutine setPBC
!
!
Subroutine NumbPara(Ttyp, Npar, Nlnr)
    integer, intent(out) :: Npar, Nlnr
    character(*) Ttyp
!set number of parameters and number of linear parameters for each interaction type
!The same infor is hard coded in FrcPairMat for performance reason    
    call upcase(Ttyp); if(Ttyp(1:3)/='BD3') Ttyp=Ttyp(1:3)

    Npar = 0; Nlnr = 0
    select case(trim(Ttyp)) !CBD type has  more than three characters
!		case('FUD'); Npar = 2; Nlnr = 0
            case('HAR'); Npar = 2; Nlnr = 2 !Harmonic Bond
            case('QUA'); Npar = 4; Nlnr = 4 !Quatic Bond
            case('COU'); Npar = 1; Nlnr = 1 !Coulombic
            case('THC'); Npar = 2; Nlnr = 1 !Thole Damped Coulombic
            case('GLJ');  Npar = 4; Nlnr = 2 !Generalized Lennard Jones 
            case('BUC'); Npar = 3; Nlnr = 2 !Buckingham.
            case('FDB'); Npar = 5; Nlnr = 2 !Fermi Damped Buckhingham. 
            case('STR'); Npar = 3; Nlnr = 1 !Shifted Truncated Power Law
            
            case('EXP'); Npar = 2; Nlnr = 1
            case('PEX'); Npar = 3; Nlnr = 1
            case('GEX'); Npar = 4; Nlnr = 3
            case('POW'); Npar = 2; Nlnr = 1 !Power Law
            case('TTP'); Npar = 3; Nlnr = 1
            case('SRD'); Npar = 3; Nlnr = 1 !short-range damped dispersion !ying
            case('FDP'); Npar = 4; Nlnr = 1 !Fermi damped power law
            case('COS'); Npar = 3; Nlnr = 2
            case('NCO'); Npar = 3; Nlnr = 1
            case('CNC'); Npar = 4; Nlnr = 1  !dihedral coupling, fixed delta
            case('CCO'); Npar = 4; Nlnr = 2  !dihedral coupling, fit delta
!            case('BND3B7');     Npar = 7; Nlnr = 7  
!            case('BND3B5');     Npar = 5; Nlnr = 5
!            case('BND3BB');     Npar = 3; Nlnr = 1
!            case('BND3BBBBA');  Npar = 6; Nlnr = 5
!            case('BND3B7QUA');  Npar = 7; Nlnr = 7
            case('MUB'); Npar = 4; Nlnr = 1 !bond-angle cross term Modifed Urey-Bradley
            case('QBB'); Npar = 5; Nlnr = 5 !bond-bond cross term, Quartic coupled
    end select
End Subroutine NumbPara
!
!
!  i------->j
!  E = E(Rij), Rij = sqrt[ (Xj-Xi)^2 + (Yj-Yi)^2 + (Zj-Zi)^2 ]
!  Fi = -gradi E(Rij) = -{ dE/dXi, dE/dYi, dE/dZi }
!     = -dE/dRij { dRij/dXi, dRij/dYi, dRij/dZi }
!	  = -dE/dRij { -(Xj-Xi)/Rij, -(Yj-Yi)/Rij, -(Zj-Zi)/Rij }
!	  =  dE/dRij { dX/Rij, dY/Dij, dZ/Rij }
!	  =  dE/dRij ^r0
!  Fj = -Fi
!
! {Fx, Fy, Fz} = dE/dR {dX, dY, dZ}/R
!
Subroutine BondFrcMat(YesFrc, Lab, dXYZ, Rij, P, Ai, Aj, Fi, Fj)
	real*8  dXYZ(3), P(*), Ai(3, *), Aj(3, *), Fi(3), Fj(3), Rij, Rsq, Rcub, Rtmp
	logical YesFrc
	character(*) Lab

	dXYZ = dXYZ/Rij
	select case(Lab)
		case('HAR')				! E = P2/2*(R-P1)^2
			if(YesFrc) then		! dE/dR = k*R-k*R0 = [P1]*R + [P2]; k=[P1], R0=-[P2]/[P1]
				Fi = P(2)*(Rij-P(1)) *dXYZ
				Fj = -Fi
			else
				Ai(1, 1:2) = dXYZ(1)*[ Rij, 1.d0 ]	! 1 [ dX*R dX ] [ P1 ]   [ Fx ]
				Ai(2, 1:2) = dXYZ(2)*[ Rij, 1.d0 ]  ! - [ dY*R dY ] |    | = [ Fy ]
				Ai(3, 1:2) = dXYZ(3)*[ Rij, 1.d0 ]  ! R [ dZ*R dZ ] [ P2 ]   [ Fz ]
				Aj(:, 1:2) = -Ai(:, 1:2)
			end if
		case('QUA') 			! E = P2/2*(R-P1)^2 + P3/3*(R-P1)^3 + P4/4*(R-P1)^4
			if(YesFrc) then		! dE/dR = k2*(R-R0) + k3*(R-R0)^2 + k4*(R-R0)^3
				Rtmp = Rij-P(1) ! dR
				Fi = ( P(2) + (P(3)+P(4)*Rtmp)*Rtmp )*Rtmp *dXYZ
				Fj = -Fi
			else				! dE/dR =   k4*R^3 + (k3-3*k4*R0)*R^2 + (k2-2*k3*R0+3*k4*R0^2)*R + (-k2*R0+k3*R0^2-k4*R0^3)
				Rsq  = Rij*Rij	! 		= [P1]*R^3 +      [P2]   *R^2 +           [P3]        *R +           [P4]
				Rcub = Rsq*Rij	! k4=P1, k3=P2+3*P1*R0, k2=P3+R0*(2*k3-3*k4*R0)=P3+R0*(2*P2+3*P1*R0)
								! P4 = -k2*R0+k3*R0^2-k4*R0^3 = -P1*R0^3 -P2*R0^2 -P3*R0
				Ai(1, 1:4) = dXYZ(1)*[ Rcub, Rsq, Rij, 1.d0 ]	! 1 [ dX*R^3 dX*R^2 dX*R dX ] [ P1 ]   [ Fx ]
				Ai(2, 1:4) = dXYZ(2)*[ Rcub, Rsq, Rij, 1.d0 ]   ! - [ dY*R^3 dY*R^2 dY*R dY ] [ P2 ] = [ Fy ]
				Ai(3, 1:4) = dXYZ(3)*[ Rcub, Rsq, Rij, 1.d0 ]   ! R [ dZ*R^3 dZ*R^2 dZ*R dZ ] [ P3 ]   [ Fz ]
				Aj(:, 1:4) = -Ai(:, 1:4)                        !                             [ P4 ]
			end if
	end select
End Subroutine BondFrcMat
!
!  i<---j--->k    Tht = Ang(i j k)
!  E = E(Tht), Rji*Rjk*cosTht = Dot(rji, rjk) = Xi*Xk+Yi*Yk+Zi*Zk
!  Rjk*gradi[Rji*cosTht] = gradi[Xi*Xk+Yi*Yk+Zi*Zk]
!  cosTht*gradi[Rji]-Rji*sinTht*gradi[Tht] = ^rjk
! -Rji*sinTht*gradi[Tht] = ^rjk - cosTht ^rji
!  gradi[Tht] = -(^rjk - cosTht ^rji)/(Rji*sinTht)

!  Fi = -gradi E(Tht) = -{ dE/dXi, dE/dYi, dE/dZi }
!     = -dE/dTht gradi[Tht]
!     = dE/dTht (^rjk - cosTht ^rji)/(Rji*sinTht)
!  Fk = dE/dTht (^rji - cosTht ^rjk)/(Rjk*sinTht)
!  Fj = -Fi-Fk
!
Subroutine AngleFrcMat(YesFrc, Lab, dXYZi, dXYZk, Rji, Rjk, P, Ai, Aj, Ak, Fi, Fj, Fk)
	real*8  dXYZi(3), dXYZk(3), P(*), Ai(3, *), Aj(3, *), Ak(3, *), Fi(3), Fj(3), Fk(3), &
	&		 Rji, Rjk, Dot, Tht, Tsq, Tcub, Rtmp
	logical YesFrc
	character(*) Lab

	dXYZi = dXYZi/Rji
	dXYZk = dXYZk/Rjk
	Dot = dot_product(dXYZi, dXYZk)
	if(Dot>1.d0) call ErrStop('Internal consistency problem in AngleFrcMat', 1)

	if (Dot<=-1.d0+10*epsilon(1.d0)) then  !if statement  considers 180 degree. --by ying
		Fi = 0.0
		Fk = 0.0
		Tht = 4.d0*atan(1.d0)
	else
		Tht = acos(Dot)
		Rtmp = 1.d0/sin(Tht)
		Fi = Rtmp*( dXYZk-Dot*dXYZi )/Rji
		Fk = Rtmp*( dXYZi-Dot*dXYZk )/Rjk
	endif

	select case(Lab)
	case('HAR') ! E = P2/2*(Tht-P1)^2
		if(YesFrc) then
			Rtmp = P(2)*(Tht-P(1))
			Fi = Rtmp*Fi
			Fk = Rtmp*Fk
			Fj = -(Fi+Fk)
		else
			Ai(1, 1:2) = Fi(1)*[ Tht, 1.d0 ]
			Ai(2, 1:2) = Fi(2)*[ Tht, 1.d0 ]
			Ai(3, 1:2) = Fi(3)*[ Tht, 1.d0 ]

			Ak(1, 1:2) = Fk(1)*[ Tht, 1.d0 ]
			Ak(2, 1:2) = Fk(2)*[ Tht, 1.d0 ]
			Ak(3, 1:2) = Fk(3)*[ Tht, 1.d0 ]

			Aj(:, 1:2) = -(Ai(:, 1:2)+Ak(:, 1:2))
		end if
	case('QUA') ! E = P2/2*(A-P1)^2 + P3/3*(A-P1)^3 + P4/4*(A-P1)^4
		if(YesFrc) then
			Rtmp = Tht-P(1)
			Rtmp = ( P(2) + (P(3)+P(4)*Rtmp)*Rtmp )*Rtmp
			Fi = Rtmp*Fi
			Fk = Rtmp*Fk
			Fj = -(Fi+Fk)
		else
			Tsq  = Tht*Tht
			Tcub = Tsq*Tht
			Ai(1, 1:4) = Fi(1)*[ Tcub, Tsq, Tht, 1.d0 ]
			Ai(2, 1:4) = Fi(2)*[ Tcub, Tsq, Tht, 1.d0 ]
			Ai(3, 1:4) = Fi(3)*[ Tcub, Tsq, Tht, 1.d0 ]

			Ak(1, 1:4) = Fk(1)*[ Tcub, Tsq, Tht, 1.d0 ]
			Ak(2, 1:4) = Fk(2)*[ Tcub, Tsq, Tht, 1.d0 ]
			Ak(3, 1:4) = Fk(3)*[ Tcub, Tsq, Tht, 1.d0 ]

			Aj(:, 1:4) = -(Ai(:, 1:4)+Ak(:, 1:4))
		end if
	end select
End Subroutine AngleFrcMat
!
!
Subroutine DihdrFrcMat(YesFrc, Lab, dXYZj, dXYZk, dXYZl, Rk, P, Ai, Aj, Ak, Al, Fi, Fj, Fk, Fl)
	real*8 Rtmp, Rk, Phi, Vijk(3), Vjkl(3), &
		&  dXYZj(3), dXYZk(3), dXYZl(3), P(*), &
		&  Fi(3), Fj(3), Fk(3), Fl(3), &
		&  Ai(3, *), Aj(3, *), Ak(3, *), Al(3, *), P12(2)
	logical YesFrc
	character(*) Lab

	call cross(dXYZj, dXYZk, Vijk)
	call cross(dXYZk, dXYZl, Vjkl)
	Phi=atan2(Rk*dot_product(dXYZj, Vjkl), dot_product(Vijk, Vjkl))

	Rtmp=Rk*Rk
	Fi= Rk*Vijk/dot_product(Vijk, Vijk)
	Fl=-Rk*Vjkl/dot_product(Vjkl, Vjkl)
	Fj = -( (Rtmp+dot_product(dXYZj, dXYZk))*Fi-dot_product(dXYZk, dXYZl)*Fl )/Rtmp

	select case(Lab(1:3))
	case('NCO') ! E = P(1)*(1+cos(P(2)*Phi-P(3)))
		if(YesFrc) then
			Rtmp = -P(1)*P(2)*sin(P(2)*Phi-P(3))
			Fi = Rtmp*Fi
			Fj = Rtmp*Fj
			Fl = Rtmp*Fl
			Fk = -(Fi+Fj+Fl)
		else
			Rtmp = -P(2)*sin(P(2)*Phi-P(3))
			Ai(:, 1) = Fi(:)*Rtmp
			Al(:, 1) = Fl(:)*Rtmp
			Aj(:, 1) = Fj(:)*Rtmp
			Ak(:, 1) = -(Ai(:, 1)+Aj(:, 1)+Al(:, 1))
		end if
	case('COS')          ! E = P2*(1+cos(P3*Phi-P1))
		if(YesFrc) then  ! dE/dPhi=-P2*P3*sin(P3*Phi-P1) = [-P2*cos(P1)]*P3*sin(P3*Phi) + [P2*sin(P1)]*P3*cos(P3*Phi)
			Rtmp = -P(2)*P(3)*sin(P(3)*Phi-P(1))
			Fi = Rtmp*Fi
			Fj = Rtmp*Fj
			Fl = Rtmp*Fl
			Fk = -(Fi+Fj+Fl)
		else
			P12=[P(3)*sin(P(3)*Phi), P(3)*cos(P(3)*Phi)]

			Ai(1, 1:2) = Fi(1)*P12
			Ai(2, 1:2) = Fi(2)*P12
			Ai(3, 1:2) = Fi(3)*P12

			Al(1, 1:2) = Fl(1)*P12
			Al(2, 1:2) = Fl(2)*P12
			Al(3, 1:2) = Fl(3)*P12

			Aj(1, 1:2) = Fj(1)*P12
			Aj(2, 1:2) = Fj(2)*P12
			Aj(3, 1:2) = Fj(3)*P12

			Ak(:, 1:2) = -(Ai(:, 1:2)+Aj(:, 1:2)+Al(:, 1:2))
		end if
	case('HAR') ! E = P2/2*(Phi-P1)^2
		if(YesFrc) then
			Rtmp = P(2)*(Phi-P(1))
			Fi = Rtmp*Fi
			Fk = Rtmp*Fk
			Fl = Rtmp*Fl
			Fj = -(Fi+Fk+Fl)
		else
			Ai(1, 1:2) = Fi(1)*[ Phi, 1.d0 ]
			Ai(2, 1:2) = Fi(2)*[ Phi, 1.d0 ]
			Ai(3, 1:2) = Fi(3)*[ Phi, 1.d0 ]

			Al(1, 1:2) = Fl(1)*[ Phi, 1.d0 ]
			Al(2, 1:2) = Fl(2)*[ Phi, 1.d0 ]
			Al(3, 1:2) = Fl(3)*[ Phi, 1.d0 ]

			Ak(1, 1:2) = Fk(1)*[ Phi, 1.d0 ]
			Ak(2, 1:2) = Fk(2)*[ Phi, 1.d0 ]
			Ak(3, 1:2) = Fk(3)*[ Phi, 1.d0 ]

			Aj(:, 1:2) = -(Ai(:, 1:2)+Al(:, 1:2)+Ak(:, 1:2))
		end if
	end select
End Subroutine DihdrFrcMat
!
!
Subroutine CDihdrFrcMat(YesFrc, Lab, dXYZj, dXYZk, dXYZl, dXYZr, Rk1, Rk2, P, Ai, Aj, Ak, Al, Ar, Fi, Fj, Fk, Fl, Fr)
	real*8 Rtmp1, Rtmp2, Rtemp3, Rtemp4, Rk1, Rk2, Phi, Psi, Vijk(3), Vjkl(3), Vklr(3), &
          &  dXYZj(3), dXYZk(3), dXYZl(3), dXYZr(3), P(*), &
          &  Fi1(3), Fj1(3), Fk1(3), Fl1(3), Fj2(3), Fk2(3), Fl2(3), Fr2(3), Fk11(3), &
          &  Fi(3), Fj(3), Fk(3), Fl(3), Fr(3), &
          &  Ai(3, *), Aj(3, *), Ak(3, *), Al(3, *), Ar(3, *), P12(2)
	logical YesFrc
        character(*) Lab

        call cross(dXYZj, dXYZk, Vijk)
        call cross(dXYZk, dXYZl, Vjkl)
        call cross(dXYZl, dXYZr, Vklr)
    
        Phi=atan2(Rk1*dot_product(dXYZj, Vjkl), dot_product(Vijk, Vjkl))
        Psi=atan2(Rk2*dot_product(dXYZk, Vklr), dot_product(Vjkl, Vklr))
    
        Rtmp1=Rk1*Rk1 ! dphi/dri + dphi/drj + dphi/drk + dphi/drl = 0 (i-j-k-l)
        Fi1= Rk1*Vijk/dot_product(Vijk, Vijk) !-dphi/dri
        Fl1=-Rk1*Vjkl/dot_product(Vjkl, Vjkl) 
        Fj1=-( (Rtmp1+dot_product(dXYZj, dXYZk))*Fi1-dot_product(dXYZk, dXYZl)*Fl1 )/Rtmp1
        Fk1=-( (Rtmp1+dot_product(dXYZk, dXYZl))*Fl1-dot_product(dXYZj, dXYZk)*Fi1 )/Rtmp1
!        Fk11=-(Fi1+Fj1+Fl1) !for debug
    
        Rtmp2=Rk2*Rk2 ! dpsi/drj + dpsi/drk + dpsi/drl + dpsi/drr = 0 (j-k-l-r)
        Fj2= Rk2*Vjkl/dot_product(Vjkl, Vjkl)
        Fr2=-Rk2*Vklr/dot_product(Vklr, Vklr) 
        Fk2=-( (Rtmp2+dot_product(dXYZk, dXYZl))*Fj2-dot_product(dXYZl, dXYZr)*Fr2 )/Rtmp2
        Fl2=-( (Rtmp2+dot_product(dXYZl, dXYZr))*Fr2-dot_product(dXYZk, dXYZl)*Fj2 )/Rtmp2
!        Fk11=-(Fj2+Fr2+Fl2) !for debug
    
        select case(Lab(1:3))
        case('CNC') ! E = P(1)*(1+cos(P(2)*Phi+P(3)*Psi-P(4)))
!            YesFrc=.True. 
            if(YesFrc) then !dE/dr=-P(1)*sin(P(2)*Phi+P(3)*Psi-P(4))*(P(2)*dphi/dr+P(3)*dPsi/dr)
                Rtmp1 = -P(1)*P(2)*sin(P(2)*Phi+P(3)*Psi-P(4))
                Rtmp2 = -P(1)*P(3)*sin(P(2)*Phi+P(3)*Psi-P(4))
                Fi = Rtmp1*Fi1
                Fj = Rtmp1*Fj1+Rtmp2*Fj2
                Fl = Rtmp1*Fl1+Rtmp2*Fl2
                Fk = Rtmp1*Fk1+Rtmp2*Fk2
                Fr = Rtmp2*Fr2
            else
                Rtmp1 = -P(2)*sin(P(2)*Phi+P(3)*Psi-P(4))
                Rtmp2 = -P(3)*sin(P(2)*Phi+P(3)*Psi-P(4))
                Ai(:, 1) = Rtmp1*Fi1(:)
                Al(:, 1) = Rtmp1*Fl1(:)+Rtmp2*Fl2(:)
                Aj(:, 1) = Rtmp1*Fj1(:)+Rtmp2*Fj2(:)
                Ak(:, 1) = Rtmp1*Fk1(:)+Rtmp2*Fk2(:)
                Ar(:, 1) = Rtmp2*Fr2(:)
            end if
        case('CCO') ! E = P2*(1+cos(P3*Phi+P4*Psi-P1))
            if(YesFrc) then !dE/dr=-P2*sin(P3*Phi+P4*Psi-P1)*(P3*dphi/dr+P4*dPsi/dr)
                            !     =-P2*[sin(P3*Phi+P4*Psi)cos(P1)-cos(P3*Phi+P4*Psi)sin(P1)]*(P3*dphi/dr+P4*dPsi/dr)
                            !     =-P2*cos(P1)*sin(P3*Phi+P4*Psi)*(P3*dphi/dr+P4*dPsi/dr)+P2*sin(P1)*cos(P3*Phi+P4*Psi)*(P3*dphi/dr+P4*dPsi/dr)
                            !     =P(1)'*sin(P3*Phi+P4*Psi)*(P3*dphi/dr+P4*dPsi/dr)
                            !     +P(2)'*cos(P3*Phi+P4*Psi)*(P3*dphi/dr+P4*dPsi/dr)
							!P(1)'=-P2*cos(P1)  
							!P(2)'=P2*sin(P1)
                Rtmp1 = -P(2)*sin(P(3)*Phi+P(4)*Psi-P(1))*P(3)
                Rtmp2 = -P(2)*sin(P(3)*Phi+P(4)*Psi-P(1))*P(4)
                Fi = Rtmp1*Fi1
                Fj = Rtmp1*Fj1+Rtmp2*Fj2
                Fl = Rtmp1*Fl1+Rtmp2*Fl2
                Fk = Rtmp1*Fk1+Rtmp2*Fk2
                Fr = Rtmp2*Fr2
            else
                Rtmp1 = sin(P(3)*Phi+P(4)*Psi)*P(3)
                Rtmp2 = sin(P(3)*Phi+P(4)*Psi)*P(4)
                Rtmp3 = cos(P(3)*Phi+P(4)*Psi)*P(3)
                Rtmp4 = cos(P(3)*Phi+P(4)*Psi)*P(4)

			    Ai(1, 1:2) = [Rtmp1*Fi1(1), Rtmp3*Fi1(1)]
			    Ai(2, 1:2) = [Rtmp1*Fi1(2), Rtmp3*Fi1(2)]
			    Ai(3, 1:2) = [Rtmp1*Fi1(3), Rtmp3*Fi1(3)]
			    
                Aj(1, 1:2) = [Rtmp1*Fj1(1)+Rtmp2*Fj2(1), Rtmp3*Fj1(1)+Rtmp4*Fj2(1)]
			    Aj(2, 1:2) = [Rtmp1*Fj1(2)+Rtmp2*Fj2(2), Rtmp3*Fj1(2)+Rtmp4*Fj2(2)]
			    Aj(3, 1:2) = [Rtmp1*Fj1(3)+Rtmp2*Fj2(3), Rtmp3*Fj1(3)+Rtmp4*Fj2(3)]
			    
			    Ak(1, 1:2) = [Rtmp1*Fk1(1)+Rtmp2*Fk2(1), Rtmp3*Fk1(1)+Rtmp4*Fk2(1)]
			    Ak(2, 1:2) = [Rtmp1*Fk1(2)+Rtmp2*Fk2(2), Rtmp3*Fk1(2)+Rtmp4*Fk2(2)]
			    Ak(3, 1:2) = [Rtmp1*Fk1(3)+Rtmp2*Fk2(3), Rtmp3*Fk1(3)+Rtmp4*Fk2(3)]
                
			    Al(1, 1:2) = [Rtmp1*Fl1(1)+Rtmp2*Fl2(1), Rtmp3*Fl1(1)+Rtmp4*Fl2(1)]
			    Al(2, 1:2) = [Rtmp1*Fl1(2)+Rtmp2*Fl2(2), Rtmp3*Fl1(2)+Rtmp4*Fl2(2)]
			    Al(3, 1:2) = [Rtmp1*Fl1(3)+Rtmp2*Fl2(3), Rtmp3*Fl1(3)+Rtmp4*Fl2(3)]
			    
			    Ar(1, 1:2) = [Rtmp2*Fr2(1), Rtmp4*Fr2(1)]
			    Ar(2, 1:2) = [Rtmp2*Fr2(2), Rtmp4*Fr2(2)]
			    Ar(3, 1:2) = [Rtmp2*Fr2(3), Rtmp4*Fr2(3)]            
            end if
            
        end select
End Subroutine CDihdrFrcMat    
!    
!    
Subroutine PairFrcMat(YesFrc, Lab, dXYZ, R, P, Ai, Aj, Fi, Fj)
    implicit none
	integer i, Npar, Nlnr
    real*8, intent(in) :: dXYZ(3) !unit vector in the direction of R
	real*8, parameter:: FcovCou = 332.063713741257D0
	real*8  P(*), E, R, Rtmp, Ptmp, igamma, &
	&		Ai(3,*), Aj(3,*), Fi(3), Fj(3)
	logical YesFrc
	character*3 Lab

!	P1=P(1); P2=P(2); P3=P(3); P4=P(4)

	select case(Lab)
    case('COU')
        Npar = 1; Nlnr = 1
! 		E = P1/R
		Ai(:, 1) = -FcovCou/(R*R)
    case('THC') ! Thole
        Npar = 2; Nlnr = 1
		Rtmp=P(2)*R**3
! 		E = P(1)/R * (1.d0-exp(-P(2)*R**3)+P(2)**(1.d0/3.d0)*igamma(2.d0/3.d0, P(2)*R**3)*R)
		Ai(:, 1) = -FcovCou/(R*R) * (1.d0-exp(-P(2)*R*R*R))
    case('GLJ')
        Npar = 4; Nlnr = 2
! 		E = P(1)*R**(-P(3)) + P(2)*R**(-P(4))
		Ai(:, 1) = -P(3)*R**(-P(3)-1.d0)
		Ai(:, 2) = -P(4)*R**(-P(4)-1.d0)
    case('BUC')
        Npar = 3; Nlnr = 2
! 		E = P(1)*exp(-P(3)*R) + P(2)/R**6
		Ai(:, 1) = -P(3)*exp(-P(3)*R)
		Ai(:, 2) = -6.d0*R**(-7.d0)
    case('FDB')
        Npar = 5; Nlnr = 2
		Rtmp = 1.d0/( 1.d0+exp(-P(4)*(R-P(5))) )
! 		E = P(1)*exp(-P(3)*R) + P(2)*Rtmp/R**6
		Ai(:, 1) = -P(3)*exp(-P(3)*R)
		Ai(:, 2) =  Rtmp*( P(4)*R*(1.d0-Rtmp)-6.d0 )*R**(-7.d0)
    case('STR')
        Npar = 3; Nlnr = 1
! 		E = P(1)*( 1.d0/R**P(2) - 1.d0/P(3)**P2 + P(2)*(R-P(3))/P(3)**(P(2)+1) )
		Ai(:, 1) = 0.d0
		if(R<P(3)) then
			Ptmp = -P(2)-1.d0
			Ai(:, 1) = P(2)*(P(3)**Ptmp-R**Ptmp)
		end if
    case('EXP')
        Npar = 2; Nlnr = 1
! 		E = P1 * exp(-P2*R)
		Ai(:, 1) = -P(2)*exp(-P(2)*R)
    case('PEX')
        Npar = 3; Nlnr = 1
! 		E = P1 * R**P2 * exp(-P3*R)
		Ai(:, 1) = (P(2)-P(3)*R)*R**(P(2)-1.d0)*exp(-P(3)*R)

    case('GEX')
        Npar = 4; Nlnr = 3
!  		E = (P1 + P2*R +P3*R**2) * exp(-P4*R)
		Rtmp=exp(-P(4)*R)
		Ai(:, 1) = -P(4)*Rtmp
		Ai(:, 2) = (1.d0-P(4)*R)  *Rtmp
		Ai(:, 3) = (2.d0-P(4)*R)*R*Rtmp

    case('POW')
        Npar = 2; Nlnr = 1
! 		E = P1*R**P2
		Ai(:, 1) = P(2)*R**(P(2)-1.d0)
    case('TTP') ! Tang-Toennies
        Npar = 3; Nlnr = 1
		Rtmp=P(3)*R
		Ptmp=-P(2)+1.d0
! 		E = P1*R**P2*(1.d0-igamma(Ptmp, Rtmp)/gamma(Ptmp))
		Ai(:, 1) = R**(-Ptmp) /gamma(Ptmp) &
		&        * ( P(2)*(gamma(Ptmp)-igamma(Ptmp, Rtmp)) + Rtmp**Ptmp*exp(-Rtmp) )
    case('SRD') !by ying
!       E = P1/(R**P2+P3**P2)
        Npar = 3; Nlnr = 1
        Rtmp=P(3)**abs(P(2))
        Ai(:,1) = -(R**abs(P(2))+Rtmp)**(-2)*abs(P(2))*R**(abs(P(2))-1)
    case('FDP')
        Npar = 4; Nlnr = 1
		Rtmp = 1.d0/( 1.d0+exp(-P(3)*(R-P(4))) )
! 		E = P1*R**P2*Rtmp
		Ai(:, 1) = Rtmp*( P(2)+P(3)*R*(1.d0-Rtmp) )*R**(P(2)-1.d0)
    case default
        call ErrStop("unsupported interaction in PairFrcMat",1)
    
    end select
    
    !    call NumbPara(Lab, Npar, Nlnr) uncomment to debug number of parameter mismatch

	 ! F = ^R * dE/dR = ^dXYZ * dE/dR
	Ai(1, 1:Nlnr) =  Ai(1, 1:Nlnr)*dXYZ(1)
	Ai(2, 1:Nlnr) =  Ai(2, 1:Nlnr)*dXYZ(2)
	Ai(3, 1:Nlnr) =  Ai(3, 1:Nlnr)*dXYZ(3)
	Aj(:, 1:Nlnr) = -Ai(:, 1:Nlnr)

	if(YesFrc) then
		Fi=0.d0
		do i=1, Nlnr
			Fi = Fi+P(i)*Ai(:, i)
		end do
		Fj = -Fi
	end if

End Subroutine PairFrcMat
!
!corss term
Subroutine Bond3bFrcMat(YesFrc, Lab, dXYZij, dXYZkj, dXYZik, Rji, Rjk, Rik, Ppot, Ai, Aj, Ak, Fi, Fj, Fk)
    implicit none
    real*8  Rsq, Rcub
    real*8  dXYZij(3), dXYZkj(3), dXYZik(3), Ppot(*), Ai(3, *), Aj(3, *), Ak(3, *), Fi(3), Fj(3), Fk(3), &
    	&	Rji, Rjk, Rik, Dot, Tht, Tsq, Tcub, Rtmp
    real*8 Rji2, Rji3, Rjk2, Rjk3
    real*8 P(20)
    logical YesFrc
    character(*) Lab
    
    !note dxyzij=\vrj-\vri
    dXYZij = dXYZij/Rji
    dXYZkj = dXYZkj/Rjk
    dXYZik = dXYZik/Rik
    Rji2=Rji*Rji
    Rji3=Rji2*Rji
    Rjk2=Rjk*Rjk
    Rjk3=Rjk2*Rjk
    
    select case(Lab)
    !nonlinear fitting: Bond-Angle cross term
    !P(1)=-krtheta
    !P(2)=2krtheta r1e^\prime
    !P(3)=krtheta r3e^\prime
    case('MUB') !bond-angle cross term E=P(1)*(Rik-P(4))*(Rij-P(2)+Rjk-P(3))
        P(1)=Ppot(1) !Ppot(1) krtheta
        P(2)=Ppot(2) !Ppot(2) r1e for rij
        P(3)=Ppot(3) !Ppot(3) r2e for rjk
        P(4)=Ppot(4) !Ppot(4) r3e for rik
        if(YesFrc) then
            Fi = -P(1)*((Rik-P(4))*dXYZij+(Rji+Rjk-P(2)-P(3))*dXYZik)
            Fk = -P(1)*((Rik-P(4))*dXYZkj-(Rji+Rjk-P(2)-P(3))*dXYZik)
            Fj = -Fi-Fk
        else
            Ai(:,1) = [ -(Rik-P(4))*dXYZij-(Rji+Rjk-P(2)-P(3))*dXYZik ]
            Ak(:,1) = [ -(Rik-P(4))*dXYZkj+(Rji+Rjk-P(2)-P(3))*dXYZik ]
            Aj(:,1) = -Ai(:,1) - Ak(:,1)
        end if
        
        case('QBB') !bond-bond cross term and two identical quartic bonds 
                       !E= Ppot3/2*(r1-Ppot1)^2 + Ppot3/3*(r1-Ppot1)^3 + Ppot5/4*(r1-Ppot1)^4
                       !  +Ppot3/2*(r2-Ppot1)^2 + Ppot3/3*(r2-Ppot1)^3 + Ppot5/4*(r2-Ppot1)^4
                       !  +Ppot2*(r1-Ppot1)*(r2-Ppot1)
        P(1)=-Ppot(3)+(2.*Ppot(4)-3.*Ppot(5)*Ppot(1))*Ppot(1)
        P(2)=-Ppot(2)
        P(3)=(Ppot(3)-(Ppot(4)-Ppot(5)*Ppot(1))*Ppot(1)+Ppot(2))*Ppot(1)
        P(4)=-Ppot(4)+3.*Ppot(5)*Ppot(1)
        P(5)=-Ppot(5)
        if(YesFrc) then
            Fi = (((P(5)*Rji+P(4))*Rji+P(1))*Rji+P(2)*Rjk+P(3)) *dXYZij
            Fk = (((P(5)*Rjk+P(4))*Rjk+P(1))*Rjk+P(2)*Rji+P(3))*dXYZkj
            Fj = -Fi-Fk
        else
            Ai(1, 1:5) = [Rji*dXYZij(1), Rjk*dXYZij(1), dXYZij(1), Rji2*dXYZij(1), Rji3*dXYZij(1)]
            Ai(2, 1:5) = [Rji*dXYZij(2), Rjk*dXYZij(2), dXYZij(2), Rji2*dXYZij(2), Rji3*dXYZij(2)]
            Ai(3, 1:5) = [Rji*dXYZij(3), Rjk*dXYZij(3), dXYZij(3), Rji2*dXYZij(3), Rji3*dXYZij(3)]
            Ak(1, 1:5) =[Rjk*dXYZkj(1), Rji*dXYZkj(1), dXYZkj(1), Rjk2*dXYZkj(1), Rjk3*dXYZkj(1)]
            Ak(2, 1:5) =[Rjk*dXYZkj(2), Rji*dXYZkj(2), dXYZkj(2), Rjk2*dXYZkj(2), Rjk3*dXYZkj(2)]
            Ak(3, 1:5) =[Rjk*dXYZkj(3), Rji*dXYZkj(3), dXYZkj(3), Rjk2*dXYZkj(3), Rjk3*dXYZkj(3)]
            Aj(:, 1:5) = -Ai(:, 1:5) -Ak(:, 1:5)
        end if
        
!    case('BND3BB')
!        P(1)=Ppot(1) !Ppot(1) krtheta
!        P(2)=Ppot(2) !Ppot(2) r1e
!        P(3)=Ppot(3)
!        if(YesFrc) then
!            Fi = -P(1)*(Rjk-P(3))*dXYZij
!            Fk = -P(1)*(Rji-P(2))*dXYZkj
!            Fj = -Fi-Fk
!        else
!            Ai(:,1) = [ -(Rjk-P(3))*dXYZij ]
!            Ak(:,1) = [ -(Rji-P(2))*dXYZkj ]
!            Aj(:,1) = -Ai(:,1) - Ak(:,1)
!        end if
!    case('BND3BBBBA')
!        P(1)=Ppot(1)
!        P(2)=Ppot(2)
!        P(3)=Ppot(3)
!        P(4)=Ppot(4)
!        P(5)=Ppot(5)
!        P(6)=Ppot(6)
!        P(7)=Rji-P(6)
!        P(8)=Rjk-P(6)
!        if(YesFrc) then
!            Fi =-(P(1)*P(7)+P(2)*P(7)**2+P(3)*P(7)**3+P(4)*(Rjk-P(6))+P(5)*(Rik-2*P(6)))*dXYZij-P(5)*(Rji+Rjk-2*P(6))*dXYZik
!            Fk =-(P(1)*P(8)+P(2)*P(8)**2+P(3)*P(8)**3+P(4)*(Rji-P(6))+P(5)*(Rik-2*P(6)))*dXYZkj+P(5)*(Rji+Rjk-2*P(6))*dXYZik
!            Fj = -Fi-Fk
!        else
!            Ai(1, 1:5) = [-P(7)*dXYZij(1),-P(7)**2*dXYZij(1), -P(7)**3*dXYZij(1),-(Rjk-P(6))*dXYZij(1), -(Rik-2*P(6))*dXYZij(1)-(Rji+Rjk-2*P(6))*dXYZik(1)]
!            Ai(2, 1:5) = [-P(7)*dXYZij(2),-P(7)**2*dXYZij(2), -P(7)**3*dXYZij(2),-(Rjk-P(6))*dXYZij(2), -(Rik-2*P(6))*dXYZij(2)-(Rji+Rjk-2*P(6))*dXYZik(2)]
!            Ai(3, 1:5) = [-P(7)*dXYZij(3),-P(7)**2*dXYZij(3), -P(7)**3*dXYZij(3),-(Rjk-P(6))*dXYZij(3), -(Rik-2*P(6))*dXYZij(3)-(Rji+Rjk-2*P(6))*dXYZik(3)]
!            Ak(1, 1:5) =[-P(8)*dXYZkj(1),-P(8)**2*dXYZkj(1), -P(8)**3*dXYZkj(1),-(Rji-P(6))*dXYZkj(1), -(Rik-2*P(6))*dXYZkj(1)+(Rji+Rjk-2*P(6))*dXYZik(1)]
!            Ak(2, 1:5) =[-P(8)*dXYZkj(2),-P(8)**2*dXYZkj(2), -P(8)**3*dXYZkj(2),-(Rji-P(6))*dXYZkj(2), -(Rik-2*P(6))*dXYZkj(2)+(Rji+Rjk-2*P(6))*dXYZik(2)]
!            Ak(3, 1:5) =[-P(8)*dXYZkj(3),-P(8)**2*dXYZkj(3), -P(8)**3*dXYZkj(3),-(Rji-P(6))*dXYZkj(3), -(Rik-2*P(6))*dXYZkj(3)+(Rji+Rjk-2*P(6))*dXYZik(3)]
!            Aj(:, 1:5) = -Ai(:, 1:5) -Ak(:, 1:5)
!        end if

!    case('BND3B7QUA')
!        P(1)=-Ppot(5)+(2.*Ppot(6)-3.*Ppot(7)*Ppot(1))*Ppot(1)
!        P(2)=-Ppot(2)
!        P(3)=-Ppot(3)
!        P(4)=2.*Ppot(3)*Ppot(1)
!        P(5)=(Ppot(5)-(Ppot(6)-Ppot(7)*Ppot(1))*Ppot(1)+Ppot(2))*Ppot(1)+Ppot(3)*Ppot(4)
!        P(6)=-Ppot(6)+3.*Ppot(7)*Ppot(1)
!        P(7)=-Ppot(7)
!        if(YesFrc) then
!            Fi = (((P(7)*Rji+P(6))*Rji+P(1))*Rji+P(2)*Rjk+P(3)*Rik+P(5)) *dXYZij
!            Fi = Fi+(P(3)*(Rji+Rjk)+P(4))*dXYZik
!            Fk = (((P(7)*Rjk+P(6))*Rjk+P(1))*Rjk+P(2)*Rji+P(3)*Rik+P(5))*dXYZkj
!            Fk = Fk-(P(3)*(Rjk+Rji)+P(4))*dXYZik
!            Fj = -Fi-Fk
!        else
!            Ai(1, 1:7) = [Rji*dXYZij(1), Rjk*dXYZij(1), Rik*dXYZij(1)+(Rji+Rjk)*dXYZik(1), dXYZik(1), dXYZij(1), Rji2*dXYZij(1), Rji3*dXYZij(1)]
!            Ai(2, 1:7) = [Rji*dXYZij(2), Rjk*dXYZij(2), Rik*dXYZij(2)+(Rji+Rjk)*dXYZik(2), dXYZik(2), dXYZij(2), Rji2*dXYZij(2), Rji3*dXYZij(2)]
!            Ai(3, 1:7) = [Rji*dXYZij(3), Rjk*dXYZij(3), Rik*dXYZij(3)+(Rji+Rjk)*dXYZik(3), dXYZik(3), dXYZij(3), Rji2*dXYZij(3), Rji3*dXYZij(3)]
!            Ak(1, 1:7) =[Rjk*dXYZkj(1), Rji*dXYZkj(1), Rik*dXYZkj(1)-(Rjk+Rji)*dXYZik(1), -dXYZik(1), dXYZkj(1), Rjk2*dXYZkj(1), Rjk3*dXYZkj(1)]
!            Ak(2, 1:7) =[Rjk*dXYZkj(2), Rji*dXYZkj(2), Rik*dXYZkj(2)-(Rjk+Rji)*dXYZik(2), -dXYZik(2), dXYZkj(2), Rjk2*dXYZkj(2), Rjk3*dXYZkj(2)]
!            Ak(3, 1:7) =[Rjk*dXYZkj(3), Rji*dXYZkj(3), Rik*dXYZkj(3)-(Rjk+Rji)*dXYZik(3), -dXYZik(3), dXYZkj(3), Rjk2*dXYZkj(3), Rjk3*dXYZkj(3)]
!            Aj(:, 1:7) = -Ai(:, 1:7) -Ak(:, 1:7)
!        end if
!    case('BND3B7') !two different harmonic bonds and two cross terms
!        P(1)=-Ppot(1) !-k1b
!        P(2)=-Ppot(2) !-k2b
!        P(3)=-Ppot(3) !-krr
!        P(4)=-Ppot(4) !-krtheta
!        P(5)=Ppot(4)*(Ppot(5)+Ppot(6)) !krtheta*(r1e+r2e)
!        P(6)=Ppot(1)*Ppot(5)+Ppot(3)*Ppot(6)+Ppot(4)*Ppot(7) !Ppot(7) is r3e
!        P(7)=Ppot(2)*Ppot(6)+Ppot(3)*Ppot(5)+Ppot(4)*Ppot(7)
!        if(YesFrc) then
!            Fi = (P(1)*Rji+P(3)*Rjk+P(4)*Rik+P(6)) *dXYZij
!           Fi = Fi+(P(4)*(Rji+Rjk)+P(5))*dXYZik
!           Fk = (P(2)*Rjk+P(3)*Rji+P(4)*Rik+P(7))*dXYZkj
!            Fk = Fk-(P(4)*(Rjk+Rji)+P(5))*dXYZik
!            Fj = -Fi-Fk
!        else
!            Ai(1, 1:7) = [Rji*dXYZij(1), 0.D0, Rjk*dXYZij(1), Rik*dXYZij(1)+(Rji+Rjk)*dXYZik(1), dXYZik(1), dXYZij(1), 0.D0]
!            Ai(2, 1:7) = [Rji*dXYZij(2), 0.D0, Rjk*dXYZij(2), Rik*dXYZij(2)+(Rji+Rjk)*dXYZik(2), dXYZik(2), dXYZij(2), 0.D0]
!            Ai(3, 1:7) = [Rji*dXYZij(3), 0.D0, Rjk*dXYZij(3), Rik*dXYZij(3)+(Rji+Rjk)*dXYZik(3), dXYZik(3),  dXYZij(3), 0.D0]
!            Ak(1, 1:7) =[0.D0, Rjk*dXYZkj(1), Rji*dXYZkj(1), Rik*dXYZkj(1)-(Rjk+Rji)*dXYZik(1), -dXYZik(1), 0.D0, dXYZkj(1)]
!            Ak(2, 1:7) =[0.D0, Rjk*dXYZkj(2), Rji*dXYZkj(2), Rik*dXYZkj(2)-(Rjk+Rji)*dXYZik(2), -dXYZik(2), 0.D0, dXYZkj(2)]
!            Ak(3, 1:7) =[0.D0, Rjk*dXYZkj(3), Rji*dXYZkj(3), Rik*dXYZkj(3)-(Rjk+Rji)*dXYZik(3), -dXYZik(3), 0.D0, dXYZkj(3)]
!            Aj(:, 1:7) = -Ai(:, 1:7) -Ak(:, 1:7)
!        end if
!    case('BND3B5') !two identical harmonic bonds and two cross terms
!        P(1)=-Ppot(1) !-k1b
!        P(2)=-Ppot(2) !-krr
!        P(3)=-Ppot(3) !-krtheta
!        P(4)=2.*Ppot(3)*Ppot(4) !2krtheta*r1e
!        P(5)=Ppot(1)*Ppot(4)+Ppot(2)*Ppot(4)+Ppot(3)*Ppot(5) !Ppot(5) is r3e
!        if(YesFrc) then
!            Fi = (P(1)*Rji+P(2)*Rjk+P(3)*Rik+P(5)) *dXYZij
!            Fi = Fi+(P(3)*(Rji+Rjk)+P(4))*dXYZik
!            Fk = (P(1)*Rjk+P(2)*Rji+P(3)*Rik+P(5))*dXYZkj
!            Fk = Fk-(P(3)*(Rjk+Rji)+P(4))*dXYZik
!            Fj = -Fi-Fk
!        else
!            Ai(1, 1:5) = [Rji*dXYZij(1), Rjk*dXYZij(1), Rik*dXYZij(1)+(Rji+Rjk)*dXYZik(1), dXYZik(1), dXYZij(1)]
!            Ai(2, 1:5) = [Rji*dXYZij(2), Rjk*dXYZij(2), Rik*dXYZij(2)+(Rji+Rjk)*dXYZik(2), dXYZik(2), dXYZij(2)]
!            Ai(3, 1:5) = [Rji*dXYZij(3), Rjk*dXYZij(3), Rik*dXYZij(3)+(Rji+Rjk)*dXYZik(3), dXYZik(3),  dXYZij(3)]
!            Ak(1, 1:5) =[Rjk*dXYZkj(1), Rji*dXYZkj(1), Rik*dXYZkj(1)-(Rjk+Rji)*dXYZik(1), -dXYZik(1), dXYZkj(1)]
!            Ak(2, 1:5) =[Rjk*dXYZkj(2), Rji*dXYZkj(2), Rik*dXYZkj(2)-(Rjk+Rji)*dXYZik(2), -dXYZik(2), dXYZkj(2)]
!            Ak(3, 1:5) =[Rjk*dXYZkj(3), Rji*dXYZkj(3), Rik*dXYZkj(3)-(Rjk+Rji)*dXYZik(3), -dXYZik(3), dXYZkj(3)]
!            Aj(:, 1:5) = -Ai(:, 1:5) -Ak(:, 1:5)
!        end if
        
    end select
End Subroutine Bond3bFrcMat
!
!
Real*8 Function getBetaKmax(Etol, RcouBeta, F)
	integer	i, n
	real*8	Etol, RcouBeta, x, low, high, F
! This is a service routine for ewald
	i=0; x=5.d0;
	do while(F(x, RcouBeta)>Etol)
		i = i+1
		x = x*2.d0
	end do

	n = i+60 ! search tolerance is 2^-60=8E-19
	low = 0.d0; high = x
	do i=0, n
		x = (low+high)*0.5d0
		if( F(x, RcouBeta)>Etol ) then
			low = x
		else
			high = x
		end if
	end do
	getBetaKmax=x
End
!
Real*8 Function getEdir(beta, r)
! Ewald Service Routine
	real*8	beta, r
	getEdir=erfc(beta*r)
End
!
Real*8 Function getErec(k, beta)
! Ewald Service Routine
	real*8	beta, k
	getErec=exp( -k*k/(4.d0*beta*beta) )/(k*k)
End
!
!
Subroutine cross(A, B, AxB)
! do a cross product with a 3D vector. 
	real*8	A(3), B(3), AxB(3)
	AxB(1) = A(2)*B(3) - A(3)*B(2)
	AxB(2) = A(3)*B(1) - A(1)*B(3)
	AxB(3) = A(1)*B(2) - A(2)*B(1)
End
!
!
Subroutine invA(A, Ainv)
! Invert a 3x3 matrix
	real*8	A(3,3), Ainv(3,3), Tmp
	Tmp =  A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
		& -A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
		& +A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
	if(abs(Tmp)<epsilon(1.d0)*1000) call ErrStop('Trying to invert a singular 3x3 Matrix', 1)
	Tmp = 1.d0/Tmp
	Ainv(1, 1) = (  A(2,2)*A(3,3)-A(2,3)*A(3,2)  )*Tmp
	Ainv(2, 1) = (-(A(2,1)*A(3,3)-A(2,3)*A(3,1)) )*Tmp
	Ainv(3, 1) = (  A(2,1)*A(3,2)-A(2,2)*A(3,1)  )*Tmp
	Ainv(1, 2) = (-(A(1,2)*A(3,3)-A(1,3)*A(3,2)) )*Tmp
	Ainv(2, 2) = (  A(1,1)*A(3,3)-A(1,3)*A(3,1)  )*Tmp
	Ainv(3, 2) = (-(A(1,1)*A(3,2)-A(1,2)*A(3,1)) )*Tmp
	Ainv(1, 3) = (  A(1,2)*A(2,3)-A(1,3)*A(2,2)  )*Tmp
	Ainv(2, 3) = (-(A(1,1)*A(2,3)-A(1,3)*A(2,1)) )*Tmp
	Ainv(3, 3) = (  A(1,1)*A(2,2)-A(1,2)*A(2,1)  )*Tmp
End Subroutine
!
!
Subroutine getPotPara(Ttypin, P, Ppot)
! translate raw parameters from SVD in (P) to physical parameters in Ppot
	real*8, parameter:: Reps=epsilon(1.d0)
	integer i
	real*8, intent(in) ::  P(*) 
    real*8, intent(out) :: Ppot(*)
    real*8 X(3)
	character(*),intent(in) :: Ttypin
    character*20 :: Ttyp
	real*8, parameter:: eps_3b = 1.D-4
	real*8  Pcub(4), Xcub(3)

	Ppot(:4) = P(:4)

!    Ttyp(1:3)=Ttypin(1:3)
    Ttyp=Ttypin
	call upcase(Ttyp); if(Ttyp(1:3)/='BD3') Ttyp=Ttyp(1:3)

	select case(Ttyp)
! 	case('COU', 'EWA', 'PEW')
! 		Ppot(1) = P(1)*FcovCou
	case('HAR')
		if(abs(P(1))>Reps) then
			Ppot(1) = -P(2)/P(1)
			Ppot(2) = P(1)
		end if
	case('QUA')
		!write(Lout, '(/, A)') 'Parameter for Quartic Potential'
		!write(Lout, *) '  #          Root                 Parameter'
		call getCubicRoot(P, X)
		do i=1, 3
			if(X(i)>0.d0 .and. X(i)<4.d0) then
				Ppot(1) = X(i)
				Ppot(4) = P(1)
				Ppot(3) = P(2)+3.d0*Ppot(4)*Ppot(1)
				Ppot(2) = P(3)+Ppot(1)*(2.d0*Ppot(3)-3.d0*Ppot(4)*Ppot(1))
				!write(Lout, '(I4, F20.9, \)') i, X(i)
				!write(Lout, '(4F20.9)') , Ppot(1:4)
			end if
		end do
	case('HARD')
		!Ppot(1) = Ppot(1)/atan(1.D0)*45.D0
		Ppot(1) = -P(2)/P(1)/atan(1.D0)*45.D0
		Ppot(2) = P(1)
	case('COS')
		Ppot(2) = sign(sqrt(P(1)**2+P(2)**2), P(2))
		Ppot(1) = acos(-P(1)/Ppot(2))
    case('CCO') !CCOS
        Ppot(2)= sign(sqrt(P(1)**2+P(2)**2), P(2))
        Ppot(1)= acos(-P(1)/Ppot(2))
	case('MUB')
		Ppot(1) = P(1) !krtheta
		Ppot(2) = P(2) !r1e^\prime
		Ppot(3) = P(3) !r2e^\prime
		Ppot(4) = P(4) !r3e
	case('QBB')
		Ppot(2) = -P(2) !krr
		Ppot(5) = -P(5) !k1b4
		Pcub(1)=-P(5)
		Pcub(2)=-P(4)
		Pcub(3)=-P(1)-P(2)
		Pcub(4)=-P(3)
		call getCubicRoot(Pcub, Xcub)
		do i=1, 3
			!write(Lout, '(I4, F20.9, \)') i, X(i)
			Ppot(1) = Xcub(i) !r1e
			Ppot(4) = -P(4)+3.*Ppot(5)*Ppot(1) !k1b3
			Ppot(3) = -P(1)+Ppot(1)*(2.D0*Ppot(4)-3.D0*Ppot(5)*Ppot(1)) !k1b2
			!  write(Lout, '(5F20.9)') , Ppot(1:5)
		end do
		do i=1, 3
			if(Xcub(i)>0.D0 .and. Xcub(i)<4.D0) then
				Ppot(1) = Xcub(i)
				Ppot(4) = -P(4)+3.*Ppot(5)*Ppot(1) !k1b3
				Ppot(3) = -P(1)+Ppot(1)*(2.D0*Ppot(4)-3.D0*Ppot(5)*Ppot(1)) !k1b2
				!    write(Lout, '(5F20.9)') , Ppot(1:5)
			end if
        end do
!    case('BND3BBBBA')
!		Ppot(1) = P(1) !k2
!		Ppot(2) = P(2) !k3
!		Ppot(3) = P(3) !k4
!		Ppot(4) = P(4) !krr
!		Ppot(5) = P(5) !krthetha
!		Ppot(6) = P(6) !re
!	case('BND3B7QUA')
!		Ppot(2) = -P(2) !krr
!		Ppot(3) = -P(3) !kr\theta
!		Ppot(7) = -P(7) !k1b4
!		!To Do
!		if(abs(P(3))>eps_3b) then
!		Ppot(1)=0.5*P(4)/P(3) !r1e
!		Ppot(6)=3.*Ppot(7)*Ppot(1)-P(6) !k1b3
!		Ppot(5)=-P(1)+(2.*Ppot(6)-3.*Ppot(7)*Ppot(1))*Ppot(1)!k1b2
!		Ppot(4) =(P(5)-Ppot(1)*(Ppot(5)-Ppot(1)*(Ppot(6)-Ppot(7)*Ppot(1)))-Ppot(2)*Ppot(1))/Ppot(3) ! r3e
!		else
!		Ppot(4)=0 !r3e is unknow, set to zero
!		Pcub(1)=-P(7)
!		Pcub(2)=-P(6)
!		Pcub(3)=-P(1)-P(2)
!		Pcub(4)=-P(5)
!		call getCubicRoot(Pcub, Xcub)

!		do i=1, 3
!			!write(Lout, '(I4, F20.9, \)') i, X(i)
!			Ppot(1) = Xcub(i) !r14
!			Ppot(6) = -P(6)+3.*Ppot(7)*Ppot(1) !k1b3
!			Ppot(5) = -P(1)+Ppot(1)*(2.D0*Ppot(6)-3.D0*Ppot(7)*Ppot(1)) !k1b2
!			!Ppot(7) = -P(7) !k1b4
!			write(Lout, '(4F20.9)') , Ppot(1:7)
!		end do
!		do i=1, 3
!			if(Xcub(i)>0.D0 .and. Xcub(i)<4.D0) then
!				Ppot(1) = Xcub(i)
!				Ppot(6) = -P(6)+3.*Ppot(7)*Ppot(1) !k1b3
!				Ppot(5) = -P(1)+Ppot(1)*(2.D0*Ppot(6)-3.D0*Ppot(7)*Ppot(1)) !k1b2
!				!Ppot(7) = -P(7) !k1b4
!				write(Lout, '(4F20.9)') , Ppot(1:7)
!			end if
!		end do
!		!=================================================
!		end if
!	case('BND3B7')
!		Ppot(1) = -P(1)
!		Ppot(2) = -P(2)
!		Ppot(3) = -P(3)
!		Ppot(4) = -P(4)
!		if(abs(P(4))>eps_3b) then
!			Ppot(5) = ((Ppot(2)-Ppot(3))/Ppot(4)*P(5)+P(6)-P(7))/(Ppot(1)+Ppot(2)-2.*Ppot(3))
!			Ppot(6) = ((Ppot(1)-Ppot(3))/Ppot(4)*P(5)-P(6)+P(7))/(Ppot(1)+Ppot(2)-2.*P(3))
!			Ppot(7) = (P(6)-Ppot(1)*Ppot(5)-Ppot(3)*Ppot(6))/Ppot(4)
!		else
!			Ppot(5)=0
!			Ppot(6)=(P(6)*Ppot(2)-P(7)*Ppot(3))/(Ppot(1)*Ppot(2)-Ppot(3)*Ppot(3))
!			Ppot(7)=(P(7)*Ppot(1)-P(6)*Ppot(3))/(Ppot(1)*Ppot(2)-Ppot(3)*Ppot(3))
!		end if
!	case('BND3B5')
!		Ppot(1) = -P(1)
!		Ppot(2) = -P(2)
!		Ppot(3) = -P(3)
!		if(abs(P(3))>eps_3b) then
!			Ppot(4) = P(4)/(2.*Ppot(3))
!		else
!			Ppot(4)=P(5)/(Ppot(1)+Ppot(2))
!			Ppot(5)=0
!		end if
        
	end select
    End Subroutine 
!
!   
Subroutine getParaRMS(Ttypin,P, KCVrms, Prms)
    use CRYOFFmod, only: Reps
    real*8, intent(in) :: P(*), KCVrms(*) 
    real*8, intent(out) :: Prms(*)
    character(*),intent(in) :: Ttypin
    character*20 :: Ttyp

    Prms(:4) = KCVrms(:4)

    Ttyp=Ttypin
    call upcase(Ttyp); if(Ttyp(1:3)/='BD3') Ttyp=Ttyp(1:3)
    
    select case(Ttyp)
    case('HAR')
    if(abs(P(1))>Reps .and. abs(P(2))>Reps) then
            Prms(1) = abs(P(2)/P(1))*sqrt(KCVrms(1)*KCVrms(1)/(P(1)*P(1))+KCVrms(2)*KCVrms(2)/(P(2)*P(2)))
            Prms(2) = KCVrms(1)
        end if
    case('QUA')
        Prms(1)=-1;Prms(2)=-1;Prms(3)=-1;Prms(4)=-1
    case('HARD')
        Prms(1)=-1;Prms(2)=-1
    case('COS')
        Prms(1)=-1;Prms(2)=-1
    case('CCO') !CCOS
        Prms(1)=-1;Prms(2)=-1
    case('MUB')
        Prms(1)=-1;Prms(2)=-1;Prms(3)=-1;Prms(4)=-1
    case('QBB')
        Prms(1)=-1;Prms(2)=-1;Prms(3)=-1;Prms(4)=-1;Prms(5)=-1      
    end select
End Subroutine getParaRMS
!
!
    
Subroutine getCubicRoot(P, X)
	real*8, parameter:: TwoPi = 8.d0*atan(1.d0), eps=1.d1*epsilon(1.d0)
	real*8 P(*), a, b, c, d, Alph, Beta, Delt, R1, R2, tht, X(3)

	X = 0.d0
	a = P(1)
	b = P(2)/(3.d0*a)
	c = P(3)/(6.d0*a)
	d = P(4)/(2.d0*a)

	Alph = -b*b*b + 3.d0*b*c - d
	Beta =  b*b - 2.d0*c
	Delt = Alph*Alph-Beta*Beta*Beta

	if(Delt>eps) then
		tht = Alph+sqrt(Delt); R1 = sign(abs(tht)**(1.d0/3.d0), tht)
		tht = Alph-sqrt(Delt); R2 = sign(abs(tht)**(1.d0/3.d0), tht)
		X(1) = -b+R1+R2
	else if(abs(Delt)<eps) then
		R1 = sign(abs(Alph)**(1.d0/3.d0), Alph)
		if(abs(R1)<eps) then
			X(1) = -b
		else
			X(1) = -b+2.d0*R1
			X(2) = -b-R1
		end if
	else if(Delt<-eps) then
		tht = acos(Alph/(sqrt(Beta)*Beta))
		X(1)  = -b+2.d0*sqrt(Beta)*cos(tht/3.d0)
		X(2)  = -b+2.d0*sqrt(Beta)*cos((tht+TwoPi)/3.d0)
		X(3)  = -b+2.d0*sqrt(Beta)*cos((tht-TwoPi)/3.d0)
	end if
End Subroutine getCubicRoot
!
!
Subroutine upcase(txt)
	integer, parameter:: IlowA=ICHAR('a'), IlowZ=ICHAR('z'), Iup=ICHAR('A')-IlowA
	integer	i, Ichr, Nchr
	character(*)	txt

	Nchr = Len_trim(txt)
	do i = 1, Nchr
		Ichr = ICHAR(txt(i:i))
		if(Ichr>=IlowA .and. Ichr<=IlowZ) txt(i:i) = CHAR(Ichr+Iup)
	end do
	txt=trim(adjustl(txt))
End Subroutine upcase
!
!
Subroutine getOption(txt, tag, option)
    implicit none
	integer i, j
	character(*) txt,tag,option
      
   
! Get option get a text option 
! It is assume the keyword(tag) and = has no space in between. 
! the option include everything after the = sign unless parentheses are used.
    
	i=index(txt, tag)
    if (i==0) option = " "
	if(i/=0) then
		j=index(txt(i:), '=')
               
		option=adjustl(txt(i+j:)) 
        
        if (j==0) then 
            option=" "
        endif 
        
        if (j>2) then
            if (index(txt(i+1:i+j-1),' ')/=0) option=" "
        endif
    
		if(option(1:1)=='(') option=option(2:index(option,')')-1)
	end if

	do i=1, len_trim(option)
		if(option(i:i)=='(' .or. option(i:i)==')' &
		&  .or.option(i:i)==',' .or.option(i:i)==':') option(i:i)=' '
	end do
End Subroutine getOption
!
!
Real*8 Function getValue(txt, tag)
	integer Ierr
	character(*) txt, tag
	character*100 option
! Look for tag in the string txt and obtain the numerical values as option of tag. 
! getvalue get a numerial value (as real number)
! it has the same limitation as getOption in misidentifying an option.  
    
	option=''; Ierr=0
	getValue=0.d0
	call getOption(txt, tag, option)
	if(len_trim(option)>0) read(option, *, IOstat=Ierr) getValue
	call ErrStop(trim(tag)//'  NOT a Number!', Ierr)
End Function getValue
!
!
!
Subroutine ErrStop(msg, Ierr)
	integer Ierr
	character(*) msg
	if(Ierr/=0) then
		print*, '!!!! ERROR !!!! '//trim(msg), Ierr
		stop
	end if
End Subroutine ErrStop
!
!
Subroutine Powell(Lout, N, P, Pinf, Psup, Pvec, Feps, Iter, MaxIter, Fnow, MyID)
! Minimization of function dErr of N variables using routin LinMin
! (dErr is not an argument, it is a fixed function name.)
! Input
! P(1:N): an initial starting point
! Pvec(1:N,1:N): an initial set of direction, usually the N unit vectors
! Ftol: (not used) the fractional tolerance in the function value such that
!       failure to decrease by more than this amount on one iteration signals doneness.
! Feps: tolerance in the function value such that
!       failure to decrease by more than this amount on one iteration signals doneness
! MaxIter: maximum allowed iterations
! Output
! P: set to the best point found
! Pvec: the then-current direction set
! Iter: the number of iterations taken.

! Put the routine in test mode. 
	integer, parameter:: TEST=0
	integer i, MyID, Lout, N, Ibig, Iter, MaxIter
	real*8	Feps, Reps, P(N), Pinf(N), Psup(N), Pvec(N, N), &
	&		Fpre, Fnow, Fxtr, dFbig, P0(N), Pxtr(N), Vect(N), dErr, dEtst

	P0   = P ! Save the initial point
	if(TEST==1) Fnow = dEtst(P0)
	if(TEST/=1) Fnow = dErr(P0)

	do Iter=1, MaxIter
		Fpre = Fnow
		if(MyID==0) then
			write(*, *)
			write(*, *)               '>> Iteration: ', Iter
			write(*, '(A, 600F12.6)') '    Y/P Initial:  ', Fpre, P
			write(Lout, *)
			write(Lout, *) ' Iteration: ', Iter
			write(Lout, '(A, 600F12.6)') '    Y/P Initial:  ', Fpre, P
		end if

		Ibig = 0; dFbig = 0.d0  ! Will be the biggest function decrease.
		do i=1, N     ! In each iteration, loop over all directions in the set.
			Vect(:) = Pvec(:, i) ! Copy the direction
			if(MyID==0) then
				write(*,    '(/, A, I3, A, 100F8.3)') '    >> Direction', i, ':', Vect
				write(Lout, '(/, A, I3, A, 100F8.3)') '      Direction', i, ':', Vect
			end if

			call LinMin(MyID, N, P, Pinf, Psup, Vect, Fnow) ! minimize along direction

			if(Fpre-Fnow>dFbig) then ! record it if it is the largest decrease so far.
				Ibig = i; dFbig = Fpre-Fnow
			end if
			if(MyID==0) then
				write(*,    '(A, 600F12.6)') '    Y/P Direction:', Fnow, P
				write(Lout, '(A, 600F12.6)') '      Y/P Direction:', Fnow, P
			end if
			if(abs(Fnow-Fpre)<Feps) return
			! if (2.*(Fpre-Fnow)<=Ftol*(abs(Fpre)+abs(Fnow))+epsilon(1.d0)) return	! Termination criterion.
		end do

		Pxtr = 2.d0*P-P0  ! extrapolated point
		Vect = P-P0       ! average direction moved
		P0   = P          ! Save the old starting point
		if(TEST==1) Fxtr = dEtst(Pxtr) ! Function value at extrapolated point
		if(TEST/=1) Fxtr = dErr(Pxtr)
		if(MyID==0) then
			write(*,    '(/, A, 100F8.3)') '    >> Direction Extrapolation:', Vect
			write(Lout, '(/, A, 100F8.3)') '      Direction Extrapolation:', Vect
		end if

		if(Fxtr<Fpre .and. 2.d0*(Fpre-2.d0*Fnow+Fxtr)*(Fpre-Fnow-dFbig)**2 &
		&				 < dFbig*(Fpre-Fxtr)**2) then ! use new direction
			call LinMin(MyID, N, P, Pinf, Psup, Vect, Fnow) ! Move to the minimum of the new direction,
			Pvec(:, Ibig) = Pvec(:, N) ! and save the new direction.
			Pvec(:, N) = Vect(:)
			if(MyID==0) then
				write(*,    '(A, 100F8.3)')  '    Direction Extrapolated:  ', Vect
				write(*,    '(A, 600F12.6)') '    Y/P Direction:', Fnow, P
				write(Lout, '(A, 100F8.3)')  '      Direction Extrapolated:  ', Vect
				write(Lout, '(A, 600F12.6)') '      Y/P Direction:', Fnow, P
			end if
		end if

		Reps = abs(Fnow-Fpre)
		if(MyID==0) then
			write(*,    '(/, A, 600F12.6)') '    Y/P Optimized:', Fnow, P
			write(*,    *)                  '   |dY|:  ', Reps
			write(Lout, '(/, A, 600F12.6)') '    Y/P Optimized:', Fnow, P
			write(Lout, *)                  '   |dY|:  ', Reps
		end if
		if(Fnow<Ftol .or. abs(Fnow-Fpre)<Feps) return
	end do
End Subroutine Powell
!
!
Real*8 Function Func1D(x)
	integer, parameter:: Nmax=100, TEST=1
	integer	N
	real*8	dErr, dEtst, x, P0(Nmax), Pvec(Nmax), P(Nmax)
	common /Func1DCom/ N, P0, Pvec

	P(1:N) = P0(1:N)+x*Pvec(1:N)
	if(TEST==1) Func1D = dEtst(P)
	if(TEST/=1) Func1D = dErr(P)
End Function Func1D
!
!
Real*8 Function dEtst(P)
	real*8	P(2), x, y

	x=P(1); y=P(2)
	dEtst= (abs(x)-5.)**2 + (abs(y)-5.)**2      ! [-10,10] Becker and Lago
	dEtst=x*x + y*y + 25.*(sin(x)**2+sin(y)**2) ! [-2pi,2pi] Eggcrate
	dEtst=100.*(y-x*x)**2 + (6.4*(y-0.5)**2-x-0.6)**2 ! [-5,5] Modified Rosenbrock
	dEtst= (1-x)**2 + 100.*(y-x*x)**2           ! [-2,2] Rosenbrock
	dEtst=(x*x+y-11.)**2 + (x+y*y-7)**2         ! [-5,5] Himmelblau
	dEtst=100.*sqrt(abs(y-.01*x*x))+.01*abs(x+10.) ![-15,-5][-3,3] Bukin
End Function dEtst
!
!
Subroutine LinMin(MyID, N, P, Pinf, Psup, Vect, Fret)
!	Given an N-dimensional point P(1:N) and an N-dimensional direction Vect(1:N),
!	moves and resets P to where the function Func1D(P) takes on a minimum
!	along the direction Vect from P
!	and replaces Vect by the actual vector displacement that P was moved.
!	Also returns as Fret the value of Func1D at the returned location P.
!	This is all accomplished by calling the routines MinBrak and BrentMin.

	integer, parameter:: Nmax=100
	real*8,  parameter:: Ptol = 1.d-3 !Tol passed to brent

	integer MyID, N, Ncom
	real*8	Fret, P(N), Pinf(N), Psup(N), Vect(N), &
	&		Xa, Xb, Xc, Xmin, &
	&		Pcom(Nmax), Vcom(Nmax), BrentMin, Func1D

	common /Func1DCom/ Ncom, Pcom, Vcom
	external Func1D

	Ncom = N  ! Set up the common block.
	Pcom = P
	Vcom = Vect

	Xa=0.d0; Xb=1.d0 ! Initial guess for brackets
	call MinBrak(MyID, Xa, Xb, Xc, Func1D)
	Fret = BrentMin(MyID, Xa, Xb, Xc, Func1D, Ptol, Xmin)

	Vect = Xmin*Vect	! Construct the vector results to return.
	P = P+Vect
End Subroutine LinMin
!
!
Subroutine MinBrak(MyID, Xa, Xb, Xc, Func)
!	Given a function Func and distinct initial points Xa and Xb,
!	this routine searches in the downhill direction
!	(defined by the function as evaluated at the initial points)
!	and returns new points Xa, Xb, Xc that bracket a minimum of the function.
!	Also returned are the function values at the three points, Fa, Fb and Fc.
!	Gold is the default ratio by which successive intervals are magnified
!	Glim is the maximum magnification allowed for a parabolic-fit step.
	integer MyID
	real*8, parameter:: Gold=1.618034D0, Glim=2.d0, Tiny=1.d-10
	real*8  Xa, Xb, Xc, Fa, Fb, Fc, Func, Fu, q, r, u, ulim
	external Func

	Fa = Func(Xa)
	Fb = Func(Xb)
	if(Fb>Fa) then ! Switch roles of a and b so that we can go downhill in the direction from a to b.
		u = Xa; Xa = Xb; Xb = u
		u = Fa; Fa = Fb; Fb = u
	end if

	Xc = Xb+Gold*(Xb-Xa)		! First guess for c.
	Fc = Func(Xc)

	do while(Fb>=Fc)					! do while": keep returning here until we bracket.
		r = (Xb-Xa)*(Fb-Fc)				! Compute u by parabolic extrapolation from a b c.
		q = (Xb-Xc)*(Fb-Fa)				! TINY is used to prevent any possible division by zero.
		u = Xb-( (Xb-Xc)*q-(Xb-Xa)*r )/(2.d0*sign(max(abs(q-r), TINY), q-r))
		ulim = Xb+Glim*(Xc-Xb)			! We won't go farther than this. Test various possibilities:
		if((Xb-u)*(u-Xc)>0.d0) then		! Parabolic u is between b and c: try it.
			Fu = Func(u)
			if(Fu<Fc) then 		! Got a minimum between b and c.
				Xa = Xb; Fa = Fb
				Xb = u;  Fb = Fu
				return
			else if(fu>Fb) then	!Got a minimum between between a and u.
				Xc = u; Fc = Fu
				return
			end if
			u = Xc+Gold*(Xc-Xb) 	! Parabolic fit was no use. Use default magnification.
			Fu = Func(u)
		else if((Xc-u)*(u-ulim)>0.d0) then 	! Parabolic fit is between c and its allowed limit.
			Fu = Func(u)
			if(Fu<Fc) then
				Xb = Xc; Xc = u
				u = Xc+Gold*(Xc-Xb)
				Fb = Fc; Fc = Fu
				Fu = Func(u)
			end if
		else if((u-ulim)*(ulim-Xc)>0.d0) then 	! Limit parabolic u to maximum allowed value.
			u = ulim; Fu = Func(u)
		else 									! Reject parabolic u, use default magnification.
			u = Xc+Gold*(Xc-Xb); Fu = Func(u)
		end if
		Xa = Xb; Xb = Xc; Xc = u		!Eliminate oldest point and continue.
		Fa = Fb; Fb = Fc; Fc = Fu
	end do

	if(Xa>Xc) then
		u = Xa; Xa = Xc; Xc = u
		u = Fa; Fa = Fc; Fc = u
	end if
	if(MyID==0) write(*, '(2(A, 3F9.3))') '       lambda Bracketed: ', Xa, Xb, Xc, '    Function Value: ', Fa, Fb, Fc
End Subroutine MinBrak
!
!
Real*8 function BrentMin(MyID, Xa, Xb, Xc, func, Rtol, Xmin)
! func: external function
! Xa, Xb, Xc: bracketing triplet of abscissas, Xa<=Xb<=Xc, func(Xb) < func(Xa)||func(Xc)
! Xmin: isolated minimum point to precision Rtol using Brent's method
! BrentMin: minimum function value
! MaxIter: Maximum allowed number of iterations
! Cgold: complance of golden ratio
! Eps: a small number protects against trying to achieve fractional
! accuracy for a minimum that happens to be exactly zero.

	integer, parameter:: MaxIter=100
	real*8,  parameter:: Cgold=0.381966011250105097d0, Eps=1.d-6

	integer MyID, Iter
	real*8	Xa, Xb, Xc, Xmin, Rtol, tol, func, &
	&	a, b, m, x0, x1, x2, u, f0, f1, f2, fu, d, e, p, q, r

	! a, b: interval within the minimum should be bracketed
	!       no function evaluation will be requested outside that range
	a = min(Xa, Xc) ! a and b must be in a scending order,
	b = max(Xa, Xc) ! though the input abscissas need not be.
	if(.not. (a<=Xb .and. Xb<=b)) call ErrStop('NOT Xa<=Xb<=Xb !', 1)

	! Initialization
	! (x0, x1, x2): Parabolic interpolation through this Memory triple, updated in each interation
	x0 = Xb; f0 = func(x0) ! x0: point with least value found so far(or most recent one for a tie)
	x1 = x0; f1 = f0       ! x1: point with second least value
	x2 = x1; f2 = f1       ! x2: previous value of x1
	d  =0.d0 ! d: step size and direction, the distance moved on the step before last
	u  =0.d0 ! u: point at which the function was evaluated most recently
	e = 0.d0 ! e: memorizes the step size (and direction) taken two iterations ago
			 !    and it is used to (definitively) fall-back to golden-section steps
			 !    when its value is too small (indicating that the polynomial fitting
			 !    is not helping to speedup the convergence.)

	do Iter=1, MaxIter
		m = 0.5d0*(a+b) ! m: midpoint of current interval (a, b), function not evaluated
		tol = Rtol*abs(x0)+Eps

		if(MyID==0) write(*, '(A, I4, 2X, 6F9.3)') '       Brent  Min.:', Iter, a, b, x0, e, d, u

		! typical ending configuration: a and b are 2×x×tol apart
		! with x (the best abscissa) at the midpoint of a and b
		! and therefore fractionally accurate to ±tol.
		! tol = rtol * abs(x0) + eps
		! x0 is the best guess found so far
		! It converges if evaluating a next guess would
		! imply evaluating f at a point that
		! is closer than tol to a previously evaluated one
		if( abs(x0-m)+0.5d0*(b-a) <= 2.d0*tol) then ! exit with best values
			Xmin = x0
			BrentMin = f0
			return
		end if

		! To be acceptable, the parabolic step must
		! 1. fall within the bounding interval (a, b), and
		! 2. imply a movement from the best current value x
		!    that is less than half the movement of the step before last
		r=0.d0; p=r; q=r
		if(abs(e)>tol) then
			! Constructa trial parabolic (Lagrange polynomial)
			! that goes through (x0, f0), (x1, f1), (x2, f2)
			r = (x0-x1)*(f0-f2)
			q = (x0-x2)*(f0-f1)
			p = (x0-x2)*q-(x0-x1)*r
			q = 2.d0*(q-r)
			if(q>0.d0) p = -p
			q = abs(q)
			r = e
			e = d
		end if

		! conditions determine the acceptability of the parabolic fit
		! function must not be evaluated too close to a or b.
		if(abs(p)<abs(0.5d0*q*r) .and. q*(a-x0)<p .and. p<q*(b-x0)) then
			d = p/q            ! Take the parabolic step.
			u = x0+d
			if (u-a<2.d0*tol .or. b-u<2.d0*tol) d = sign(tol, m-x0)
		else
			if(x0>=m) e = a-x0 ! for golden section step, take into the larger segments
			if(x0< m) e = b-x0
			d = Cgold*e        ! Take the golden section step
		end if

		! function must not be evaluated too close to x
		u = x0+d ! with d computed from parabolic fit or golden section
		if(abs(d)<tol) u = x0+sign(tol, d)
		! Notice that we have u in [a+tol, x0-tol] or
		!                     u in [x0+tol, b-tol]
		! (if one ignores rounding errors.)

		!  Is the most recently evaluated point better (or equal) than the best so far?
		fu = func(u)         ! one function evaluation per iterati[p
		if(fu<=f0) then      ! find new lower point at u, Decrease interval size.
			if(u>=x0) a = x0 ! u on the right side of x
			if(u< x0) b = x0 ! u on the left  side
			x2 = x1; f2 = f1 ! Shift: drop the previous third best point out and
			x1 = x0; f1 = f0 ! include the newest point (found to be the best so far)
			x0 = u;  f0 = fu
		else
			if(u< x0) a = u
			if(u>=x0) b = u
			! Is the most recently evaluated point at better (or equal) than the second best one?
			if(fu<=f1 .or. x1==x0) then
				x2 = x1; f2 = f1 ! Insert u between (rank-wise) x0 and x1 in the triple # (x0, x1, x2).
				x1 = u;  f1 = fu
			else if(fu<=f2 .or. x2==x0 .or. x2==x1) then
				x2 = u;  f2 = fu ! Insert u in the last position of the triple (x0, x1, x2)
			end if
		end if
	end do
	call ErrStop('BrentMin Exceed Maximum Iterations!', 1)
End function BrentMin
!
!
subroutine DE_Fortran90(obj, Dim_XC, XCmin, XCmax, VTR, NP, itermax, F_XC, &
           CR_XC, strategy, refresh, iwrite, bestmem_XC, bestval, nfeval, &
		   F_CR, method, restart,nrestart,MyID, np_world,comm_world, &
		   MyID_pop, np_pop,comm_pop,MyID_force, np_force,comm_force)
!.......................................................................
!
! Differential Evolution for Optimal Control Problems
!
!.......................................................................
!  This Fortran 90 program translates from the original MATLAB
!  version of differential evolution (DE). This FORTRAN 90 code
!  has been tested on Compaq Visual Fortran v6.1.
!  Any users new to the DE are encouraged to read the article of Storn and Price.
!
!  Refences:
!  Storn, R., and Price, K.V., (1996). Minimizing the real function of the
!    ICEC'96 contest by differential evolution. IEEE conf. on Evolutionary
!    Comutation, 842-844.
!-------------------------------------------------------------------------
!  This subroutine parallizes the original F90 code of Dr. Feng-Sheng Wang.
!  Department of Chemical Engineering, National Chung Cheng University,
!  Chia-Yi 621, Taiwan, e-mail: chmfsw@ccunix.ccu.edu.tw
!.........................................................................
!                obj : The user provided file for evlauting the objective function.
!                      subroutine obj(xc,fitness)
!                      where "xc" is the real decision parameter vector.(input)
!                            "fitness" is the fitness value.(output)
!             Dim_XC : Dimension of the real decision parameters.
!      XCmin(Dim_XC) : The lower bound of the real decision parameters.
!      XCmax(Dim_XC) : The upper bound of the real decision parameters.
!                VTR : The expected fitness value to reach.
!                 NP : Population size.
!            itermax : The maximum number of iteration.
!               F_XC : Mutation scaling factor for real decision parameters.
!              CR_XC : Crossover factor for real decision parameters.
!           strategy : The strategy of the mutation operations is used in HDE.
!            refresh : The intermediate output will be produced after "refresh"
!                      iterations. No intermediate output will be produced if
!                      "refresh < 1".
!             iwrite : The unit specfier for writing to an external data file.
! bestmen_XC(Dim_XC) : The best real decision parameters.
!              bestval : The best objective function.
!             nfeval : The number of function call.
!         method(1) = 0, Fixed mutation scaling factors (F_XC)
!                   = 1, Random mutation scaling factors F_XC=[0, 1]
!                   = 2, Random mutation scaling factors F_XC=[-1, 1]
!         method(2) = 1, Random combined factor (F_CR) used for strategy = 6
!                        in the mutation operation
!                   = other, fixed combined factor provided by the user
!         method(3) = 1, Saving results in a data file.
!                   = other, displaying results only.
!         restart   = .true., reading NP XC from restart.dat;
!                   = .false.,randomly initializing; saving NP XC to restart.dat
!         nrestart  = an integer, saving NP XC every nrestart iterations.
  implicit none
#ifdef MPI
  include "mpif.h"
#endif
  integer(kind=4), parameter :: IB=4, RP=8
  integer(kind=IB), intent(in) :: NP, Dim_XC, itermax, strategy,   &
       iwrite, refresh
  integer(kind=IB) MyID,np_world,comm_world
  integer(kind=IB) MyID_pop,np_pop,comm_pop
  integer(kind=IB) MyID_force,np_force,comm_force
  real(kind=RP), intent(in) :: VTR, CR_XC
  real(kind=RP) :: F_XC, F_CR
  real(kind=RP), dimension(Dim_XC), intent(in) :: XCmin, XCmax
  real(kind=RP), dimension(Dim_XC), intent(inout) :: bestmem_XC
  real(kind=RP), intent(out) :: bestval
  integer(kind=IB), intent(out) :: nfeval
  real(kind=RP), dimension(NP,Dim_XC) :: pop_XC, bm_XC, mui_XC, mpo_XC,   &
       popold_XC, rand_XC, ui_XC
  integer(kind=IB) :: i, ibest, jter
  integer(kind=IB), dimension(NP) :: rot, a1, a2, a3, a4, a5, rt
  integer(kind=IB), dimension(4) :: ind
  real(kind=RP) :: tempval
  real(kind=RP), dimension(NP) :: val
  real(kind=RP), dimension(Dim_XC) :: bestmemit_XC
  real(kind=RP), dimension(Dim_XC) :: rand_C1
  integer(kind=IB), dimension(3), intent(in) :: method
  real(kind=RP) :: obj
  logical:: restart
  integer irestart,nrestart
  integer Ierr, j
#ifdef MPI
  integer idiv,irem
  integer nvec_proc, myvec_st, myct, my_nvec
  integer,allocatable:: displs(:),scounts(:),rdispls(:),rcounts(:)
  real*8, allocatable::mypop(:,:)
  real*8 pop_XC_T(Dim_XC,NP),val_all(NP),val_it(NP)
  real*8, allocatable::mypop_T(:,:), val_proc(:)
#endif
  external  obj
  intrinsic max, min, random_number, mod, abs, any, all, maxloc
  interface
     function randperm(num)
       !use data_type, only : IB
       implicit none
       integer(kind=4), parameter :: IB=4, RP=8
       integer(kind=IB), intent(in) :: num
       integer(kind=IB), dimension(num) :: randperm
     end function randperm
  end interface

  !-----Initialize a population --------------------------------------------!
  if(MyID==0) then
     write(*,*) "Initializing DE population"
     bestval=HUGE(1.d0)
     if (restart==.true.) then
        open(unit=irestart, file='restart.dat', status='old', IOstat=Ierr)
        if(Ierr/=0) then
           print*, '>>Input file  restart.dat  NOT Exist !'
#ifdef MPI
           print*, '>>Normal Exit of Program CRYOFF.'
           call MPI_Abort(comm_world, Ierr)
#else
           stop '>>Normal Exit of Program CRYOFF.'
#endif
        endif
        do i=1,NP
           do j=1,Dim_XC
              read(irestart,*) pop_XC(i,j)
           enddo
        enddo
        do j=1,Dim_XC
           read(irestart,*) bestmemit_XC(j)
        enddo
        do j=1,Dim_XC
           read(irestart,*) bestmem_XC(j)
        enddo
        read(irestart,*) bestval
     else
        pop_XC=0.0_RP
        pop_XC(1,:)=bestmem_XC
        do i=2,NP
           call random_number(rand_C1)
           pop_XC(i,:)=XCmin+rand_C1*(XCmax-XCmin)
        end do
        open(unit=irestart, file='restart.dat', IOstat=Ierr)
        if(Ierr/=0) then
           print*, '>>Can not creat restart.dat  !'
#ifdef MPI
           print*, '>>Normal Exit of Program CRYOFF.'
           call MPI_Abort(comm_world, Ierr)
#else
           stop '>>Normal Exit of Program CRYOFF.'
#endif
        endif
     end if
#ifdef MPI
     !trans pose
     do i=1,NP
        do j=1,Dim_XC
           pop_XC_T(j,i)=pop_XC(i,j)
        end do
     end do
#endif
  end if

!-----------------------------------------------------------------------
#ifdef MPI
  allocate(displs(np_pop),scounts(np_pop),stat=Ierr)
  allocate(rdispls(np_pop),rcounts(np_pop),stat=Ierr)
  idiv=int(NP/np_pop)
  irem=mod(NP,np_pop)
  my_nvec=idiv
  if(MyID_pop<irem) then
     my_nvec=idiv+1
  end if
  myct=my_nvec*Dim_XC

  allocate(mypop(my_nvec,Dim_XC),stat=Ierr)
  allocate(mypop_T(Dim_XC,my_nvec),stat=Ierr)

  do i=0, np_pop-1
     nvec_proc=idiv
     if(i<irem) then
        nvec_proc=idiv+1
     end if
     if(i<=irem) then
        myvec_st=i*(idiv+1)+1
     else
        myvec_st=irem*(idiv+1)+(i-irem)*idiv+1
     end if
     displs(i+1)=(myvec_st-1)*Dim_XC
     scounts(i+1)=nvec_proc*Dim_XC
     rdispls(i+1)=myvec_st-1
     rcounts(i+1)=nvec_proc
  end do
  call MPI_Bcast(pop_XC_T, Dim_XC*NP, MPI_REAL8, 0, comm_force, Ierr)
  if (np_pop > 1 ) then
  call MPI_ScatterV(pop_XC_T,scounts,displs, MPI_REAL8, mypop_T, myct, MPI_REAL8, 0, comm_pop, Ierr)
  endif
  !transpose back
  do i=1,my_nvec
     do j=1, Dim_XC
        mypop(i,j)=mypop_T(j,i)
     end do
  end do
#endif
!-----------------------------------------------------------------

!------Evaluate fitness functions and find the best member-----------------!
  val=0.0_RP
  nfeval=0
#ifdef MPI
  allocate(val_proc(my_nvec),stat=Ierr)
  do i=1,my_nvec
     val_proc(i)=obj(mypop(i,:))
  end do
  if (np_pop > 1 ) then
  call MPI_GATHERV(val_proc,my_nvec,MPI_REAL8,val_all,rcounts,rdispls,MPI_REAL8, 0, comm_pop,Ierr)
  endif
  if(MyID==0) then
     ibest=1
     bestval=val_all(1)
     nfeval=nfeval+1
     do i=2,NP
        nfeval=nfeval+1
        if (val_all(i) < bestval) then
           ibest=i
           bestval=val_all(i)
        end if
     end do
     bestmemit_XC=pop_XC(ibest,:)
     bestmem_XC=bestmemit_XC
  end if
#else
  val(1)=obj(pop_XC(1,:))
  ibest=1
  bestval=val(1)
  nfeval=nfeval+1
  do i=2,NP
     !call obj(pop_XC(i,:), val(i))
     val(i)=obj(pop_XC(i,:))

     nfeval=nfeval+1
     if (val(i) < bestval) then
        ibest=i
        bestval=val(i)
     end if
  end do
     bestmemit_XC=pop_XC(ibest,:)
     bestmem_XC=bestmemit_XC
#endif

!!--------------------------------------------------------------------------!!

     bm_XC=0.0_RP
     rot=(/(i,i=0,NP-1)/)
     jter=1

     if(MyID==0) then
        write(unit=*,FMT=206) jter, nfeval

        do i=1,Dim_XC
           write(*,FMT=202) i,bestmem_XC(i)
        end do
        write(unit=*, FMT=201) bestval
!---Perform evolutionary computation------------------------------------!!
        write(*,*)
        write(*,*) "Perfrom DE evolution"
        write(*,*)
     end if

  do while (jter <= itermax)
     popold_XC=pop_XC
     !dump out pop_XC
     if (MyID==0 .and. mod(jter,nrestart)==1) then
        rewind(irestart)
        do i=1,NP
           do j=1,Dim_XC
              write(irestart,FMT=207) pop_XC(i,j)
           enddo
        enddo
        do j=1,Dim_XC
           write(irestart,FMT=207) bestmemit_XC(j)
        enddo
        do j=1,Dim_XC
           write(irestart,FMT=207) bestmem_XC(j)
        enddo
        write(irestart,FMT=207) bestval
     endif
!!-----Mutation operation--------------------------------------------------!!
!      ind=randperm(4)
!      a1=randperm(NP)
     rt=mod(rot+ind(1),NP)
     a2=a1(rt+1)
     rt=mod(rot+ind(2),NP)
     a3=a2(rt+1)
     rt=mod(rot+ind(3),NP)
     a4=a3(rt+1)
     rt=mod(rot+ind(4),NP)
     a5=a4(rt+1)
     bm_XC=spread(bestmemit_XC, DIM=1, NCOPIES=NP)

!----- Generating a random sacling factor--------------------------------!
     select case (method(1))
     case (1)
        call random_number(F_XC)
     case(2)
        call random_number(F_XC)
        F_XC=2.0_RP*F_XC-1.0_RP
     end select

     !---- select a mutation strategy-----------------------------------------!
     select case (strategy)
     case (1)
        ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

     case default
        ui_XC=popold_XC(a3,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

     case (3)
        ui_XC=popold_XC+F_XC*(bm_XC-popold_XC+popold_XC(a1,:)-popold_XC(a2,:))

     case (4)
        ui_XC=bm_XC+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:)-popold_XC(a4,:))

     case (5)
        ui_XC=popold_XC(a5,:)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:)+popold_XC(a3,:) &
             -popold_XC(a4,:))
     case (6) ! A linear crossover combination of bm_XC and popold_XC
        if (method(2) == 1) call random_number(F_CR)
        ui_XC=popold_XC+F_CR*(bm_XC-popold_XC)+F_XC*(popold_XC(a1,:)-popold_XC(a2,:))

     end select

!------Crossover operation-------------------------------------------------!
     call random_number(rand_XC)
     mui_XC=0.0_RP
     mpo_XC=0.0_RP
     where (rand_XC < CR_XC)
        mui_XC=1.0_RP
        !           mpo_XC=0.0_RP
     elsewhere
        !           mui_XC=0.0_RP
        mpo_XC=1.0_RP
     end where

     ui_XC=popold_XC*mpo_XC+ui_XC*mui_XC

!------Evaluate fitness functions and find the best member-----------------

#ifdef MPI
     if(MyID==0) then
        do i=1,NP
           ui_XC(i,:)=max(min(ui_XC(i,:),XCmax),XCmin)
        end do
        !trans pose
        do i=1,NP
           do j=1,Dim_XC
              pop_XC_T(j,i)=ui_XC(i,j)
           end do
        end do
     end if
     call MPI_Bcast(pop_XC_T, Dim_XC*NP, MPI_REAL8, 0, comm_force, Ierr)
     if (np_pop > 1) call MPI_ScatterV(pop_XC_T,scounts,displs, MPI_REAL8, mypop_T, myct, MPI_REAL8, 0, comm_pop, Ierr)
     !transpose back
     do i=1,my_nvec
        do j=1, Dim_XC
           mypop(i,j)=mypop_T(j,i)
        end do
     end do

     do i=1,my_nvec
        val_proc(i)=obj(mypop(i,:))
     end do
     if (np_pop > 1) call MPI_GATHERV(val_proc,my_nvec,MPI_REAL8,val_it,rcounts,rdispls,MPI_REAL8, 0, comm_pop,Ierr)
     if(MyID==0) then
        do i=1,NP
           tempval=val_it(i)
           nfeval=nfeval+1
           if (tempval < val_all(i)) then
              pop_XC(i,:)=ui_XC(i,:)
              val_all(i)=tempval
              if (tempval < bestval) then
                 bestval=tempval
                 bestmem_XC=ui_XC(i,:)
              end if
           end if
        end do
     end if
#else
     do i=1,NP
!------Confine each of feasible individuals in the lower-upper bound-------
        ui_XC(i,:)=max(min(ui_XC(i,:),XCmax),XCmin)
        tempval=obj(ui_XC(i,:))
        nfeval=nfeval+1
        if (tempval < val(i)) then
           pop_XC(i,:)=ui_XC(i,:)
           val(i)=tempval
           if (tempval < bestval) then
              bestval=tempval
              bestmem_XC=ui_XC(i,:)
           end if
        end if
     end do
#endif

     if(MyID==0) then
        bestmemit_XC=bestmem_XC

        if( (refresh > 0) .and. (mod(jter,refresh)==0)) then
           if (method(3)==1) write(unit=iwrite,FMT=203) jter
           write(unit=*, FMT=203) jter
           do i=1,Dim_XC
              if (method(3)==1) write(unit=iwrite, FMT=202) i, bestmem_XC(i)
              write(*,FMT=202) i,bestmem_XC(i)
           end do
           if (method(3)==1) write(unit=iwrite, FMT=201) bestval
           write(unit=*, FMT=201) bestval
        end if
     end if
#ifdef MPI
     if (np_pop > 1) call MPI_BARRIER(comm_pop,Ierr) !????
#endif
     jter=jter+1
     if(MyID==0) then
        if ( bestval <= VTR .and. refresh > 0) then
           write(unit=iwrite, FMT=*) ' The best fitness is smaller than VTR'
           write(unit=*, FMT=*) 'The best fitness is smaller than VTR'
           exit
        endif
     end if

     if(MyID ==0) then
        write(unit=*,FMT=206) jter, nfeval
        do i=1,Dim_XC
           write(*,FMT=202) i,bestmem_XC(i)
        end do
        write(unit=*, FMT=201) bestval
        write(unit=*, FMT=205) VTR

        write(unit=iwrite,FMT=206) jter, nfeval
        do i=1,Dim_XC
           write(iwrite,FMT=202) i,bestmem_XC(i)
        end do
        write(unit=iwrite, FMT=201) bestval
        write(unit=iwrite, FMT=205) VTR
     end if
     !-----------------------------------
  end do

  if(MyID ==0) then
     write(*,*) "end the evolutionary computation"
  end if

#ifdef MPI
  deallocate(displs,scounts,stat=Ierr)
  deallocate(rdispls,rcounts,stat=Ierr)
  deallocate(mypop,stat=Ierr)
  deallocate(mypop_T,stat=Ierr)
  deallocate(val_proc,stat=Ierr)
#endif

!!------end the evolutionary computation------------------------------!!
201 format(2x, 'bestval =', ES14.7)
202 format(5x, 'bestmem_XC(', I3, ') =', ES12.5)
203 format(2x, 'No. of iteration  =', I8)
204 format(2x, 'No. of function evaluation  =', I8)
205 format(2x, 'targetval =', ES14.7,/)
206 format(2x, 'No. of iter. =', I8, '  No. of fun. eval.  =', I8)
207 format(ES14.7)
end subroutine DE_Fortran90
!
!
function randperm(num)
  !use data_type, only : IB, RP
  implicit none
  integer(kind=4), parameter :: IB=4, RP=8
  integer(kind=IB), intent(in) :: num
  integer(kind=IB) :: number, i, j, k
  integer(kind=IB), dimension(num) :: randperm
  real(kind=RP), dimension(num) :: rand2
  intrinsic random_number
  call random_number(rand2)
  do i=1,num
     number=1
     do j=1,num
        if (rand2(i) > rand2(j)) then
	       number=number+1
        end if
     end do
     do k=1,i-1
        if (rand2(i) <= rand2(k) .and. rand2(i) >= rand2(k)) then
	       number=number+1
        end if
     end do
     randperm(i)=number
  end do
  return
end function randperm
!
!
Subroutine Simplex(Lout, Mopt, Nopt, P, Y, funk, deps, Iter, MaxIter, MyID)
    use CRYOFFmod, only : YesOut
	integer i, MyID, Lout, Mopt, Nopt, Neva, Iter, MaxIter, Ihig, Ilow, Inhi
	real*8  P(Mopt, Nopt), Y(Mopt), Psum(Nopt)
	real*8 deps, Reps, Ysav, Ytry
	real*8, external:: funk, Simptry

	Neva = 0
	Iter = 0

100	do i=1, Nopt ! starting or overall contracted, Recompute psum
		Psum(i) = sum(P(1:Mopt, i))
	end do

200	Ihig = maxloc(Y, 1) ! highest (worst)
	Ilow = minloc(Y, 1) ! lowest (best)
	Inhi = maxloc(Y, 1, Y<Y(Ihig)) ! next-highest

	Reps = Y(Ihig)-Y(Ilow)
	if(MyID==0 .and. (mod(Iter,5)==0 .or. YesOut(3))) then
		write(*, *) '>Iteration: ', Iter
		write(*, *) ' Ymin = ', Y(Ilow)
		write(*, *) ' Ymax = ', Y(Ihig)
		write(*, '(A7, 600F12.6)') 'Pmin:', P(Ilow,:)
		write(*, '(A7, 600F12.6)') 'Pmax:', P(Ihig,:)
		write(*, *) ' dY   = ', Reps
		write(*, *)
		write(Lout, *) ' Iteration:', Iter
        write(Lout, *) ' Vertices of the simplex followed by Y of each vertex'
		do  i = 1, Mopt
			write(Lout, '(4X, 600F12.6)')  P(i,:), Y(i)
		end do
		write(Lout, '(A, 600F12.6)') '  min Y/P:', Y(Ilow), P(Ilow,:)
		write(Lout, '(A, 600F12.6)') '  max Y/P:', Y(Ihig), P(Ihig,:)
		write(Lout, *) ' Ymax-Ymin:', Reps, achar(10)
	end if

	if (Reps<deps .or. Iter>=MaxIter) then
		Y( [1,Ilow]   ) = Y( [Ilow,1] )   ! returning, swap best value
		P( [1,Ilow], :) = P( [Ilow,1], :) ! and point to slot 1
		return
	end if

	Neva = Neva+2
	Iter = Iter+1

	! reflect the simplex from the high point
	Ytry = Simptry(Mopt, Nopt, P, Y, funk, Psum, Ihig, -1.d0)

	if (Ytry<=Y(Ilow)) then ! try extrapolation by a factor 2
		Ytry = Simptry(Mopt, Nopt, P, Y, funk, Psum, Ihig, 2.d0)
	else if (Ytry>=Y(Inhi)) then ! worse than the 2nd-highest, contraction
		Ysav = Y(Ihig)
		Ytry = Simptry(Mopt, Nopt, P, Y, funk, Psum, Ihig, 0.5d0)
		if (Ytry>=Ysav) then ! Can't get rid of high point, contract around the lowest(best) point
			do i=1, Mopt
				if(i/=Ilow) then
					Psum(:) = 0.5d0*(P(i,:)+P(Ilow,:))
					P(i, :) = Psum(:)
					Y(i) = funk(Psum)
				end if
			end do
			Neva = Neva+Nopt ! Keep track of function evaluations
			goto 100         ! back for the test of doneness and the next iteration
		end if
	else
		Neva = Neva-1        ! Correct the evaluation count
	end if

	goto 200
End Subroutine Simplex
!
!
Real*8  Function Simptry(Mopt, Nopt, P, Y, funk, Psum, Ihig, fac)
	! Extrapolates by factor fac through the face across from the high point
	! tries it, and replaces the high point if the new point is better
	integer Mopt, Nopt, Ihig
	real*8  fac, Rfac, Y(Mopt), P(Mopt, Nopt), Psum(Nopt), Ptry(Nopt)
	real*8, external:: funk

	Rfac = (1.d0-fac)/dble(Nopt)
	Ptry(:) = Psum(:)*Rfac-P(Ihig, :)*(Rfac-fac)
	Simptry = funk(Ptry)
	if (Simptry<Y(Ihig)) then ! better than the highest, replace the highest
		Y(Ihig) = Simptry
		Psum(:) = Psum(:)-P(Ihig, :)+Ptry(:)
		P(Ihig, :) = Ptry(:)
	end if
End Function Simptry
!
!
!
Integer Function Level(txt)
	integer i, j, k
    character(*) txt
	character*200 key
! Determine the level of the force field definition in the .ff file.
! sections start with [# will have negative levels to allow whole section to be skipped.
! 3 levels are defined for the keywords.
! Data lines have levels of 1000. 
    
   	k=0; j=min(200, len_trim(txt))
    if (len_trim(txt) > 200) then 
        print*, 'warning: line over 200 characters long'
        print*, txt
    endif 
    
    !remove blanks from input and put txt in key
	do i=1, j
		if(txt(i:i)/=' ' .and. txt(i:i)/='	') then
			k=k+1
			key(k:k)=txt(i:i)
		end if
    end do

    ! k variable reused. It can only be 1 or -1 in this part. 
	k=1; i=index(key, ']')-1
	if(key(1:2)=='[#') then
		k=-1
		key=key(3:i)
	else
		key=key(2:i)
	end if

	Level=1000
	select case(trim(key(1:3)))
		case('HAR', 'QUA', 'COS', 'NCO', 'CNC', 'CCO'); Level=3
		case('ATO', 'BON', 'ANG', 'DIH', 'BD3', 'EXC', 'PAI','FUD'); Level=2
		case('MOL', 'FIL', 'OPT', 'KEY', &
		&	'COU',  'THC',  'GLJ',  'BUC', 'FDB', &
		&	'EXP', 'PEX', 'GEX', 'POW', 'TTP', 'SRD', 'FDP', 'STR', 'CHA', 'CST', 'EQV'); Level=1 !ying  EQV
    end select
        
    if (trim(key(1:4)) .eq. 'CDIH') then
     Level=2
    endif 

	Level =k*Level
end
!
!
Subroutine uniq(n, m, arr)
    integer, intent(in) :: n
    integer, intent(out) :: m
	integer i, j
	character(*) arr(*)
!  arr is an array of strings. uniq remove duplicates
! and create a new array with only uniq names
! the number of uniq names m is returned. 
	do i=1, n
		if(arr(i)/='') then
			do j=i+1, n
				if(arr(j)==arr(i)) arr(j)=''
			end do
		end if
	end do

	do i=1, n
		if(arr(i)=='') then
			do j=i+1, n
				if(arr(j)/='') then
					arr(i)=arr(j)
					arr(j)=''
					exit
				end if
			end do
		end if
	end do

	do i=1, n
		if(arr(i)=='') exit
	end do
	m=i-1
End Subroutine uniq
!
!
Real*8 Function gammai(x) ! Γ(x) for x>0
	real*8 :: x, C(0:6)=[ 1.000000000190015d0, &
	&  76.18009172947146d0, -86.50532032941677d0,    24.01409824083091d0, &
	&  -1.231739572450155d0,  0.1208650973866179d-2, -0.5395239384953d-5 ]

	gammai = (x+5.5d0)**(x+0.5d0) * exp(-x-5.5d0) * sqrt(8.d0*atan2(1.d0,1.d0)) / x &
	&      * ( C(0) + C(1)/(x+1.d0) + C(2)/(x+2.d0) + C(3)/(x+3.d0)    &
	&        + C(4)/(x+4.d0) + C(5)/(x+5.d0) + C(6)/(x+6.d0) )
End
!
!
Real*8 Function igamma(a, x) ! Γ(a, x)=Γ(a) if x<=0 or a<=0
! 	ITMAX: the maximum allowed number of iterations
! 	EPS  : the relative accuracy
! 	FPMIN: a number near the smallest representable floating-point number
	integer:: n, ITMAX=100
	real*8 :: a, x, an, Tn, b, c, d, EPS=1.d3*epsilon(1.d0), FPMIN=1.d-30

	if(x<=0.d0 .or. a<=0.d0) then
		igamma=gamma(a)
		return
	end if

	igamma=0.d0
	if(x<a+1.d0) then ! use series representation
		an=a; Tn=1.d0/a; igamma=Tn
		do n=1, ITMAX
			an = an+1.d0
			Tn = Tn*x/an
			igamma = igamma+Tn
			if(abs(Tn)<abs(igamma)*EPS) then
				igamma = gamma(a)-igamma * exp(-x) * x**a
				return
			end if
		end do
	else ! use continued fraction representation
		b=x+1.d0-a  ! modified Lentz's method with b0=0
		c=1.d0/FPMIN
		d=1.d0/b
		igamma=d
		do n=1, ITMAX
			an=dble(-n)*(dble(n)-a)
			b=b+2.d0
			d=an*d+b; if(abs(d)<FPMIN) d=FPMIN
			c=b+an/c; if(abs(c)<FPMIN) c=FPMIN
			d=1.d0/d
			Tn=d*c
			igamma = igamma*Tn
			if(abs(Tn-1.d0)<EPS) then
				igamma = igamma * exp(-x) * x**a
				return
			end if
		end do
	end if
End
!
!
Subroutine Shuffle(n, a, s)
	integer n, i, Ipos, Itmp, a(*), s
	real  r

	if(s<0) call random_seed()
	if(s>0) call random_seed(put=[s])

	do i=n, 2, -1
		call random_number(r)
		Ipos = int(r*i) + 1
		Itmp    = a(Ipos)
		a(Ipos) = a(i)
		a(i)    = Itmp
	end do
End Subroutine Shuffle

   

Subroutine readFF(Ftop)
! read input file into txtFF.
! determine number of none comments lines in the input file Nff 
! the lines in txtFF has been cleaned up to make future processing simplier.

use CRYOFFmod
implicit none
character*80, intent(in) :: Ftop
integer :: Level
character*1024 :: txt, tag
integer :: ipass, il, ierr, i, j, k

external Level

open(unit=Ltop, file=trim(adjustl(Ftop))//'.ff', status='old', IOstat=Ierr)
		call ErrStop('Input file  '//trim(adjustl(Ftop))//'.ff  cannot be found!', Ierr)

		print*, '>> Reading Input File'
        
        do ipass=1,2
        ! Figure out number of non commented lines in the .ff file in the first pass then read the file in the second pass.

        if (ipass == 2) then 
        rewind(Ltop)   
        allocate(txtFF(Nff))
        endif
        
		Nff=0 ! line counter for .ff file
		read(Ltop, '(A)', IOstat=Ierr) txt

		do while(Ierr==0)
			txt=trim(adjustl(txt))

			if( txt(1:1)/=';' .and. txt(1:1)/='#' .and. txt(1:1)/='!' &
			&	.and. Len_trim(txt)>0) then
				tag=txt; call upcase(tag)
				i=Level(tag); j=MaxLev
				if(i<0) then
                ! skip commented sections
					do
						read(Ltop, '(A)', IOstat=Ierr) txt
						tag=txt; call upcase(tag); j=Level(tag)
						if(abs(j)<=abs(i)) exit
					end do
					backspace(Ltop)
                else
					i=index(txt, ';'); if(i/=0) txt = txt(1:i-1)
					i=index(txt, '#'); if(i/=0) txt = txt(1:i-1)
					i=index(txt, '!'); if(i/=0) txt = txt(1:i-1)
					i=index(txt, '	')
					do while(i/=0)
						txt(i:i)=' '
						i=index(txt, '	')
					end do

					txt=trim(adjustl(txt))
					il=len_trim(txt)
                    
					if(il>0) then
						if(txt(1:1)=='[') then
							k=index(txt, ']')
							txt(1:1)=' '; txt(k:k)=' '
                        end if
                        
                    txt=trim(adjustl(txt))
 					il=len_trim(txt)
                    
					if(il>0) then                   
                    Nff = Nff+1
					if (ipass == 2) then 
                        tag=txt; call upcase(tag)
                        if(tag(1:3)/='FIL')  call upcase(txt) !convert case except the FILe line
                        txtFF(Nff)=txt
                    endif 
                    endif
                        
                    endif
				end if
			end if
			read(Ltop, '(A)', IOstat=Ierr) txt
        end do
        end do !ipass
		close(Ltop)
End subroutine readFF

Subroutine SetSystemOptions
! Go through txtFF(input .ff file) to set general system options
! This include the keyword (key) lines and the optimization (opt) line. 
use CRYOFFmod
implicit none

integer :: i, ierr, ii, j
character*1024 :: txt
character*1024 :: tag
real*8  :: Rtmp
real*8, external:: getValue,  getBetaKmax, getEdir

!------------Going through txtFF for a second time (read "KEY" directive)
	IniSnp=1; LstSnp=1E9;      ! Number of Snapshots Nsnp 
    ! IniSnp is number of snapshot to skip at the begining of the ref file.
    
	np_pop=1; Iwgt=3; wgtfac=1.0
	Iopt=0; Pgen=2
	deps=0; MaxIter=100
	YesInt=.false.; YesPBC=.false.; YesBug=.false.; YesHyb=.false.
!	YesRest=.false.; Nrest=1
	Rcou=1.d9; Rvdw=1.d9; Etol=1.d-6; Rcon=1.d4*Reps
    ! increased default Rcon. 
    
	Nkcv=1; ItypKCV=0
    

	do i=1, Nff
		read(txtFF(i), '(A)', IOstat=Ierr) txt

       ! Parsing the KEY keyword     
		if(txt(1:3)=='KEY') then
            
        if (MyID==0) print *, "reading keyword directives:"
        if (MyID==0) print *, trim(txt)
        
        Rtmp=getValue(txt, 'NPP')      ; if(Rtmp>0.d0) np_pop =int(Rtmp)
            
        text=' '
        if(index(txt, 'INTE')/=0) then 
            
           YesInt = .true.
            ! When INTE is not found, txt here is left over from the previous call. 
            ! Should probably fix get Option to make txt empty
            call getOption(txt, 'INTE', tag)
			if(len_trim(tag)>0)    text=tag//' '           
        endif
        
        if(index(txt, 'HYBR')/=0) then 
            
           YesHyb = .true.
 
        endif
            
! Two level parallelization. np_pop is for number of CPUs used for non-linear optimization. 
! np_force control the number of CPUs to be used to distribute conformations over number of CPUs.
! NPpop can be used to run differential evoluation over np_pop CPUs. 
! Most likely This part of the code is broken. only np_pop=1 works.
            

			if(index(txt, 'PBC') /=0) YesPBC = .true.            
			Rtmp=getValue(txt, 'RCOU')     ; if(Rtmp>0.d0) Rcou=Rtmp
			Rtmp=getValue(txt, 'RVDW')     ; if(Rtmp>0.d0) Rvdw=Rtmp
			Rtmp=getValue(txt, 'ETOL')     ; if(Rtmp>0.d0) Etol=Rtmp
            beta = getBetaKmax(Etol, Rcou, getEdir)
            
			Rtmp=getValue(txt, 'RCON')*Reps; if(Rtmp>0.d0) Rcon=Rtmp
            if (MyID==0) print*, 'Setting Condition number for SVD to ', Rcon
            
			if(index(txt,'USEFR')/=0) then
            call getOption(txt, 'USEFR', tag)
			read(tag, *, IOstat=Ierr) IniSnp, LstSnp
			if(Ierr/=0) then
				IniSnp=1;
				read(tag, *, IOstat=Ierr) LstSnp
            end if
            end if
            
            if(index(txt,'MAXFR')/=0) then
			call getOption(txt, 'MAXFR', tag)
			read(tag, *, IOstat=Ierr) LstSnp
            end if
            
            if (MyID==0) print *, "The fit will be performed using frames,", IniSnp, " to ",LstSnp,"in the .ref file"
            
! Not sure if normalization is an accurate word.

! Default to Iwgt 3, which is the published method. 
            
			call getOption(txt, 'NORM', tag)
            Rtmp=0.0
			if(len_trim(tag)>0) then
                tag=trim(adjustl(tag))
				if(tag(1:1)=='W' .or. tag(1:3)=='RMS') Iwgt = 3
				if(tag(1:1)=='U' .or.  tag(1:3)=='REL') Iwgt = 4
                read(tag, *, IOstat=Ierr) tag, Rtmp
                if(Rtmp>0.d0 .and. Ierr==0) wgtfac=Rtmp
            end if
            
            if (MyID==0) then
            if (Iwgt == 3) print *,"Weight scheme for the fitting is inverse RMSF."
            if (Iwgt == 3 .and. wgtfac /= 1.0) print *,"Weight scheme weights NETF ",wgtfac, " times more than TORQ."
            if (Iwgt == 3 .and. wgtfac /= 1.0) print *,"This only affects intermolecular fit."
            if (Iwgt == 4) print *,"Fraction weight scheme with forces larger than", wgtfac, " times RMSF has smaller weight."
            endif

            if(index(txt, 'PRI') /=0) then 
                YesBug = .true.
			call getOption(txt, 'PRI', tag)
			if(len_trim(tag)>0) then
				if(index(tag, 'MAT')/=0) YesOut(1)=.true.
				if(index(tag, 'FIT')/=0) YesOut(2)=.true.
                if(index(tag, 'DEBUG')/=0) YesOut(3)=.true.
            end if
            if (MyID==0) print *,"extra printing enabled with option ",trim(tag)
            end if
         
			if(index(txt,'KCV')/=0) then 
            call getOption(txt, 'KCV', tag)
			if(len_trim(tag)>0) then
				read(tag, *) Nkcv
				if(index(tag, 'R')/=0) then
					ItypKCV=1; Iseed=0
					read(tag, *, IOstat=Ierr) tag, tag, Iseed
				end if
            end if
            if (MyID==0) then 
            print *,"Compute Cross Validation with ",Nkcv, "folds"
            if (ItypKCV==0) print *,"folds will be performed in sequence."
            if (ItypKCV==1) print *,"folds will be performed after randomization with seed", Iseed
            if (Nkcv>1) print *,"Full data fits will not be performed."
            endif 
            endif !KCV
            
            if(index(txt,'CMD')/=0) then
                call getOption(txt, 'CMD', tag)
                CMDatmlist=trim(tag)
                YesCMD=.true.
                print*, "Charges will be determined with Charge Matrix Decomposition:", trim(CMDatmlist)
                
                if(len_trim(tag)==0) then
				print*, "Warning: List of charges to be determined with CMD must be provided."
                endif
            endif
            
            !print potential index    EQV
            if(index(txt, 'PIND')/=0) then
                YesPotIdx=.true.
            endif            
            
        end if   !endif key
            
           
        if(txt(1:3)=='OPT' .and. Nopt>0) then
            if (MyID==0) print *, "reading optimziation directives:"
            if (MyID==0) print *, trim(txt)
            
			read(txt, *) tag, tag
			if(index(tag, 'SIM')/=0) Iopt = 1
			if(index(tag, 'POW')/=0) Iopt = 2
			if(index(tag, 'DE') /=0) Iopt = 3
			if(index(tag, '.R') /=0) Pgen = 0
            !. SIM.R means using simplex with read in vertices.
 
            if(index(txt,'INI')/=0) then
                call getOption(txt,'INI',tag) 
                if (index(tag, 'READ') /=0 ) Pgen = 0
                if (index(tag, 'UNIT') /=0 ) Pgen = 1
                if (index(tag, 'FRAC') /=0 ) Pgen = 2
            endif
            
            if (Pgen == 2) optstp=0.1 
            if (Pgen == 1) optstp=1.0
            if(index(txt,'STEP')/=0) then
                Rtmp=getValue(txt,'STEP'); if (Rtmp > 0.0) optstp=Rtmp
            endif 
 
 
            
!			read(txt, *) tag, tag, tag
!			if(tag(1:3)=='RES') then
!				YesRest=.true.
!				ii=index(tag, '=')+1
!				read(tag(ii:), *) Nrest
!           end if
! Nrest is a keyword for differential evolution.
! I do not believe this is ever used.
            
!           Control non-linear optimization convergence. 
            
! Yeps removed. It terminates nonlinear optimization when the lowest objective function is less than this value. 
! It does not make much sense for simplex and Powell. 
! For Differential Evolution, the fitness of the best population is A termination criterion. 
! In that case, we will reinterprete deps. 
            if (Nkcv>1) then 
                Iopt=0 
                if (MyID==0) print*, "nonlinear optimziation disabled for cross validation"
            endif 
            
            if (MyID==0) then
            if (Iopt==1) print*, "Simplex Optimization"
            if (Iopt==2) print*, "Powell Optimization"
            if (Iopt==3) print*, "Differential Evolution"
            if (Pgen==0) print*, "Read initial points for nonlinear optimizations"
            if (Pgen==1) print*, "Generate Initial Vectors of length", optstp
            if (Pgen==2) print*, "Generate Initial Vectors with frational scaling with step size", optstp
            endif
            
            Rtmp=getValue(txt, 'DP');    if(Rtmp>0.d0) deps=Rtmp
			Rtmp=getValue(txt, 'CONV');    if(Rtmp>0.d0) deps=Rtmp
            Rtmp=getValue(txt, 'MAXIT'); if(Rtmp>0.d0) MaxIter=int(Rtmp)

			Mopt = Nopt+1
			if(Iopt==1) allocate(P(Mopt, Nopt), Y(Mopt), stat=Ierr)
			if(Iopt==2) allocate(Xi(Nopt, Nopt), stat=Ierr)
			call ErrStop('P/Y/Xi Memory Allocation Error', Ierr)

! Read in initial values for various type of non-linear optimizations. 
			if(Pgen==0) then
				if(Iopt==1) then
                    print*, "reading", Mopt, "vectors of length ", Nopt
					do ii=1, Mopt
						read(txtFF(i+ii), *) (P(ii, j), j=1, Nopt)
					end do
				else if(Iopt==2) then
					do ii=1, Nopt
						read(txtFF(i+ii), *) (Xi(j, ii), j=1, Nopt)
					end do
				else if(Iopt==3) then
					read(txtFF(i+1), *) (dei(j), j=1, 6)
					read(txtFF(i+2), *) (der(j), j=1, 3)
				end if
			end if
		end if
    end do
    
        if (MyID==0) then
        print *, "The non-linear optimization will be performed using ",np_pop," parallel processes."
        np_force=np_world/np_pop
        print *, "Each process will use ",np_force, "CPUs for force evaluations."
        print *, "The convergence is set to ",deps
        print *, "The maximum number of optimization steps is ",MaxIter
        
        if (YesInt .and. YesHyb) then 
            call ErrStop('Intermolecular fit and hybrid fit can not be requested together.', 1)
        endif 
        
        if (YesInt) then 
            print*, 'Intermolecular fit'
        else if (YesHyb) then
            print*, 'Hybrid Fit including both Intermolecular and Intramolecular forces'
        else
            print*, 'Fitting atomic forces only (ignore net force and net torque)'
        end if
         
        if(len_trim(text)>0)    then           
             print *,"selective fit directive for intermolecular terms: ",trim(text)
             print *,"This is an experimental feature, the use of solvation factor in .ref file is preferred."
        endif
        
        if(YesPBC) then             
			print *, "Periodic boundary condition turned on. Ewald summation will be used for Coulobmic interactions"
            print *, "Rcou is set to ", Rcou
            print *, "Rvdw is set to ", Rvdw
            print *, "Ewald tolerance is set to ",Etol
            print *, "Ewald beta is computed to be ",beta
        endif
        endif
       
End Subroutine SetSystemOptions

Subroutine ReadMolecularDef
use CRYOFFmod , NffTrue => Nff
implicit none
! Read molecular definition and intramolecular terms. 

integer Nff, i, j, k, ii ,Ierr, Itmp, Jtmp
integer Natm, Nvir, Ntyp, Ngrp, Ntmp, Ityp
integer Nlnr
character*1024 :: txt,tag,txt2
character*80 :: txtPara(MaxPara)
real*8  dtemp, ftmp1, ftmp2
logical NoCouTyp, isequiv

    NoCouTyp=.false.
    
    Nff=1; Nopt=0
    do
        read(txtFF(Nff), *) txt
        if(txt(1:3)=='MOL') exit
        Nff=Nff+1
    end do
    
    if(YesOut(3) .and. MyID==0) print*, trim(txt)

    NatmTyp=0

! The idea of "text" to allow only netforce on some molecules to be fit.
! The default is YesNetF=.true. unless it is overwritten with INTER=MOLNAME
! I note YesNetF here are arrays.
! I believe a cleaner way to implment such a function is to use different solvation factor.
! Maybe it is a good idea to remove this function.
! might have been better to move this to the actual reading of INTE
    
    if(len_trim(text)==0) then
        YesNetF(:)=.true.; YesTorq(:)=.true.
    else
        YesNetF(:)=.false.; YesTorq(:)=.false.
    end if

! for each MOL type readin its potential parameter
    do i=1, NmolTyp
        isequiv=.false.
        Mole(i)%EquivID=i;
        read(txtFF(Nff), *, IOstat=Ierr) txt, Mole(i)%Name, txt2, tag
        if(YesOut(3) .and. MyID==0) print*, trim(txt)
        if (index(txt2, 'EQU')==1 ) then 
            if (MyID==0) print*, Mole(i)%Name, trim(txt2), tag
            isequiv=.true.
            do j=1,i-1
                if (Mole(j)%Name .eq. tag) then 
                    Mole(i)%EquivID=j
                    exit
                endif 
            enddo
            if (Mole(i)%EquivID==i) call ErrStop("Equivalent Molecule not found!",i)
        endif 
        
        Nff=Nff+1
        read(txtFF(Nff), *) txt, Natm; Nff=Nff+1
        NatmTyp=NatmTyp+Natm
        ! NatmTyp is the total number of atoms in all molecules.
        ! It is not the number of different atom types since each molecule can have more than one atoms of the same type
        Mole(i)%NmolTot = 0
        Mole(i)%Natm = Natm
        allocate( Mole(i)%Satm(Natm), Mole(i)%VDWatm(Natm), Mole(i)%COUatm(Natm), &
                & Mole(i)%YesExc(Natm, Natm), Mole(i)%YesPair(Natm, Natm) )

        !Find MaxVir for allocate memory. 
        !MaxVir is the number of vectors needed to calculate the virtual site. (not max number of virtual sites)
        Nvir=0; MaxVir=0
        do j=1, Natm
            read(txtFF(Nff), '(A)', IOstat=Ierr) txt
            k=index(txt, '*')
            if(k/=0) then
                if (j==1) then
                    call ErrStop("The first atoms of a molecule can not be virtual site. Define virtual site after physical atom sites",1)
                endif
                txt(k:k)=''         !virtual site debug wang
                ! remove : from input txt. Might best write a funciton for it. 
                do ii=k, len_trim(txt)
                    if(txt(ii:ii)==':') txt(ii:ii)=''
                end do
                if (NoCouTyp) then
                    read(txt, *, IOstat=Ierr) tag, tag, Itmp
                else
                    read(txt, *, IOstat=Ierr) tag, tag, tag, Itmp
                endif
                call ErrStop("error reading virtual site", Ierr)
                
                Nvir = Nvir+1
                MaxVir = max(MaxVir, Itmp)
            end if
            ! The input allows COUatm and VDWatm to be different. I am not sure if this is a direction we want to take. 
            read(txt, *, IOstat=Ierr) tag, Mole(i)%VDWatm(j), Mole(i)%COUatm(j)
            if (Ierr/=0 .or. NoCouTyp) then
                read(txt, *, IOstat=Ierr) tag, Mole(i)%VDWatm(j)
                Mole(i)%COUatm(j)=Mole(i)%VDWatm(j)
                NoCouTyp=.true.
            endif 
            call ErrStop('Atom Line '//trim(txt), Ierr)
            Nff=Nff+1
        end do

		Ntyp = 0; Ngrp = 0
		do
			read(txtFF(Nff), *, IOstat=Ierr) txt, Ntmp
			if( Ierr/=0 .or. index(txt, 'BON')+index(txt, 'ANG')+index(txt, 'DIH')+index(txt, 'CDIH') &
				& +index(txt, 'BD3')+index(txt, 'EXC')+index(txt, 'PAI')+index(txt, 'FUD')==0) exit

			Nff=Nff+1
			if(txt(1:3)=='EXC' .or. txt(1:3)=='PAI' .or. txt(1:3)=='FUD') then
				Nff = Nff+Ntmp
			else
				Ntyp = Ntyp+Ntmp
				do j=1, Ntmp
					read(txtFF(Nff), *) txt, Itmp; Nff=Nff+1
					Ngrp = max(Ngrp, Itmp)
					Nff = Nff+Itmp
				end do
			end if
		end do

		Mole(i)%Ntyp = Ntyp
		Mole(i)%Nvir = Nvir
		allocate( Mole(i)%Ttyp(Ntyp), Mole(i)%Ityp(Ntyp), Mole(i)%Ngrp(Ntyp), &
			& Mole(i)%Npar(Ntyp), Mole(i)%Nlnr(Ntyp), Mole(i)%Ilnr(Ntyp), Mole(i)%ParaRMS(Ntyp, MaxPara), &
			& Mole(i)%YesFit(Ntyp), Mole(i)%Para(Ntyp, MaxPara), &
			& Mole(i)%Iatm(Ntyp, Ngrp, MaxAdj), &
			& Mole(i)%IatmVir(Nvir), Mole(i)%NatmVir(Nvir), &
			& Mole(i)%Ivir(Nvir, MaxVir), Mole(i)%Rvir(Nvir, MaxVir), stat=Ierr)
		call ErrStop('All Mole(i)', Ierr)

! This also part of the selective enable of fitting net force or torque for certain types of molecules.
! It will be easy to realize this using solvation factors. 
        
        if(    index(text, trim(Mole(i)%Name)//' ')/=0 &
		& .or. index(text, trim(Mole(i)%Name)//'.NETF')/=0) YesNetF(i)=.true.
		if(    index(text, trim(Mole(i)%Name)//' ')/=0 &
		& .or. index(text, trim(Mole(i)%Name)//'.TORQ')/=0) YesTorq(i)=.true.

		if(MyID==0) then
			print*, '   >> Molecule Type', i, ': '//trim(Mole(i)%Name)
			if(YesInt .and. YesNetF(i) .and. YesTorq(i)) then
				print*, '      Inter-Fitting using:    Net Force and Torque'
			else if(YesInt .and. YesNetF(i) ) then
				print*, '      Inter-Fitting using:    Net Force'
			else if(YesInt .and. YesTorq(i) ) then
				print*, '      Inter-Fitting using:    Torque'
			end if

			print*, '      No. of Atoms:          ', Mole(i)%Natm
			if(Mole(i)%Nvir>0) print*, '      No. of Virtual Sites:  ', Mole(i)%Nvir, '    No. of vec. that defines the vsite', MaxVir
			if(Mole(i)%Ntyp>0) print*, '      No. of Intra-Pot Types:', Mole(i)%Ntyp, '    Max of Groups:', Ngrp, LF
        end if

        ! rewind to atoms
        
		do
			read(txtFF(Nff), *) txt
			if(index(txt, 'ATO')/=0) exit
			Nff = Nff-1
        end do

        ! Read Virtual sites
		Nff=Nff+1; Nvir=0
		do j=1, Natm
			read(txtFF(Nff), '(A)', IOstat=Ierr) txt
			k=index(txt, '*')
			if(k/=0) then
				txt(k:k)=''
				do ii=k, len_trim(txt)
					if(txt(ii:ii)==':' .or. txt(ii:ii)=='+' &
					& .or. txt(ii:ii)=='(' .or. txt(ii:ii)==')') txt(ii:ii)=''
				end do

				Nvir = Nvir+1
                if (NoCouTyp) then 
                    read(txt, *) tag, tag, Itmp, &
					& (Mole(i)%Rvir(Nvir,k), Mole(i)%Ivir(Nvir,k), k=1, Itmp)
                else
				read(txt, *) tag, tag, tag, Itmp, &
					& (Mole(i)%Rvir(Nvir,k), Mole(i)%Ivir(Nvir,k), k=1, Itmp)
                endif
				Mole(i)%IatmVir(Nvir)=j
				Mole(i)%NatmVir(Nvir)=Itmp
			end if
			Nff=Nff+1
		end do

		if(MyID==0) then
			print*, '      First Atom:', 1,    '      VDW: '//trim(Mole(i)%VDWatm(1))   //'      COU: '//trim(Mole(i)%COUatm(1))
			print*, '      Last  Atom:', Natm, '      VDW: '//trim(Mole(i)%VDWatm(Natm))//'      COU: '//trim(Mole(i)%COUatm(Natm)), LF

			if(Nvir>0) then
				Itmp=Mole(i)%IatmVir(1); Ntmp=Mole(i)%NatmVir(1)
				write(*, '(7X, 2(A,I4))') 'First Vsite:', 1, ' Atom:', Itmp
				write(*, '(11X, A, 99(I4, F9.6))')  'Virtual Site Combination Rule:', (Mole(i)%Ivir(1,k), Mole(i)%Rvir(1,k), k=1, Ntmp)

				Itmp=Mole(i)%IatmVir(Nvir); Ntmp=Mole(i)%NatmVir(Nvir)
				write(*, '(7X, 2(A,I4))') 'Last  Vsite:', Nvir, ' Atom:', Itmp
				write(*, '(11X, A, 99(I4, F9.6))')  'Virtual Site Combination Rule:', (Mole(i)%Ivir(Nvir,k), Mole(i)%Rvir(Nvir,k), k=1, Ntmp)
				print*
			end if
        end if

        ! finished reading molecule definition continue with interactions. 
        
		Ntyp = 0; Ngrp = 0
		Mole(i)%YesExc(:,:)=.false.
		Mole(i)%YesPair(:,:)=.false.
        do
            read(txtFF(Nff), *, IOstat=Ierr) txt

            Ityp = -1
            if(txt(1:3)=='EXC') Ityp=0
            if(txt(1:3)=='PAI') Ityp=1
            if(txt(1:3)=='FUD') Ityp=1
            if(txt(1:3)=='BON') Ityp=2
            if(txt(1:3)=='ANG') Ityp=3
            if(txt(1:3)=='DIH') Ityp=4
            if(txt(1:4)=='CDIH') Ityp=5
            if(txt(1:3)=='BD3') Ityp=6  !bonded 3 body term including bond-bond and bond-angle coupling
            if(Ityp<0) exit
            
            read(txtFF(Nff), *, IOstat=Ierr) txt, Ntmp
            call ErrStop(trim(txtFF(Nff)), Ierr)

			if(Ityp<=1) then
            ! pair list or exclusion list
				if(Ityp==1) read(txtFF(Nff), *, IOstat=Ierr) &
				&		txt, txt, Mole(i)%FudgeVDW, Mole(i)%FudgeQQ
                
! Ngrp here is number of terms in this type of interaction. 
! Be careful, I do not think meaning of Ngrp is the same accross the code. 
                
				Ngrp = Ntmp
				do k=1, Ngrp
                    
                    if (Ityp==0) then 
					read(txtFF(Nff+k), *, IOstat=Ierr) Itmp, Jtmp
					call ErrStop(trim(txtFF(Nff+k)), Ierr)
                        if (Mole(i)%YesPair(Itmp, Jtmp)) then 
                            print *,"Warning: excluding pairs that are fudged (self contradictory ff file)", Itmp, Jtmp
                        endif
                        Mole(i)%YesExc(Itmp, Jtmp)=.true.
						Mole(i)%YesExc(Jtmp, Itmp)=.true.
                    else !Ityp == 1
                    read(txtFF(Nff+k), *, IOstat=Ierr) Itmp, Jtmp, ftmp1, ftmp2
                        if (Mole(i)%YesExc(Itmp, Jtmp)) then 
                            print *,"Warning: fudge pairs that are excluded (self contradictory ff file)", Itmp, Jtmp
                        endif
                        if (Mole(i)%FudgeVDW/=ftmp1 .or. Mole(i)%FudgeQQ/=ftmp2) then
                            call ErrStop("FudgeVDW and FudgeQQ must be the same for all pairs in each molecule",1)
                        endif
                        Mole(i)%YesPair(Itmp, Jtmp)=.true.
						Mole(i)%YesPair(Jtmp, Itmp)=.true.
                    endif 
                end do
				Nff = Nff+Ntmp+1
            else !Ityp>1
            ! read intramolecular terms. 
                
				Nff=Nff+1
				do j=1, Ntmp
					Ntyp = Ntyp+1
					read(txtFF(Nff), *, IOstat=Ierr) tag
					call ErrStop(trim(txtFF(Nff)), Ierr)

					call NumbPara(tag, Npar, Nlnr)
					if(Npar+Nlnr==0) call ErrStop('Potential Type NOT defined! >>> '//trim(tag), 1)

					read(txtFF(Nff), *, IOstat=Ierr) txt, Ngrp, txt, (txtPara(k), k=1, Npar)
					call ErrStop(trim(txtFF(Nff)), Ierr)

					Mole(i)%Ityp(Ntyp) = Ityp
					Mole(i)%Npar(Ntyp) = Npar
					Mole(i)%Nlnr(Ntyp) = Nlnr
					Mole(i)%Ngrp(Ntyp) = Ngrp
					Mole(i)%Ttyp(Ntyp) = trim(tag)
					Mole(i)%YesFit(Ntyp) = .false.
                    
                    if (IsEquiv.and.Ityp>1) then
                        ii=Mole(i)%EquivID
                        if (Mole(ii)%Ityp(Ntyp) /= Ityp) IsEquiv=.false.
                        if (Mole(ii)%Ngrp(Ntyp) /= Ngrp) IsEquiv=.false.
                        if (.not. IsEquiv) call ErrStop("Equivalent molecules must have identical bonded terms in identical order.",1)
                    endif 
                        
					if(index(txt, 'FIT')/=0) Mole(i)%YesFit(Ntyp) = .true.

					do k=1, Npar
						ii=index(txtPara(k), '_[')
						if(ii/=0) then
							txtPara(k)(ii:ii+1)='  '
							do Itmp=ii, len_trim(txtPara(k))
								if(    txtPara(k)(Itmp:Itmp)==':'  &
								& .or. txtPara(k)(Itmp:Itmp)==']') &
								&      txtPara(k)(Itmp:Itmp)=' '
							end do

							if(Mole(i)%YesFit(Ntyp)) then
								Nopt = Nopt+1
								Popt(Nopt)%P => Mole(i)%Para(Ntyp, k)
                                
                                ii=index(txtPara(k),'R')
                                if (ii/=0) then
                                    txtPara(k)(ii:ii)=' '
                                    PLR(Nopt)=.true. !Restraint
                                    read(txtPara(k), *,IOstat=Ierr) dtemp, PRC(Nopt), PRP(Nopt), dP(Nopt)
                                    if (Ierr/=0) then 
                                          if (Pgen==1) dP(Nopt)=optstp
                                          if (Pgen==2) dP(Nopt)=dtemp*optstp
                                    endif
                          
                                else
                                    PLR(Nopt)=.false.

								read(txtPara(k), *) dtemp, Pinf(Nopt), Psup(Nopt)
								if(Iopt==1 .or. Iopt==2) then 
                                    read(txtPara(k), *, IOstat=Ierr) txt, txt, txt, dP(Nopt)
                                            if (Ierr/=0) then 
                                                if (Pgen==1) dP(Nopt)=optstp
                                                if (Pgen==2) dP(Nopt)=dtemp*optstp
                                            endif
                                endif
                                endif
							end if
						end if
! 						print*, 'txtPara', ii, Nopt, k, Npar, trim(txtPara(k))
						read(txtPara(k), *) Mole(i)%Para(Ntyp, k)
                    end do

                    ! unit conversion for angles.
                    if(Ityp==3) Mole(i)%Para(Ntyp,1) = Mole(i)%Para(Ntyp,1)*Deg2Rad ! input unit of Angle is Deg
                    if(Ityp==4 .and. tag(1:3)=='NCO')  Mole(i)%Para(Ntyp,3) = Mole(i)%Para(Ntyp,3)*Deg2Rad ! input unit of Angle is Deg
                    if(Ityp==4 .and. tag(1:3)=='COS')  Mole(i)%Para(Ntyp,1) = Mole(i)%Para(Ntyp,1)*Deg2Rad ! input unit of Angle is Deg
                    if(Ityp==4 .and. tag(1:3)=='HAR')  Mole(i)%Para(Ntyp,1) = Mole(i)%Para(Ntyp,1)*Deg2Rad ! input unit of Angle is Deg
                    if(Ityp==5 .and. tag(1:3)=='CNC')  Mole(i)%Para(Ntyp,4) = Mole(i)%Para(Ntyp,4)*Deg2Rad
                    if(Ityp==5 .and. tag(1:3)=='CCO')  Mole(i)%Para(Ntyp,1) = Mole(i)%Para(Ntyp,1)*Deg2Rad
                    do k=1, Ngrp
! Always try to read 5 atoms for each term. MaxAdj = 5.
                        read(txtFF(Nff+k), *, IOstat=Ierr) (Mole(i)%Iatm(Ntyp, k, ii), ii=1,MaxAdj)
! 						if(k==1 .or. k==Ngrp) print*, '    >>', Mole(i)%Iatm(Ntyp, k, :)
                    end do
                    Nff = Nff+Ngrp+1
                    end do
                    end if
                end do

		if(MyID==0 .and. Ntyp>0) then
			write(*, '(7X, A, I4, A, I2, A, A5, A, I4)') 'First Intra-Potential:', 1, &
				& ' Type:', Mole(i)%Ityp(1), ' Function:', trim(Mole(i)%Ttyp(1)), '  Group:', Mole(i)%Ngrp(1)
			write(*, '(19X, A, 9F18.9)') 'Parameter:', Mole(i)%Para(1,1:Mole(i)%Npar(1))
			write(*, '(7X, A, I4, A, I2, A, A5, A, I4)') 'Last  Intra-Potential:', Ntyp, &
				& ' Type:', Mole(i)%Ityp(Ntyp), ' Function:', trim(Mole(i)%Ttyp(Ntyp)), '  Group:', Mole(i)%Ngrp(Ntyp)
			write(*, '(19X, A, 9F18.9)') 'Parameter:', Mole(i)%Para(Ntyp,1:Mole(i)%Npar(Ntyp))
			print*
		end if
    end do ! for NmolTyp
    
    ! Done reading Intramolecular
End subroutine ReadMolecularDef

Subroutine ReadIntermolecularTerms
use CRYOFFmod , NffTrue => Nff
use extrafuncs

implicit none
! Read nonbonded terms. 

integer Nff, i, j, k, ii, Ierr, Itmp, Jtmp
integer Ngrp, Nlnr, Ntmp
real*8 dtemp
character*1024 :: txt,tag
character*80 :: txtPara(MaxPara), tagIatm, tagJatm

	if(MyID==0) print*, '>> Parsing Inter-Potential'

	Nff=1
	do
		read(txtFF(Nff), *, IOstat=Ierr) txt
		if(Ierr/=0 .or. index(AllPairType,txt(1:3))/=0) exit
		Nff=Nff+1
    end do

	Itmp=Nff; Npair=0; NpairQQ=0
	do
        ! Count total number of intermolecular interactions of all types to be read in
		read(txtFF(Itmp), *, IOstat=Ierr) txt, Ngrp
		if(Ierr/=0 .or. index(AllPairType, txt(1:3))==0) exit
		Itmp  = Itmp+1+Ngrp
		Npair = Npair+Ngrp
		if(txt(1:3)=='COU' .or. txt(1:3)=='THC') NpairQQ=NpairQQ+Ngrp
    end do
! Is PEW particle mesh ewald? If so, it should also check for PBC, right? 
    
	if(MyID==0) print*, ' Total Number of Inter-Molecular Interactions: ', Npair
    if(MyID==0) print*, ' Out of which, the number of coulobmic interactions (counting THC) is ',NpairQQ

	if(Npair>0) then
		allocate( TtypPair(Npair), TatmPair(Npair), IcouPair(Npair,2), IvdwPair(Npair,2), ParaPair(Npair, MaxPara), ParaRMS(Npair, MaxPara),&
			&	  NparPair(Npair), NlnrPair(Npair), IlnrPair(Npair),YesFitPair(Npair), Iuniq(Npair), stat=Ierr)   !EQV
		call ErrStop('Pair', Ierr)

		Npair = 0
		YesFitPair = .false.
		do
			read(txtFF(Nff), *, IOstat=Ierr) tag, Ngrp
			Nff=Nff+1
			if(Ierr/=0 .or. index(AllPairType, tag(1:3))==0) exit
            if(MyID==0) then 
                print*, "Reading ", tag(1:3)
                if (YesOut(3)) then 
                    print*, "First Line: ", trim(txtFF(Nff+1))
                    print*, "Last Line: ", trim(txtFF(Nff+Ngrp))
                endif 
            endif 

			call NumbPara(tag(1:3), Npar, Nlnr)
			if(Npar+Nlnr==0) call ErrStop('Interaction type unknown! >>> '//trim(tag), 1)

			do i=1, Ngrp
				Ntmp = Npair+i
				read(txtFF(Nff), *) TagIatm, TagJatm, txt

				NparPair(Ntmp) = Npar
				NlnrPair(Ntmp) = Nlnr
				TtypPair(Ntmp) = trim(tag)
				TatmPair(Ntmp) = trim(TagIatm)//'~'//trim(TagJatm)

				IcouPair(Ntmp, 1)=findatm(NcouTyp, TagIatm, COUatmTyp(1:NcouTyp))
				IcouPair(Ntmp, 2)=findatm(NcouTyp, TagJatm, COUatmTyp(1:NcouTyp))

				IvdwPair(Ntmp, 1)=findatm(NvdwTyp, TagIatm, VDWatmTyp(1:NvdwTyp))
				IvdwPair(Ntmp, 2)=findatm(NvdwTyp, TagJatm, VDWatmTyp(1:NvdwTyp))

				if(index(txt, 'FIT')/=0) YesFitPair(Ntmp) = .true.

				read(txtFF(Nff), *) txt, txt, txt, (txtPara(j), j=1, Npar)

				do j=1, Npar
					k=index(txtPara(j), '_[')
					if(k/=0) then
						txtPara(j)(k:k+1)='  '
						do Itmp=k, len_trim(txtPara(j))
							if(    txtPara(j)(Itmp:Itmp)==':'  &
							& .or. txtPara(j)(Itmp:Itmp)==']') &
							&      txtPara(j)(Itmp:Itmp)=' '
						end do

						if(YesFitPair(Ntmp)) then
							Nopt = Nopt+1
							Popt(Nopt)%P => ParaPair(Ntmp, j)
                            
                         ii=index(txtPara(j),'R')
                                if (ii/=0) then
                                    txtPara(j)(ii:ii)=' '
                                    PLR(Nopt)=.true. !Restraint
                                    read(txtPara(j), *,IOstat=Ierr) dtemp, PRC(Nopt), PRP(Nopt), dP(Nopt)
                                    if (Ierr/=0) then 
                                          if (Pgen==1) dP(Nopt)=optstp
                                          if (Pgen==2) dP(Nopt)=dtemp*optstp
                                    endif
                          
                                else
                                    PLR(Nopt)=.false.
                                    
							read(txtPara(j), *) dtemp, Pinf(Nopt), Psup(Nopt)
							if(Iopt==1 .or. Iopt==2) then 
                                 read(txtPara(j), *, IOstat=Ierr) txt, txt, txt, dP(Nopt)
                                        if (Ierr/=0) then 
                                            if (Pgen==1) dP(Nopt)=optstp
                                            if (Pgen==2) dp(Nopt)=dtemp*optstp
                                         endif
                            endif
                                endif
                                
						end if
					end if
					read(txtPara(j), *) ParaPair(Ntmp, j)
				end do
				Nff  = Nff+1
			end do
			Npair = Npair+Ngrp
		end do

		if(NpairQQ>0) then
			allocate(IpairQQ(NpairQQ), tagQQ(2*NpairQQ), TatmIpairQQ(NpairQQ), stat=Ierr)
			call ErrStop('PairQQ', Ierr)

            if (YesPBC) YesEwa=.true.
			NpairQQ=0; Nqq=0
			do i=1, Npair
				txt=TtypPair(i)
				if(txt(1:3)=='COU' .or. txt(1:3)=='THC' ) then
					NpairQQ=NpairQQ+1
                    ! IpairQQ is for charge constraint
					IpairQQ(NpairQQ)=i
                    TatmIpairQQ(NpairQQ)=Tatmpair(i)

					if(YesEwa) then
						tag=TatmPair(i)
						k=index(tag, '~')
						TagIatm=tag(1:k-1)
						TagJatm=tag(k+1:)

						YesQQ=.true.
						do j=1, Nqq
							if(tagQQ(j)==TagIatm) then
								YesQQ=.false.; exit
							end if
						end do
						if(YesQQ) then
							Nqq=Nqq+1
							tagQQ(Nqq)=TagIatm
						end if

						YesQQ=.true.
						do j=1, Nqq
							if(tagQQ(j)==TagJatm) then
								YesQQ=.false.; exit
							end if
						end do
						if(YesQQ) then
							Nqq=Nqq+1
							tagQQ(Nqq)=TagJatm
						end if
                    end if !YesEwa, Compute Nqq and create tagQQ for ewald
                    ! Looks line tagQQ is all the atoms involved in ewald sum. 
                    ! Isn't it the same as number of atoms with partial charges? 
				end if
            end do
 
		end if

		if(MyID==0) then
			write(*, '(4X, A, I4, A, A5, A, A10)') 'First Inter-Potential:', 1, &
			&	' Function:', trim(TtypPair(1)), ' Atom:', trim(TatmPair(1))
			write(*, '(16X, A, 99F18.9)') 'Parameter:', (ParaPair(1,j), j=1,NparPair(1))
			write(*, '(4X, A, I4, A, A5, A, A10)') 'Last  Inter-Potential:', Npair, &
				& ' Function:', trim(TtypPair(Npair)), ' Atom:', trim(TatmPair(Npair))
			write(*, '(16X, A, 99F18.9)') 'Parameter:', (ParaPair(Npair,j), j=1,NparPair(Npair))
			print*
		end if
    end if
    
    ! done reading intermolecular interactions
End subroutine ReadIntermolecularTerms   
    
Subroutine ReadChargeConstraint
use CRYOFFmod , NffTrue => Nff
use extrafuncs

implicit none
integer :: Nff, Ierr,i,j,k
real*8 :: Wgt
character*1024 :: txt

if(MyID==0) print*, '>> Parsing Charge Constraints'

	Nff=1; Ierr=0
	do while(Ierr==0 .and. Nff<size(txtFF))
		read(txtFF(Nff), *, IOstat=Ierr) txt
		if(index(txt, 'CSTR')+index(txt, 'CHAR')/=0) exit
		Nff=Nff+1
	end do

	Ncst=0
	if(Nff<size(txtFF)) read(txtFF(Nff), *, IOstat=Ierr) txt, Ncst

	if(MyID==0) print*, '   No. of Charge Constraints:', Ncst

	if(Ncst>0) then
		allocate( NchgGrp(Ncst), Fchg(Ncst), Wchg(Ncst), Nchg(Ncst, MaxChg), Ichg(Ncst, MaxChg) )
		Nff=Nff+1
		do i=1, Ncst
			txt=trim(adjustl(txtFF(Nff)))
			do j=1, len_trim(txt)
				if(txt(j:j)==':' .or. txt(j:j)=='+' &
				& .or. txt(j:j)=='(' .or. txt(j:j)==')' &
				& .or. txt(j:j)=='=') txt(j:j)=''
			end do
			read(txt, *) j, (Nchg(i, k), Ichg(i, k), k=1,j), Fchg(i), Wgt
            if (j>MaxChg) then
                Call ErrStop('Too many terms in charge constraint',j)
            endif 
            Wchg(i)=Wgt
			! Wchg(i)=sqrt(Wgt)
			NchgGrp(i) = j
			Nff=Nff+1
		end do

		if(MyID==0) then
			do i=1, Ncst
				if(i==1 .or. i==Ncst) then
					txt='Last'; if(i==1) txt='First'
					write(*, '(14X, A, I4, 99(I3, A))', advance='no') txt(1:5)//' Equation:', i, &
					&	 ( Nchg(i,j), ' '//trim(TatmPair(IpairQQ(Ichg(i,j))))//' + ', j=1,NchgGrp(i)-1 )
					j=NchgGrp(i); txt=TatmPair(IpairQQ(Ichg(i,j)))
					write(*, '(I4, A, F6.3, A, F9.3)') Nchg(i,j), ' '//trim(txt)//' = ', Fchg(i), '  Weight: ', Wchg(i)
				end if
			end do
			print*
		end if
    end if

End subroutine ReadChargeConstraint

Subroutine ReadSetEQVPotential()			!EQV
use CRYOFFmod , NffTrue => Nff
use extrafuncs
implicit none
integer :: Nff, Ierr,i,j, Neqv, tagNum, itemp, Ntyp, NshiftTot
integer, allocatable::  Tag(:),potmap(:),parashift(:)
character*10, allocatable:: MoleName(:)
character*1024 :: txt

    if(MyID==0) print*, '>> Parsing Equivalent potentials'
    Nff=1; Ierr=0
	do while(Ierr==0 .and. Nff<size(txtFF))
		read(txtFF(Nff), *, IOstat=Ierr) txt
		if(index(txt, 'EQV')/=0) exit
		Nff=Nff+1
    end do
    Neqv=0
	if(Nff<size(txtFF)) read(txtFF(Nff), *, IOstat=Ierr) txt, Neqv
	if(MyID==0) print*, '   No. of Equivalent nonbonded potentials:', Neqv
	if(Neqv>0) then
        YesEQV=.true.
        
        !original potential map 
        allocate(potmap(NlnrTot), parashift(NlnrTot))
        do i=1,NlnrTot
            potmap(i)=i
        enddo
        
        !update potential map
		Nff=Nff+1
		do i=1, Neqv
			txt=trim(adjustl(txtFF(Nff)))
            call getNumOftag(txt,tagNum) !get number of potential index in each line
            allocate(Tag(tagNum))
            read(txt,*) Tag(1:tagNum)
            do j=2,tagNum
                if(Tag(j) > Tag(1)) then
                    if (potmap(Tag(j)) /= Tag(j)) call ErrStop ("Equivalence Redefined",Tag(j))
                    potmap(Tag(j))=potmap(Tag(1))
                else
                    call ErrStop ("The first index should be the smallest.",Tag(1))
                endif
            enddo
            deallocate(Tag)
            Nff=Nff+1
        end do
    endif

    ! Take care of EQV: 
    if(YesEQV) then
	    !first pass EQV, update NlnrTot, Ilnr and parashift
        NshiftTot=0; NlnrTot=0; parashift(1:NlnrTot)=0
	    do i=1, NmolTyp
	    	Ntyp = Mole(i)%Ntyp
	    	do j=1, Ntyp
	    		if(Mole(i)%YesFit(j)) then
                    itemp=Mole(i)%Ilnr(j)
                    if (itemp /= potmap(itemp) ) then
                        Mole(i)%Ilnr(j)=potmap(itemp)
                        NshiftTot=NshiftTot+Mole(i)%Nlnr(j)
                    else
                        NlnrTot = NlnrTot+Mole(i)%Nlnr(j)
                        parashift(itemp)=NshiftTot
                    endif
                end if
	    	end do
        end do
        
	    do i=1, Npair
	    	if(YesFitPair(i)) then
                itemp=Ilnrpair(i)
                if(itemp /=potmap(itemp)) then
	    		    IlnrPair(i) = potmap(itemp)
                    NshiftTot=NshiftTot+NlnrPair(i)
                else
                    NlnrTot = NlnrTot+NlnrPair(i)
                    parashift(itemp)=NshiftTot
                endif
	    	end if
        end do
        ! done first pass EQV
        
        ! second pass EQV, remove empty slots. 
        do i=1, NmolTyp
	    	Ntyp = Mole(i)%Ntyp
	    	do j=1, Ntyp
	    		if(Mole(i)%YesFit(j)) then
                    itemp = Mole(i)%Ilnr(j)
	    			Mole(i)%Ilnr(j) = Mole(i)%Ilnr(j)-parashift(itemp)
	    		end if
            end do
        end do
        
	    do i=1, Npair
	    	if(YesFitPair(i)) then
                itemp=Ilnrpair(i)
                IlnrPair(i) = IlnrPair(i) - parashift(itemp)
	    	end if
        end do
        ! done second pass EQV
    endif !End EQV

End subroutine ReadSetEQVPotential

Subroutine getNumOftag(tag,Tagnum)				!EQV
    use CRYOFFmod
    implicit none
	integer i, Tagnum
	character(*) tag
! Get the number of fields
     
    do i=1,len_trim(tag)
        if(tag(i:i)==':' .or. tag(i:i)=='+' &
            & .or. tag(i:i)=='(' .or. tag(i:i)==')' &
            & .or. tag(i:i)=='=') tag(i:i)=''
    enddo
    Tagnum=0
    if(len_trim(tag)>0) then
        do i=1,len_trim(tag)
            if(i>1) then
            if (tag(i:i)==' ' .and. tag(i-1:i-1)==' ') then
                tag(i:i)=''
            elseif(tag(i:i)==' ') then
                Tagnum=Tagnum+1
            endif
            if (i==len_trim(tag)) then
                Tagnum=Tagnum+1
            endif
            endif
        enddo
    endif

End Subroutine getNumOftag 
    
subroutine ReadRefFile(Fref, ipass)
! Read the ref file. (only MyID==0 called this routine)
! The first pass (ipass) get dimension and calculate statistics.
! The second pass store data

use CRYOFFmod
use extrafuncs

implicit none
integer, intent(in) :: ipass
character*80, intent(in) :: Fref
integer :: Ierr, i, j, Itmp
integer :: NwgtFatm, NwgtNetF, NwgtTorq, NatmSnp
integer :: Iatm, Iref
integer :: NrefSnp, Natm !ipass 2 only 
character*256 txt, tag, tag2, TagIatm
real *8 ::  SsqFatm,  MaxFatm, SsqNetF, MaxNetF, SsqTorq,  MaxTorq
real *8 :: wgtsca, wgtsol
real *8, external :: getValue, getBetaKmax, getErec

real*8 :: Rtmp,  XYZi(3), Fi(3), Wgt
        
        if (ipass == 1) then
		print*, '>> Reading Geometery and Reference Force'

		open(unit=Lref, file=trim(adjustl(Fref))//'.ref', status='old', IOstat=Ierr)
		call ErrStop('File  '//trim(adjustl(Fref))//'.ref  can not be opened!', Ierr)

		MaxAtm  = 0 ! Max Natm for all frm
		NsnpTot = 0; IsnpTot=0
		NatmTot = 0
		NrefTot = 0
		NatmRefTot = 0
		NwgtFatm=0;    NwgtNetF=0;    NwgtTorq=0
		SsqFatm =0.d0; SsqNetF =0.d0; SsqTorq =0.d0
		WssqFatm=0.d0; WssqNetF=0.d0; WssqTorq=0.d0
		MaxFatm =0.d0; MaxNetF =0.d0; MaxTorq =0.d0
        else !ipass == 2
        rewind(Lref)
      ! Note recycling variables.       
		NsnpTot = 0; IsnpTot=0
		NrefTot = 0  !Count Number of Reference Molecules
		NatmTot = 0
		NatmRefTot = 0

		IatmRefTot(:)=0
		IatmAtmTot(:)=0
		NatmSnpTot(:)=0
		NmolRefTot(:)=0
		ImolRefTot(:)=0
		IbgnAtmRef(:)=0
		NatmRefSnp(:)=0
		IatmRefSnp(:)=0
        Mole(:)%NmolTot = 0

		RminPair=1.d99
		RmaxPair=0.d0
        
        end if 
                    
		do while(IsnpTot<LstSnp)
			read(Lref, *, IOstat=Ierr) NatmSnp
			if(Ierr/=0) exit

			IsnpTot = IsnpTot+1
			if(IsnpTot<IniSnp) then  
				do i=0, NatmSnp
					read(Lref, *)
				end do
            else
            NsnpTot = NsnpTot+1

            if (ipass==2) then 
			NatmSnpTot(NsnpTot)=NatmSnp
			IbgnAtmSnp(NsnpTot)=NatmTot+1
            endif 
            
!           reading the comment line in the ref file.
            read(Lref, '(A)', IOstat=Ierr) txt
            if (ipass==2) then
!           scaling weight for each frame
			tag=txt; call upcase(tag);      i=index(trim(tag), 'WGTSC')
            wgtsca=1.0
            Rtmp=getValue(tag,'WGTSC'); if (Rtmp>0.d0) wgtsca=Rtmp
            
			if(YesPBC) then
                ! read in box dimension when doing PBC.
				i=index(trim(txt), 'box', back=.true.)
                if (i==0) then
                    print *, "frame number ", NsnpTot
                    call ErrStop("box dimension unspecified", 1)
                endif

				txt=txt(i+3:); Itmp=len_trim(txt)
				do i=1, Itmp ! remove non-digital chars
					j = ICHAR(txt(i:i))
					if( .NOT.(j==ICHAR('-') .or.  j==ICHAR('.') &
					&   .or. (ICHAR('0')<=j .and. j<=ICHAR('9'))) ) txt(i:i)=' '
				end do
				read(txt, *,IOstat=Ierr) Avec, Bvec, Cvec
                call ErrStop('Error reading box vectors', Ierr)
                
				Abox=norm2(Avec)    !Length of the a vector
				Bbox=norm2(Bvec)
				Cbox=norm2(Cvec)

				YesCub=.false.
				if(abs(dot_product(Avec,Bvec)) &
				& +abs(dot_product(Avec,Cvec)) &
				& +abs(dot_product(Bvec,Cvec))<Reps) YesCub=.true.

                call cross(Avec, Bvec, Uvec); Vol=dot_product(Uvec, Cvec)
                
                if(YesCub) then
					Uvec=0.d0; Uvec(1)=TwoPi/Abox
					Vvec=0.d0; Vvec(2)=TwoPi/Bbox
					Wvec=0.d0; Wvec(3)=TwoPi/Cbox
				else
					Rtmp=TwoPi/Vol
					call cross(Bvec, Cvec, Uvec); Uvec = Uvec*Rtmp
					call cross(Cvec, Avec, Vvec); Vvec = Vvec*Rtmp
					call cross(Avec, Bvec, Wvec); Wvec = Wvec*Rtmp
				end if
! 				print*, Abox, Bbox, Cbox, Vol, Beta, Etol, Rcou

				Rtmp = getBetaKmax(Etol*Vol/TwoPi, beta, getErec)

                ! This is computed and stored for all the frames. 
                ! I note these are recomputed "wrongly" everytime setPBC is called.
				VolTot( NsnpTot)=Vol
				AboxTot(NsnpTot)=Abox
				BboxTot(NsnpTot)=Bbox
				CboxTot(NsnpTot)=Cbox
				AvecTot(NsnpTot,:)=Avec; UvecTot(NsnpTot,:)=Uvec
				BvecTot(NsnpTot,:)=Bvec; VvecTot(NsnpTot,:)=Vvec
				CvecTot(NsnpTot,:)=Cvec; WvecTot(NsnpTot,:)=Wvec
				ImaxTot(NsnpTot)=int(Rtmp/norm2(Uvec))+1
				JmaxTot(NsnpTot)=int(Rtmp/norm2(Vvec))+1
				KmaxTot(NsnpTot)=int(Rtmp/norm2(Wvec))+1
            end if !YesPBC
            endif ! ipass == 2
            
            Iatm=0
            NrefSnp=0
			NrefAtmSnp=0
            
			do while(Iatm<NatmSnp)
				read(Lref, *) txt, XYZi, Fi, Wgt, tag
				call upcase(tag)
				do i=1, NmolTyp
					if(index(tag, trim(Mole(i)%Name))/=0) exit
				end do
				if(i>NmolTyp) call ErrStop('Mol. Type '//trim(adjustl(tag))//'  NOT  found in .ff File!', 1)

                
				Iref=0
				backspace(Lref)
                Natm=Mole(i)%Natm
				do j=1, Natm
					Iatm = Iatm+1
                    NatmTot=NatmTot+1
                    if (ipass==2) 		IatmAtmTot(NatmTot)=j
                    
					read(Lref, *) txt, XYZi, Fi, Wgt, tag2
					call upcase(txt)
                    call upcase(tag2)
                    if (tag/=tag2) call ErrStop ('Insufficient number of atoms read for molecule'//trim(adjustl(tag)), 1) 

                    if (ipass==2) then
                    SatmTot(NatmTot) = txt
                    TagIatm=txt !for KCV 
                    ! Tot here means the full set of equations not just this frame. 
	                XatmTot(NatmTot)=XYZi(1); YatmTot(NatmTot)=XYZi(2); ZatmTot(NatmTot)=XYZi(3)
                    NameTot(NatmTot) = tag

                    ! The idea is to find the record number of the string Mole(i)%COUatm(j) in the array COUatmTyp(1:NcouTyp)
					IcouTot(NatmTot)=findatm(NcouTyp, Mole(i)%COUatm(j), COUatmTyp(1:NcouTyp))
					IvdwTot(NatmTot)=findatm(NvdwTyp, Mole(i)%VDWatm(j), VDWatmTyp(1:NvdwTyp))
                    endif ! ipass==2
                    
                    ! set Wgt to zero when seletive fitting is used.  
                    ! I think it is clearer to simply construct different .ref file if only want of fit NETF or TORQ
                    
					if( YesInt .and. &
					&  (     txt(1:4)=='NETF' .and. .not.YesNetF(i) &
					&	.or. txt(1:4)=='TORQ'.and. .not.YesTorq(i) ) ) Wgt=0.d0

                    if(Wgt>WgtTol) then
                        ! At this stage the Wgt is the solvation factor. 
                        ! All statistics only counting weighted atoms. 
                        Iref=1
                        NatmRefTot = NatmRefTot+1
                        if (ipass==1) then 

						Rtmp=dot_product(Fi,Fi)
						if(txt(1:4)=='NETF') then
							NwgtNetF=NwgtNetF+1
							SsqNetF =SsqNetF +Rtmp
							WssqNetF=WssqNetF+Wgt*Rtmp
							Rtmp=sqrt(Rtmp)
							MaxNetF=max(Rtmp, MaxNetF)
						else if(txt(1:4)=='TORQ') then
							NwgtTorq=NwgtTorq+1
							SsqTorq =SsqTorq +Rtmp
							WssqTorq=WssqTorq+Wgt*Rtmp
							Rtmp=sqrt(Rtmp)
							MaxTorq=max(Rtmp, MaxTorq)
						else
							if(count(Mole(i)%IatmVir(:)==j)==0) then
                                ! Not counting virtual sites. 
								NwgtFatm=NwgtFatm+1
								SsqFatm =SsqFatm +Rtmp
								WssqFatm=WssqFatm+Wgt*Rtmp
								Rtmp=sqrt(Rtmp)
								MaxFatm=max(Rtmp, MaxFatm)
                            end if
                            endif
                        else !ipass == 2 
						NrefAtmSnp = NrefAtmSnp+1
						IatmRefTot(NatmTot)=NatmRefTot
						Itmp=0
                        Wgtsol=Wgt ! Save the solvation factor.
						if(TagIatm(1:4)=='NETF') then
							call setWgt(Iwgt, Wgt, Fi, WRmsNetF,.false.)
						else if(TagIatm(1:4)=='TORQ') then
							call setWgt(Iwgt, Wgt, Fi, WRmsTorq, .true.)
						else
							call setWgt(Iwgt, Wgt, Fi, WRmsFatm,.false.)
						end if
						Itmp = 3*NatmRefTot
						FatmTot(Itmp-2:Itmp) = Fi
						FwgtTot(Itmp-2:Itmp) = Wgt*wgtsca*Wgtsol
						IkcvTot(Itmp-2:Itmp) = IkcvFrm(NsnpTot)
						write(txt, '(3(I0, A))') IkcvFrm(NsnpTot), '_', NsnpTot, '.'//trim(tag)//'.', j, '-'//trim(TagIatm)
                        if (len_trim(txt) > MaxLab) call ErrStop("Label Length Too short",1)
						LabAtmTot(Itmp-2:Itmp) = trim(adjustl(txt))
                        ! This seems to be a uniqe atom name and encode more information than a simple atom name. 
                        end if !ipass
                        endif !Wgt
                    end do
				if(Iref>0) then
                    NrefTot = NrefTot+1
                    if (ipass==2) then
                    NrefSnp = NrefSnp+1
					ImolRefTot(NrefTot) = i
					IbgnAtmRef(NrefTot) = Iatm-Natm
                    Mole(i)%NmolTot = Mole(i)%NmolTot + 1
                    endif
                endif
            end do

            if (ipass==1) then 
                if(NatmSnp>MaxAtm) MaxAtm = NatmSnp
            else 
            NmolRefTot(NsnpTot)=NrefSnp
			IbgnMolRef(NsnpTot)=NrefTot-NrefSnp+1
			NatmRefSnp(NsnpTot)=NrefAtmSnp
			IatmRefSnp(NsnpTot)=NatmRefTot-NrefAtmSnp+1
                       
            !           Only print the first and last 6 frames when reading. 
            ! Should print box vectors for PBC runs. 
			if(NsnpTot<6 .or. MaxSnp-NsnpTot<6) then
				write(txt, '(A, I0)') '/', MaxSnp
				write(*, '(A, I6, 2A, I9, 3(A, I6))') '   >> Frame:', NsnpTot, trim(txt), &
				& '      Atom:',      IbgnAtmSnp(NsnpTot), '  This Frame', NatmSnp, &
				& '      Fit Mol.:', IbgnMolRef(NsnpTot), '  This Frame', NrefSnp
            end if
            endif !ipass
          
		end if !IsnpTot
        end do
        
        if (ipass==1) then
        ! Note these are weighted per force component. 
		if(NwgtFatm>0) 	RmsFatm=sqrt(SsqFatm/dble(NwgtFatm))
		if(NwgtNetF>0)  RmsNetF=sqrt(SsqNetF/dble(NwgtNetF))
		if(NwgtTorq>0) 	RmsTorq=sqrt(SsqTorq/dble(NwgtTorq))
        if(NwgtFatm>0) 	WRmsFatm=sqrt(WSsqFatm/dble(NwgtFatm))
		if(NwgtNetF>0)  WRmsNetF=sqrt(WSsqNetF/dble(NwgtNetF))
		if(NwgtTorq>0) 	WRmsTorq=sqrt(WSsqTorq/dble(NwgtTorq))
        if(NwgtTorq>0 .and. NwgtNetF>0) wgtfac=wgtfac*sqrt(dble(NwgtNetF)/dble(NwgtTorq))
        endif !ipass == 1
        
        
        if (ipass==1)  then 
        print*, '   Total Frames Read   ', NsnpTot
		print*, '   Total Atoms    ', NatmTot, '  Max.  Atom/Frame', MaxAtm
		print*, '   Total Ref. Mol.', NrefTot, '  Total Ref. Atom ', NatmRefTot
		print*, '                |Atomic Force|       |Net Force|     |Torque|'
		write(*, '(A, 5X, I16, 2X, I16, 2X, I16)') '  Number', NwgtFatm, NwgtNetF, NwgtTorq
        write(*, '(A, 5X, F16.6, 2X, F16.6, 2X, F16.6)') '    RMS   ', RmsFatm,  RmsNetF,  RmsTorq
		write(*, '(A, 5X, F16.6, 2X, F16.6, 2X, F16.6)') '    Max.  ', MaxFatm,  MaxNetF,  MaxTorq
		write(*, '(A, 5X, F16.6, 2X, F16.6, 2X, F16.6,/)')'    WRMS  ', WRmsFatm, WRmsNetF, WRmsTorq
        
		write(Lout, '(A)') 'Reference Geometry/Force File: '//trim(Fref)//'.ref'
		write(Lout, *) ' Total Frames Read    ', NsnpTot
		write(Lout, *) ' Total Atoms    ', NatmTot, '  Max.  Atom/Frame', MaxAtm
		write(Lout, *) ' Total Ref. Mol.', NrefTot, '  Total Ref. Atom ', NatmRefTot, LF
		write(Lout, *) ' item         |Atomic Force|       |Net Force|     |Torque|'
		write(Lout, '(A, 5X, I16, 2X, I16, 2X, I16)')    '  Number', NwgtFatm, NwgtNetF, NwgtTorq
        write(Lout, '(A, 5X, F16.6, 2X, F16.6, 2X, F16.6)') '    RMS   ', RmsFatm,  RmsNetF,  RmsTorq
		write(Lout, '(A, 5X, F16.6, 2X, F16.6, 2X, F16.6)') '    Max.  ', MaxFatm,  MaxNetF,  MaxTorq
		write(Lout, '(A, 5X, F16.6, 2X, F16.6, 2X, F16.6,/)')'    WRMS  ', WRmsFatm, WRmsNetF, WRmsTorq
        else 
        write(*, *) ' Finished Reading    ', NsnpTot,' frames'
        write(Lout, *) ' Finished Reading    ', NsnpTot,' frames'
        close(Lref)
        endif

End subroutine ReadRefFile 
