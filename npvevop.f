! NonParametric Variational Eigenvector OPtimization 
      implicit none
      character*100 arg,val
      integer iarg,narg
      
      integer niter,plev,seed,ny,iwrev,iwrzc,nev,idt1
      character*100 idcd
      real*8 r
      
      integer natom,nsets
      integer natom0,nsets0
      logical ft

      integer omp_get_num_threads,omp_get_thread_num
      
      idcd='in.dcd'
      seed=1234
      niter=1000000
      plev=0
      nev=1
      ny=3
      natom=0
      nsets=0
      iwrev=2000
      iwrzc=200
      idt1=4096
      
      narg=command_argument_count()
      do iarg=1,narg,2
        call get_command_argument(iarg,arg)
        arg=trim(arg)
        call get_command_argument(iarg+1,val)
        val=trim(val)
        if (arg=='--idcd')then
          idcd=val
        elseif (arg=='--nev')then
          read(val,*)nev
        elseif (arg=='--niter')then
          read(val,*)niter
        elseif (arg=='--plev')then
          read(val,*)plev
        elseif (arg=='--natom')then
          read(val,*)natom
        elseif (arg=='--nsets')then
          read(val,*)nsets
        elseif (arg=='--seed')then
          read(val,*)seed
        elseif (arg=='--ny')then
          read(val,*)ny
        elseif (arg=='--iwrev')then
          read(val,*)iwrev
        elseif (arg=='--iwrzc')then
          read(val,*)iwrzc
        elseif (arg=='--idt1')then
          read(val,*)idt1
        else
          write(*,*)'unknown command: ',trim(arg)
          stop 'exiting...' 
        endif
      enddo
      if (plev>=0) write(*,*)
     $ 'NonParametric Variational EigenVector Optimization'
      if (plev>=1)then
        write(*,*)'input parameters:'
        write(*,*)'--idcd = ', trim(idcd)
        write(*,*)'--nev =',nev
        write(*,*)'--niter =',niter
        write(*,*)'--plev =',plev
        write(*,*)'--natom =',natom
        write(*,*)'--nsets =',nsets
        write(*,*)'--seed =',seed
        write(*,*)'--ny =',ny
        write(*,*)'--iwrev =',iwrev
        write(*,*)'--iwrzc =',iwrzc
        write(*,*)'--idt1 =',idt1
      endif
      call opendcd(idcd,1,natom0,nsets0,ft)
      close(1)
      if (plev>=0)then
        write(*,*)'number of atoms in the trajectory',natom0
        write(*,*)'number of snapshots in the trajectory',nsets0
      endif
      if (nsets==0)nsets=nsets0
      if (nsets>nsets0)stop'nsets > trajectory length'
      if (natom==0)natom=natom0
      if (natom>natom0)stop'natom > number of atoms in the trajectory'
      if (nev>1)stop'nev>1 is not implemented yet'
      write(*,*)'RNG(seed)=',rand(-abs(seed))
!$OMP PARALLEL
!$    if (plev>=0 .and. omp_get_thread_num()==0)
!$   $  write(*,*)'number of threads', omp_get_num_threads()
!$OMP END PARALLEL
      call driver(idcd,natom,natom0,nsets,nev,niter,ny,plev,iwrev,iwrzc
     $  ,idt1)
      end
      
      subroutine driver
     $ (idcd,natom,natom0,nsets,nev,niter,ny,plev,iwrev,iwrzc,idt1)
      implicit none
      character*(*) idcd
      integer natom,natom0,nsets,nev,niter,plev,ny,iwrev,iwrzc,idt1

      real*8 ev(nsets,nev)
      real*8 rij(nsets),eval,compeval,eval0
      integer rcfix(nsets)

      integer iset,idt,iter
      
      real*4 x(nsets,natom),y(nsets,natom),z(nsets,natom)
      
      integer i21(natom)
      real*8 evt
      real*8 disp,cd
      
      ! initialization
      call readxyz(idcd,x,y,z,natom0,natom,nsets,i21) 
      idt=1
      call comprij(x,y,z,natom,nsets,ev(1,1))
      eval=10000
      call anev(ev(1,1),nsets,idt,eval)
      rcfix=0
      cd=disp(ev,nsets,16)
      write(*,*)eval,cd
      
      if (plev>=0)write(*,*)'  iter   eval      eval0'
      do iter=0,niter
        if (iter>0 .and. (natom/=natom0 .and. mod(iter,100000)==0))
     $    call readxyz(idcd,x,y,z,natom0,natom,nsets,i21) 
         
10      call comprij(x,y,z,natom,nsets,rij)
        evt=100
        call anev(rij,nsets,idt,evt)
c        write(*,*)evt
        if (evt<1e-4)then
          write(*,*)'eval(rij)=',evt,'<1e-4, recomputing rij'
          goto 10
        endif
        if (evt>1e-1)goto 10
        call optimdx2ev2(ev,nev,nsets,ny,rij,ny,rcfix,idt,eval,cd)
        if (mod(iter,5)==0)
     $    call anev(ev,nsets,idt,eval)
        if (mod(iter,10)==0)
     $    call optimdx2ev(ev,nev,nsets,10,rij,0,rcfix,idt,eval)
        
        eval0=compeval(ev(1,1),nsets,idt1)
        if (mod(iter,5)==0)
     $    cd=disp(ev,nsets,16)
        if (eval<eval0 .and. eval0<1)exit
        if (mod(iter,10)==0 .and. plev>=0) write(*,*)iter,eval,eval0,cd
        if (iwrev>0 .and. iter>0 .and. mod(iter,iwrev)==0)
     $    call writeev('ev',ev,nsets,nev)
        if (iwrzc>0 .and. iter>0 .and. mod(iter,iwrzc)==0)
     $    call writeevcriterion('ev.zc',ev,nsets,16,idt1)
      enddo
      call anev(ev,nsets,idt,eval)
      write(*,*)iter,eval,eval0
      call writeev('ev',ev,nsets,nev)
      call writeevcriterion('ev.zc',ev,nsets,16,idt1)
      end
      
      
      
