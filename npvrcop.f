      implicit none
      character*100 arg,val
      integer iarg,narg
      
      integer niter,plev,seed,ny,seedrcind,iwrcr,iwrzc
      character*100 idcd,seedrc
      real*8 xA, xB,dr2,r
      
      integer natom,nsets
      integer natom0,nsets0
      logical ft

      integer omp_get_num_threads,omp_get_thread_num
      
      idcd='in.dcd'
      seedrc='seedrc'
      seedrcind=1
      seed=1234
      xA=0
      xB=1
      niter=1000000
      dr2=-1
      plev=0
      ny=4
      natom=0
      nsets=0
      iwrcr=2000
      iwrzc=200
      
      narg=command_argument_count()
      do iarg=1,narg,2
        call get_command_argument(iarg,arg)
        arg=trim(arg)
        call get_command_argument(iarg+1,val)
        val=trim(val)
        if (arg=='--idcd')then
          idcd=val
        elseif (arg=='--seedrc')then
          seedrc=val
        elseif (arg=='--seedrcind')then
          read(val,*)seedrcind
        elseif (arg=='--xA')then
          read(val,*)xA
        elseif (arg=='--xB')then
          read(val,*)xB
        elseif (arg=='--niter')then
          read(val,*)niter
        elseif (arg=='--dr2')then
          read(val,*)dr2
        elseif (arg=='--printlev')then
          read(val,*)plev
        elseif (arg=='--natom')then
          read(val,*)natom
        elseif (arg=='--nsets')then
          read(val,*)nsets
        elseif (arg=='--seed')then
          read(val,*)seed
        elseif (arg=='--ny')then
          read(val,*)ny
        elseif (arg=='--iwrcr')then
          read(val,*)iwrcr
        elseif (arg=='--iwrzc')then
          read(val,*)iwrzc
        else
          write(*,*)'unknown command: ',trim(arg)
          stop 'exiting...' 
        endif
      enddo
      if (plev>=0) write(*,*)
     $ 'NonParametric Variational Reaction Coordinate Optimization'
      if (plev>=1)then
        write(*,*)'input parameters:'
        write(*,*)'--idcd = ', trim(idcd)
        write(*,*)'--seedrc = ', trim(seedrc)
        write(*,*)'--seedrcind = ', seedrcind
        write(*,*)'--xA =',xA, '--xB =',xB
        write(*,*)'--niter =',niter
        write(*,*)'--dr2 =',dr2
        write(*,*)'--printlev =',plev
        write(*,*)'--natom =',natom
        write(*,*)'--nsets =',nsets
        write(*,*)'--seed =',seed
        write(*,*)'--ny =',ny
        write(*,*)'--iwrcr =',iwrcr
        write(*,*)'--iwrzc =',iwrzc
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
      r=rand(-abs(seed))
!$OMP PARALLEL
!$    if (plev>=0 .and. omp_get_thread_num()==0)
!$   $  write(*,*)'number of threads', omp_get_num_threads()
!$OMP END PARALLEL
      call driver
     $ (idcd,natom,natom0,nsets,seedrc,seedrcind,xA,xB,niter,dr2,ny
     $,plev,iwrcr,iwrzc)
      end
      
      subroutine driver
     $ (idcd,natom,natom0,nsets,seedrc,seedrcind,xA,xB,niter,dx2m,ny
     $,plev,iwrcr,iwrzc)
      implicit none
      character*(*) idcd,seedrc
      integer natom,natom0,nsets,niter,plev,ny,seedrcind,iwrzc,iwrcr
      real*8 xA, xB,dx2m
      
      integer rcfix(nsets),rcind(nsets),idt
      real*8 rc(nsets),rij(nsets),NAB

      real*4 x(nsets,natom),y(nsets,natom),z(nsets,natom)
      integer i21(natom)
      
      real*8 compdx2ring,dx2
      
      integer iset,iter
      
      call readrc(rc,rcind,rcfix,seedrc,nsets,seedrcind,xA,xB,NAB)
      do iset=1,nsets
        if (rcfix(iset)==0)rc(iset)=2*rand()-1 ! random initialization
      enddo
      if (dx2m<0)dx2m=NAB
      idt=1
      dx2=compdx2ring(rc,rcind,nsets,idt)
      if (plev>=0)write(*,*)'  iter   dr^2/2      NAB'
      do iter=0,niter
        if (iter==0 .or. (natom/=natom0 .and. mod(iter,100000)==0))
     $    call readxyz(idcd,x,y,z,natom0,natom,nsets,i21) 
         
        call comprij(x,y,z,natom,nsets,rij)
        call optimdx2ring(ny,ny,rc,rij,rcind,rcfix,nsets,idt,dx2)
        if (mod(iter,5)==0) 
     $    call intdx2(rc,rcfix,rcind,nsets,idt,dx2)
        if (mod(iter,10)==0)
     $    call optimdx2ring(10,0,rc,rij,rcind,rcfix,nsets,idt,dx2)
     
        if (mod(iter,10)==0 .and. plev>=0)
     $    write(*,'(i7,f9.2,f9.2)')iter,dx2/8,NAB 
        if (dx2/8<dx2m)exit
        if (iwrcr>0 .and. iter>0 .and. mod(iter,iwrcr)==0)
     $    call writecoor('rc',rc,nsets)
        if (iwrzc>0 .and. iter>0 .and. mod(iter,iwrzc)==0)
     $    call writezc('rc.zc',rc,rcind,nsets,16)
      enddo
      call writecoor('rc',rc,nsets)
      end
      
      
      
