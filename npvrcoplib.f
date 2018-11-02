! NonParametric Variational Reaction Coordinate OPtimization 
! library of subroutines

      subroutine readxyz(dcd,x,y,z,natom0,natom,nsets,i21)
      implicit none 
      character*(*) dcd
      integer natom0,natom,i21(natom),nsets
      real*4 x(nsets,natom),y(nsets,natom),z(nsets,natom)
      real*4 temp(natom0),t(12)
        
      integer i,iset
      logical ft
      call opendcd(dcd,1,natom0,i,ft)
      if (natom0==natom)then
        do i=1,natom
          i21(i)=i
        enddo
      else
        call select_atoms(natom0,natom,i21)
      endif

      do iset=1,nsets
        if(ft) read(1)t
        read(1)temp
        do i=1,natom
          x(iset,i)=temp(i21(i))
        enddo  
        read(1)temp
        do i=1,natom
          y(iset,i)=temp(i21(i))
        enddo  
        read(1)temp
        do i=1,natom
          z(iset,i)=temp(i21(i))
        enddo  
      enddo
      close(1)
      end

      subroutine select_atoms(natom0,natom,i21)
      implicit none
      integer natom0,natom,i21(natom)
      
      integer i,l,taken(natom0)
      do i=1,natom0
          taken(i)=0
      enddo
      do l=1,natom
10      i=rand()*natom0+1
        if (taken(i)/=0)goto 10
        taken(i)=1
        i21(l)=i
      enddo  
      end

      
      subroutine optimdx2ring(nr,ny,rc,y,rcind,rcfix,nsets,idt,dx2)
      implicit none
      integer nr,ny,nsets,rcind(nsets),rcfix(nsets),idt
      real*8 rc(nsets),y(nsets),dx2

!!!!  minimum of [x(t)-x(t+dt)]^2 +[x(t)-x_A]^2+[x(t)-x_B]^2 where x_A=-1 and x_B=1
!!!!  for all i: \sum_j al_j[<dx_idx_j>_t+<x_ix_j>_AB] = x_AB <x_i>_AB 
      
      real*8 dxdx((ny+1)*(ny+1)+nr-ny,(ny+1)*(ny+1)+nr-ny)
      real*8 al((ny+1)*(ny+1)+nr-ny),alal
      real*8 rhs((ny+1)*(ny+1)+nr-ny)
      
     
      integer i,j,i1,nij,iset,k,isize
      real*8 r
      
      integer info,ipiv((ny+1)*(ny+1)+nr-ny)
      real*8 sval((ny+1)*(ny+1)+nr-ny),work(((ny+1)*(ny+1)+nr-ny)*1000)
      integer rank,lwork,iwork(((ny+1)*(ny+1)+nr-ny)*1000)
      
      integer i12((ny+1)*(ny+1)+nr-ny)
      integer i21((ny+1)*(ny+1)+nr-ny),nij2
      
      real*8 compdx2ring,dx2n
      real*8 rcn(nsets)

      real*8 rcp1,yp1,dfij((ny+1)*(ny+1)+nr-ny)
     
      isize=(ny+1)*(ny+1)+nr-ny
      lwork=isize*1000
      
      do i=1,isize
        al(i)=0
        do j=1,i
          dxdx(i,j)=0
        enddo
      enddo
      alal=0
      
      call compxxring(rc,y,nr,ny,rcfix,rcind,nsets,idt,alal,al,dxdx
     $   ,isize)
      do j=1,isize
        do k=1,j
          dxdx(j,k)=dxdx(j,k)/idt
        enddo
        al(j)=al(j)/idt
      enddo
      alal=alal/idt
      nij2=0
      do i=1,isize
        i21(i)=0
      enddo  
      do i=1,isize  ! select only non zero
        if (dxdx(i,i)<1e-5)cycle
        nij2=nij2+1
        i12(nij2)=i
        i21(i)=nij2
      enddo
      do i=1,nij2
        do j=1,i
          dxdx(i,j)=dxdx(i12(i),i12(j))
          dxdx(j,i)=dxdx(i,j)
        enddo
        al(i)=al(i12(i))
        rhs(i)=-al(i)/2
      enddo    


      if (.true.)then !solving SLE
        call dgetrf(nij2,nij2,dxdx,isize,ipiv,info)
        if(info/=0)then
          write(*,*)'dgetrf info=',info
          return
          stop'info/=0'
        endif
        call dgetrs('N',nij2,1,dxdx,isize,ipiv,rhs,isize,info)
        if(info/=0)then
          write(*,*)'dgetrs info=',info
          stop'info/=0'
        endif
      else           ! solving least square problem
        call DGELSD(nij2,nij2,1,dxdx,isize,rhs,isize,sval,-1d-14, rank,
     $                         work,lwork,iwork,info)   
        if(info/=0)then
          write(*,*)'DGELSD info=',info
          stop'info/=0'
        endif
      endif
      

!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,nij2,ny,nr,y,rc,rcfix,rcn,rhs,i12) 
!$OMP&   PRIVATE(nij,rcp1,yp1,dfij,i1,r)

!$OMP DO
      do iset=1,nsets
          rcn(iset)=rc(iset)
          if (rcfix(iset)/=0) cycle
          nij=0
          rcp1=1
          do i=0,ny
            yp1=rcp1
            do j=0,ny
              nij=nij+1
c              dfij(nij)=rc(iset)**i*y(iset)**j
              dfij(nij)=yp1
              yp1=yp1*y(iset)
            enddo
            rcp1=rcp1*rc(iset)
          enddo  
          do i=ny+1,nr
            nij=nij+1
c            dfij(nij)=rc(iset)**i
            dfij(nij)=rcp1
            rcp1=rcp1*rc(iset)
          enddo    
          r=0
          do i=1,nij2
            i1=i12(i)
            r=r+rhs(i)*dfij(i1)
          enddo    
          rcn(iset)=r
      enddo
!$OMP END PARALLEL
      dx2n=compdx2ring(rcn,rcind,nsets,idt)
      if (dx2n>dx2) then
c        write(*,*)'optim: rejected dx2n,dx2=',dx2n,dx2
        return
      endif
      rc=rcn
      dx2=dx2n 
      end    

      subroutine compxxring(rc,y,nr,ny,rcfix,rcind,nsets,dt,alal,al
     $ ,dxdx,isize2)
      implicit none
      integer nsets,dt,isize2,nr,ny
      real*8 rc(nsets),y(nsets),alal,al(isize2),dxdx(isize2,isize2)
      integer rcfix(nsets),rcind(nsets)
      
      integer iset,nij,i,j,ix
      real*8 rcp1,rcp2,yp1,yp2,dfij(isize2)
      real*8 dind
      integer bi,b1,b2
      
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dt,rcind,rcfix,ny,nr,y,rc) 
!$OMP&   PRIVATE(bi,b1,b2,dind,ix,nij,rcp2,rcp1,yp2,yp1,dfij)
!$OMP&   reduction(+:alal,al,dxdx) 

!$OMP DO
        do iset=1,nsets-dt
  !       f(t)=2rcind(t)+(-1)**rcind(t)*rc(t)=k+b*rc
          dind=rcind(iset+dt)-rcind(iset)
          b1=1
          if (mod(rcind(iset),2)/=0)b1=-1 ! (-1)**rcind
          b2=1
          if (mod(rcind(iset+dt),2)/=0)b2=-1 ! (-1)**rcind
          if (rcfix(iset)==0 .and. rcfix(iset+dt)==0)then
            nij=0
            rcp2=b2
            rcp1=b1
            do i=0,ny
              yp2=rcp2
              yp1=rcp1
              do j=0,ny
                nij=nij+1
c                dfij(nij)=b2*rc(iset+dt)**i*y(iset+dt)**j
c     $                     -b1*rc(iset)**i*y(iset)**j
                 dfij(nij)=yp2-yp1
                 yp2=yp2*y(iset+dt)
                 yp1=yp1*y(iset)
              enddo
              rcp2=rcp2*rc(iset+dt)
              rcp1=rcp1*rc(iset)
            enddo
            do i=ny+1,nr
              nij=nij+1
c              dfij(nij)=b2*rc(iset+dt)**i-b1*rc(iset)**i
              dfij(nij)=rcp2-rcp1
              rcp2=rcp2*rc(iset+dt)
              rcp1=rcp1*rc(iset)
            enddo   
            alal=alal+4*dind*dind
            do j=1,nij
              al(j)=al(j)+4*dind*dfij(j)
            enddo
            do i=1,nij
              do j=1,i
                dxdx(i,j)=dxdx(i,j)+dfij(i)*dfij(j)
              enddo
            enddo
          else
            if (rcfix(iset)/=0 .and. rcfix(iset+dt)/=0)then
              alal=alal+(2*dind+b2*rc(iset+dt)-b1*rc(iset))**2
              cycle
            endif
            if (rcfix(iset)/=0)then 
              ix=iset+dt ! point in the middle
              dind=2*dind-b1*rc(iset)     ! distance offset
              bi=b2
            endif
            if (rcfix(iset+dt)/=0)then
              ix=iset ! point in the middle
              dind=-2*dind-b2*rc(iset+dt)    ! distance offset
              bi=b1
            endif
            nij=0
            rcp1=bi
            do i=0,ny
              yp1=rcp1
              do j=0,ny
                nij=nij+1
c                dfij(nij)=bi*rc(ix)**i*y(ix)**j
                dfij(nij)=yp1
                yp1=yp1*y(ix)
              enddo
              rcp1=rcp1*rc(ix)
            enddo  
            do i=ny+1,nr
              nij=nij+1
c              dfij(nij)=bi*rc(ix)**i
              dfij(nij)=rcp1
              rcp1=rcp1*rc(ix)
            enddo    
            alal=alal+dind*dind
            do j=1,nij
              al(j)=al(j)+2*dind*dfij(j)
            enddo
            do i=1,nij
              do j=1,i
                dxdx(i,j)=dxdx(i,j)+dfij(i)*dfij(j)
              enddo
            enddo
          endif  
        enddo
!$OMP END PARALLEL
      end

      subroutine comprij(x,y,z,natom,nsets,rij)
      implicit none
      integer natom,nsets
      real*4 x(nsets,natom),y(nsets,natom),z(nsets,natom)
      real*8 rij(nsets)
      
      integer i1,j1,iset
      real*8 r
       
      i1=rand()*natom+1
10    j1=rand()*natom+1
      if (i1==j1)goto 10
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,x,y,z,rij,i1,j1) 
!$OMP&   PRIVATE(r)

!$OMP DO
      do iset=1,nsets
        r=(x(iset,i1)-x(iset,j1))**2
     $   +(y(iset,i1)-y(iset,j1))**2
     $   +(z(iset,i1)-z(iset,j1))**2
        r=sqrt(r)
        rij(iset)=r
      enddo
!$OMP END PARALLEL
      end
      

      subroutine intdx2(rc,rcfix,rcind,nsets,idt,dx2)
      implicit none
      integer nsets,rcind(nsets),idt,rcfix(nsets)
      real*8 rc(nsets),dx2
      
      real*8 cii(10000),q(10000),dqdx(10000)

      integer jdt,i,i1,i2,iset,j
      real*8 dx2n,rcn(nsets),s,dx,compdx2ring
      
      
      dx=2.2/10000
      call zc1ring(rc,rcind,nsets,cii,10000,idt)
      
      i1=(-1+1.1)/dx
      i2=(1+1.1)/dx
      s=0
      do i=i1,i2
        q(i)=s
        if (cii(i)>0)dqdx(i)=1/cii(i)
        if (cii(i+1)>0)dqdx(i)=1/cii(i+1)
        if (cii(i)>0 .and. cii(i+1)>0)dqdx(i)=(1/cii(i)+1/cii(i+1))/2
        s=s+dqdx(i)
      enddo
      do i=i1,i2
        q(i)=q(i)/s*2-1
        dqdx(i)=dqdx(i)*2/s
      enddo
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(rcfix,nsets,rc,dx,i1,i2,q,dqdx,rcn) 
!$OMP&   PRIVATE(j)

!$OMP DO
      do iset=1,nsets
        rcn(iset)=rc(iset)
        if (rcfix(iset)/=0)cycle
        j=(rc(iset)+1.1)/dx
        if (j<i1)j=i1
        if (j>i2)j=i2
        rcn(iset)=q(j)+dqdx(j)*(rc(iset)-j*dx-1.1)
      enddo  
!$OMP END PARALLEL
      dx2n=compdx2ring(rcn,rcind,nsets,idt)
      if (dx2n<dx2) then
        rc=rcn
        dx2=dx2n
      endif
      end

      subroutine opendcd(name,iunit,natom,nsets,ft)
      implicit none
      integer natom,nsets 
      character*4 hdr
      integer cntrl(20), ntitle,i
      character*80 title(20)
      character*(*)name
      integer iunit
      logical ft

      open(unit=iunit, file=name,form='unformatted',status='old',err=10)
      read(iunit)hdr,cntrl
      nsets=cntrl(1)
      ft=cntrl(11)>0
      read(iunit)ntitle,(title(i),i=1,ntitle)
      read(iunit)natom
      return
10    write(*,*)'can not open the file'
      write(*,*)name
      stop
      end 
      
      subroutine writecoor(name,rc,nsets)
      implicit none
      integer nsets
      real*8 rc(nsets)
      character*(*) name
      
      integer iset
      
      open(unit=1,file=name)
      do iset=1,nsets
        write(1,*)rc(iset)/2+0.5
      enddo
      close(1)
      end

      subroutine writeind(name,rcind,nsets)
      implicit none
      integer nsets
      integer rcind(nsets)
      character*(*) name
      
      open(unit=1,file=name,form='unformatted')
      write(1)rcind
      close(1)
      end
 
      real*8 function compdx2ring(rc,rcind,nsets,idt)
      implicit none
      INCLUDE "omp_lib.h"
      integer nsets,rcind(nsets),idt
      real*8 rc(nsets)
      
      integer iset
      real*8 r
      
      integer dind
      real*8 b1,b2
      r=0
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,idt,rcind,rc) 
!$OMP&   PRIVATE(b1,b2,dind)
!$OMP&   reduction(+:r) 

!$OMP DO
      do iset=1,nsets-idt
        b1=1
        if (mod(rcind(iset),2)/=0)b1=-1 ! (-1)**rcind
        b2=1
        if (mod(rcind(iset+idt),2)/=0)b2=-1 ! (-1)**rcind
        dind=rcind(iset+idt)-rcind(iset)
        r=r+(2*dind+b2*rc(iset+idt)-b1*rc(iset))**2
      enddo    
!$OMP END PARALLEL 
      compdx2ring=r/idt
      end


      subroutine readrc(rc0,rcind,rcfix,rc0file,nsets,ind,xA,xB,JAB)
      implicit none
      integer nsets,ind,rcind(nsets),rcfix(nsets)
      real*8 rc0(nsets),temp(ind),xA,xB,JAB
      character *(*)rc0file
      
      real*8 r
      integer i,j,curind,lastb
      
      open(unit=1,file=rc0file,status='old',err=10)
      
      curind=0
      lastb=0
      JAB=0
      do i=1,nsets
        read(1,*,end=20)(temp(j),j=1,ind)
        r=temp(ind)
        r=2*(r-xA)/(xB-xA)-1
        if (r>=1) r=1
        if (r<=-1)r=-1
        rc0(i)=r
        rcfix(i)=0
        if (r>=1) rcfix(i)=1
        if (r<=-1) rcfix(i)=-1
        if (rcfix(i)/=0)then
          if (lastb/=rcfix(i))then
            if (lastb/=0)JAB=JAB+1
            lastb=rcfix(i)
          endif
        endif      
        if (rcfix(i)/=0 .and. rand()<0.5 )then !on the boundary node, switch to another branch of the multivalued RC
            if (mod(curind,2)==0 .and. rcfix(i)==-1)then
              curind=curind-1 
            elseif(mod(curind,2)/=0 .and. rcfix(i)==-1)then
              curind=curind+1   
            elseif(mod(curind,2)==0 .and. rcfix(i)==1)then
              curind=curind+1   
            elseif(mod(curind,2)/=0 .and. rcfix(i)==1)then
              curind=curind-1   
            endif
        endif   
        rcind(i)=curind
      enddo
      JAB=JAB/2
      close(1)
      if (JAB<0.5)then
        write(*,*)'small NAB ', JAB
        write(*,*)'Possible reasons are the wrong positions of the'
        write(*,*)' boundaries xA and xB, or wrong format of'
        write(*,*)' the seedrc file'
        stop
      endif
      return 
10    write(*,*)'cannot open file ',rc0file
      stop
20    write(*,*)'end of file ',rc0file
      stop
      end  
      
      subroutine writezc(name,rc,rcind,nsets,ldt)
      implicit none
      integer nsets,ldt,rcind(nsets)
      real*8 rc(nsets)
      character*(*) name
      
      real*8 cii(1000,ldt),dx,x,r
      integer idt,iset,j
      
      real*8 zh(1000),azh2(1000)


      cii=0
      zh=0
      azh2=0

      dx=2.2/1000
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dx,rc) 
!$OMP&   PRIVATE(j,r)
!$OMP&   reduction(+:zh) 

!$OMP DO
      do iset=1,nsets
        r=rc(iset)
        if (r<-1)r=-1
        if (r>1)r=1
        j=(r+1.1)/dx
        if (j>1000)j=1000
        if (j<1)j=1
        zh(j)=zh(j)+1
      enddo
!$OMP END PARALLEL 
      
      call zh2(rc,nsets,azh2,1000)

      call zc1TPsldt(rc,nsets,cii,1000,ldt)
c      call zc1ringldt(rc,rcind,nsets,cii,1000,ldt)

      open(unit=11,file=name)
      dx=2.2/1000
      do j=1,1000
        x=(dx*j-0.1)/2
        if (x<0 .or. x>1) cycle
        write(11,'(f6.4, 20f9.3 )')x,-log(max(zh(j),1d0)/dx*2)
     $   ,-log(max(azh2(j),1d0)*2*2)
     $   ,(-log(max(cii(j,idt)/2.0,1.0)),idt=1,ldt)
      enddo    
      close(11)
      end  

      subroutine zc1ring(rc,rcind,msets,ci,ni,dt)
      implicit none
      integer msets,rcind(*),ni,dt
      real*8 rc(*)
      real*8 ci(ni)
      
      integer j,jl,iset,jp1,jm1,nloops

      real*8 s,dx,d,absd,b,bl,dind
      
      do j=1,ni
        ci(j)=0
      enddo
      dx=2.2d0/ni

      jp1=(1+1.1)/dx
      if (jp1>ni)jp1=ni
      jm1=(-1+1.1)/dx
      if (jm1<1)jm1=1


!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(msets,dt,rcind,rc,jp1,jm1,dx) 
!$OMP&   PRIVATE(b,bl,dind,nloops,absd,d,jl,j)
!$OMP&   reduction(+:ci) 

!$OMP DO
      do iset=1,msets-dt
        j=(rc(iset+dt)+1.1)/dx
        if (j<jm1)j=jm1
        if (j>jp1)j=jp1
        jl=(rc(iset)+1.1)/dx
        if (jl<jm1)jl=jm1
        if (jl>jp1)jl=jp1
        b=1
        if (mod(rcind(iset+dt),2)/=0)b=-1 ! (-1)**rcind
        bl=1
        if (mod(rcind(iset),2)/=0)bl=-1 ! (-1)**rcind
        dind=rcind(iset+dt)-rcind(iset)
        d=2*dind+b*rc(iset+dt)-bl*rc(iset)
        absd=abs(d)
        nloops=abs(dind)
        if (rcind(iset)==rcind(iset+dt))then
          if (j>jl) then
            ci(j)=ci(j)-absd
            ci(jl)=ci(jl)+absd
          endif
          if (j<jl) then
            ci(j)=ci(j)+absd
            ci(jl)=ci(jl)-absd
          endif
        else
          if (b==bl)then
            ci(jm1)=ci(jm1)+absd*nloops
            ci(jp1)=ci(jp1)-absd*nloops
          else  
            ci(jm1)=ci(jm1)+absd*(nloops-1)
            ci(jp1)=ci(jp1)-absd*(nloops-1)
          endif  
          if     (dind>0 .and. bl>0 .and. b<0)then !going through 1
            ci(jp1)=ci(jp1)-2*absd
            ci(j)=ci(j)+absd
            ci(jl)=ci(jl)+absd
          elseif (dind<0 .and. bl<0 .and. b>0)then !going through 1
            ci(jp1)=ci(jp1)-2*absd
            ci(j)=ci(j)+absd
            ci(jl)=ci(jl)+absd
          elseif (dind<0 .and. bl>0 .and. b<0)then !going through -1
            ci(jm1)=ci(jm1)+2*absd
            ci(j)=ci(j)-absd
            ci(jl)=ci(jl)-absd
          elseif (dind>0 .and. bl<0 .and. b>0)then !going through -1
            ci(jm1)=ci(jm1)+2*absd
            ci(j)=ci(j)-absd
            ci(jl)=ci(jl)-absd
          elseif (dind<0 .and. bl<0 .and. b<0)then !going through -1
            ci(j)=ci(j)-absd
            ci(jl)=ci(jl)+absd
          elseif (dind<0 .and. bl>0 .and. b>0)then !going through 1
            ci(j)=ci(j)+absd
            ci(jl)=ci(jl)-absd
          elseif (dind>0 .and. bl<0 .and. b<0)then !going through -1
            ci(j)=ci(j)+absd
            ci(jl)=ci(jl)-absd
          elseif (dind>0 .and. bl>0 .and. b>0)then !going through 1
            ci(j)=ci(j)-absd
            ci(jl)=ci(jl)+absd
          endif
        endif
      enddo      
!$OMP END PARALLEL
      s=0
      do j=1,ni
        s=s+ci(j)
        ci(j)=s/2/dt
      enddo
      ci(jp1)=ci(jp1-1)
      end      

      subroutine zc1ringldt(rc,rcind,nsets,cildt,ni,ldt)
      implicit none
      integer nsets,rcind(nsets),ni,dt,ldt
      real*8 rc(nsets)
      real*8 cildt(ni,ldt),ci(ni)
      
      integer j,jl,iset,jp1,jm1,nloops,idt

      real*8 s,dx,d,absd,b,bl,dind,absdnloops,r
      
      integer jset(nsets),bset(nsets),isetdt
      real*8 brcset(nsets)
      
      integer dt0,mt
      common /params2/dt0,mt
      
      dt0=1
      mt=2

      dx=2.2d0/ni

      jp1=(1+1.1)/dx
      if (jp1>ni)jp1=ni
      jm1=(-1+1.1)/dx
      if (jm1<1)jm1=1
      
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,jm1,jp1,rcind,bset,brcset,jset,dx,rc) 
!$OMP&   PRIVATE(b,dind,j)

!$OMP DO
      do iset=1,nsets
        j=(rc(iset)+1.1)/dx
        if (j<jm1)j=jm1
        if (j>jp1)j=jp1
        jset(iset)=j
        b=1
        if (mod(rcind(iset),2)/=0)b=-1 ! (-1)**rcind
        bset(iset)=b
        brcset(iset)=b*rc(iset)
      enddo  
!$OMP END PARALLEL
        
      
      do idt=1,ldt
        dt=dt0*mt**(idt-1)
      do j=1,ni
        ci(j)=0
      enddo
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dt,rcind,rc,jp1,jm1,brcset,bset,jset) 
!$OMP&   PRIVATE(b,bl,dind,absdnloops,nloops,absd,d,jl,j,isetdt)
!$OMP&   reduction(+:ci) 

!$OMP DO
      do iset=1,nsets-dt
        isetdt=iset+dt
        j=jset(isetdt)
        jl=jset(iset)
        b=bset(isetdt)
        bl=bset(iset)
        dind=rcind(isetdt)-rcind(iset)
        d=dind+dind+brcset(isetdt)-brcset(iset)
        absd=abs(d)
        nloops=abs(dind)
        if (rcind(iset)==rcind(iset+dt))then
          if (j>jl) then
            ci(j)=ci(j)-absd
            ci(jl)=ci(jl)+absd
          endif
          if (j<jl) then
            ci(j)=ci(j)+absd
            ci(jl)=ci(jl)-absd
          endif
        else
          if (b==bl)then
            absdnloops=absd*nloops
            ci(jm1)=ci(jm1)+absdnloops
            ci(jp1)=ci(jp1)-absdnloops
          else  
            absdnloops=absd*nloops-absd
            ci(jm1)=ci(jm1)+absdnloops
            ci(jp1)=ci(jp1)-absdnloops
          endif  
          if (dind>0) then
            if (bl>0) then
              if (b>0) then
                ci(j)=ci(j)-absd
                ci(jl)=ci(jl)+absd
              else
                ci(jp1)=ci(jp1)-2*absd
                ci(j)=ci(j)+absd
                ci(jl)=ci(jl)+absd
              endif
            else
              if (b>0) then            
                ci(jm1)=ci(jm1)+2*absd
                ci(j)=ci(j)-absd
                ci(jl)=ci(jl)-absd
              else
                ci(j)=ci(j)+absd
                ci(jl)=ci(jl)-absd
              endif
            endif
          else      
            if (bl>0) then
              if (b>0) then
                ci(j)=ci(j)+absd
                ci(jl)=ci(jl)-absd
              else
                ci(jm1)=ci(jm1)+2*absd
                ci(j)=ci(j)-absd
                ci(jl)=ci(jl)-absd
              endif
            else
              if (b>0) then            
                ci(jp1)=ci(jp1)-2*absd
                ci(j)=ci(j)+absd
                ci(jl)=ci(jl)+absd
              else
                ci(j)=ci(j)-absd
                ci(jl)=ci(jl)+absd
              endif
            endif
          endif          
        endif
      enddo      
!$OMP END PARALLEL
      s=0
      r=1./2/dt
      do j=1,ni
        s=s+ci(j)
        cildt(j,idt)=s*r
      enddo
      cildt(jp1,idt)=cildt(jp1-1,idt)
      enddo
      end   
      
      
         
      subroutine zc1TPsldt(rc,nsets,cildt,ni,ldt)
      implicit none
      integer nsets,ni,ldt
      real*8 rc(nsets)
      real*8 cildt(ni,ldt)
      
      integer j,jl,iset,iset1,iset2,i1,i2,jp1,jm1,idt,dt
      integer jset(nsets)
      logical isbdr(nsets)
      real*8 s,dx,absd,ci(ni),r
      
      integer dt0,mt
      common /params2/dt0,mt
      
      dt0=1
      mt=2

      dx=2.2d0/ni

      jp1=(1+1.1)/dx
      if (jp1>ni)jp1=ni
      jm1=(-1+1.1)/dx
      if (jm1<1)jm1=1
      
      
      
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,jm1,jp1,jset,dx,rc,isbdr) 
!$OMP&   PRIVATE(j)

!$OMP DO
      do iset=1,nsets
        isbdr(iset)=(rc(iset)>=1.or.rc(iset)<=-1)
        j=(rc(iset)+1.1)/dx
        if (j<jm1)j=jm1
        if (j>jp1)j=jp1
        jset(iset)=j
      enddo  
!$OMP END PARALLEL
      
      do idt=1,ldt
        dt=dt0*mt**(idt-1)
      do j=1,ni
        ci(j)=0
      enddo
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dt,isbdr,rc,jp1,jm1,jset) 
!$OMP&   PRIVATE(absd,jl,j,i1,i2,iset1,iset2)
!$OMP&   reduction(+:ci) 

!$OMP DO
      do iset1=1,nsets-1
        if ((.not. isbdr(iset1)) .and. iset1/=1)cycle ! start of a TP
        do iset2=iset1+1,nsets-1
          if (isbdr(iset2))exit ! end of a TP
        enddo
        ! padding traj at the start 
        do i1=iset1+1,iset2  !i1 if the right end of dt interval
          i2=max(iset1,i1-dt) ! i2 is the left end of dt interval
          j=jset(i1)
          jl=jset(i2)
          if (j==jl)cycle            
          absd=abs(rc(i1)-rc(i2))
          if (j>jl) then
            ci(j)=ci(j)-absd
            ci(jl)=ci(jl)+absd
          endif
          if (j<jl) then
            ci(j)=ci(j)+absd
            ci(jl)=ci(jl)-absd
          endif
        enddo
        ! padding traj at the end
        j=jset(iset2)
        do i1=max(iset1+1,iset2-dt+1),iset2-1  
          jl=jset(i1)
          if (j==jl)cycle
          absd=abs(rc(iset2)-rc(i1))
          if (j>jl) then
            ci(j)=ci(j)-absd
            ci(jl)=ci(jl)+absd
          endif
          if (j<jl) then
            ci(j)=ci(j)+absd
            ci(jl)=ci(jl)-absd
          endif
        enddo  
        if (iset2-iset1<dt) then !transitions between the end points
          j=jset(iset2)   
          jl=jset(iset1)
          if (j/=jl)then
          absd=abs(rc(iset2)-rc(iset1))*(dt-(iset2-iset1))
          if (j>jl) then
            ci(j)=ci(j)-absd
            ci(jl)=ci(jl)+absd
          endif
          if (j<jl) then
            ci(j)=ci(j)+absd
            ci(jl)=ci(jl)-absd
          endif
          endif
        endif
      enddo
!$OMP END PARALLEL
      s=0
      r=1./2/dt
      do j=1,ni
        s=s+ci(j)
        cildt(j,idt)=s*r
      enddo
      cildt(jp1,idt)=cildt(jp1-1,idt)
      enddo
      end  

      subroutine zh2(rc,msets,ci,ni)
      implicit none
      integer msets,ni
      real*8 rc(msets)
      real*8 ci(ni)
      
      integer j,jl,iset,jp1,jm1

      real*8 s,dx,absd
      
      ci=0
      dx=2.2d0/ni

      jp1=(1+1.1)/dx
      if (jp1>ni)jp1=ni
      jm1=(-1+1.1)/dx
      if (jm1<1)jm1=1


!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(msets,rc,jp1,jm1,dx) 
!$OMP&   PRIVATE(absd,jl,j)
!$OMP&   reduction(+:ci) 

!$OMP DO
      do iset=1,msets-1
        j=(rc(iset+1)+1.1)/dx
        if (j<jm1)j=jm1
        if (j>jp1)j=jp1
        jl=(rc(iset)+1.1)/dx
        if (jl<jm1)jl=jm1
        if (jl>jp1)jl=jp1
        if (j==jl)cycle
        absd=1/abs(rc(iset+1)-rc(iset))
        if (j>jl) then
            ci(j)=ci(j)-absd
            ci(jl)=ci(jl)+absd
        endif
        if (j<jl) then
          ci(j)=ci(j)+absd
          ci(jl)=ci(jl)-absd
        endif
      enddo
!$OMP END PARALLEL
      s=0
      do j=jm1,jp1
        s=s+ci(j)
        ci(j)=s/2
      enddo
      ci(jp1)=ci(jp1-1)
      end      


