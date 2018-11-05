! NonParametric Variational Eigenvector OPtimization 
! library of subroutines
      subroutine filter1(ev,nsets,ra,rb,changed)
      implicit none
      integer nsets
      real*8 ev(nsets),ra,rb,eval
      logical changed
      
      integer iset
      
      changed=.False.
      do iset=1,nsets
        if (ev(iset)<ra)then
          ev(iset)=ra
          changed=.True.
        endif
        if (ev(iset)>rb)then
          ev(iset)=rb
          changed=.True.
        endif
      enddo
      end

      subroutine optimdx2ev(ev,nev,nsets,npev,y,npy,rcfix,idt,eval)
      implicit none
      integer nev,npev,npy,nsets,rcfix(nsets),idt
      real*8 ev(nsets,nev),y(nsets),eval

      real*8 dxdx((npy+1)*(npy+1)+npev-npy,(npy+1)*(npy+1)+npev-npy)
      real*8 xx((npy+1)*(npy+1)+npev-npy,(npy+1)*(npy+1)+npev-npy)
      real*8 al((npy+1)*(npy+1)+npev-npy,(npy+1)*(npy+1)+npev-npy)

      integer info,lwork
      real*8 work(10*((npy+1)*(npy+1)+npev-npy))
      real*8 w((npy+1)*(npy+1)+npev-npy)

      real*8 yp1,rcp1,r,fij((npy+1)*(npy+1)+npev-npy)
      real*8 evalt,evt(nsets,nev)
      integer nij2,i12((npy+1)*(npy+1)+npev-npy)
      integer i21((npy+1)*(npy+1)+npev-npy),nij
      integer i,j,k,isize,iset,iev
      
      real*8 compeval
     
      isize=(npy+1)*(npy+1)+npev-npy
      lwork=isize*10
      
      dxdx=0
      xx=0
      
      call compxxev(nev,ev,npev,npy,y,rcfix,nsets,idt,xx,dxdx,isize)
      nij2=0
      i21=0
      do i=1,isize  ! select only non zero
        if (i>1 .and. dxdx(i,i)<1e-5)cycle
        nij2=nij2+1
        i12(nij2)=i
        i21(i)=nij2
      enddo
      do i=1,nij2
        do j=1,i
          dxdx(i,j)=dxdx(i12(i),i12(j))/nsets
          dxdx(j,i)=dxdx(i,j)
          xx(i,j)=xx(i12(i),i12(j))/nsets
          xx(j,i)=xx(i,j)
        enddo
      enddo
      call dsygv(1,'V','U',nij2,dxdx,isize,xx,isize,w,work,lwork,info)
      if (info/=0)then
c        write(*,*)'info=',info,nij2,npev,npy,i12(info-nij2)
        return 
        stop
      endif

      if (w(2)<1e-7)then
        write(*,*)'small w',w(2)
        return
      endif

      do j=1,nev
        do i=1,nij2
          al(i,j)=dxdx(i,1+j)  ! linear term should be positive
        enddo  
      enddo
      
!$OMP PARALLEL DO default(none) 
!$OMP&   SHARED(nsets,rcfix,npy,npev,y,ev,evt,i12,i21,nij2,al,nev) 
!$OMP&   PRIVATE(nij,rcp1,yp1,fij,r)

      do iset=1,nsets
        nij=0
        rcp1=1
        do i=0,npy
          yp1=rcp1
          do j=0,npy
            nij=nij+1
            fij(nij)=yp1
            if (rcfix(iset)/=0)then !leave only powers of ev
              yp1=0
            else
              yp1=yp1*y(iset)
            endif
          enddo
          rcp1=rcp1*ev(iset,1)
        enddo
        do i=npy+1,npev
          nij=nij+1
          fij(nij)=rcp1
          rcp1=rcp1*ev(iset,1)
        enddo
        do iev=1,nev
          r=0
          do i=1,nij2
            r=r+fij(i12(i))*al(i,iev)
          enddo
          evt(iset,iev)=r
        enddo
      enddo 
!$OMP END PARALLEL DO
      
      evalt=compeval(evt(1,1),nsets,idt)
      if (evalt>eval)then
        write(*,*)'skipping result ', eval,evalt,w(2)/2
        return
      endif
      r=0
      do iset=1,nsets,1000
        r=r+ev(iset,1)*evt(iset,1)
      enddo
      eval=evalt
      ev=evt
      if (r<0)ev=-evt
      end    


      subroutine compxxev(nev,ev,npev,npy,y,rcfix,nsets,dt,xx
     $ ,dxdx,isize)
      implicit none
      integer nev,npev,npy,nsets,rcfix(nsets),dt,isize
      real*8 ev(nsets,nev),y(nsets)

      real*8 dxdx(isize, isize)
      real*8 xx(isize, isize)
      real*8 fij(isize)
      

      real*8 yp1,yp2,rcp1,rcp2,r
      integer i,j,k,nij,ix,iset



      
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dt,rcfix,npy,npev,y,ev) 
!$OMP&   PRIVATE(ix,nij,rcp1,yp1,fij)
!$OMP&   reduction(+:xx) 

!$OMP DO
        do iset=1,nsets !computing xx
          nij=0
          rcp1=1
          do i=0,npy
            yp1=rcp1
            do j=0,npy
              nij=nij+1
              fij(nij)=yp1
              if (rcfix(iset)/=0)then !leave only powers of ev
                yp1=0
              else
                yp1=yp1*y(iset)
              endif
            enddo
            rcp1=rcp1*ev(iset,1)
          enddo
          do i=npy+1,npev
            nij=nij+1
            fij(nij)=rcp1
            rcp1=rcp1*ev(iset,1)
          enddo
          do j=1,nij
            do k=1,j
              xx(j,k)=xx(j,k)+fij(j)*fij(k)
            enddo
          enddo
        enddo
!$OMP END PARALLEL
          
          
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dt,rcfix,npy,npev,y,ev) 
!$OMP&   PRIVATE(ix,nij,rcp2,rcp1,yp2,yp1,fij,r)
!$OMP&   reduction(+:dxdx) 

!$OMP DO
        do iset=1,nsets-dt !computing dxdx
          nij=0
          rcp1=1
          rcp2=1
          do i=0,npy
            yp1=rcp1
            yp2=rcp2
            do j=0,npy
              nij=nij+1
              fij(nij)=yp2-yp1
              if (rcfix(iset)/=0)then !leave only powers of ev
                yp1=0
              else
                yp1=yp1*y(iset)
              endif
              if (rcfix(iset+dt)/=0)then !leave only powers of ev
                yp2=0
              else
                yp2=yp2*y(iset+dt)
              endif
            enddo
            rcp1=rcp1*ev(iset,1)
            rcp2=rcp2*ev(iset+dt,1)
          enddo
          do i=npy+1,npev
            nij=nij+1
            fij(nij)=rcp2-rcp1
            rcp1=rcp1*ev(iset,1)
            rcp2=rcp2*ev(iset+dt,1)
          enddo
          do j=1,nij
            do k=1,j
              dxdx(j,k)=dxdx(j,k)+fij(j)*fij(k)
            enddo
          enddo
        enddo
!$OMP END PARALLEL
      end

      subroutine anev(rc,nsets,idt,eval)
      implicit none
      integer nsets,idt
      real*8 rc(nsets),eval
      
      integer ni,nev
      parameter (ni=10000,nev=3)
      real*8 zc1(ni),zh(ni),u(ni),dudx(ni)

      integer jdt,i,i1,i2,iset,j,mev
      real*8 rcn(nsets),s,xmin,dx,evaln,compeval,s2
      real*8 dzc1,d(ni),e(ni-1),w(ni),z(ni,nev)

      real*8 ud(ni),ld(ni),rho(ni)
      
      real*8 work(5*ni)
      integer iwork(5*ni),ifail(ni),info
      
      
      call compzc1(rc,nsets,zc1,ni,xmin,dx,idt)
      call compzh2(rc,nsets,zh,ni,xmin,dx)
      
      do i=2,ni
        if (zc1(i)<1e-5)zc1(i)=zc1(i-1)
        if (zh(i)<1e-5)zh(i)=zh(i-1)
      enddo
      
      do i=1,ni
        if (i==1) then
          ud(i)=sqrt(abs(zc1(i+1)*zc1(i)))/abs(zh(i))
          d(i)=-ud(i)
        elseif(i==ni) then
          ld(i-1)=sqrt(abs(zc1(i-1)*zc1(i)))/abs(zh(i))
          d(i)=-ld(i-1)
        else
          ud(i)=sqrt(abs(zc1(i+1)*zc1(i)))/abs(zh(i))
          ld(i-1)=sqrt(abs(zc1(i-1)*zc1(i)))/abs(zh(i))
          d(i)=-ud(i)-ld(i-1)
        endif
      enddo
      !symmetrize
      do i=1,ni-1
        ud(i)=ud(i)*sqrt(abs(zh(i)/zh(i+1)))
        ld(i)=ld(i)*sqrt(abs(zh(i+1)/zh(i)))
      enddo
c          dstevx (JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL, M, W,
c           Z, LDZ, WORK, IWORK, IFAIL, INFO)
      call dstevx('V','I',ni,d,ld,0d0,0d0,ni-nev+1,ni,1d-15,mev,w,z,ni
     $            ,work,iwork,ifail,info)
      if (mev==0)return
c      write(*,*)mev,w(1:nev)/dx**2
      do j=1,nev
        do i=1,ni
          z(i,j)=z(i,j)/sqrt(zh(i))
        enddo
      enddo
      do i=2,ni
        u(i)=z(i-1,nev-1)
        dudx(i)=(z(i,nev-1)-z(i-1,nev-1))/dx
      enddo
      u(1)=u(2)
      dudx(1)=dudx(2)
          
      s=0
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,rc,dx,i1,i2,u,dudx,rcn,xmin) 
!$OMP&   PRIVATE(j)
!$OMP&   reduction(+:s) 

!$OMP DO
      do iset=1,nsets
        j=(rc(iset)-xmin)/dx
        if (j<1)j=1
        if (j>ni)j=ni
        rcn(iset)=u(j)+dudx(j)*(rc(iset)-j*dx-xmin)
        s=s+rcn(iset)
      enddo  
!$OMP END PARALLEL
      s=s/nsets
      rcn=rcn-s
      s=0
      do iset=1,nsets
        s=s+rcn(iset)**2
      enddo
      s=sqrt(s/nsets)
      rcn=rcn/s
      s=0
      do i=1,nsets,1000
        s=s+rc(i)*rcn(i)
      enddo
      if(s<0)rcn=-rcn
      
      evaln=compeval(rcn,nsets,1)
      if (evaln<eval) then
        rc=rcn
        eval=evaln
      endif
      end

      subroutine compzc1(rc,nsets,ci,ni,xmin,dx,dt)
      implicit none
      integer nsets,ni,dt
      real*8 rc(nsets),xmin,xmax,dx
      real*8 ci(ni)
      
      integer j,jl,iset
      integer jset(nsets)
      real*8 s,absd,r
      

      xmin=1e5
      xmax=-1e5
      do iset=1,nsets
        xmin=min(xmin,rc(iset))
        xmax=max(xmax,rc(iset))
      enddo
      dx=(xmax-xmin)/ni

      
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,iset,dx,xmin,ni,jset,rc) 
!$OMP&   PRIVATE(j)

!$OMP DO
      do iset=1,nsets
        j=(rc(iset)-xmin)/dx
        if (j<1)j=1
        if (j>ni)j=ni
        jset(iset)=j
      enddo  
!$OMP END PARALLEL
      
      ci=0
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dt,rc,jset) 
!$OMP&   PRIVATE(absd,jl,j,iset)
!$OMP&   reduction(+:ci) 

!$OMP DO
      do iset=1,nsets-dt
        j=jset(iset+dt)
        jl=jset(iset)
        absd=abs(rc(iset+dt)-rc(iset))
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
      r=1./2/dt
      do j=1,ni
        s=s+ci(j)
        ci(j)=s*r
      enddo
c      ci(ni)=ci(ni-1)
      end  

      subroutine writeev(name,ev,nsets,nev)
      implicit none
      character*(*) name
      integer nsets,nev
      real*8 ev(nsets,nev)
      
      character*1 dig
      integer iev,iset
      
      write(*,*)'writing eigenvectors'
      open(unit=1,file=name)
      do iset=1,nsets
        write(1,*)(ev(iset,iev),iev=1,nev)
      enddo
      close(1)
      end
      

      real*8 function compeval(ev,nsets,dt)
      implicit none
      integer nsets,dt
      real*8 ev(nsets)
      
      integer iset
      real*8 dx2,x2

      dx2=0
      x2=0
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dt,ev) 
!$OMP&   reduction(+:dx2,x2) 

!$OMP DO
      do iset=1,nsets-dt
        dx2=dx2+(ev(iset+dt)-ev(iset))**2
        x2=x2+ev(iset)**2
      enddo    
!$OMP END PARALLEL 
      if (dx2>x2*2)then
        compeval=777
      else
        compeval=-log(1-dx2/x2/2)/dt
      endif
      end

      subroutine writeevcriterion(name,ev,nsets,ldt,dt0)
      implicit none
      character*(*) name
      integer nsets,ldt,dt0
      real*8 ev(nsets)
      
      character*2 pref
      integer ni
      parameter(ni=1000)
      
      real*8 zc1i(ni,ldt),zc10i(ni),zh(ni),zh2(ni),dx,xmin
      real*8 sc(ldt),eval,compeval

      integer idt,i,dt,iset,j
      real*8 r

      

      call zc1ldt(ev,nsets,zc1i,ni,xmin,dx,ldt)
      call zc10(ev,nsets,zc10i,ni,xmin,dx)
      
      zh=0
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dx,ev,xmin) 
!$OMP&   PRIVATE(j)
!$OMP&   reduction(+:zh) 

!$OMP DO
      do iset=1,nsets
        j=(ev(iset)-xmin)/dx
        if (j>1000)j=1000
        if (j<1)j=1
        zh(j)=zh(j)+1
      enddo
!$OMP END PARALLEL 
      call compzh2(ev,nsets,zh2,ni,xmin,dx)
      
      eval=compeval(ev,nsets,dt0)
      do idt=1,ldt
        dt=2**(idt-1)
        if (eval*dt<10)then
          sc(idt)=(1.-exp(-eval*dt))/dt
        else
          sc(idt)=1./dt
        endif
      enddo
      open(unit=11,file=name)
      do i=1,ni-1
        do idt=1,ldt
        if (i>ni-10 .and. zc1i(i-1,idt) >1e-5 .and.
     $   abs(zc1i(i,idt)/zc10i(i)/(zc1i(i-1,idt)/zc10i(i-1))-1)>0.5)
     $   goto 10 ! to avoid boundary spikes
        enddo
        write(11,'(f8.4, 20f9.3)')dx*i+xmin
     $   ,-log(max(zh(i),1d0)/dx)
     $   ,-log(max(zh2(i),1d0)*2)
     $   ,(-log(zc1i(i,idt)/zc10i(i)/sc(idt)/2),idt=1,ldt)
        enddo
10    close(11)
      end  

      subroutine zc1ldt(rc,nsets,cildt,ni,xmin,dx,ldt)
      implicit none
      integer nsets,ni,ldt
      real*8 rc(nsets),xmin,xmax,dx
      real*8 cildt(ni,ldt)
      
      integer j,jl,iset,idt,dt
      integer jset(nsets)
      real*8 s,absd,ci(ni),r
      

      xmin=1e5
      xmax=-1e5
      do iset=1,nsets
        xmin=min(xmin,rc(iset))
        xmax=max(xmax,rc(iset))
      enddo
      dx=(xmax-xmin)/ni

      
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,iset,dx,xmin,ni,jset,rc) 
!$OMP&   PRIVATE(j)

!$OMP DO
      do iset=1,nsets
        j=(rc(iset)-xmin)/dx
        if (j<1)j=1
        if (j>ni)j=ni
        jset(iset)=j
      enddo  
!$OMP END PARALLEL
      
      do idt=1,ldt
        dt=2**(idt-1)
        ci=0
!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,dt,rc,jset) 
!$OMP&   PRIVATE(absd,jl,j,iset)
!$OMP&   reduction(+:ci) 

!$OMP DO
      do iset=1,nsets-dt
        j=jset(iset+dt)
        jl=jset(iset)
        absd=abs(rc(iset+dt)-rc(iset))
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
      r=1./2/dt
      do j=1,ni
        s=s+ci(j)
        cildt(j,idt)=s*r
      enddo
      enddo
      end  

      subroutine zc10(rc,nsets,ci,ni,xmin,dx)
      implicit none
      integer nsets,ni
      real*8 rc(nsets),xmin,dx
      real*8 ci(ni)
      
      integer j,j0,iset
      real*8 s,absd
      

      j0=(0-xmin)/dx
      ci=0

!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(nsets,rc,j0,xmin,dx,ni) 
!$OMP&   PRIVATE(absd,j,iset)
!$OMP&   reduction(+:ci) 

!$OMP DO
      do iset=1,nsets
        j=(rc(iset)-xmin)/dx
        if (j<1)j=1
        if (j>ni)j=ni
        absd=abs(rc(iset))
        if (j>j0) then
          ci(j)=ci(j)-absd
          ci(j0)=ci(j0)+absd
        endif
        if (j<j0) then
          ci(j)=ci(j)+absd
          ci(j0)=ci(j0)-absd
        endif
      enddo
!$OMP END PARALLEL
      s=0
      do j=1,ni
        s=s+ci(j)
        ci(j)=s/2
      enddo
      end  

      subroutine compzh2(rc,msets,ci,ni,xmin,dx)
      implicit none
      integer msets,ni
      real*8 rc(msets),xmin,dx
      real*8 ci(ni)
      
      integer j,jl,iset,jp1,jm1

      real*8 s,absd
      
      ci=0
      jm1=1
      jp1=ni

!$OMP PARALLEL  default(none) 
!$OMP&   SHARED(msets,rc,jp1,jm1,dx,xmin) 
!$OMP&   PRIVATE(absd,jl,j)
!$OMP&   reduction(+:ci) 

!$OMP DO
      do iset=1,msets-1
        j=(rc(iset+1)-xmin)/dx
        if (j<jm1)j=jm1
        if (j>jp1)j=jp1
        jl=(rc(iset)-xmin)/dx
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


      subroutine readxyz(dcd,x,y,z,natom,natom2,nsets,i21)
      implicit none 
      character*(*) dcd
      integer natom,natom2,i21(natom2),nsets
      real*4 x(nsets,natom2),y(nsets,natom2),z(nsets,natom2)
      real*4 temp(natom),t(12)
        
      integer i,iset
      logical ft
      call opendcd(dcd,1,natom,i,ft)
      if (natom==natom2)then
        do i=1,natom
          i21(i)=i
        enddo
      else
        call select_atoms(natom,natom2,i21)
      endif
      do iset=1,nsets
        if(ft) read(1)t
        read(1)temp
        do i=1,natom2
          x(iset,i)=temp(i21(i))
        enddo  
        read(1)temp
        do i=1,natom2
          y(iset,i)=temp(i21(i))
        enddo  
        read(1)temp
        do i=1,natom2
          z(iset,i)=temp(i21(i))
        enddo  
      enddo
      close(1)
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
        if (r>1000)then
          write(*,*)'r',r
          write(*,*)x(iset,i1),y(iset,i1),z(iset,i1)
          write(*,*)x(iset,j1),y(iset,j1),z(iset,j1)
        endif
        rij(iset)=r
      enddo
!$OMP END PARALLEL
      end

      subroutine select_atoms(natom,natom2,i21)
      implicit none
      integer natom,natom2,i21(natom2)
      
      integer i,l,taken(natom)
      do i=1,natom
          taken(i)=0
      enddo
      do l=1,natom2
10      i=rand()*natom+1
        if (taken(i)/=0)goto 10
        taken(i)=1
        i21(l)=i
      enddo  
      end

