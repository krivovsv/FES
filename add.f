      real*8 function disp(ev,nsets,ldt)
      implicit none
      integer nsets,ldt,dt0
      real*8 ev(nsets)
      
      character*2 pref
      integer ni
      parameter(ni=1000)
      
      real*8 zc1i(ni,ldt),dx,xmin
      real*8 s

      integer idt,i,dt,iset,j
      real*8 zz,zzm,zz2m

      

      call zc1ldt(ev,nsets,zc1i,ni,xmin,dx,ldt)
      
      disp=0
      do idt=2,ldt
        zzm=0
        zz2m=0
        do i=51,ni-50
          zz=log(zc1i(i,idt)/zc1i(i,1))
          zzm=zzm+zz
          zz2m=zz2m+zz**2
        enddo
        zzm=zzm/(ni-100)
        zz2m=zz2m/(ni-100)-zzm**2
        disp=max(disp,sqrt(zz2m))
      enddo
      end
          

      subroutine optimdx2ev2(ev,nev,nsets,npev,y,npy,rcfix,idt,eval,cd)
      implicit none
      integer nev,npev,npy,nsets,rcfix(nsets),idt
      real*8 ev(nsets,nev),y(nsets),eval,cd

      real*8 dxdx((npy+1)*(npy+1)+npev-npy,(npy+1)*(npy+1)+npev-npy)
      real*8 dxdx0((npy+1)*(npy+1)+npev-npy,(npy+1)*(npy+1)+npev-npy)
      real*8 xx((npy+1)*(npy+1)+npev-npy,(npy+1)*(npy+1)+npev-npy)
      real*8 xx0((npy+1)*(npy+1)+npev-npy,(npy+1)*(npy+1)+npev-npy)
      real*8 al((npy+1)*(npy+1)+npev-npy,(npy+1)*(npy+1)+npev-npy)

      integer info,lwork
      real*8 work(10*((npy+1)*(npy+1)+npev-npy))
      real*8 w((npy+1)*(npy+1)+npev-npy)

      real*8 yp1,rcp1,r,fij((npy+1)*(npy+1)+npev-npy)
      real*8 evalt,evt(nsets,nev)
      integer nij2,i12((npy+1)*(npy+1)+npev-npy)
      integer i21((npy+1)*(npy+1)+npev-npy),nij
      integer i,j,k,isize,iset,iev
      
      real*8 compeval,disp,dn
      real*8 ym
      integer iskip
     
      isize=(npy+1)*(npy+1)+npev-npy
      lwork=isize*10
      
      dxdx0=0
      xx0=0
      
      
      call compxxev(nev,ev,npev,npy,y,rcfix,nsets,idt,xx0,dxdx0,isize)
      iskip=0
10    nij2=0
      i21=0
      do i=1,isize-iskip  ! select only non zero
        if (i>1 .and. dxdx0(i,i)<1e-5)cycle
        nij2=nij2+1
        i12(nij2)=i
        i21(i)=nij2
      enddo
      do i=1,nij2
        do j=1,i
          dxdx(i,j)=dxdx0(i12(i),i12(j))/nsets
          dxdx(j,i)=dxdx(i,j)
          xx(i,j)=xx0(i12(i),i12(j))/nsets
          xx(j,i)=xx(i,j)
        enddo
      enddo
      call dsygv(1,'V','U',nij2,dxdx,isize,xx,isize,w,work,lwork,info)
      if (info/=0)then
        iskip=iskip+1
        if (iskip<isize-2)goto 10
        write(*,*)'info=',info,nij2,npev,npy,i12(info-nij2)
        return
      endif

      if (w(2)<1e-7)then
        write(*,*)'small w',w(2)
        return
      endif

      do j=1,nev
        do i=1,nij2
          al(i,j)=dxdx(i,1+j) 
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
      if (evalt<1e-7)then
        write(*,*) 'w<0',w(1:3)/2
        return
      endif
c      write(*,*)'result ', eval,evalt,-log(1-w(2)/2)/idt,iskip
c     $ ,isize-iskip
      if (evalt>eval)then
        iskip=iskip+1
        if (iskip<isize-2)goto 10
        write(*,*)'skipping result ', eval,evalt,-log(1-w(2)/2)/idt
     $   ,iskip,isize-iskip
        return
      endif
      dn=disp(evt,nsets,16)
c      write(*,*)'disp',dn,cd,evalt
      if (eval<1e-3.and.dn>cd .and. exp(-(dn-cd)/0.0002)<rand())return
      r=0
      do iset=1,nsets,1000
        r=r+ev(iset,1)*evt(iset,1)
      enddo
      eval=evalt
      ev=evt
      if (r<0)ev=-evt
      cd=dn
      end    

