  module comm
    use omp_lib
    implicit none

!++ contant
    integer,parameter   :: r_size=8
    real(r_size),parameter   :: pi=atan(1.d0)*4.d0

!++ model param
    integer,parameter   :: IMAX=96,JMAX=48,KMAX=8   ! X,Y,Z
    integer,parameter   :: edim=<EDIM>                  ! ensemble size
    integer,parameter   :: mdim=(IMAX*JMAX*KMAX)*9 &! UVTQ
                               +(IMAX*JMAX)*32      ! many other 2D vars !
                                                    ! model state dim.
    integer,parameter   :: hdim=KMAX*9+32           ! history
    integer,parameter   :: NGP=IMAX*JMAX
    integer,parameter   :: mdim_small=NGP*KMAX*4+NGP*3  ! TUVQ,SP,LST,SST

!++ DA vars
    real(r_size)        :: xa(mdim,edim)    ! analysis
    real(r_size)        :: xf(mdim,edim)    ! forecast
    real(r_size),allocatable    :: obs(:)
    real(r_size),allocatable    :: err(:)
    real(r_size),allocatable    :: hx(:,:)
    integer,allocatable :: ijobs_1(:),ijobs_2(:)
    integer             :: nobs             ! number of observation

    real(r_size),parameter  :: undef=9.99e19
    real(r_size),parameter  :: c_local=<LOCAL>    ! km for 1dy
    real(r_size),parameter  :: roi=c_local*2.d0   ! radius of influenc (km)
    real(r_size),parameter  :: delta=0.d0     ! covariance inflation (MULT)
    real(r_size),parameter  :: alpha=<INFL>   ! covariance inflation 
    character*4         :: c_infl='<INFLTYPE>'    ! covariance inflation

!++ report
    real(r_size)        :: xf_save(mdim,edim)

!++
    character           :: cobsnet*256      ! observation network
    character           :: cobserr*256      ! observation error
    character           :: cobs*256         ! observation
    character           :: cobsope*256      ! Hxf
    character           :: cgues*256        ! xf
    character           :: cdate*8,cyyyy*4,cmm*2,cdd*2,chh*2
    


  end module comm

!>-----------------------------------------------------------------------
!> main 
  program main
    use comm
    implicit none

!++ read settings
    call get_args
!++ read forecast
    call read_fcst
!!++ read obs.
    call read_yobs
!!++ convert X to HX
    call trans_XtoY     ! in reality, this just reads file
!!++ assimilation!
    call letkf
!++ write analysis
    call write_anl
!++ flash report 
    call report
!++ check
    !call check

  contains

!>-----------------------------------------------------------------------
!> read args
  subroutine get_args
    use comm
    implicit none

    call getarg(1,cdate)
    call getarg(2,cobsnet)
    call getarg(3,cobs)
    call getarg(4,cobserr)
    call getarg(5,cgues)
    call getarg(6,cobsope)

    cyyyy = cdate(1:4)
    cmm   = cdate(5:6)
    cdd   = cdate(7:8)
    chh   = "00"

    return

  end subroutine get_args

  !>-----------------------------------------------------------------------
  !> read forecast data
  subroutine read_fcst
    use comm
    implicit none

    real(r_size) :: mean(mdim),sprd(mdim)
    real(4) :: r4var(NGP)
    integer :: mem,ij,n
    character :: cmem*3,cif*256,cof*256
    
    do mem=1,edim
      write(cmem,'(I3.3)') mem
      cif='gues/'//cmem//'/'//trim(adjustl(cgues))
      open(33,file=cif,form='unformatted',access='sequential')
      do n=1,hdim
        read(33) r4var
        xf((n-1)*NGP+1:n*NGP,mem)=dble(r4var(:))
      enddo
    enddo

    xf_save = xf

    ! Stats
    call stats(xf,mean,sprd)

    ! Write fcst ensemble mean
    cof='gues/mean/'//trim(adjustl(cgues))
    open(78,file=cof,form='unformatted',access='sequential')
    do n=1,hdim
      r4var(:)=real(mean((n-1)*NGP+1:n*NGP))
      write(78) r4var
    enddo
    close(78)

    ! Write fcst ensemble spread
    cof='gues/sprd/'//trim(adjustl(cgues))
    open(79,file=cof,form='unformatted',access='sequential')
    do n=1,hdim
      r4var(:)=real(sprd((n-1)*NGP+1:n*NGP))
      write(79) r4var
    enddo
    close(79)

    return

  end subroutine read_fcst

  !>-----------------------------------------------------------------------
  !> read obs. for OSSE
  !> read full fields of all the variables
  subroutine read_yobs
    use comm
    implicit none

    real(r_size) :: flag(NGP,hdim),ytmp(NGP,hdim),yerr(NGP,hdim)
    real(4) :: r4var(NGP)
    integer :: p,ij,n

!++ read observation network
    open(22,file=cobsnet,form='unformatted',access='sequential')
    do n=1,hdim
      read(22) r4var
      flag(:,n)=dble(r4var(:))
    enddo
    close(22)

!++ read observation 
    open(22,file=cobs   ,form='unformatted',access='sequential')
    do n=1,hdim
      read(22) r4var
      ytmp(:,n)=dble(r4var(:))
    enddo
    close(22)

!++ read observation error
    open(22,file=cobserr,form='unformatted',access='sequential')
    do n=1,hdim
      read(22) r4var
      yerr(:,n)=dble(r4var(:))
    enddo
    close(22)

!++ count Nobs
    p=0
    do n=1,hdim ; do ij=1,IMAX*JMAX
      if( flag(ij,n).eq.1.d0 ) p=p+1
    enddo ; enddo
    nobs=p

    allocate(obs(nobs),err(nobs),ijobs_1(nobs),ijobs_2(nobs))

!++ creat yobs
    p=0
    do n=1,hdim ; do ij=1,IMAX*JMAX
      if( flag(ij,n).eq.1.d0 )then
        p=p+1
        obs(p)=dble(ytmp(ij,n))
        err(p)=dble(yerr(ij,n))
        ijobs_1(p)=ij
        ijobs_2(p)=n
      endif
    enddo ; enddo
    
    return

  end subroutine read_yobs

!>-----------------------------------------------------------------------
!> write analysis data
  subroutine write_anl
    use comm
    implicit none

    real(r_size) :: mean(mdim),sprd(mdim)
    real(4) :: r4var(NGP)
    integer :: mem,ij,n
    character :: cmem*3,cof*256,canal*256

    canal = cgues

    ! Analysis (member)
    do mem=1,edim
      write(cmem,'(I3.3)') mem
      cof='anal/'//cmem//'/'//trim(adjustl(canal))
      open(77,file=cof,form='unformatted',access='sequential')
      do n=1,hdim
        r4var(:)=real(xa((n-1)*NGP+1:n*NGP,mem))
        write(77) r4var
      enddo
      close(77)
    enddo

    ! Stats
    call stats(xa,mean,sprd)

    ! Analysis (ensemble mean)
    cof='anal/mean/'//trim(adjustl(canal))
    open(78,file=cof,form='unformatted',access='sequential')
    do n=1,hdim
      r4var(:)=real(mean((n-1)*NGP+1:n*NGP))
      write(78) r4var
    enddo
    close(78)

    ! Analysis (spread)
    cof='anal/sprd/'//trim(adjustl(canal))
    open(79,file=cof,form='unformatted',access='sequential')
    do n=1,hdim
      r4var(:)=real(sprd((n-1)*NGP+1:n*NGP))
      write(79) r4var
    enddo
    close(79)

    return

  end subroutine write_anl

!>-----------------------------------------------------------------------
!> convert X to Y
!> make annual using previous 1yr
  subroutine trans_XtoY
    use comm
    implicit none

    real(4) :: r4var(NGP),hist(NGP,hdim)
    integer :: cnt,mem,yyyy,mm,dd,hh,days,p,ij,n
    character :: cyyyy2*4,cmem*3

    allocate(hx(nobs,edim))

    read(cyyyy,'(I4)') yyyy

    do mem=1,edim
      write(cmem,'(I3.3)') mem
  
      open(33,FILE='gues/'//cmem//'/'//trim(adjustl(cobsope)),FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      do n=1,hdim
        read (33) R4VAR
        hist (:,n) = R4VAR(:)
      enddo
      close(33)
  
      do p=1,nobs
        hx(p,mem)=dble(hist(ijobs_1(p),ijobs_2(p)))
      enddo


    enddo ! mem
    
    return

  end subroutine trans_XtoY
!>-----------------------------------------------------------------------
!> return the number of days in the month specified
  subroutine endday(yy,mm,dd)
    implicit none
    
    integer,intent(in) :: yy,mm
    integer,intent(out) :: dd
    integer :: days(12)

    data days / 31,28,31,30,31,30,31,31,30,31,30,31 /

    dd=days(mm)
    if( mm.eq.2.and.mod(yy,4).eq.0 )then
      dd=29
      if(mod(yy,100).eq.0) dd=28
      if(mod(yy,400).eq.0) dd=29
    endif

    return

  end subroutine endday

!>-----------------------------------------------------------------------
!> local ensemble transform Kalman filter
!> Hunt et al. (2007)
  subroutine letkf
    use comm
    implicit none

    !--- 
    integer,parameter :: nxy=imax*jmax
    real(r_size) :: Xfm(mdim)
    real(r_size) :: dXf(mdim,edim)
    real(r_size) :: mI (edim,edim)
    real(r_size) :: Pa_tilda(edim,edim) ! Pa in ensemble space
    real(r_size) :: V(edim,edim)    ! Eigen vector
    real(r_size) :: D(edim,edim)    ! Eigen value
    ! LAPACK
    integer,parameter :: lda = edim
    integer,parameter :: lwork = (edim*edim+2)*edim
    integer,parameter :: iwork = (edim*edim+2)*edim
    integer,parameter :: liwork = (edim*edim+2)*edim
    !integer,parameter :: lwork = (edim*edim+2)*edim*edim
    integer :: info
    real(r_size) :: eval(edim)
    real(r_size) :: work(lwork)
    real(r_size) :: Vtmp(edim,edim)
    ! global
    real(r_size),allocatable :: HXm (:)
    real(r_size),allocatable :: dY  (:,:)
    ! local 
    real(r_size) :: yo_l_tmp(nxy)
    real(r_size) :: dy_l_tmp(nxy,edim)
    real(r_size) :: HX_l_tmp(nxy)
    real(r_size) :: R_inv_tmp(nxy,nxy)
    real(r_size) :: Xfm_l
    real(r_size) :: dXf_l(edim)
    real(r_size) :: Xam_l
    real(r_size) :: dXa_l(edim)
    real(r_size),allocatable :: yo_l(:)
    real(r_size),allocatable :: HX_l(:)
    real(r_size),allocatable :: dY_l(:,:)
    real(r_size),allocatable :: R_inv(:,:)
    real(r_size),allocatable :: K(:)
    !
    real(r_size) :: dist,wgt,errobs
    integer :: nobs_l

    integer :: p,ij,i,j,cnt, ens,ii
    integer :: ij_small(mdim_small)

    !--- update simultaneously

    allocate(HXm(nobs),dY(nobs,edim))
    call decompose_x(xf,Xfm,dXf)
    call decompose_y(hx,HXm,dY)
    call unitmat(edim,mI)

    !-- make list
    cnt=0
    do ij=NGP*KMAX+1,NGP*KMAX*5 ! TUVQ
      cnt=cnt+1
      ij_small(cnt)=ij
    enddo
    do ij=NGP*KMAX*9+1,NGP*KMAX*9+NGP ! SP
      cnt=cnt+1
      ij_small(cnt)=ij
    enddo
    do ij=NGP*KMAX*9+NGP*14+1,NGP*KMAX*9+NGP*16 ! LST,SST
      if( all(xf(ij,:) /= undef) )then
        cnt=cnt+1
        ij_small(cnt)=ij
      endif
    enddo


    !--- initial value
    xa = xf

!$omp parallel default(private),shared(xa)
!$omp do
    do i=1,mdim_small
      ij=ij_small(i)

      !--- covariance inflation (multiplicative)
      select case(c_infl)
        case('MULT')
          call mult(dXf(ij,:))
      end select

      !---
      xfm_l = xfm(ij)
      dxf_l = dXf(ij,:)

      !--- count num. of obs
      cnt=0
      do p=1,nobs
        call getdistance_t30(dist,ij,ijobs_1(p))
        if( dist.lt.roi )then
          errobs=err(p) 
          cnt=cnt+1
          call gc99(dist,wgt)
          yo_l_tmp(cnt)=obs(p)
          HX_l_tmp(cnt)=HXm(p)
          dY_l_tmp(cnt,:)=dY(p,:)
          R_inv_tmp(cnt,cnt)=wgt/errobs
        endif
      enddo

      nobs_l = cnt  ! num obs in local patch
      if( nobs_l.eq.0 )then
        cycle
      endif

      !--- allocate for local patch
      allocate( &
                yo_l (nobs_l), &
                HX_l (nobs_l), &
                K    (nobs_l), &
                dY_l (nobs_l,edim), &
                R_inv(nobs_l,nobs_l) &
              )

      R_inv = 0.d0
      do p=1,nobs_l
        yo_l(p) = yo_l_tmp(p)
        HX_l(p) = HX_l_tmp(p)
        dY_l(p,:) = dY_l_tmp(p,:)
        R_inv(p,p) = R_inv_tmp(p,p)
      enddo

      !--- eigen decomposition
      Pa_tilda = dble(edim-1)/sqrt(1.d0+delta)*mI + matmul(matmul(transpose(dY_l),R_inv),dY_l)
      V = 0.d0
      do ii=1,edim
        do j=ii,edim
          V(ii,j) = Pa_tilda(ii,j)
        enddo
      enddo

      info = 0
      call dsyev('V','U',edim,V,lda,eval,work,lwork,info)
      if(eval(1).le.0.d0) cycle

      !--- mean update
      D = 0.d0
      do ens=1,edim
        D(ens,ens) = 1.d0/eval(ens)  ! inv D
      enddo
      Pa_tilda = matmul(matmul(V,D),transpose(V))
      K = matmul(matmul(matmul(dXf_l,Pa_tilda),transpose(dY_l)),R_inv)
      Xam_l = Xfm_l
      do p=1,nobs_l
        Xam_l = Xam_l + K(p)*(yo_l(p)-HX_l(p))
      enddo
  
      !--- ensemble update
      D = 0.d0
      do ens=1,edim
        D(ens,ens) = sqrt(dble(edim-1)/dble(eval(ens)))    ! (inv D)**1/2
      enddo
      Pa_tilda = matmul(matmul(V,D),transpose(V))
      dXa_l = matmul(dXf_l,Pa_tilda)
  
      ! update (only center of patch)
      do ens=1,edim
        xa(ij,ens) = Xam_l + dXa_l(ens)
      enddo

      !--- covariance inflation
      select case(c_infl)
        case('RTPS')
          call rtps(xf_save(ij,:),xa(ij,:))
        case('RTPP')
          call rtpp(xf_save(ij,:),xa(ij,:))
      end select

      deallocate(yo_l,dY_l,HX_l,R_inv,K)

    enddo ! ij-loop
!$omp end do
!$omp end parallel

    !call rtpp(xf_save,xa)
    !call rtps(xf_save,xa)
    !call mult(xa)

    deallocate(HXm,dY)

    return

  end subroutine letkf


!>-----------------------------------------------------------------------
  subroutine getdistance_t30(dist,i,j)
    use comm
    implicit none

    integer,intent(in) :: i,j
    real(r_size),intent(out) :: dist
    real(r_size),parameter ::  dxs = 3.75      !! delta longitude
    real(r_size),parameter ::  x0  = 1.875     !!  
    real(r_size),parameter ::  req = 6378.137  !! equatorial radius
    real(r_size),parameter ::  rpo = 6356.752  !! polar radius
    real(r_size),parameter :: d2r=pi/180.d0
    real(r_size) :: dx,dy,x1,y1,x2,y2
    real(r_size) :: myuy,e2,w,m,n
    integer      ::  ix,iy,jx,jy
    real(r_size) :: lat0(48)
    data      lat0    /   &
        -87.159,-83.479,-79.777,-76.070,-72.362,-68.652,-64.942,-61.232, & 
        -57.521,-53.810,-50.099,-46.389,-42.678,-38.967,-35.256,-31.545, & 
        -27.833,-24.122,-20.411,-16.700,-12.989, -9.278, -5.567, -1.856, & 
          1.856,  5.567,  9.278, 12.989, 16.700, 20.411, 24.122, 27.833, & 
         31.545, 35.256, 38.967, 42.678, 46.389, 50.099, 53.810, 57.521, & 
         61.232, 64.942, 68.652, 72.362, 76.070, 79.777, 83.479, 87.159/

    ix=mod(i,IMAX)
    if (ix.eq.0) ix=IMAX
    iy=mod((i-1)/IMAX+1,JMAX)
    if (iy.eq.0) iy=JMAX
    jx=mod(j,IMAX)
    if (jx.eq.0) jx=IMAX
    jy=mod((j-1)/IMAX+1,JMAX)
    if (jy.eq.0) jy=JMAX

    x1=(x0+(ix-1)*dxs)*d2r 
    y1=(lat0(iy))*d2r
    x2=(x0+(jx-1)*dxs)*d2r 
    y2=(lat0(jy))*d2r

    dx=x1-x2
    if(abs(dx).gt.pi) dx=2*pi-abs(dx)
    dy=y1-y2
    myuy=(y1+y2)/2

    e2=(req**2-rpo**2)/(req**2)
    w=sqrt(1-e2*(sin(myuy)**2))
    m=req*(1-e2)/(w**3)
    n=req/w

    dist=sqrt((dy*m)**2+(dx*n*cos(myuy))**2)

    return

  end subroutine getdistance_t30
!>-----------------------------------------------------------------------
  subroutine gc99(dist,l)
    use comm
    implicit none

    real(r_size),intent(in) :: dist ! physical distance
    real(r_size),intent(out) ::   l ! weight
    real(r_size) :: r

    r = dist/c_local
    if (r.ge.0.and.r.le.1) then
      l = 1.d0 - 0.25*(r**5) + 0.5*(r**4) + 0.625*(r**3) -5.d0/3.d0*(r**2)
    else if (r.gt.1.and.r.le.2) then
      l = (r**5)/12.d0 - 0.5d0*(r**4) + 0.625*(r**3) + 5.d0/3.d0*(r**2) - 5.d0*r + 4.d0 - 2.d0/3.d0/r
    else
      l = 0
    endif

    return

    end subroutine gc99
!>-----------------------------------------------------------------------
  subroutine decompose_x(xorg,xmean,xanom)
    use comm
    implicit none

    real(r_size),intent(in) :: xorg(:,:)
    real(r_size),intent(out) :: xmean(:) ! ensemble mean
    real(r_size),intent(out) :: xanom(:,:) !

    integer :: mem

    xmean=0.d0
    do mem=1,edim
      xmean(:)=xmean(:)+xorg(:,mem)
    enddo
    xmean=xmean/edim

!$omp parallel
!$omp do private(mem)
    do mem=1,edim
      xanom(:,mem)=xorg(:,mem)-xmean(:)
    enddo
!$omp end do
!$omp end parallel

    return

  end subroutine decompose_x

!>-----------------------------------------------------------------------
  subroutine decompose_y(yorg,ymean,yanom)
    use comm
    implicit none

    real(r_size),intent(in) :: yorg(:,:)
    real(r_size),intent(out) :: ymean(:) ! ensemble mean
    real(r_size),intent(out) :: yanom(:,:) !

    integer :: mem

    ymean=0.d0
    do mem=1,edim
      ymean(:)=ymean(:)+yorg(:,mem)
    enddo
    ymean=ymean/edim

!$omp parallel
!$omp do private(mem)
    do mem=1,edim
      yanom(:,mem)=yorg(:,mem)-ymean(:)
    enddo
!$omp end do
!$omp end parallel

    return

  end subroutine decompose_y
!>-----------------------------------------------------------------------
  subroutine combine_x(xmean,xanom,xorg)
    use comm
    implicit none

    real(r_size),intent(in) :: xmean(:) ! ensemble mean
    real(r_size),intent(in) :: xanom(:,:) !
    real(r_size),intent(out) :: xorg(:,:)

    integer :: mem
!$omp parallel
!$omp do private(mem)
    do mem=1,edim
      xorg(:,mem)=xmean(:)+xanom(:,mem)
    enddo
!$omp end do
!$omp end parallel

    return

  end subroutine combine_x

!>-----------------------------------------------------------------------
  !> return m x m dimension unit matrix
  subroutine unitmat( &
        m, &
        I )
    implicit none

    integer,intent(in) :: m
    real(r_size),intent(out):: I(m,m)

    integer :: ii

    I = 0.d0
    do ii=1,m
      I(ii,ii) = 1.d0
    enddo

    return

  end subroutine unitmat

!>-----------------------------------------------------------------------
  !> Calc. ensemble mean and spread
  subroutine stats( &
        x, &
        xmean, &
        xsprd )
    implicit none

    real(r_size),intent(in)  :: x(mdim,edim)
    real(r_size),intent(out) :: xmean(mdim)
    real(r_size),intent(out) :: xsprd(mdim)

    integer :: i,mem

    ! Ensemble mean
    xmean=0.d0
    do mem=1,edim
      xmean(:) = xmean(:) + x(:,mem)
    enddo
    xmean=xmean/edim

    ! Ensemble spread
    xsprd=0.d0
    do mem=1,edim
      xsprd(:) = xsprd(:) + (x(:,mem)-xmean(:))**2
    enddo
    xsprd=sqrt(xsprd/edim)

    return

  end subroutine stats
  !>-----------------------------------------------------------------------
  !> multipricative
  subroutine mult( dxb )
    use comm, only : r_size, edim, alpha
    real(r_size),intent(inout) :: dxb(1,edim)
    integer :: i,e

    dxb = dxb * alpha

    return
  end subroutine mult
  !>-----------------------------------------------------------------------
  !> Relaxation to prior pertabation 
  subroutine rtpp( xb,xa )
    use comm, only : r_size, edim, alpha
    implicit none
    real(r_size),intent(in) :: xb(1,edim)
    real(r_size),intent(inout) :: xa(1,edim)
    real(r_size) :: mxa(1)
    real(r_size) :: mxb(1)
    real(r_size) :: dxa(1,edim)
    real(r_size) :: dxb(1,edim)


    call decompose_x(xb,mxb,dxb)
    call decompose_x(xa,mxa,dxa)
    xa(1,:)=mxa(1)+(1.d0-alpha)*dxa(1,:)+alpha*dxb(1,:)

    return
  end subroutine rtpp
  !>-----------------------------------------------------------------------
  !> Relaxation to prior spreada (Whitaker and Hamill, 2012)
  subroutine rtps( xb,xa )
    use comm, only : r_size, edim, alpha
    implicit none
    real(r_size),intent(in) :: xb(1,edim)
    real(r_size),intent(inout) :: xa(1,edim)
    real(r_size) :: mxa(1)
    real(r_size) :: mxb(1)
    real(r_size) :: dxa(1,edim)
    real(r_size) :: dxb(1,edim)
    real(r_size) :: a,b,inffct
    integer :: i,e


    call decompose_x(xb,mxb,dxb)
    call decompose_x(xa,mxa,dxa)
    a=0.d0
    b=0.d0
    do e=1,edim
      a=a+dxa(1,e)*dxa(1,e)
      b=b+dxb(1,e)*dxb(1,e)
    enddo
    a=sqrt(a)
    b=sqrt(b)
    if(a<=0.d0.or.b<=0.d0) return
    inffct=(b-a)/a*alpha+1.d0
    if(inffct.gt.1.1) inffct=1.1
    dxa=inffct*dxa
    xa(1,:)=mxa(1)+dxa(1,:)

    return
  end subroutine rtps
  !>-----------------------------------------------------------------------
  !>
  subroutine check
    use comm
    implicit none

    integer :: i,e

    do i=ngp*kmax*3+1,ngp*kmax*4 ! q
      do e=1,edim
        if (xa(i,e).le.0.d0 ) xa(i,e)=0.d0
      enddo
    enddo
    do i=ngp*kmax*4+ngp*5+1,ngp*kmax*4+ngp*6
      do e=1,edim
        if (xa(i,e).lt.0.d0 ) xa(i,e)=0.d0
        if (xa(i,e).gt.0.d0 ) xa(i,e)=1.d0
      enddo
    enddo

    return

    end subroutine check

  !>-----------------------------------------------------------------------
  !>
  subroutine report
    use comm
    implicit none

    real(r_size) :: rmse1(2), rmse2(2)
    real(r_size) :: xm(mdim)
    real(r_size) :: dx(mdim,edim)
    integer :: cnt(2)
    integer :: p

    call decompose_x(xf_save,xm,dx)
    cnt(:)=0
    rmse1(:)=0.d0
    do p=1,nobs
      if( ijobs_2(p).eq.9 )then
        cnt(1)=cnt(1)+1
        rmse1(1)=rmse1(1)+(xm(ngp*8+ijobs_1(p))-obs(p))**2
      elseif( ijobs_2(p).eq.88 )then
        cnt(2)=cnt(2)+1
        rmse1(2)=rmse1(2)+(xm(ngp*87+ijobs_1(p))-obs(p))**2
      endif
    enddo
    rmse1=sqrt(rmse1/cnt)

    call decompose_x(xa,xm,dx)
    cnt(:)=0
    rmse2(:)=0.d0
    do p=1,nobs
      if( ijobs_2(p).eq.9 )then
        cnt(1)=cnt(1)+1
        rmse2(1)=rmse2(1)+(xm(ngp*8+ijobs_1(p))-obs(p))**2
      elseif( ijobs_2(p).eq.88 )then
        cnt(2)=cnt(2)+1
        rmse2(2)=rmse2(2)+(xm(ngp*87+ijobs_1(p))-obs(p))**2
      endif
    enddo
    rmse2=sqrt(rmse2/cnt)

    write(*,*) '===================================================='
    write(*,*) 'RMSE against obs. (',nobs,')'
    write(*,*) '----------------------------------------------------'
    write(*,*) '            temp.(',cnt(1),')         sst(',cnt(2),') '
    write(*,*) '----------------------------------------------------'
    write(*,*) 'gues      |',rmse1(1),rmse1(2)
    write(*,*) 'anal      |',rmse2(1),rmse2(2)
    write(*,*) '===================================================='

    return

  end subroutine report

!>-----------------------------------------------------------------------

  end program
