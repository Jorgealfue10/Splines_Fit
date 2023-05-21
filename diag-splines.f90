program diagsp

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Program to perform the splines fitting and the posterior diagonalization of the SO couplings.!
  ! Reads all couplings automatically (with format s??-??.dat) based on the number of SO states. !
  ! The outputs are: splines of the couplings in "splines-coupling" directory and the PECs       !
  ! after diagonalization. NO FIELD NO THETA APPLIED.                                            !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! To do: add subroutine for checking the continuity of 2 couplings and change it if needed.    !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  implicit none
  integer,parameter :: ns=34,nr=500,nrr=56
  integer :: i,j,k,l
  integer :: counti,countj,deci,decj,uni
  integer :: ntheta,nfield,aux
  character(len=10) :: input
  character(len=27) :: outsp
  character(len=1)  :: s11,s12,s21,s22
  real(kind=8),dimension(nrr,ns,ns) :: matr,auxma
  complex(kind=8),dimension(nrr,ns,ns) :: matdiag
  real(kind=8),dimension(nrr) :: dists,b,c,d
  real(kind=8) :: theta,field,dist,v1,v2,eval,r,dr
  real(kind=8) :: diff
  integer,dimension(nrr) :: wut
  real(8),dimension(nrr-1) :: der

  complex(kind=8),dimension(nr,ns,ns) :: matcomp
  real(kind=8),dimension(nr,ns,ns) :: u
  real(kind=8),dimension(nr) :: t

  real(kind=8),dimension(:) :: w(ns),rwork(3*ns-2),wr(nr,ns)
  complex(kind=8),dimension(:) :: work(2*ns-1)
  integer :: info
  complex(kind=8),dimension(nr,ns,ns) :: dia
  real(kind=8),dimension(nr,ns) :: adia

  real(kind=8),external :: ispline

  counti=0
  deci=0
  uni=15
  do i=1,ns
    counti=counti + 1
    decj=0
    countj=0
    do j=1,ns
      countj=countj+1  
      if (counti.eq.10) then
        deci=deci+1
        s11=char(48+deci)
        s12=char(48)
        counti=0
      else
        s11=char(48+deci)
        s12=char(48+counti)
      endif

      if (countj.eq.10) then
        decj=decj+1
        s21=char(48+decj)
        s22=char(48)
        countj=0
      else
        s21=char(48+decj)
        s22=char(48+countj)
      endif
      input="s"//s11//s12//"-"//s21//s22//".dat"
      open(uni,file=input,action="read")      

      do k=1,nrr
        read(uni,*) dists(k),matr(k,i,j)
      enddo

      close(uni) 
     
      uni=uni+1
!      write(*,*) input
    enddo
  enddo

  auxma=matr

!  do i=1,ns
!    do j=1,ns
!      l=1
!      do k=1,nrr
!        if (k.ne.nrr) then
!          diff=abs(matr(k,i,j)-matr(k+1,i,j))
!          !der(k)=abs(abs((matr(k+1,i,j)-matr(k,i,j))/(dists(k+1)-dists(k)))-abs((matr(k,i,j)-matr(k-1,i,j))/(dists(k)-dists(k-1))))
!          if ((diff.gt.4.d0)) then
!            call reordering(matr,dists,auxma,ns,nrr,i,j,k)!,wut(l))
!          endif
!        endif
!      enddo
!    enddo
!  enddo

!  counti=0
!  deci=0
!  uni=15
!  do i=1,ns
!    counti=counti + 1
!    decj=0
!    countj=0
!    do j=1,ns
!      countj=countj+1  
!      if (counti.eq.10) then
!        deci=deci+1
!        s11=char(48+deci)
!        s12=char(48)
!        counti=0
!      else
!        s11=char(48+deci)
!        s12=char(48+counti)
!      endif

!      if (countj.eq.10) then
!        decj=decj+1
!        s21=char(48+decj)
!        s22=char(48)
!        countj=0
!      else
!        s21=char(48+decj)
!        s22=char(48+countj)
!      endif
!      outsp="splines-reorders/s"//s11//s12//"-"//s21//s22//".dat"
!      open(uni,file=outsp,action="write")   

!      do k=2,nrr-1
!        !der(k)=abs(abs((matr(k+1,i,j)-matr(k,i,j))/(dists(k+1)-dists(k)))-abs((matr(k,i,j)-matr(k-1,i,j))/(dists(k)-dists(k-1))))
!        write(uni,*) dists(k),i,j,matr(k,i,j),der(k)
!      enddo

!      close(uni)
!    enddo
!  enddo

  dr=(dists(nrr)-dists(1))/(nr)

  t=0.d0
  u=0.d0
  counti=0
  deci=0
  uni=15
  do i=1,ns
    counti=counti + 1
    decj=0
    countj=0
    do j=1,ns
      countj=countj+1  
      if (counti.eq.10) then
        deci=deci+1
        s11=char(48+deci)
        s12=char(48)
        counti=0
      else
        s11=char(48+deci)
        s12=char(48+counti)
      endif

      if (countj.eq.10) then
        decj=decj+1
        s21=char(48+decj)
        s22=char(48)
        countj=0
      else
        s21=char(48+decj)
        s22=char(48+countj)
      endif
      outsp="splines-coupling/s"//s11//s12//"-"//s21//s22//".dat"
      open(uni,file=outsp,action="write")      

      call spline(dists(:),matr(:,i,j),b,c,d,nrr)
      r=dists(1)
      do k=1,nr
        u(k,i,j)=ispline(r,dists,matr(:,i,j),b,c,d,nrr)
        t(k)=r
        write(uni,*) r,u(k,i,j)
        r=r+dr
      enddo
      close(uni)
    enddo
  enddo


  do k=1,nr
    do i=1,ns
      do j=1,ns
        if (u(k,i,j).ne.0.d0) then
          v1=u(k,i,j)
          v2=u(k,j,i)
          if ((v1+v2).eq.0.d0) then
            matcomp(k,i,j)=cmplx(0.d0,v1,8)
            matcomp(k,j,i)=cmplx(0.d0,v2,8)
          else
            matcomp(k,i,j)=cmplx(v1,0.d0,8)
            matcomp(k,j,i)=cmplx(v2,0.d0,8)
          endif
        else
          matcomp(k,i,j)=cmplx(0.d0,0.d0,8)
          matcomp(k,j,i)=cmplx(0.d0,0.d0,8)
        endif
      enddo
    enddo
  enddo

  open(15,file="PECs.out",action="write")

  do k=1,nr
    call zheev('V','U',ns,matcomp(k,:,:),ns,w,work,2*ns-1,rwork,info)
    wr(k,:)=w(:)
    if (info.ne.0) print*,"nope",k
  enddo

  do k=1,nr
    write(15,*) t(k),wr(k,:)-minval(wr(:,1))
  enddo

  close(15)

!  t=0.d0
!  u=0.d0
!  counti=0
!  deci=0
!  uni=15
!  do i=1,ns
!    counti=counti + 1
!    decj=0
!    countj=0
!    do j=1,ns
!      countj=countj+1  
!      if (counti.eq.10) then
!        deci=deci+1
!        s11=char(48+deci)
!        s12=char(48)
!        counti=0
!      else
!        s11=char(48+deci)
!        s12=char(48+counti)
!      endif

!      if (countj.eq.10) then
!        decj=decj+1
!        s21=char(48+decj)
!        s22=char(48)
!        countj=0
!      else
!        s21=char(48+decj)
!        s22=char(48+countj)
!      endif
!      outsp="eigenvect/s"//s11//s12//"-"//s21//s22//".dat"
!      open(uni,file=outsp,action="write")      

!      do k=1,nr
!        write(uni,*) t(k),matcomp(k,i,j)
!      enddo

!      close(uni)
  
!    enddo
!  enddo

end program diagsp

subroutine spline (x, y, b, c, d, n)
  !======================================================================
  !  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
  !  for cubic spline interpolation
  !  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
  !  for  x(i) <= x <= x(i+1)
  !  Alex G: January 2010
  !----------------------------------------------------------------------
  !  input..
  !  x = the arrays of data abscissas (in strictly increasing order)
  !  y = the arrays of data ordinates
  !  n = size of the arrays xi() and yi() (n>=2)
  !  output..
  !  b, c, d  = arrays of spline coefficients
  !  comments ...
  !  spline.f90 program is based on fortran version of program spline.f
  !  the accompanying function fspline can be used for interpolation
  !======================================================================
  implicit none
  integer n
  double precision x(n), y(n), b(n), c(n), d(n)
  integer i, j, gap
  double precision h
  
  gap = n-1
  ! check input
  if ( n < 2 ) return
  if ( n < 3 ) then
    b(1) = (y(2)-y(1))/(x(2)-x(1))   ! linear interpolation
    c(1) = 0.
    d(1) = 0.
    b(2) = b(1)
    c(2) = 0.
    d(2) = 0.
    return
  end if
  !
  ! step 1: preparation
  !
  d(1) = x(2) - x(1)
  c(2) = (y(2) - y(1))/d(1)
  do i = 2, gap
    d(i) = x(i+1) - x(i)
    b(i) = 2.0*(d(i-1) + d(i))
    c(i+1) = (y(i+1) - y(i))/d(i)
    c(i) = c(i+1) - c(i)
  end do
  !
  ! step 2: end conditions 
  !
  b(1) = -d(1)
  b(n) = -d(n-1)
  c(1) = 0.0
  c(n) = 0.0
  if(n /= 3) then
    c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
    c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
    c(1) = c(1)*d(1)**2/(x(4)-x(1))
    c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
  end if
  !
  ! step 3: forward elimination 
  !
  do i = 2, n
    h = d(i-1)/b(i-1)
    b(i) = b(i) - h*d(i-1)
    c(i) = c(i) - h*c(i-1)
  end do
  !
  ! step 4: back substitution
  !
  c(n) = c(n)/b(n)
  do j = 1, gap
    i = n-j
    c(i) = (c(i) - d(i)*c(i+1))/b(i)
  end do
  !
  ! step 5: compute spline coefficients
  !
  b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
  do i = 1, gap
    b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
    d(i) = (c(i+1) - c(i))/d(i)
    c(i) = 3.*c(i)
  end do
  c(n) = 3.0*c(n)
  d(n) = d(n-1)
  end subroutine spline
  
    function ispline(u, x, y, b, c, d, n)
  !======================================================================
  ! function ispline evaluates the cubic spline interpolation at point z
  ! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
  ! where  x(i) <= u <= x(i+1)
  !----------------------------------------------------------------------
  ! input..
  ! u       = the abscissa at which the spline is to be evaluated
  ! x, y    = the arrays of given data points
  ! b, c, d = arrays of spline coefficients computed by spline
  ! n       = the number of data points
  ! output:
  ! ispline = interpolated value at point u
  !=======================================================================
  implicit none
  double precision ispline
  integer n
  double precision  u, x(n), y(n), b(n), c(n), d(n)
  integer i, j, k
  double precision dx
  
  ! if u is ouside the x() interval take a boundary value (left or right)
  if(u <= x(1)) then
    ispline = y(1)
    return
  end if
  if(u >= x(n)) then
    ispline = y(n)
    return
  end if
  
  !*
  !  binary search for for i, such that x(i) <= u <= x(i+1)
  !*
  i = 1
  j = n+1
  do while (j > i+1)
    k = (i+j)/2
    if(u < x(k)) then
      j=k
      else
      i=k
     end if
  end do
  !*
  !  evaluate spline interpolation
  !*
  dx = u - x(i)
  ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
  end function ispline

  subroutine reordering(mat,dists,matref,ns,nrr,ni,nj,nk)
    implicit none
    integer :: i,j,pk,l
    integer,intent(in) :: ns,nrr,ni,nj,nk
    real(kind=8),dimension(nrr,ns,ns),intent(inout) :: mat
    real(kind=8),dimension(nrr),intent(in) :: dists
    real(kind=8),dimension(nrr,ns,ns),intent(out) :: matref
    real(8) :: diff
    real(kind=8),dimension(nrr) :: aux

    !print*,"ee"
    !do l=1,5
    do i=1,ns
      do j=1,ns
        diff=abs(mat(nk+1,ni,nj)-mat(nk,ni,nj))
        !diff=abs(abs((mat(nk+1,ni,nj)-mat(nk,ni,nj))/(dists(nk+1)-dists(nk)))-abs((mat(nk,i,j)&
        !-mat(nk-1,i,j))/(dists(nk)-dists(nk-1))))
        if (((diff.le.2.d0)))then!.and.(mat(nk,ni,nj).ne.0.d0)) then
          !--------------------------!
          pk=nk+1
          !--------------------------!
          mat(pk,ni,nj) = matref(pk,i,j)
          !--------------------------!
        !  diff=abs(abs((mat(pk+1,ni,nj)-mat(pk,ni,nj))/(dists(pk+1)-dists(pk)))-abs((mat(pk,i,j)&
        !  -mat(pk-1,i,j))/(dists(pk)-dists(pk-1))))
          !--------------------------!
          diff=abs(mat(nk+1,ni,nj)-mat(nk,ni,nj))
          !--------------------------!
          do while ((diff.lt.0.8d0).and.(pk.ne.nrr))
            !-------------------------!
            diff=abs(mat(nk+1,ni,nj)-mat(nk,ni,nj))
            !-------------------------!
        !    diff=abs(abs((mat(pk+1,ni,nj)-mat(pk,ni,nj))/(dists(pk+1)-dists(pk)))-abs((mat(pk,i,j)&
        !    -mat(pk-1,i,j))/(dists(pk)-dists(pk-1))))
            !-----------------------------!
            mat(pk+1,ni,nj) = matref(pk+1,i,j)
            !-----------------------------!
            pk=pk+1
          enddo
        else
          mat(nk,ni,nj)=mat(nk,ni,nj)
        endif
      enddo
    enddo
    !enddo
end subroutine reordering
