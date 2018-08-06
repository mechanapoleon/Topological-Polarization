program matrix_diagonalization
implicit none
 CHARACTER :: JOBZ, UPLO
 INTEGER   :: N, LDA
 COMPLEX (kind=8), allocatable, dimension(:,:) :: A, CarPosE, ZX  !sotto
 COMPLEX (kind=8), allocatable, dimension (:) :: WORK
 INTEGER ::  LWORK, INFO
 REAL (kind=8), allocatable, dimension(:) :: RWORK
 REAL (kind=8), allocatable, dimension(:) :: W
 integer :: i !dummy index ciclo do
 real :: DELTA, t , l, dipole !parametri dell'hamiltoniana, lunghezza della cella primitiva, dipolo
 real (kind=8), allocatable, dimension(:) :: groundstate !groundstate dell'hamiltoniana

DELTA = 0.25
t = 1
l=2



do N=4,400,2



JOBZ = 'V'
UPLO = 'U' 
LDA  = N
LWORK = max(1,3*(N-1))

allocate(A(LDA,N))
allocate(W(N))
allocate(WORK(max(1,LWORK)))
allocate(RWORK(max(1,3*N-2)))
allocate(groundstate(N))
allocate(CarPosE(N,N))
allocate(ZX(N,N))
A=0
do i=1,N
   A(i,i)= ((-1)**i ) * Delta

   if ( i /= N ) then
      A(i+1,i)= -t
      A(i,i+1)= -t
   else
   end if

end do

write(unit=11,fmt=*) A

call zheev(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO)

groundstate=0                  !definiamo il groundstate su cui vado a calcolare il dipolo
do i=1,N
   groundstate(i)= A(i,1)
end do


!andiamo a definire gli operatori posizione del reticolo X e posizione degli elettroni

ZX=0   !operatore CARICAxPOSIZIONE dei cationi-anioni

do i=1,N,2
   ZX(i,i)= l*(i+0.25)
   ZX(i+1,i+1)=-l*(i+0.75)
end do

CarPosE = 0 !operatore CaricaXPosizione dei cationi-anioni

do i=1,N
   CarPosE(i,i)=-groundstate(i) * l /2
end do

dipole=0

do i=1,N
   dipole= dipole + CarPosE(i,i) + ZX(i,i)
end do

dipole = dipole/ (N * l)
!dipole = dipole / l


write(unit=20,fmt=*),i, dipole, mod(dipole,1.0)



deallocate(A)
deallocate(W)
deallocate(WORK)
deallocate(RWORK)
deallocate(groundstate)
deallocate(CarPosE)
deallocate(ZX)


end do






!print*, INFO
!print*, W
!do i=1,n
!write(unit=12,fmt=*) i, W(i)   !autovalori
!end do

end program 
