program matrix_diagonalization
implicit none
 CHARACTER :: JOBZ
 INTEGER   :: N, LDZ
 REAL (kind=8), allocatable, dimension(:,:) :: CarPosE, ZX, Z  !sotto
 REAL (kind=8), allocatable, dimension (:) :: WORK, D, E
 INTEGER, allocatable, dimension(:) :: IWORK
 INTEGER ::  LWORK, LIWORK, INFO
 integer :: i, k !dummy index cicli do
 real :: DELTA, t , l, dipole !parametri dell'hamiltoniana, lunghezza della cella primitiva, dipolo
 real (kind=8), allocatable, dimension(:) :: groundstate !groundstate dell'hamiltoniana

DELTA = 0.25
t = 1
l=2


N=100


JOBZ = 'V'
LDZ  = N
LWORK  = 1 + 4*N + N*N   ! dstevd, autovettori richiesti (JOBZ='V')
LIWORK = 3 + 5*N

allocate(D(N))
allocate(E(max(1,N-1)))
allocate(Z(LDZ,N))
allocate(WORK(max(1,LWORK)))
allocate(IWORK(max(1,LIWORK)))
allocate(groundstate(N))
allocate(CarPosE(N,N))
allocate(ZX(N,N))

! l'Hamiltoniana e' gia' tridiagonale (on-site +-Delta, hopping -t primi vicini):
! la passiamo direttamente a dstevd invece di costruire la matrice piena e usare
! dsyev, che dovrebbe prima ridurla a tridiagonale internamente (lavoro sprecato)
do i=1,N
   D(i) = ((-1)**i) * Delta
end do
do i=1,N-1
   E(i) = -t
end do

call dstevd(JOBZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO)

groundstate=0                  !densita' elettronica: somma su tutti gli N/2 stati occupati (energia negativa), non solo il piu' basso
do k=1,N/2
   do i=1,N
      groundstate(i)= groundstate(i) + Z(i,k)**2
   end do
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

print*, dipole


print*, INFO
!print*, D
do i=1,n
write(unit=12,fmt=*) i, D(i)   !autovalori
end do

end program
