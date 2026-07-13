program matrix_diagonalization_bench
implicit none
! Confronta la diagonalizzazione densa (dsyev, usata finora) con quella
! specializzata per matrici tridiagonali (dstevd), sfruttando il fatto che
! l'Hamiltoniana della catena (on-site +-Delta, hopping -t primi vicini) e'
! gia' tridiagonale: costruirla come matrice piena e ridurla a tridiagonale
! dentro dsyev e' lavoro sprecato.

integer, parameter :: nsizes = 6
integer, parameter :: nreps = 3
integer :: sizes(nsizes) = (/100, 200, 400, 800, 1600, 3200/)
integer :: isize, N, i, LWORK, LIWORK, INFO, M, irep
real (kind=8), allocatable :: A(:,:), W(:), WORK(:)
real (kind=8), allocatable :: D(:), E(:), Z(:,:), Wr(:), Dr(:), Er(:), Zr(:,:)
integer, allocatable :: IWORK(:), ISUPPZ(:)
real (kind=8) :: DELTA, t, ABSTOL
real :: t0, t1, time_dense, time_tri, time_partial

DELTA = 0.25d0
t = 1.0d0
ABSTOL = 0.0d0

print *, "     N   dense_dsyev[s]   tridiag_dstevd[s]   occupied_dstevr[s]   speedup_tri   speedup_partial"

do isize = 1, nsizes
   N = sizes(isize)
   time_dense = huge(1.0)
   time_tri = huge(1.0)
   time_partial = huge(1.0)

   do irep = 1, nreps
      ! ---- metodo attuale: costruisci la matrice piena N x N e chiama dsyev ----
      allocate(A(N,N), W(N))
      A = 0.0d0
      do i=1,N
         A(i,i) = ((-1)**i) * DELTA
         if (i /= N) then
            A(i+1,i) = -t
            A(i,i+1) = -t
         end if
      end do
      LWORK = max(1,3*N-1)
      allocate(WORK(LWORK))
      call cpu_time(t0)
      call dsyev('V','U', N, A, N, W, WORK, LWORK, INFO)
      call cpu_time(t1)
      time_dense = min(time_dense, t1 - t0)
      deallocate(A, W, WORK)

      ! ---- metodo ottimizzato: matrice gia' in forma tridiagonale, dstevd ----
      allocate(D(N), E(max(1,N-1)), Z(N,N))
      do i=1,N
         D(i) = ((-1)**i) * DELTA
      end do
      do i=1,N-1
         E(i) = -t
      end do
      LWORK  = 1 + 4*N + N*N
      LIWORK = 3 + 5*N
      allocate(WORK(LWORK), IWORK(LIWORK))
      call cpu_time(t0)
      call dstevd('V', N, D, E, Z, N, WORK, LWORK, IWORK, LIWORK, INFO)
      call cpu_time(t1)
      time_tri = min(time_tri, t1 - t0)
      deallocate(D, E, Z, WORK, IWORK)

      ! ---- solo gli N/2 stati occupati (energia piu' negativa) con dstevr ----
      allocate(Dr(N), Er(max(1,N-1)), Wr(N/2), Zr(N,N/2), ISUPPZ(2*max(1,N/2)))
      do i=1,N
         Dr(i) = ((-1)**i) * DELTA
      end do
      do i=1,N-1
         Er(i) = -t
      end do
      LWORK  = 20*N
      LIWORK = 10*N
      allocate(WORK(LWORK), IWORK(LIWORK))
      call cpu_time(t0)
      call dstevr('V', 'I', N, Dr, Er, 0.0d0, 0.0d0, 1, N/2, ABSTOL, M, Wr, Zr, N, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO)
      call cpu_time(t1)
      time_partial = min(time_partial, t1 - t0)
      deallocate(Dr, Er, Wr, Zr, ISUPPZ, WORK, IWORK)
   end do

   print '(I6, F18.6, F18.6, F18.6, F14.2, F16.2)', N, time_dense, time_tri, time_partial, &
         time_dense/time_tri, time_dense/time_partial

end do

end program matrix_diagonalization_bench
