program main
    use numbertypes
    use customla_mod
    use asvqrd_mod
    implicit none

    integer, parameter :: N = 5, L = 25, north = 5, lwork = N * N * N
    real(dp) :: A(N, N), Q(N, N), R(N, N), matwork(N, N), T(N, N), D(N), F(N), tau(N), work(lwork)
    integer :: P(N)
    type(ASvQRDdata) :: S

    S%N = L
    S%l = 12
    allocate(S%work1(N, N), S%work2(N, N))



    call ASvQRD(A, S, N, L, north, getj, make_B, left_Bmult, &
                Q, P, tau, work, lwork, D, F, T, matwork, R)

    call print_matrix(A)

    contains

        function getj(i, S) result(j)
            integer         , intent(in) :: i
            type(ASvQRDdata), intent(in) :: S

            integer :: j

            if (i .le. S%N - S%l) then
                j = S%l + i
            else
                j = S%l - S%L + i
            endif


        endfunction getj


        subroutine make_B(A, j, N, S)
            real(dp)        , intent(out)    :: A(N, N)
            integer         , intent(in)     :: j
            integer         , intent(in)     :: N
            type(ASvQRDdata), intent(inout)  :: S

            call random_number(A)

            
        endsubroutine make_B


        subroutine left_Bmult(A, j, N, S)
            real(dp)        , intent(out)    :: A(N, N)
            integer         , intent(in)     :: j
            integer         , intent(in)     :: N
            type(ASvQRDdata), intent(inout)  :: S

            call make_B(S%work1, j, N, S)
            call left_matmul(A, S%work1, N, S%work2)


        endsubroutine left_Bmult


        subroutine print_matrix(A)
            real(dp), intent(in) :: A(:, :)
            
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            do j = 1, n
                do i = 1, m
                    write(*, "(F12.6)", advance="no") A(i, j)
                enddo
                write(*, *) ""
            enddo


        endsubroutine print_matrix

endprogram main