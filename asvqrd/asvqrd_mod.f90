module asvqrd_mod
    use numbertypes
    use customla_mod
    implicit none

    type :: ASvQRDdata
        integer               :: N
        integer               :: l
        real(dp), allocatable :: work1(:, :)
        real(dp), allocatable :: work2(:, :)
    endtype ASvQRDdata


    contains

    
        subroutine ASvQRD(A, S, N, L, north, getj, make_B, left_Bmult, &
                          Q, P, tau, work, lwork, D, F, T, matwork, R)
            real(dp)        , intent(out)   :: A(N, N)
            type(ASvQRDdata), intent(inout) :: S
            integer         , intent(in)    :: N
            integer         , intent(in)    :: L
            integer         , intent(in)    :: north
            real(dp)        , intent(out)   :: Q(N, N)
            integer         , intent(out)   :: P(N, N)
            real(dp)        , intent(out)   :: tau(N)
            real(dp)        , intent(out)   :: work(lwork)
            integer         , intent(in)    :: lwork
            real(dp)        , intent(out)   :: D(N)
            real(dp)        , intent(out)   :: F(N)
            real(dp)        , intent(out)   :: T(N, N)
            real(dp)        , intent(out)   :: matwork(N, N)
            real(dp)        , intent(out)   :: R(N, N)
            interface
                function getj(j, S)
                    import                       :: ASvQRDdata
                    integer         , intent(in) :: j
                    type(ASvQRDdata), intent(in) :: S
                    integer                      :: getj
                endfunction getj
                subroutine make_B(A, j, N, S)
                    import                          :: ASvQRDdata, dp
                    real(dp)        , intent(out)   :: A(N, N)
                    integer         , intent(in)    :: j
                    integer         , intent(in)    :: N
                    type(ASvQRDdata), intent(inout) :: S
                endsubroutine make_B
                subroutine left_Bmult(A, j, N, S)
                    import                          :: ASvQRDdata, dp
                    real(dp)        , intent(out)   :: A(N, N)
                    integer         , intent(in)    :: j
                    integer         , intent(in)    :: N
                    type(ASvQRDdata), intent(inout) :: S
                endsubroutine left_Bmult
            endinterface

            integer :: j ! B matrix counter
            integer :: i ! north counter
            integer :: info

            ! Iteration j = 1
            j = 1

            ! Q = B(l)
            call make_B(Q, getj(j, S), N, S)

            ! north-1 more B multiplications before QRP factorisation
            do j = 2, north
                call left_Bmult(Q, getj(j, S), N, S)
            enddo

            ! QRP factorise Q
            P = 0
            call dgeqp3(N, N, Q, N, P, tau, work, lwork, info)

            ! D = diag(Q)
            call diag(Q, D, N)

            ! T = uppertri(Q)
            call uppertri(Q, T, N)
            ! T = inv(D) * T
            call left_diaginvmult(T, D, N)
            ! T = T * inv(P)
            call colpivswap(T, P, N, matwork)

            ! north B matrices have been multiplied and QRP factorised by now
            ! Multiply the rest of the B matrices
            j = north + 1
            do while (j .lt. L)
                ! Q = full Q (from previous iteration QRP factorisation)
                call dorgqr(N, N, N, Q, N, tau, work, lwork, info)

                ! Multiply north B matrices into Q (unless j reaches L along the way)
                i = 0
                do while ((i .lt. north) .and. (j .lt. L))
                    ! Q = B(getj(j)) * Q
                    call left_Bmult(Q, getj(j, S), N, S)
                    i = i + 1
                    j = j + 1
                enddo

                ! Q = Q * D (D from previous iteration)
                call right_diagmult(Q, D, N)

                ! QRP factorise Q
                P = 0
                call dgeqp3(N, N, Q, N, P, tau, work, lwork, info)

                ! D = diag(Q)
                call diag(Q, D, N)
                ! R = uppertri(Q)
                call uppertri(Q, R, N)

                ! R = inv(D) * R
                call left_diaginvmult(R, D, N)
                ! R = R * inv(P)
                call colpivswap(R, P, N, matwork)

                ! T = R * T
                call left_matmul(T, R, N, matwork) ! Should try to get rid of calling this subroutine
            enddo

            ! Now to finish calculating G

            ! D(L) = Db * Ds decomposition, see ASvQRD algorithm
            ! D = Db, F = Ds
            call DbDs(D, F, N)

            ! Q = full Q
            call dorgqr(N, N, N, Q, N, tau, work, lwork, info)

            ! matwork = trans(Q)
            call trans(matwork, Q, N)
            ! matwork = inv(Db) * B
            call left_diaginvmult(matwork, D, N)

            ! T = Ds * T                  
            call left_diagmult(T, F, N)    
            ! T = T + matwork                                              
            call add_matrix(T, matwork, N)
            ! T = inv(T)
            call invert(T, N, P, work, lwork, info)   

            ! A = T * matwork
            call dgemm('n', 'n', N, N, N, 1.0_dp, T, N, matwork, N, 0.0_dp, A, N)

        endsubroutine ASvQRD

        
        subroutine DbDs(D, F, n)
            !
            ! The ASvQRD algorithm asks for the last diagonal D matrix (stored as a vector)
            ! to be decomposed as follows:
            !
            !         | D(i) if abs(D(i)) > 1
            ! Db(i) = | 1    otherwise
            !
            !         | D(i) if abs(D(i)) <= 1
            ! Ds(i) = | 1    otherwise
            !
            ! This subroutine takes the last diagonal D matrix and sets it to Db,
            ! and takes another (unset) diagonal matrix F (stored as a vector) and sets it to Ds
            !
            real(dp), intent(inout) :: D(n)
            real(dp), intent(out)   :: F(n)
            integer , intent(in)    :: n

            integer :: i

            do i = 1, N
                if (abs(D(i)) .gt. 1) then
                    F(i) = 1
                else
                    F(i) = D(i)
                    D(i) = 1
                endif
            enddo


        endsubroutine DbDs


endmodule asvqrd_mod