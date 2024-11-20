# ASvQRD
This folder contains a generic implementation of the ASvQRD algorithm.
See the paper Stable solutions of linear systems involving long chain of matrix multiplications https://doi.org/10.1016/j.laa.2010.06.023.

The ASvQRD algorithm is the workhorse behind stably computing equal time Green's functions in DQMC.
This implementation does the following computation:

$$(\text{id} + B_LB_{L-1}\cdots B_1)^{-1}$$

For use specified matrices $B_l$ (see the usage section).

The code can easily be adapted to compute just $B_LB_{L-1}\cdots B_1$ or $(B_LB_{L-1}\cdots B_1)^{-1}$ in a stable manner if desired, and is planned to be added here in the future.

# Usage
The ASvQRD algorithm resides in the ```asvqrd_mod.f90``` module, and depends on the ```numbertypes.f90``` and ```custom_la.f90``` modules.
The subroutine that should be called is named ```ASvQRD```. It has a calling signature as follows:

```
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
```
Most arguments just hold workspaces, and should be input according to their required sizes.

The arguments of importance are as follows:

```A```: the ```NxN``` matrix in which the result will be stored in.

```S```: a datatype defined in the ```asvqrd_mod.f90``` module that the user should alter for their specific problem, see more later.

```N```: the dimension of ```A```.

```L```: how many matrices to multiply in the chain.

```north```: how many ```B_l``` matrices to multiply together before performing a QRP factorisation. A lower number is more stable, but slower.

```getj, make_B, left_Bmult``` are a function, subroutine, and subroutine respectfully, and all work together.
Typically, as is the case in DQMC simulations, the $B_l$ matrices aren't quite in such a linear order like $B_LB_{L-1}\cdots B_1$.
In DQMC, computing the $l$th equal time Green's function requires the product $B_l\cdots B_1B_L\cdot B_{l+1}$.
The ASvQRD algorithm iterates a variable ```j=1,...,L```  that counts what number matrix on the right is up for multiplication.
Given a ```j```, ```getj``` should return the index of what ```B_l``` matrix is up, which is input into ```left_Bmult``` (this subroutine's ```j``` argument) so the proper ordering of multiplications will take place.
The subroutine ```left_Bmult``` should implement left multiplication by the properly indexed $B$ matrix into its ```A``` matrix input.
The subroutine ```make_B``` is only called once, when ```j=1```, and should explicitly construct the very first $B_l$ matrix on the right.

In ```asvqrd_mod.f90``` near the top, before the ```contains``` statement, there is an ```ASvQRDdata``` data type definition.
In this repository this data type is filled with variables required by the example program ```test.f90```.
This data type should contain all extra data required to call ```getj, make_B, left_Bmult``` as illustrated in their interface.

An example program ```test.f90``` is provided illustrating calling the ASvQRD algorithm.
There are two important parts to keep in mind:

1. the ```ASvQRDdata``` data type has been appropriately defined in the ```asvqrd_mod.f90``` module
2. ```getj, make_B, left_Bmult``` have all been properly defined for this example problem (in this case, random $B_l$ matrices are multiplied together)
