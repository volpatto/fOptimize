        !> Contains subroutine to compute numerical solution of linear
        !! systems Ax=b. 
        !!
        !! The present module has a general purpose such
        !! that the routines here intend to be independent of others
        !! module.
        !!
        !! @author Diego T. Volpatto
        module msolver

            implicit none

            contains
            
            !> Applies Gauss reduction in A(n,n) to obtain a superior
            !! triangular equivalent form.
            !! @param A     [in/out]A matrix A(n,n)
            !! @param n     [in]Number of rows/columns of matrix A
            !! @author Diego T. Volpatto
            subroutine tri(A, n)

                implicit none

                integer :: n
                real*8 :: A(n,n)

                integer :: n1, j1, i, j, k
                real*8 :: fac
                real*8, parameter :: tiny_val = 1.d-30

                n1 = n-1
                ! Eliminate degree of freedom i
                do i=1,n1

                ! Check for excessively small pivot
                if (dabs(A(i,i)).lt.tiny_val) then
                    print*, i, A(i,i)
                    print*, "Reduction failed due to small pivot"
                    stop
                endif

                j1 = i+1
                ! Modify row j
                do j=j1,n
                if (A(j,i).eq.0.0d0) cycle
                fac = A(j,i)/A(i,i)
                do k=j1,n
                A(j,k)=A(j,k)-A(i,k)*fac
                enddo
                enddo
                enddo

            endsubroutine

            !> Does the forward substitution on the right-side-hand.
            !! @param A     [in]A matrix A(n,n)
            !! @param x     [out]Solution vector
            !! @param b     [in/out]RHS-vector
            !! @param n     [in]Number of solution points
            !! @author Diego T. Volpatto
            subroutine rhsub(A, x, b, n)

                implicit none

                real*8 :: A(n,n), x(n), b(n)
                integer :: n
                integer :: n1, j1, i, j, ib

                n1 = n-1
                ! Begin forward reduction of right hand side
                do i=1,n1
                j1=i+1
                do j=j1,n
                b(j)=b(j)-b(i)*A(j,i)/A(i,i)
                enddo
                enddo
                ! Begin back substitution
                x(n)=b(n)/A(n,n)
                do i=1,n1
                ib=n-i
                j1=ib+1
                do j=j1,n
                b(ib)=b(ib)-A(ib,j)*x(j)
                enddo
                x(ib)=b(ib)/A(ib,ib)
                enddo

            endsubroutine

            subroutine solverJacobi(A, b, x, itol, ikmax)

                implicit none

                real*8, intent(inout) :: A(:,:), b(:)
                real*8, allocatable, intent(out) :: x(:)
                real*8, intent(in), optional :: itol
                integer, optional :: ikmax
                real*8, allocatable :: xprev(:)
                real*8 :: tol = 1.d-10, summation
                integer :: i, j, k, n, kmax = 2000

                if (present(itol)) tol = itol
                if (present(ikmax)) kmax = ikmax

                if ((size(A(1,:)).ne.(size(A(:,1)))).or.(size(b).ne.size(A(1,:)))) then
                    print*, "Invalid system Ax = b: Incompatibility of dimensions."
                else
                    n = size(b)
                endif

                allocate(x(n)); x = 0.d0
                allocate(xprev(n)); xprev = b

                !write(*,'(/(a)/)') "****************** Begin solver Jacobi ******************"
                do k=1,kmax
                    print*, "Iteration =", k
                    do i=1,n
                        summation = 0.d0
                        do j=1,n
                            if (i.ne.j) then 
                                summation = summation + A(i,j)*xprev(j)
                            endif
                        enddo
                        x(i) = (b(i)-summation)/A(i,i)
                    enddo
                    !if (norm2(x-xprev)/norm2(x).le.tol) then
                    if (dabs(maxval(x-xprev))/dabs(maxval(x)).le.tol) then
                        exit
                    else
                        xprev = x
                    endif
                enddo

                write(*,'(1x,(a),2x,i0,2x,(a),2x,es10.2)') "Solver Jacobi iterations:",k,"tol =",tol
                !write(*,'(/(a)/)') "****************** End solver Jacobi ******************"

            endsubroutine

            subroutine solverGS(A, b, x, itol, ikmax)

                implicit none

                real*8, intent(inout) :: A(:,:), b(:)
                real*8, allocatable, intent(out) :: x(:)
                real*8, intent(in), optional :: itol
                integer, optional :: ikmax
                real*8, allocatable :: xprev(:)
                real*8 :: tol = 1.d-10, summation
                integer :: i, j, k, n, kmax = 2000

                if (present(itol)) tol = itol
                if (present(ikmax)) kmax = ikmax

                if ((size(A(1,:)).ne.(size(A(:,1)))).or.(size(b).ne.size(A(1,:)))) then
                    print*, "Invalid system Ax = b: Incompatibility of dimensions."
                else
                    n = size(b)
                endif

                allocate(xprev(n)); xprev = b
                allocate(x(n)); x = xprev

                !write(*,'(/(a)/)') "****************** Begin solver Gauss-Seidel ******************"
                do k=1,kmax
                    !print*, "Iteration =", k
                    do i=1,n
                        summation = 0.d0
                        do j=1,n
                            if (i.ne.j) then 
                                summation = summation + A(i,j)*x(j)
                            endif
                        enddo
                        x(i) = (b(i)-summation)/A(i,i)
                    enddo
                    !if (norm2(x-xprev)/norm2(x).le.tol) then
                    if (dabs(maxval(x-xprev))/dabs(maxval(x)).le.tol) then
                        exit
                    else
                        xprev = x
                    endif
                enddo

                write(*,'(1x,(a),2x,i0,2x,(a),2x,es10.2)') "Solver GS iterations:",k,"tol =",tol
                !write(*,'(/(a)/)') "****************** End solver Gauss-Seidel ******************"

            endsubroutine

            subroutine solverSOR(A, b, x, iomega, itol, ikmax)

                implicit none

                real*8, intent(inout) :: A(:,:), b(:)
                real*8, allocatable, intent(out) :: x(:)
                real*8, intent(in), optional :: itol, iomega
                integer, optional :: ikmax
                real*8, allocatable :: xprev(:)
                real*8 :: tol = 1.d-10, summation, omega=1.2d0
                integer :: i, j, k, n, kmax = 2000

                if (present(itol)) tol = itol
                if (present(ikmax)) kmax = ikmax
                if (present(iomega)) omega = iomega

                if ((size(A(1,:)).ne.(size(A(:,1)))).or.(size(b).ne.size(A(1,:)))) then
                    print*, "Invalid system Ax = b: Incompatibility of dimensions."
                else
                    n = size(b)
                endif

                allocate(xprev(n)); xprev = b
                allocate(x(n)); x = 0.d0

                !write(*,'(/(a)/)') "****************** Begin solver SOR ******************"
                do k=1,kmax
                    !print*, "Iteration =", k
                    do i=1,n
                        summation = 0.d0
                        do j=1,n
                            if (i.ne.j) then 
                                summation = summation + A(i,j)*x(j)
                            endif
                        enddo
                        x(i) = (b(i)-summation)/A(i,i)
                        x(i) = omega*x(i) + (1.d0-omega)*xprev(i)
                    enddo
                    !if (norm2(x-xprev)/norm2(x).le.tol) then
                    if (dabs(maxval(x-xprev))/dabs(maxval(x)).le.tol) then
                        exit
                    else
                        xprev = x
                    endif
                enddo
                write(*,'(1x,(a),2x,i0,2x,(a),2x,es10.2)') "Solver SOR iterations:",k,"tol =",tol

                !write(*,'(/(a)/)') "****************** End solver SOR ******************"

            endsubroutine
        endmodule
