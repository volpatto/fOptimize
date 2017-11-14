module mUtilities

    implicit none

    contains

        function matId(mdim)

            implicit none

            integer, intent(in) :: mdim
            real*8, allocatable :: matId(:,:)
            integer :: i

            allocate(matId(mdim,mdim)); matId = 0.d0

            do i=1,mdim
                matId(i,i) = 1.d0
            enddo

        endfunction

        function matvec(A,x)

            implicit none

            real*8, dimension(:,:), intent(in) :: A
            real*8, dimension(:), intent(in) :: x
            real*8, allocatable :: matvec(:)
            integer :: i, j

            if (size(A(1,:)).ne.size(x).or.size(A(:,1)).ne.size(x)) stop "Dimensional incompatibility"

            allocate(matvec(size(x))); matvec = 0.d0

            do i=1,size(x)
                do j=1,size(x)
                   matvec(i) = matvec(i) + A(i,j)*x(j) 
                enddo
            enddo

        endfunction

        function dyadic_prod(u,v)

            implicit none

            real*8, dimension(:), intent(in) :: u, v
            real*8, allocatable :: dyadic_prod(:,:)
            integer :: i, j

            if (size(u).ne.size(v)) stop "Dimensional incompatibility in dyadic product"
            allocate(dyadic_prod(size(u),size(u))); dyadic_prod = 0.d0

            do i=1,size(u)
                do j=1,size(u)
                    dyadic_prod(i,j) = u(i)*v(j)
                enddo
            enddo

        endfunction

endmodule
