module mObjective

    implicit none

    contains

        function f1(x)

            implicit none

            real*8, intent(in) :: x
            real*8 :: f1

            f1 = 3.d0*x**4.d0 - 4.d0*x**3.d0 + 1.d0

        endfunction

        function df1(x)

            implicit none

            real*8, intent(in) :: x
            real*8 :: df1

            df1 = 12.0d0*x**3.d0 - 12.0d0*x**2.d0

        endfunction

        function ddf1(x)

            implicit none

            real*8, intent(in) :: x
            real*8 :: ddf1

            ddf1 = 36.0d0*x**2.d0 - 24.0d0*x

        endfunction

        function f2(x)

            implicit none

            real*8, intent(in) :: x
            real*8 :: f2

            f2 = (1.d0/4.d0)*x**4.d0 - (5.d0/3.d0)*x**3.d0 - 6.d0*x**2.d0 + 19.d0*x - 7.d0

        endfunction

        function df2(x)

            implicit none

            real*8, intent(in) :: x
            real*8 :: df2

            df2 = x**3.d0 - 5.d0*x**2.d0 - 12.d0*x + 19.d0

        endfunction

        function ddf2(x)

            implicit none

            real*8, intent(in) :: x
            real*8 :: ddf2

            ddf2 = 3.d0*x**2.d0 - 10.d0*x - 12.d0

        endfunction

        function fRosenbrock(x)

            implicit none

            real*8, dimension(:), intent(in) :: x
            real*8 :: fRosenbrock

            fRosenbrock = 100.d0*(x(2)-x(1)**2.d0)**2.d0 + (x(1)-1.d0)**2.d0

        endfunction

        function gRosenbrock(x)

            implicit none

            real*8, dimension(:), intent(in) :: x
            real*8, allocatable :: gRosenbrock(:)
            integer :: xdim

            xdim = size(x)
            allocate(gRosenbrock(xdim)); gRosenbrock = 0.d0

            gRosenbrock(1) = 400.d0*(x(1)**3.d0-x(1)*x(2)) + 2.d0*(x(1)-1.d0)
            gRosenbrock(2) = 200.d0*(x(2)-x(1)**2.d0)

        endfunction

        function hRosenbrock(x)

            implicit none

            real*8, dimension(:), intent(in) :: x
            real*8, allocatable :: hRosenbrock(:,:)
            integer :: xdim

            xdim = size(x)
            allocate(hRosenbrock(xdim,xdim)); !hRosenbrock = 0.d0

            hRosenbrock(1,1) = 400.d0*(3.d0*x(1)**2.d0-x(2))+2.d0
            hRosenbrock(1,2) = -400.d0*x(1)
            hRosenbrock(2,1) = hRosenbrock(1,2)
            hRosenbrock(2,2) = 200.d0

        endfunction

        function phiRosenbrock(x,d,lamb)

            implicit none

            real*8, dimension(:), intent(in) :: x, d
            real*8, intent(in) :: lamb
            real*8 :: phiRosenbrock

            phiRosenbrock = fRosenbrock(x+lamb*d)

        endfunction

        function dphiRosenbrock(x,d,lamb)

            implicit none

            real*8, dimension(:), intent(in) :: x, d
            real*8, intent(in) :: lamb
            real*8 :: dphiRosenbrock

            dphiRosenbrock = dot_product(gRosenbrock(x+lamb*d),d)

        endfunction

endmodule
