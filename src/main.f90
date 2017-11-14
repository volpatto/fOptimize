program main

    use mObjective
    use mOptimize
    use mUtilities

    implicit none

    real*8 :: xin(2)
    real*8, allocatable :: mId(:,:)

    xin(1) = -1.2d0; xin(2) = 1.d0
    mId = matId(2)

    call unconstMultidim(xin,mId,fRosenbrock,gRosenbrock,hRosenbrock,phiRosenbrock,imethod=4,iomega=0.5d0,ifac=.true.)
    !call unconstMultidim(xin,mId,fRosenbrock,gRosenbrock,hRosenbrock,phiRosenbrock,imethod=4,iomega=0.5d0)

endprogram
