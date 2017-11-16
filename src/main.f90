program main

    use mObjective
    use mOptimize
    use mUtilities

    implicit none

    real*8 :: xin(2), delta(2)
    real*8, allocatable :: mId(:,:)

    xin(1) = -1.d0; xin(2) = 1.d0
    delta(1) = 0.5d0; delta(2) = 0.5d0
    mId = matId(2)

    !call unconstBox(xin,delta,fRosenbrock,ikmax=10000,itol=1.d-5)
    call unconstHookeJeeves2(xin,delta,fRosenbrock,ikmax=10000,itol=1.d-5)
    !call unconstMultidim(xin,mId,fRosenbrock,gRosenbrock,hRosenbrock,phiRosenbrock,imethod=4,iomega=0.5d0,ifac=.true.)
    !call unconstMultidim(xin,mId,fRosenbrock,gRosenbrock,hRosenbrock,phiRosenbrock,imethod=4,iomega=0.5d0)

endprogram
