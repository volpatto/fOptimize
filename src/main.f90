program main

    use mObjective
    use mOptimize
    use mUtilities

    implicit none

    real*8 :: xin(2), delta(2)
    real*8, allocatable :: mId(:,:)

    xin(1) = -1.2d0; xin(2) = 1.d0
    delta(1) = 0.2d0; delta(2) = 0.2d0
    mId = matId(2)

    call unconstBox(xin,delta,fRosenbrock,ikmax=10000,itol=1.d-5,idfile=1)

    xin(1) = -1.2d0; xin(2) = 1.d0
    delta(1) = 0.2d0; delta(2) = 0.2d0

    call unconstHookeJeeves(xin,delta,fRosenbrock,idfile=1)
    call unconstMultidim(xin,mId,fRosenbrock,gRosenbrock,hRosenbrock,phiRosenbrock,imethod=4,iomega=0.5d0,ifac=.true.,idfile=1)

endprogram
