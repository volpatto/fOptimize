program main

    use mObjective
    use mOptimize
    use mUtilities

    implicit none

    real*8 :: xin(2), delta(2)
    real*8, allocatable :: mId(:,:)

    mId = matId(2)

    xin(1) = 0.d0; xin(2) = 0.d0
    delta(1) = 1.d0; delta(2) = 1.0d0
    call unconstBox(xin,delta,fHimmelblau,ikmax=10000,itol=1.d-5,idfile=1)

    delta(1) = 0.2d0; delta(2) = 0.2d0
    call unconstHookeJeeves(xin,delta,fHimmelblau,idfile=1)

    xin(1) = -1.2d0; xin(2) = 1.d0
    call unconstMultidim(xin,mId,fRosenbrock,gRosenbrock,hRosenbrock,phiRosenbrock,imethod=4,iomega=0.5d0,ifac=.true.,idfile=2)

endprogram
