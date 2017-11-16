module mOptimize

    implicit none

    contains

        subroutine biss_1d(f,a,b,itol,ikmax,imode,id)

            implicit none

            real*8, intent(in) :: a,b
            real*8, intent(in), optional :: itol
            integer, intent(in), optional :: ikmax, imode, id
            interface
                function f(x)
                    real*8, intent(in) :: x
                    real*8 :: f
                endfunction
            endinterface
            real*8 :: tol = 1.d-5, lambda_k, mu_k, a_k, b_k, eps
            integer :: k, kmax = 500, mode = 1
            character(len=50) :: filename, fname, char_id

            if (present(itol)) tol = itol
            if (present(ikmax)) kmax = ikmax
            if (present(imode)) mode = imode
            if (present(id)) then
                ! Opening solution output file
                ! Making an id char variable
                write(char_id, '(i0)') id
                write(filename, '(a)') "bissdf"
                ! Concatenate strings filename, char_id and ".dat"
                fname = trim(filename)//trim(char_id)//".dat"
                ! Opening file with unit as id and file name as attributed in fname
                open(unit=id,file=fname)
            endif

            a_k = a; b_k = b

            k = 0
            do while(((b_k-a_k).ge.tol).and.(k.lt.kmax))

                eps = 0.1d0*(b_k - a_k)
                lambda_k = (a_k+b_k)/2.d0 - eps
                mu_k = (a_k+b_k)/2.d0 + eps

                if (mode.eq.1) then
                    if (f(lambda_k).lt.f(mu_k)) then
                        b_k = mu_k
                    else
                        a_k = lambda_k
                    endif
                else
                    if (-f(lambda_k).lt.-f(mu_k)) then
                        b_k = mu_k
                    else
                        a_k = lambda_k
                    endif
                endif
                k = k + 1
                if (present(id)) write(id,2000) k, (a_k+b_k)/2.d0, f((a_k+b_k)/2.d0)

            enddo

            if (k.ge.kmax) stop "biss_df nao convergiu"
            print*, "k =", k, "x_opt =", (a_k+b_k)/2.d0, "f(x_opt) =", f((a_k+b_k)/2.d0)

            if (present(id)) close(id)

2000        format((i0),2(4x,(es20.8)))

        endsubroutine

        subroutine biss_phi(f,x_k,d_k,a,b,lamb,itol,ikmax,imode,id)

            implicit none

            real*8, dimension(:), intent(in) :: x_k, d_k
            real*8, intent(in) :: a,b
            real*8, intent(out) :: lamb
            real*8, intent(in), optional :: itol
            integer, intent(in), optional :: ikmax, imode, id
            interface
                function f(x)
                    real*8, dimension(:), intent(in) :: x
                    real*8 :: f
                endfunction
            endinterface
            !real*, allocatable :: xnext
            real*8 :: tol = 1.d-5, lambda_k, mu_k, a_k, b_k, eps
            integer :: k, kmax = 500, mode = 1
            character(len=50) :: filename, fname, char_id

            if (present(itol)) tol = itol
            if (present(ikmax)) kmax = ikmax
            if (present(imode)) mode = imode
            if (present(id)) then
                ! Opening solution output file
                ! Making an id char variable
                write(char_id, '(i0)') id
                write(filename, '(a)') "bissdf"
                ! Concatenate strings filename, char_id and ".dat"
                fname = trim(filename)//trim(char_id)//".dat"
                ! Opening file with unit as id and file name as attributed in fname
                open(unit=id,file=fname)
            endif

            !allocate(xnext(size(x))); xnext = 0.d0

            a_k = a; b_k = b

            k = 0
            do while(((b_k-a_k).ge.tol).and.(k.lt.kmax))

                eps = 0.1d0*(b_k - a_k)
                lambda_k = (a_k+b_k)/2.d0 - eps
                mu_k = (a_k+b_k)/2.d0 + eps

                if (mode.eq.1) then
                    if (f(x_k+lambda_k*d_k).lt.f(x_k+mu_k*d_k)) then
                        b_k = mu_k
                    else
                        a_k = lambda_k
                    endif
                else
                    if (-f(x_k+lambda_k*d_k).lt.-f(x_k+mu_k*d_k)) then
                        b_k = mu_k
                    else
                        a_k = lambda_k
                    endif
                endif
                k = k + 1
                if (present(id)) write(id,2000) k, (a_k+b_k)/2.d0, f(x_k+(a_k+b_k)/2.d0*d_k)

            enddo

            if (k.ge.kmax) stop "biss_phi nao convergiu"
            !print*, "k =", k, "x_opt =", (a_k+b_k)/2.d0, "f(x_opt) =", f(x_k+(a_k+b_k)/2.d0*d_k)
            lamb = (a_k+b_k)/2.d0

            if (present(id)) close(id)

2000        format((i0),2(4x,(es20.8)))

        endsubroutine

        subroutine golden_1d(f,a,b,itol,ikmax,imode,id)

            implicit none

            real*8, intent(in) :: a,b
            real*8, intent(in), optional :: itol
            integer, intent(in), optional :: ikmax, imode, id
            interface
                function f(x)
                    real*8, intent(in) :: x
                    real*8 :: f
                endfunction
            endinterface
            real*8 :: tol = 1.d-5, lambda_k, mu_k, a_k, b_k, eps
            integer :: k, kmax = 500, mode = 1
            character(len=50) :: filename, fname, char_id

            if (present(itol)) tol = itol
            if (present(ikmax)) kmax = ikmax
            if (present(imode)) mode = imode
            if (present(id)) then
                ! Opening solution output file
                ! Making an id char variable
                write(char_id, '(i0)') id
                write(filename, '(a)') "goldendf"
                ! Concatenate strings filename, char_id and ".dat"
                fname = trim(filename)//trim(char_id)//".dat"
                ! Opening file with unit as id and file name as attributed in fname
                open(unit=id,file=fname)
            endif

            a_k = a; b_k = b; eps = 6.18d-1
            lambda_k = a_k + (1.d0-eps)*(b_k-a_k)
            mu_k = a_k + eps*(b_k-a_k)

            k = 0
            do while(((b_k-a_k).ge.tol).and.(k.lt.kmax))
                if (mode.eq.1) then
                    if (f(lambda_k).lt.f(mu_k)) then
                        b_k = mu_k; mu_k = lambda_k
                        lambda_k = a_k + (1.d0-eps)*(b_k-a_k)
                    else
                        a_k = lambda_k; lambda_k = mu_k
                        mu_k = a_k + eps*(b_k-a_k)
                    endif
                else
                    if (-f(lambda_k).lt.-f(mu_k)) then
                        b_k = mu_k; mu_k = lambda_k
                        lambda_k = a_k + (1.d0-eps)*(b_k-a_k)
                    else
                        a_k = lambda_k; lambda_k = mu_k
                        mu_k = a_k + eps*(b_k-a_k)
                    endif
                endif
                k = k + 1
                write(id,2000) k, (a_k+b_k)/2.d0, f((a_k+b_k)/2.d0)

            enddo

            if (k.ge.kmax) stop "golden_df nao convergiu"
            print*, "k =", k, "x_opt =", (a_k+b_k)/2.d0, "f(x_opt) =", f((a_k+b_k)/2.d0)

            if (present(id)) close(id)

2000        format((i0),2(4x,(es20.8)))

        endsubroutine

        subroutine newton_1d(f0,f,fprime,x0,itol,ikmax,id)

            implicit none

            real*8, intent(in) :: x0
            real*8, intent(in), optional :: itol
            integer, intent(in), optional :: ikmax, id
            interface
                function f0(x)
                    real*8, intent(in) :: x
                    real*8 :: f0
                endfunction
                function f(x)
                    real*8, intent(in) :: x
                    real*8 :: f
                endfunction
                function fprime(x)
                    real*8, intent(in) :: x
                    real*8 :: fprime
                endfunction
            endinterface
            real*8 :: tol = 1.d-5, x_k, x_prev
            integer :: k, kmax = 500
            character(len=50) :: filename, fname, char_id

            if (present(itol)) tol = itol
            if (present(ikmax)) kmax = ikmax
            if (present(id)) then
                ! Opening solution output file
                ! Making an id char variable
                write(char_id, '(i0)') id
                write(filename, '(a)') "newton"
                ! Concatenate strings filename, char_id and ".dat"
                fname = trim(filename)//trim(char_id)//".dat"
                ! Opening file with unit as id and file name as attributed in fname
                open(unit=id,file=fname)
            endif

            x_k = x0
            x_prev = x0
            k = 0
            do while ((dabs(f(x_k)).ge.tol).and.(k.lt.kmax))
                if (fprime(x_prev).ne.0.d0) then
                    x_k = x_prev - f(x_prev)/fprime(x_prev)
                    x_prev = x_k 
                else
                    print*, "Divisao por 0 em newton"; stop
                endif
                k = k + 1
                write(id,2000) k, x_k, f0(x_k), f(x_k)
            enddo

            if (k.ge.kmax) stop "Newton nao convergiu"
            print*, "k =", k, "x_opt =", x_k, "f'(x_opt) =", f(x_k)

            if (present(id)) close(id)

2000        format((i0),3(4x,(es20.8)))

        endsubroutine

        subroutine secant_1d(f,fprime,x0,x1,itol,ikmax,id)

            implicit none

            real*8, intent(in) :: x0, x1
            real*8, intent(in), optional :: itol
            integer, intent(in), optional :: ikmax, id
            interface
                function f(x)
                    real*8, intent(in) :: x
                    real*8 :: f
                endfunction
                function fprime(x)
                    real*8, intent(in) :: x
                    real*8 :: fprime
                endfunction
            endinterface
            real*8 :: tol = 1.d-5, x_k, x_prev1, x_prev2
            integer :: k, kmax = 500, mode = 1
            character(len=50) :: filename, fname, char_id

            if (present(itol)) tol = itol
            if (present(ikmax)) kmax = ikmax
            if (present(id)) then
                ! Opening solution output file
                ! Making an id char variable
                write(char_id, '(i0)') id
                write(filename, '(a)') "secant"
                ! Concatenate strings filename, char_id and ".dat"
                fname = trim(filename)//trim(char_id)//".dat"
                ! Opening file with unit as id and file name as attributed in fname
                open(unit=id,file=fname)
            endif

            x_k = x0
            x_prev1 = x0; x_prev2 = x1
            x_k = (fprime(x_prev1)*x_prev2 - fprime(x_prev2)*x_prev1)/(fprime(x_prev1)-fprime(x_prev2))
            x_prev2 = x_prev1; x_prev1 = x_k; k = 1
            write(id,2000) k, x_k, f(x_k), fprime(x_k)
            do while ((dabs(fprime(x_k)).ge.tol).and.(k.lt.kmax))
                if ((fprime(x_prev1)-fprime(x_prev2)).ne.0.d0) then
                    x_k = (fprime(x_prev1)*x_prev2 - fprime(x_prev2)*x_prev1)/(fprime(x_prev1)-fprime(x_prev2))
                    x_prev2 = x_prev1; x_prev1 = x_k
                else
                    print*, "Divisao por 0 em secante"; stop
                endif
                k = k + 1
                write(id,2000) k, x_k, f(x_k), fprime(x_k)
            enddo

            if (k.ge.kmax) stop "Secante nao convergiu"
            print*, "k =", k, "x_opt =", x_k, "f'(x_opt) =", fprime(x_k)

            if (present(id)) close(id)

2000        format((i0),3(4x,(es20.8)))

        endsubroutine

        subroutine unconstMultidim(x0,B0,f,gf,hf,phif,imethod,iomega,ieps,ifac)

            use mUtilities
            use msolver

            implicit none

            real*8, dimension(:), intent(in) :: x0
            real*8, dimension(:,:), intent(in) :: B0
            integer, intent(in), optional :: imethod
            real*8, intent(in), optional :: iomega, ieps
            logical, intent(in), optional :: ifac
            interface
                function f(x)
                    real*8, dimension(:), intent(in) :: x
                    real*8 :: f
                endfunction
                function gf(x)
                    real*8, dimension(:), intent(in) :: x
                    real*8, allocatable :: gf(:)
                endfunction
                function hf(x)
                    real*8, dimension(:), intent(in) :: x
                    real*8, allocatable :: hf(:,:)
                endfunction
                function phif(x,d,lamb)
                    real*8, dimension(:), intent(in) :: x, d
                    real*8, intent(in) :: lamb
                    real*8 :: phif
                endfunction
            endinterface
            real*8, allocatable :: d_k(:), xprev(:), xnext(:), Bmat(:,:), p_k(:), q_k(:), Id(:,:)
            real*8, allocatable :: Bdfp(:,:), Bbfgs(:,:), Hnewton(:,:), gnewton(:)
            real*8 :: lambda_k, eps, rho_bfgs, omega, fac_scale
            real*8, parameter :: llambda = 0.d0, rlambda=2.d1
            integer :: kmax = 1000, k, method

            ! Standard optional input values
            eps = 1.d-5; method = 0; omega = 0.5d0; fac_scale = 1.d0

            if (present(imethod)) method = imethod
            if (present(iomega)) omega = iomega
            if (present(ieps)) eps = ieps

            if (present(iomega).and.present(imethod).and.(method.eq.4).and.(omega.gt.1.d0 .or. omega.lt.0.d0)) then
                print*, "Invalid Broyden combination: Combination must be convex."; stop
            endif
        
            allocate(d_k(size(x0))); allocate(xprev(size(x0))); xprev = x0
            allocate(xnext(size(x0))); allocate(Bmat(size(x0),size(x0)))
            allocate(p_k(size(x0))); allocate(q_k(size(x0)))
            if (method.eq.1) then 
                allocate(Hnewton(size(x0),size(x0))); allocate(gnewton(size(x0))) 
            endif
            if (method .eq. 4) then
                allocate(Bdfp(size(x0),size(x0))); allocate(Bbfgs(size(x0),size(x0)))
            endif
            Id = matId(size(x0));
            Bmat = B0
            k = 0

            ! Newton
            if (method.eq.1) then
            do while((norm2(gf(xprev)).gt.eps).and.(k.lt.kmax))
                gnewton = -gf(xprev); Hnewton = hf(xprev)
                call solverSOR(Hnewton,gnewton,d_k,itol=1.d-7,ikmax=5000) 
                call biss_phi(f,xprev,d_k,llambda,rlambda,lambda_k)
                xnext = xprev + lambda_k*d_k
                xprev = xnext
                k = k + 1
            enddo
            ! DFP
            else if (method.eq.2) then
            do while((norm2(gf(xprev)).gt.eps).and.(k.lt.kmax))
                d_k = -matmul(Bmat,gf(xprev))
                call biss_phi(f,xprev,d_k,llambda,rlambda,lambda_k)
                xnext = xprev + lambda_k*d_k
                p_k = xnext-xprev; 
                q_k = gf(xnext)-gf(xprev)
                if (present(ifac).and.(ifac.eqv..True.)) fac_scale = dot_product(p_k,q_k)/dot_product(q_k,matmul(Bmat,q_k))
                Bmat = fac_scale*Bmat + (1.d0/dot_product(p_k,q_k))*dyadic_prod(p_k,p_k) &
                    - fac_scale*matmul(matmul(Bmat,dyadic_prod(q_k,q_k)),Bmat)*(1.d0/dot_product(q_k,matmul(Bmat,q_k)))
                xprev = xnext
                k = k + 1
            enddo
            ! BFGS
            else if (method.eq.3) then
            do while((norm2(gf(xprev)).gt.eps).and.(k.lt.kmax))
                d_k = -matmul(Bmat,gf(xprev))
                call biss_phi(f,xprev,d_k,llambda,rlambda,lambda_k)
                xnext = xprev + lambda_k*d_k
                p_k = xnext-xprev; 
                q_k = gf(xnext)-gf(xprev)
                rho_bfgs = (1.d0/dot_product(q_k,p_k))
                if (present(ifac).and.(ifac.eqv..True.)) fac_scale = dot_product(p_k,q_k)/dot_product(q_k,matmul(Bmat,q_k))
                Bmat = fac_scale*matmul((Id-rho_bfgs*dyadic_prod(p_k,q_k)),matmul(Bmat,(Id-rho_bfgs*dyadic_prod(q_k,p_k)))) &
                    + rho_bfgs*dyadic_prod(p_k,p_k)
                xprev = xnext
                k = k + 1
            enddo
            ! DFP + BFGS
            else if (method.eq.4) then
            do while((norm2(gf(xprev)).gt.eps).and.(k.lt.kmax))
                d_k = -matmul(Bmat,gf(xprev))
                call biss_phi(f,xprev,d_k,llambda,rlambda,lambda_k)
                xnext = xprev + lambda_k*d_k
                p_k = xnext-xprev; 
                q_k = gf(xnext)-gf(xprev)
                if (present(ifac).and.(ifac.eqv..True.)) fac_scale = dot_product(p_k,q_k)/dot_product(q_k,matmul(Bmat,q_k))
                Bdfp = fac_scale*Bmat + (1.d0/dot_product(p_k,q_k))*dyadic_prod(p_k,p_k) &
                    - fac_scale*matmul(matmul(Bmat,dyadic_prod(q_k,q_k)),Bmat)*(1.d0/dot_product(q_k,matmul(Bmat,q_k)))
                rho_bfgs = (1.d0/dot_product(q_k,p_k))
                Bbfgs = fac_scale*matmul((Id-rho_bfgs*dyadic_prod(p_k,q_k)),matmul(Bmat,(Id-rho_bfgs*dyadic_prod(q_k,p_k)))) &
                    + rho_bfgs*dyadic_prod(p_k,p_k)
                Bmat = omega*Bbfgs + (1.d0-omega)*Bdfp
                xprev = xnext
                k = k + 1
            enddo
            ! Gradient
            else
            do while((norm2(gf(xprev)).gt.eps).and.(k.lt.kmax))
                d_k = -matmul(Bmat,gf(xprev))
                call biss_phi(f,xprev,d_k,llambda,rlambda,lambda_k)
                xnext = xprev + lambda_k*d_k
                xprev = xnext
                k = k + 1
            enddo
            endif

            print*, "****** End unconstrained optimization procedure ******"
            print*, "Iterations =", k 
            print*, "x_min =", xprev 
            print*, "f(x_min) =", f(xprev)

        endsubroutine

        subroutine unconstBox(x0,delta,f,itol,ifac,ikmax)

            implicit none

            real*8, dimension(:), intent(in) :: x0
            real*8, dimension(:), intent(inout) :: delta
            real*8, intent(in), optional :: itol, ifac
            integer, intent(in), optional :: ikmax
            interface
                function f(x)
                    real*8, dimension(:), intent(in) :: x
                    real*8 :: f
                endfunction
            endinterface
            real*8 :: eps, f_min, f_temp, fac
            real*8, allocatable :: x_k(:), x_temp(:), x_next(:)
            integer :: i, j, k, kmax

            eps = 1.d-5; fac = 0.8d0; kmax = 5000

            if (present(itol)) eps = itol
            if (present(ifac)) fac = ifac 
            if (present(ikmax)) kmax = ikmax

            if (size(x0).ne.size(delta)) stop "Error: Reductor vector must be same size than x0"
            if (minval(delta).lt.0.d0) stop "Error: All reductor values must be positive"
            if (norm2(delta).lt.eps) stop "Error: Provided reductor is too small"
            if (fac.ge.1.d0) stop "Error: Contraction factor must be lower than 1.0"

            allocate(x_k(size(x0))); x_k = x0
            allocate(x_temp(size(x0))); allocate(x_next(size(x0)));

            k = 0
            f_min = f(x_k)
            do while ((norm2(delta).gt.eps).and.(k.le.kmax))
                do i=1,size(delta)
                    do j=1,2
                        x_temp = x_k
                        x_temp(i) = x_k(i) + ((-1.d0)**j)*fac*delta(i)
                        f_temp = f(x_temp)
                        if (f_temp .lt. f_min) then 
                            f_min = f_temp
                            x_next = x_temp
                        endif
                    enddo
                enddo
                if (all(x_k.eq.x_next)) then
                    delta = fac*delta
                else
                    x_k = x_next
                endif
                k = k + 1
            enddo

            print*, "****** End unconstrained Box method ******"
            print*, "Iterations =", k 
            print*, "x_min =", x_k 
            print*, "f(x_min) =", f(x_k)

        endsubroutine

        subroutine exploreHJ(x,delta,f,flag)

            implicit none

            real*8, dimension(:), intent(inout) :: x, delta
            logical, intent(inout) :: flag
            interface
                function f(x)
                    real*8, dimension(:), intent(in) :: x
                    real*8 :: f
                endfunction
            endinterface
            integer :: i, j
            real*8 :: eps, f_min, f_temp
            real*8, allocatable :: x_k(:), x_temp(:)

            flag = .false.

            allocate(x_k(size(x))); x_k = x
            allocate(x_temp(size(x)))
            do i=1,size(x)
                f_min = f(x_k)
                do j=1,2
                    x_temp = x_k
                    x_temp(i) = x_k(i) + ((-1.d0)**j)*delta(i)
                    f_temp = f(x_temp)
                    if (f_temp .lt. f_min) then 
                        f_min = f_temp
                        x_k = x_temp
                    endif
                enddo
            enddo
            if (any(x.ne.x_k)) then 
                flag = .true.
                x = x_k
            endif

        endsubroutine

        subroutine unconstHookeJeeves(x0,delta,f,itol,ifac,ikmax)

            implicit none

            real*8, dimension(:), intent(in) :: x0
            real*8, dimension(:), intent(inout) :: delta
            real*8, intent(in), optional :: itol, ifac
            integer, intent(in), optional :: ikmax
            interface
                function f(x)
                    real*8, dimension(:), intent(in) :: x
                    real*8 :: f
                endfunction
            endinterface
            real*8 :: eps, fac
            real*8, allocatable :: x_k(:), x_prev(:), x_next(:)
            integer :: k, kmax
            logical :: flag = .true.

            eps = 1.d-5; fac = 0.5d0; kmax = 5000

            if (present(itol)) eps = itol
            if (present(ifac)) fac = ifac 
            if (present(ikmax)) kmax = ikmax

            if (size(x0).ne.size(delta)) stop "Error: Reductor vector must be same size than x0"
            if (minval(delta).lt.0.d0) stop "Error: All reductor values must be positive"
            if (norm2(delta).lt.eps) stop "Error: Provided reductor is too small"
            if (fac.ge.1.d0) stop "Error: Contraction factor must be lower than 1.0"

            allocate(x_k(size(x0))); x_k = x0
            allocate(x_prev(size(x0))); allocate(x_next(size(x0)));

            k = 0
            200 x_prev = x_k
            call exploreHJ(x_k,delta,f,flag)
            if (flag .eqv. .true.) then 
                goto 400 
            else
                goto 300
            endif
            300 if (norm2(delta).lt.eps) then
                    goto 700
                else
                    delta = fac*delta
                    goto 200
                endif
            400 k = k + 1
            if (k.gt.kmax) goto 700
            x_next = x_k + (x_k - x_prev)
            call exploreHJ(x_next,delta,f,flag)
            if (f(x_next).lt.f(x_k)) then
                x_prev = x_k
                x_k = x_next
                goto 400
            else
                goto 300
            endif

700         print*, "****** End unconstrained Hooke and Jeeves method ******"
            print*, "Iterations =", k 
            print*, "x_min =", x_k 
            print*, "f(x_min) =", f(x_k)

        endsubroutine
endmodule
