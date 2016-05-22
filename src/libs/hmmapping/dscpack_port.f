C Initialization of mapping. Returns conformal module.
        real*8 function dscpack_init(
C Input
     &           n1, n2,       !number of outer/inner polygon points
     &           z1, z2,       !array of cmplx[n1/n2] outer/inner polygon pts
                               !in anti-clockwise. z1(0) is close to z2(0)
     &           prec,         !precision of integrating (int > 3)
C Output. Arrays should be allocated
     &           alfa1, alfa2, !array of doubles [N1/N2] - angles
     &           w1, w2,       !array of dcmplx [N1/N2] - canonic coords
     &           phi1, phi2,   !double array[N1/N2] - args of the prevertices
     &           u, c,         !double, dcmplx
     &           qwork         !double array[prec*(2(N1+N2)+3)]
     &  )
C input data declaration
        integer n1, n2, prec
        complex*16 z1(n1), z2(n2)
C output data declaration
        real*8 alfa1(n1), alfa2(n2), phi1(n1), phi2(n2), u, 
     &         qwork((2*(n1+n2)+3)*prec)
        complex*16 w1(n1), w2(n2), c
C other
        integer ishape, iguess, linearc
        real*8 tol

C Initializing angles
        !external polygon
        call dsc_angles(n1,z1,alfa1,0)
        !internal polygon
        call dsc_angles(n2,z2,alfa2,1)

C generate the gauss-jacobi weights & nodes and check the input
        !0 for finite regions, 1 for infinite regions
        ishape = 0
        call dsc_qinit(n1,n2,alfa1,alfa2,prec,qwork)
C        call dsc_check(alfa1,alfa2,n1,n2,ishape)

C   SOLVE THE ACCESSORY PARAMETER PROBLEM:
        !(=0) a non-equally spaced initial guess
        !(=1) the other  equally-spaced  initial guess provided
        !(=2) user-supplied initial guess which is required
        !     to be the arguments of the initial prevertices.
        iguess = 1
        !(=0) integrating along line segment whenever possible
        !(=1) integrating along circular arc whenever possible
        linearc = 1
        !tolerance to control the convergence
        tol = 1.D-10
        ! on return
        ! u, c, w0, w1 will contain computed parameters.
        ! phi0, phi1 will contain the arguments of the prevertices
        call dscsolv(tol,iguess,n1,n2,
     &               u,c,w1,w2,phi1,phi2,
     &               z1,z2,alfa1,alfa2,prec,qwork,ishape,linearc)

C compute module and return
        dscpack_init = abs(w2(1))
        return
        end


C From canonic to original
        complex*16 function dscpack_backward(
     &           ww1,          ! input point
     &           n1, n2,       ! see description at dscpack_init
     &           z1, z2,       ! ...
     &           prec,         ! ...
     &           alfa1, alfa2, ! ...
     &           w1, w2,       ! ...
     &           phi1, phi2,   ! ...
     &           u, c,         ! ...
     &           qwork         ! ...
     &  )
C input
        integer n1, n2, prec
        complex*16 ww1
        complex*16 z1(n1), z2(n2)
        real*8 alfa1(n1), alfa2(n2), phi1(n1), phi2(n2), u, 
     &         qwork((2*(n1+n2)+3)*prec)
        complex*16 w1(n1), w2(n2), c
C external 
        complex*16 zdsc
        external zdsc
C other
        integer iopt, kww, ic
        !=1 IS NORMALLY ASSUMED FOR ANY GEOMETRY.
        !=# ASSUMES THAT ONLY LINE SEGMENT PATH IS USED.
        iopt = 1
        !KWW=K, IC=0 IF WW=W0(K), OR 
        !KWW=K, IC=1 IF WW=W1(K), OR
        !KWW=0 AND IC=2 OTHERWISE.
        kww = 0
        ic = 2
        dscpack_backward = zdsc(ww1,kww,ic,n1,n2,u,c,w1,w2,z1,z2,
     &                          alfa1,alfa2,phi1,
     &                          phi2,prec,qwork,iopt)
        return
        end


C From original to canonic
        complex*16 function dscpack_forward(
     &           zz1,          ! input point
     &           n1, n2,       ! see description at dscpack_init
     &           z1, z2,       ! ...
     &           prec,         ! ...
     &           alfa1, alfa2, ! ...
     &           w1, w2,       ! ...
     &           phi1, phi2,   ! ...
     &           u, c,         ! ...
     &           qwork         ! ...
     &  )
C input
        integer n1, n2, prec
        complex*16 zz1
        complex*16 z1(n1), z2(n2)
        real*8 alfa1(n1), alfa2(n2), phi1(n1), phi2(n2), u, 
     &         qwork((2*(n1+n2)+3)*prec)
        complex*16 w1(n1), w2(n2), c
C external 
        complex*16 wdsc
        external wdsc
C other
        integer iopt
        real*8 eps

        ! eps is the required accuracy supplied by the user
        eps = 1.0d-10
        !=1 IS NORMALLY ASSUMED FOR ANY GEOMETRY.
        !=# ASSUMES THAT ONLY LINE SEGMENT PATH IS USED.
        iopt = 1
        dscpack_forward = wdsc(zz1,n1,n2,u,c,w1,w2,z1,z2,
     &                         alfa1,alfa2,phi1,phi2,
     &                         prec,qwork,eps,iopt)

        return
        end


