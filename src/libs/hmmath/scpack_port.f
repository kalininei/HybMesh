        real*8 function scpack_init(
     &           n,                   !number of polygon points
     &           wcoords, w0,         !polygon points and center
     &           corners,             !indexes of corners (from 1)
     &           prec,                !precision of integrating (int > 3)
     &           zcoords, 
     &           factor, factor2,
     &           wcoords2, w02,
     &           betam, qwork, qwork2
     &  )

C Builds a map from arbitrary polygon to rectangle
C Returns conformal module of input polygon
C input
        integer n, corners(4), prec
        complex*16 wcoords(n), w0
        complex*16 wcoords2(4), w02
C output
        real*8 betam(n), qwork(prec*(2*n+3)), qwork2(prec*11)
        complex*16 zcoords(n), factor, factor2
C internal
        real*8 typical_size, tol, eest, betam2(4)
        complex*16 z2(4), zero, wsc

C dimensions
        data zero /(0.d0,0.d0)/
        data betam2 /-5d-1, -5d-1, -5d-1, -5d-1/

C calculate tolerance
        typical_size = abs(wcoords(2)-wcoords(1))
        do 100 i=1,n-1
                typical_size = typical_size
     &            + abs(wcoords(i+1)-wcoords(i))
100             continue
        typical_size = typical_size/(n-1);
        tol = 10.d0**(-prec-1) * typical_size

C transform from disk to polygon
        call angles(n, wcoords, betam)
        call qinit(n, betam, prec, qwork)
        call scsolv(-2, 0, tol, eest, n, factor, zcoords, w0, wcoords,
     &          betam, prec, qwork)

c map from disk to rectangle
        do 200 k = 1, 4
                z2(k) = zcoords(corners(k))
200             continue
        call qinit(4, betam2, prec, qwork2)

        factor2 = (1.0d0, 0.0d0)
        do 300 k = 1,4
                wcoords2(k) = wsc(z2(k), k, zero, zero, 0, 4,
     &                      factor2, z2, betam2,
     &                      prec, qwork2)
300             continue

        factor2 = (0.d0,1.d0) / (wcoords2(1)-wcoords2(2))
        w02 = -factor2*wcoords2(2)
        do 400 k = 1,4
                wcoords2(k) = factor2*wcoords2(k) + w02
400             continue
        scpack_init = abs(wcoords2(3)-wcoords2(2))

        return
        end

c ******************************* mapping function
        complex*16 function scpack_forward(
     &          point,                !point in polygon
     &          n,                    !number of points in polygon
     &          wcoords, w0,          !polygon points and center
     &          wcoords2, w02,        !rectangle points and center
     &          zcoords, zcoords2,    !poly/rect coordinates in disk
     &          factor, factor2,      !complex scaling constants
     &          betam, qwork, qwork2, !working arrays
     &          prec                  !precision
     &  )
c input
        complex*16 point, zcoords(n), wcoords(n), w0, factor, factor2
        complex*16 zcoords2(4), wcoords2(4), w02
        integer n, prec
        real*8 betam(n), qwork(prec*(2*n+3)), qwork2(prec*11)
c internal
        complex*16 zn1, wn1, zn2, wn2, zz, zsc, wsc
        real*8 ztolerance, betam2(4)
        integer kn1, kn2, ier

        data ztolerance/1d-8/
        data betam2/-5d-1, -5d-1, -5d-1, -5d-1/

c polygon -> disk
        call nearw(point, zn1, wn1, kn1, n, zcoords,
     &    w0, wcoords, betam)
        !error flag
        ier = 0
        !no initial value (iguess=0), always use ode
        zz = zsc(point, 0, (0d0, 0d0), zn1, wn1, kn1, ztolerance,
     &           ier, n, factor,
     &           zcoords, w0, wcoords, betam, prec, qwork)

c disk -> rectangle
        call nearz(zz, zn2, wn2, kn2, 4, zcoords2,
     &             w02, wcoords2, betam2)

        scpack_forward = wsc(zz, 0, zn2, wn2, kn2, 4, factor2, 
     &                       zcoords2, betam2, prec, qwork2)

        return
        end



c ******************************* mapping function
c From rectangle to Polygon
        complex*16 function scpack_backward(
     &          point,                !point in rectangle
     &          n,                    !number of points in polygon
     &          wcoords, w0,          !polygon points and center
     &          wcoords2, w02,        !rectangle points and center
     &          zcoords, zcoords2,    !poly/rect coordinates in disk
     &          factor, factor2,      !complex scaling constants
     &          betam, qwork, qwork2, !working arrays
     &          prec                  !precision
     &  )
c input
        complex*16 point, zcoords(n), wcoords(n), w0, factor, factor2
        complex*16 zcoords2(4), wcoords2(4), w02
        integer n, prec
        real*8 betam(n), qwork(prec*(2*n+3)), qwork2(prec*11)
c internal
        complex*16 zn1, wn1, zn2, wn2, zz, zsc, wsc
        real*8 ztolerance, betam2(4)
        integer kn1, kn2, ier

        data ztolerance/1d-8/
        data betam2/-5d-1, -5d-1, -5d-1, -5d-1/

c rectangle -> disk
        call nearw(point, zn1, wn1, kn1, 4, zcoords2,
     &    w02, wcoords2, betam2)
        !error flag
        ier = 0
        !no initial value (iguess=0), always use ode
        zz = zsc(point, 0, (0d0, 0d0), zn1, wn1, kn1, ztolerance,
     &           ier, 4, factor2,
     &           zcoords2, w02, wcoords2, betam2, prec, qwork2)

c disk -> polygon
        call nearz(zz, zn2, wn2, kn2, n, zcoords,
     &             w0, wcoords, betam)

        scpack_backward = wsc(zz, 0, zn2, wn2, kn2, n, factor, 
     &                       zcoords, betam, prec, qwork)

        return
        end
