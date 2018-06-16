! Get the code from [Bremer 2016] and run the following commands:
! make
! gfortran -O3 -c laguerre_quad_Adj.f90 
! gfortran -O3 -o comparison explExpa.f90 utils.o chebyshev.o odesolve.o kummer.o laguerre_quad_Adj.o
! To run the timings into Timings.out: ./comparison
! To run the memory tests into the other *.out: ./comparison asdf
program explExpa
    use kummer_laguerre_routines
    double precision, parameter :: pi = acos(-1.d0)
    double precision, allocatable :: x(:), w(:), xt(:), wt(:)
    double precision :: alpha = 0.0
    integer :: n = 100, mode = 1, ni
    character(len=80) :: tmp

    mode = command_argument_count()
    if (mode == 0) then
	! kummer only seems correct for alpha = 0
	call testAll(0.0d0);
	return
    elseif (mode == 1) then
	open(7, file='Ns.out')
	call SYSTEM('echo " " > NbAlloc.out')
	call SYSTEM('echo " " > BytesAlloc.out')
	call SYSTEM('echo " " > MrssTime.out')
	do ni = 1, 15
	    n = nint(exp( log(100*1.0)*(15 -ni)/(15 -1.0) + (ni-1.0)*log(10**9*1.0)/(15-1.0) ) )
	    do mode = 1,4
		write (tmp, '(a, i10, a, i1, a)'), 'valgrind ./comparison ', n, ' 0.0 ', mode, ' 2> valg'
		call SYSTEM(tmp)
		call SYSTEM('grep -o -P  "(?<=heap usage:).*(?=allocs,)" valg | tr -d ''\n'' | tr -d '','' >> NbAlloc.out')
		call SYSTEM('grep -o -P  "(?<=frees,).*(?=bytes al)" valg | tr -d ''\n'' | tr -d '','' >> BytesAlloc.out')
		write (tmp, '(a, i10, a, i1, a)'), '/usr/bin/time -v ./comparison ', n, ' 0.0 ', mode, ' 2> ubtime'
		call SYSTEM(tmp)
		call SYSTEM('grep -Po  "(?<=Maximum resident set size \(kbytes\):).*(?=)" ubtime | tr -d ''\n'' >> MrssTime.out')
	    enddo
	    call SYSTEM('echo " " >> NbAlloc.out')
	    call SYSTEM('echo " " >> BytesAlloc.out')
	    call SYSTEM('echo " " >> MrssTime.out')
	    write(7,*) n
	enddo
	close(7)
	return
    endif

    call get_command_argument(1,tmp)
    read(tmp, '(i10)' ) n
    call get_command_argument(2,tmp)
    read(tmp, '(f15.10)' ) alpha
    call get_command_argument(3,tmp)
    read(tmp, '(i10)' ) mode

    if (mode == 1) then ! Compressed representation with bes(zer)
	call laguerreExp(n, alpha, x, w, .true., .true., .true.)
    else
	allocate(x(n), w(n))
    endif

    if (mode == 2) then ! Uncompressed representation without bes(zer)
	call laguerreExp(n, alpha, x, w, .false., .false., .false.)
    elseif (mode == 3) then ! Uncompressed representation with bes(zer)
	call laguerreExp(n, alpha, x, w, .true., .true., .false.)
    elseif (mode == 4) then ! Nonoscillatory phase
	call kummer_laguerre(n, alpha, x, w)
    endif
    ! Use the result such that the compiler does not optimize out the computation
    if (abs(sum(w)/GAMMA(alpha+1) -1) > 1e-8) then
	print *, sum(w), "error", GAMMA(alpha+1), w(1:10)
    endif
    deallocate(x,w)
contains

subroutine testAll(alpha)
    use kummer_laguerre_routines
    implicit none
    integer, parameter :: minN = 100, maxN = 10**9, nbRep = 5, lN = 15
    integer :: ni, n, nTest, t1, t2, clm, dNt(8), rate, ns(lN), mode, row, whTim
    double precision, allocatable :: xn(:), wn(:), xe(:), we(:), xcpr(:), wcpr(:)
    double precision :: alpha, ers(lN,12)
    real :: ts, tm, te, td, tarr(2), elTimes(lN,nbRep,4,4), start, prevToc = 0.0
    elTimes = -1.d0
    ers = -1.d0
    ns = (/ (nint(exp( log(minN*1.0d0)*((lN -ni)/(lN -1.0d0)) + &
log(maxN*1.0d0)*((ni-1.0d0)/(lN-1.0d0)) ) ), ni=1,lN) /)
    do ni =1,lN
	n = ns(ni)
	allocate(xn(n), wn(n), xe(n), we(n))
	do mode = 1, 4
          do nTest = 1, nbRep
	    call etime(tarr, ts) 
	    call cpu_time(td)
	    call date_and_time(values=dNt)
	    start = 0.001d0*dNt(8) + dNt(7) + 60*(dNt(6) + 60*(dNt(5) + 24*dNt(3) ))
	    call system_clock(t1, rate, clm)

	    if (mode == 1) then ! Compressed representation with bes(zer)
		call laguerreExp(n, alpha, xcpr, wcpr, .true., .true., .true.)
		deallocate(xcpr, wcpr)
	    end if
	    if (mode == 2) then ! Uncompressed representation without bes(zer)
		call laguerreExp(n, alpha, xe, we, .false., .false., .false.)
	    elseif (mode == 3) then ! Uncompressed representation with bes(zer)
		call laguerreExp(n, alpha, xe, we, .true., .true., .false.)
	    elseif (mode == 4) then ! Nonoscillatory phase
		call kummer_laguerre(n, alpha, xn, wn)
	    endif

	    call etime(tarr, tm)
	    elTimes(ni,nTest,1, mode) = tm-ts
	    call date_and_time(values=dNt)
	    elTimes(ni,nTest,2, mode) = 0.001d0*dNt(8) + dNt(7) + 60*(dNt(6) + 60*(dNt(5) + 24*dNt(3) )) - start
	    call cpu_time(te)
	    elTimes(ni,nTest,3, mode) = te-td
	    call system_clock(t2, rate, clm)
	    elTimes(ni,nTest,4, mode) = (t2 -t1 +0.0)/rate
	    if (tm - prevToc > 3.0) then ! print progress info every three seconds
		print '(i2,a,i1,a, i2, a, es12.3, a, es12.3, a)', ni, "=ni, mode=", mode, ", nTest=", nTest, ", etime=", tm, & 
", ", tm*(sum(ns)*4.0*nbRep/( sum(ns(1:ni-1))*4.0*nbRep +n*(mode-1.0)*nbRep +n*nTest) -1.0)/3600.0, " hr left."
		prevToc = tm
	    endif
	  enddo
	enddo
	! xe and we are now the last simulations (ntest=nbRep) with the uncompressed representation with computation of Bessel function and zeros (mode = 3)
	nTest = floor(0.9*n)-1
	ers(ni,1) = norm2(xn -xe)/norm2(xn)
	ers(ni,2) = norm2(wn -we)/norm2(wn)
	ers(ni,3) = norm2((xn -xe)/xn)
	ers(ni,4) = norm2((wn -we)/wn)
	ers(ni,5) = norm2(xn -xe)
	ers(ni,6) = norm2(wn -we)
	ers(ni,7) = norm2(xn(1:nTest) -xe(1:nTest))/norm2(xn(1:nTest))
	ers(ni,8) = norm2(wn(1:nTest) -we(1:nTest))/norm2(wn(1:nTest))
	ers(ni,9) = norm2((xn(1:nTest) -xe(1:nTest))/xn(1:nTest))
	ers(ni,10) = norm2((wn(1:nTest) -we(1:nTest))/wn(1:nTest))
	ers(ni,11) = norm2(xn(1:nTest) -xe(1:nTest))
	ers(ni,12) = norm2(wn(1:nTest) -we(1:nTest))
	deallocate(xn, wn, xe, we)
	open(8, file="Timings.out")
	write(8, *), ns, "=ns, elTimes = "
	do mode =1,size(elTimes,4)
	    do whTim =1,size(elTimes,3)
		write(8,*), "Mode", mode, " time nr ", whTim, " = "
		do row = 1,size(elTimes,1)
		    write (8,*), elTimes(row,:,whTim,mode)
		enddo
	    enddo
	enddo
	write (8, *), ", ers = "
	do row = 1,size(ers,1)
	    write (8,*), ers(row,:)
	enddo
	close(8)
    end do
end subroutine 











subroutine laguerreExp(n, alpha, x, w, compZer, compBes, compRepr)
    implicit none
    integer :: n, T, ibes, iair, it, k, mn, mxb
    double precision :: alpha, d, pt, term, co, so, cf, zeta, prev, jta
    double precision, parameter :: ak(11)  = (/ -13.69148903521072d0, &
 -12.828776752865757d0, -11.93601556323626d0,  -11.00852430373326d0, & 
 -10.04017434155809d0,  -9.02265085340981d0, -7.944133587120853d0, &
  -6.786708090071759d0, -5.520559828095551d0, -4.08794944413097d0, &
  -2.338107410459767d0 /) ! First 11 roots of the Airy function in reverse order
    double precision, parameter :: dak(10)  = (/  &
 -1.0677938592d0, 1.0487206486d0,  -1.0277386888d0, & 
 1.0043701227d0,  -0.9779228086d0, 0.9473357094d0, &
  -0.9108507370d0, 0.8652040259d0, -0.8031113697d0, &
  0.7012108227d0 /) ! Derivatives at first 11 roots of the Airy function in reverse order [DLMF tab 9.9.1]
    double precision, allocatable :: x(:), w(:)
    logical compZer, compBes, compRepr

    double precision mu, b, a1, a3, a5, a7, a9, a11, a13, TT(30)
    double precision, parameter :: tab(20) = (/ 2.4048255576957728d0, &
 5.5200781102863106d0,  8.6537279129110122d0,   11.791534439014281d0, &
 14.930917708487785d0,  18.071063967910922d0,   21.211636629879258d0, &
 24.352471530749302d0,  27.493479132040254d0,   30.634606468431975d0, &
 33.775820213573568d0,  36.917098353664044d0,   40.058425764628239d0, &
 43.199791713176730d0,  46.341188371661814d0,   49.482609897397817d0, &
 52.624051841114996d0,  55.765510755019979d0,   58.906983926080942d0, &
 62.048469190227170d0  /)  ! First 20 roots of J0(x) are precomputed (using Wolfram Alpha)
    double precision, parameter :: C(6,30) = reshape( &
(/2.883975316228d0,  8.263194332307d0, 11.493871452173d0, 14.689036505931d0, &
 17.866882871378d0, 21.034784308088d0, &
  0.767665211539d0,  4.209200330779d0,  4.317988625384d0,  4.387437455306d0, &
  4.435717974422d0,  4.471319438161d0, &
 -0.086538804759d0, -0.164644722483d0,  -0.130667664397d0, -0.109469595763d0, &
 -0.094492317231d0, -0.083234240394d0, &
  0.020433979038d0,  0.039764618826d0,  0.023009510531d0,  0.015359574754d0, &
  0.011070071951d0,  0.008388073020d0, &
 -0.006103761347d0, -0.011799527177d0, -0.004987164201d0, -0.002655024938d0, &
 -0.001598668225d0, -0.001042443435d0, &
  0.002046841322d0,  0.003893555229d0,  0.001204453026d0,  0.000511852711d0, &
  0.000257620149d0,  0.000144611721d0, &
 -0.000734476579d0, -0.001369989689d0, -0.000310786051d0, -0.000105522473d0, &
 -0.000044416219d0, -0.000021469973d0, &
  0.000275336751d0,  0.000503054700d0,  0.000083834770d0,  0.000022761626d0, &
  0.000008016197d0,  0.000003337753d0, &
 -0.000106375704d0, -0.000190381770d0, -0.000023343325d0, -0.000005071979d0, &
 -0.000001495224d0, -0.000000536428d0, &
  0.000042003336d0,  0.000073681222d0,  0.000006655551d0,  0.000001158094d0, &
  0.000000285903d0,  0.000000088402d0, &
 -0.000016858623d0, -0.000029010830d0, -0.000001932603d0, -0.000000269480d0, &
 -0.000000055734d0, -0.000000014856d0, &
  0.000006852440d0,  0.000011579131d0,  0.000000569367d0,  0.000000063657d0, &
  0.000000011033d0,  0.000000002536d0, &
 -0.000002813300d0, -0.000004672877d0, -0.000000169722d0, -0.000000015222d0, &
 -0.000000002212d0, -0.000000000438d0, &
  0.000001164419d0,  0.000001903082d0,  0.000000051084d0,  0.000000003677d0, &
  0.000000000448d0,  0.000000000077d0, &
 -0.000000485189d0, -0.000000781030d0, -0.000000015501d0, -0.000000000896d0, &
 -0.000000000092d0, -0.000000000014d0, &
  0.000000203309d0,  0.000000322648d0,  0.000000004736d0,  0.000000000220d0, &
  0.000000000019d0,  0.000000000002d0, &
 -0.000000085602d0, -0.000000134047d0, -0.000000001456d0, -0.000000000054d0, &
 -0.000000000004d0,             0.d0, &
  0.000000036192d0,  0.000000055969d0,  0.000000000450d0,  0.000000000013d0, &
  0.d0, 0.d0, &
 -0.000000015357d0, -0.000000023472d0, -0.000000000140d0, -0.000000000003d0, &
  0.d0, 0.d0, &
  0.000000006537d0,  0.000000009882d0,  0.000000000043d0,  0.000000000001d0, &
  0.d0, 0.d0, &
 -0.000000002791d0, -0.000000004175d0, -0.000000000014d0, 0.d0, 0.d0, 0.d0, &
  0.000000001194d0,  0.000000001770d0,  0.000000000004d0, 0.d0, 0.d0, 0.d0, &
 -0.000000000512d0, -0.000000000752d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.000000000220d0,  0.000000000321d0, 0.d0, 0.d0, 0.d0, 0.d0, &
 -0.000000000095d0, -0.000000000137d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.000000000041d0,  0.000000000059d0, 0.d0, 0.d0, 0.d0, 0.d0, &
 -0.000000000018d0, -0.000000000025d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.000000000008d0,  0.000000000011d0, 0.d0, 0.d0, 0.d0, 0.d0, &
 -0.000000000003d0, -0.000000000005d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.000000000001d0,  0.000000000002d0, 0.d0, 0.d0, 0.d0, 0.d0 /), &
  (/ 6, 30 /) ) 

! McMahon's expansion. This expansion gives very accurate approximation
! for the k-th zero (k >= 7) for extreme nu, and moderate approximation 
! in other cases for low k or medium-sized Bessel parameters.
    mu = 4*alpha**2
    a1 = 1.d0/8
    a3 = (7*mu-31)/384
    a5 = 4*(3779 +mu*(-982 +83*mu))/61440
    a7 = 6*(-6277237 +mu*(1585743 +mu*(-153855 +6949*mu)))/20643840;
    a9 = 144*(2092163573 +mu*(-512062548 +mu*(48010494 +mu*(-2479316 + & 
70197*mu)) )) / 11890851840.d0
    a11 = 720.d0*(-8249725736393.d0 +mu*(1982611456181.d0 +mu*(-179289628602.d0 &
+mu*(8903961290.d0 +mu*(-287149133.d0 +5592657.d0*mu) ) ) ) ) / 10463949619200.d0
    a13 = (576.d0*(423748443625564327.d0 + mu*(-100847472093088506.d0 &
+mu*(8929489333108377.d0 + mu*(-426353946885548.d0 +mu*(13172003634537.d0 &
+mu*(-291245357370.d0 + mu*4148944183.d0))) ))) / 13059009124761600.d0)

    if (compRepr) then
	mn = min(nint(17*sqrt(n*1.0d0)), n)
	allocate(x(mn), w(mn))
    else
	mn = n
    endif

    ! This is a heuristic for the number of terms in the expansions that follow.
    T = ceiling(34/log(n+0.0) )
    d = 1/(4*n+2*alpha+2)

    ibes = max(ceiling(sqrt(n*1.0) ), 7)
    ! Heuristics to switch between Bessel, bulk and Airy initial.
    iair = floor(0.9*n)
    mxb = min(mn, iair-1)
!print *, ibes, iair, "=iair, mxb=", mxb
    if (compBes) then
	b = (alpha-2)/3
	TT(1) = 1.d0
	TT(2) = b
	do k = 2,29
	        TT(k+1) = 2*b*TT(k) - TT(k-1)
	enddo
    endif
    jta = 1.0d0
    x = 0.d0
    w = 0.d0
    do k=1,ibes !Bessel region
	! Roots of the function J_v(x) from chebfun by L. L. Peixoto, 2015, later modified by A. Townsend to work in Julia
	if ((compBes) .and. (alpha .eq. 0.d0) .and. (k < 21)) then
	   jta = tab(k)
	elseif ((compBes) .and. (alpha < 5) .and. (k < 7)) then
	   ! Piessens's Chebyshev series approximations (1984). Calculates the 6 first
	   ! zeros to at least 12 decimal figures in region -1 <= V <= 5:
	   
	   jta = dot_product(C(k,:),TT) 
	   if (k == 1) then
		jta = jta*sqrt(alpha+1) ! Scale the first root.
	   endif
	elseif (compBes) then
	   b = 0.25d0*(2*alpha+4*k-1)*pi
	   ! Evaluate using Horner's scheme: Possibly inaccurate results
	   jta = (((a13/b**2 + a11)/b**2 + a9)/b**2 + a7)/b**2 + a5
	   jta = b - (mu-1)*( ((jta/b**2 + a3)/b**2 + a1)/b)
	endif

	!Computed j_{alpha,k}, now add higher order terms
        if (T >= 7) then
        ! These higher order terms in the left and bulk region are derived in [Huybrechs and Opsomer 2018, in preparation]
	   x(k) = x(k) + (657*jta**6 +36*jta**4*( &
73*alpha**2-181) +2*jta**2*(2459*alpha**4 -10750*alpha**2 &
 +14051) + 4*(1493*alpha**6 -9303*alpha**4 +19887*alpha**2 - 12077))*d**6/2835
	   w(k) = w(k) + (11944*alpha**6 + 5256*jta**6 &
- (5061*alpha**5 + 5085*alpha**4 + 4830*alpha**3 -22724*alpha**2 -22932*alpha &
+ 39164)*jta**4 - 74424*alpha**4 + 8*(2459*alpha**4 -10750*alpha**2 &
+ 14051)*jta**2 + 159096*alpha**2 - 96616)/2835/2*d**6
       endif
       if (T >= 5) then
	   x(k) = x(k) + (11*jta**4 +3*jta**2*( &
11*alpha**2-19) +46*alpha**4 -140*alpha**2 +94)*d**4/45
	   w(k) = w(k) + (46*alpha**4 + 33*jta**4 &
+6*jta**2*(11*alpha**2 -19) -140*alpha**2 +94)/45*d**4
       endif
       if (T >= 3) then
	   x(k) = x(k) + (jta**2 + 2*alpha**2 - 2)*d**2/3
	   w(k) = w(k) + (alpha**2 + jta**2 -1)*2/3*d**2
       endif
       x(k) = jta**2*d*(1 + x(k))
       if ((compBes) .and. (k < 6+10*alpha**2/pi)) then
	    ! j_{alpha,k} \sim k*pi + O(k^0) with j_{alpha,k} > (alpha+1)^2 
            pt = 1/GAMMA(alpha+2)
	    term = pt
	    prev = pt*jta**2/(alpha+2)
            do it = 1, 100
		term = -term*jta**2/4/it/(alpha+it+1)
		if (abs(term)  < abs(pt)*10.d0**(-12) ) then 
		    ! Stop when no significant increase any more
		    exit
		endif
		pt = term + pt
		prev = term
	    enddo
	    ! Computed J_{alpha+1}(j_{alpha,k}) but = -J_{alpha-1} and square
w(k)=4*d*x(k)**alpha*exp(-x(k))/(pt*(jta/2)**(alpha+1))**2*(1+w(k))
	elseif (compBes) then
	    ! Approximate J_{\alpha-1}(j_{\alpha,k}) for large arguments
	    pt = 0.d0
	    co = cos(jta -(alpha-1)*pi/2 -pi/4)
	    so = sin(jta -(alpha-1)*pi/2 -pi/4)
            cf = 1.d0
	    prev = 2.0d0
            do it = 0, 50
		term = co*(-1)**it*cf
		cf = cf*(4*(alpha-1)**2 -(4*it +1)**2)/(2*it +1)/8/jta
		term = term - so*(-1)**it*cf
		if ((abs(term) + abs(prev) < abs(pt)*10.d0**(-12) ) .or. &
(abs(term) > abs(prev))) then
		    ! Stop when no significant increase any more taking into account [1; 0; 0.5; 0; 0.25 ...], or when the series starts diverging
		    exit
		endif
		pt = term + pt
		prev = term
		cf = cf*(4*(alpha-1)**2 -(4*it +3)**2)/(2*it+2)/8/jta
	    enddo
	    w(k) = 2*d*x(k)**alpha*exp(-x(k))/pt**2*jta*pi*(1+w(k)) 
	else
	    w(k) = 4*d*x(k)**alpha*exp(-x(k))*(1+w(k))
	endif
    enddo ! End loop bessel region

    do k=ibes+1,mxb !Loop lens
	if (compBes) then
	    pt = (4*n -4*k +3)*d
	    jta = pi**2/16*(pt -1)**2 ! This t is not a very good initial approximation of the inverse function of f(y) = (4*n -4*k +3)*d*pi +2*sqrt(y)*sqrt(1-y) -acos(2*y-1); but better approximations are much more time-consuming
	    do it = 1,6
jta = jta - (pt*pi +2*sqrt(jta -jta**2) -acos(2*jta -1) )*sqrt(jta/(1-jta))/2
	    end do
	endif
        if (T >= 7) then
        ! These higher order terms in the left and bulk region are derived in [Huybrechs and Opsomer 2018, in preparation]
	   x(k) = x(k) -d**5/181440*(9216*(21*alpha**6 &
- 105*alpha**4 + 147*alpha**2 - 31)*jta**10 &
-69120*(21*alpha**6 - 105*alpha**4 + 147*alpha**2 - 31)*jta**9 &
+384*(12285*alpha**6 -61320*alpha**4 +85785*alpha**2 &
-18086)*jta**8 - 64*(136080*alpha**6 - 675675*alpha**4 &
+943110*alpha**2 - 198743)*jta**7 &
+ 144*(70560*alpha**6 - 345765*alpha**4 +479850*alpha**2 &
- 101293)*jta**6 + 72576*alpha**6 - (8128512*alpha**6 &
- 38656800*alpha**4+ 52928064*alpha**2 - 13067711)*jta**5 &
+ 5*(1016064*alpha**6 - 4581360*alpha**4 +6114528*alpha**2 &
+ 113401)*jta**4 - 317520*alpha**4 - 10*(290304*alpha**6 &
-1245888*alpha**4 + 1620864*alpha**2 - 528065)*jta**3 &
+ 5*(290304*alpha**6 -1234800*alpha**4 + 1598688*alpha**2 &
- 327031)*jta**2 + 417312*alpha**2 -5*(96768*alpha**6 &
- 417312*alpha**4 + 544320*alpha**2 - 111509)*jta &
-85616)/(jta-1)**8/jta**2
	   w(k) = w(k) + d**6/362880*(9216*(21*alpha**6 &
- 105*alpha**4 + 147*alpha**2 - 31)*jta**10 &
-1536*(945*alpha**6-4830*alpha**4 +6825*alpha**2 -1444)*jta**9 &
+ 384*(11340*alpha**6 -60165*alpha**4 + 86310*alpha**2 &
- 18289)*jta**8 - 2*(2903040*alpha**6 - 17055360*alpha**4 &
+ 25401600*alpha**2 - 5*alpha - 5252997)*jta**7 &
- (11753280*alpha**4 - 23506560*alpha**2+ 67*alpha &
- 13987519)*jta**6 - 290304*alpha**6 + 12*(1016064*alpha**6 &
-3578400*alpha**4 +4108608*alpha**2 +16*alpha +7100871)*jta**5 &
 - 5*(4064256*alpha**6 -16559424*alpha**4 + 20926080*alpha**2 + 61*alpha &
- 15239393)*jta**4 + 1270080*alpha**4 +10*(1741824*alpha**6 &
-7386624*alpha**4 +9547776*alpha**2 +29*alpha -1560107)*jta**3 &
- 15*(580608*alpha**6 - 2503872*alpha**4 + 3265920*alpha**2 + 11*alpha &
- 669051)*jta**2- 1669248*alpha**2 + 4*(604800*alpha**6 &
- 2630880*alpha**4 + 3447360*alpha**2 + 13*alpha- 706850)*jta &
- 7*alpha + 342463)/(jta -1)**9/jta**3
	endif
	if (T >= 5) then
	   x(k) = x(k) - d**3/720*(32*(15*alpha**4 &
- 30*alpha**2 + 7)*jta**6 -144*(15*alpha**4 - 30*alpha**2 &
+ 7)*jta**5 + 16*(225*alpha**4 - 450*alpha**2 &
+104)*jta**4 - 240*alpha**4 - 480*(5*alpha**4 - 10*alpha**2 &
+ 1)*jta**3 + 480*alpha**2 +45*(16*alpha**4 - 32*alpha**2 &
+ 7)*jta +990*jta**2 -105)/(jta &
-1)**5/jta
	   w(k) = w(k) + d**4/720*(16*(15*alpha**4 &
- 30*alpha**2 + 7)*jta**6 - 32*(45*alpha**4 - 90*alpha**2 &
+22)*jta**5 + 48*(75*alpha**4 - 150*alpha**2 + &
74)*jta**4 + 240*alpha**4 - 600*(8*alpha**4- 16*alpha**2 &
- 5)*jta**3 + 45*(80*alpha**4 - 160*alpha**2 + &
57)*jta**2 - 480*alpha**2 -90*(16*alpha**4 - 32*alpha**2 &
+ 7)*jta + 105)/(jta -1)**6/jta**2
	endif
	if (T >= 3) then
	   x(k) = x(k) - d/12*(4*(3*alpha**2 &
-1)*jta**2 +12*alpha**2 -12*(2*alpha**2 -1)*jta - 3)/(jta-1)**2
	   w(k) = w(k)  + d**2/6*(2*jta + 3)/(jta-1)**3
	endif
	x(k) = x(k) + jta/d;
	w(k) = x(k)**alpha*exp(-x(k) )*2*pi*sqrt(jta/(1-jta))*(1 +w(k))
    enddo !Loop lens


    do k =iair,mn !Loop Airy
	if ((compBes) .and. (k < n - 10)) then
	    pt = 3*pi/2*( n -k +0.75d0) ! [DLMF (9.9.6)]
	    jta = -pt**(2.d0/3)*(1 + 5/48.d0/pt**2 - 5/36.d0/pt**4 + &
77125/82944.d0/pt**6 -10856875/6967296.d0/pt**8)
	elseif (compBes) then
	    jta = ak(k -n +11)
	endif
	! For the Airy region, only O(n**{-4}) relative error for x and O(n**{-2/3}) for w as the latter are extremely small or even underflow
	if (T >= 5) then
	   x(k) = x(k) -(15152.d0/3031875*jta**5 &
+1088.d0/121275*jta**2)*2**(1.d0/3)*d**(7.0d0/3) ! [Gatteshi 2002 (4.9)], Gives an O(n**{-4}) relative error
	endif
	if (T >= 3) then
	   x(k) = x(k) + jta**2*(d*16)**(1.d0/3)/5 &
+ (11.d0/35 -alpha**2 -12.d0/175*jta**3)*d + (16.d0/1575*jta &
+ 92.d0/7875*jta**4)*2**(2.d0/3)*d**(5.0d0/3)
	endif
	x(k) = x(k) + 1/d +jta*(d/4)**(-1.d0/3)
	if ((compBes) .and. (k <= n-10)) then
	    zeta = -2*jta*sqrt(-jta)/3 ! = zeta in [DLMF 9.7.10]
	    co = cos(zeta -pi/4)
	    so = sin(zeta -pi/4)
! Asymptotic expansion of Airy(z) = exp(-zeta)/2/sqrt(pi)/z^(1/4) sum_{n=0}^\infty (-1)^n Gamma(n+5/6) Gamma(n+1/6) (3/4)^n /2/pi/ n!/z^(3n/2) but for pos z and
! (6*k+1)/(1-6*k)/216^k/factorial(k)*product  (2*j+1) for j in range(k,3*k-1)  = -(2^(-5 k) 27^(-k) (6 k + 1) Γ(6 k - 1))/(Γ(2 k) Γ(k + 1))
	    pt = 0.d0
	    cf = 1.d0 ! = GAMMA(5.d0/6)*GAMMA(7.d0/6)*3/pi
	    prev = 2.0d0
	    do it = 0, 50
		term = so*cf/(1-12*it)*(-1)**it
		cf = cf*(2*it +5.d0/6)*(2*it+7.d0/6)/2/(2*it+1)/zeta
		term = term - co*cf/(1-6*(2*it+1))*(-1)**it 
		cf = cf*(2*it +11.d0/6)*(2*it+13.d0/6)/2/(2*it+2)/zeta
		if ((abs(term) + abs(prev) < abs(pt)*10.d0**(-12) ) .or. (abs(term) > abs(prev))) then
		    ! Stop when no significant increase any more taking into account [1; 0; 0.5; 0; 0.25 ...], or when the series starts diverging
		    exit
		 endif
		pt = pt + term
		prev = term
	    enddo
w(k) = 4**(1.d0/3)*x(k)**(alpha +1.d0/3)*exp(-x(k) )/( pt**2*sqrt(-jta)/pi )
	elseif (compBes) then
	    ! Could also use a series expansion of the Airy function as Ai"(z)=z Ai(z)
	    w(k) = 4**(1.d0/3)*x(k)**(alpha+1.d0/3)*exp(-x(k))/(dak(k-n+10))**2
	else
	    w(k) = 4**(1.d0/3)*x(k)**(alpha+1.d0/3)*exp(-x(k))
	endif
    enddo !Loop Airy

    if (((minval(x) < 0.d0) .or. (maxval(x) > 4*n +2*alpha +2) .or. &
(minval(w) < 0.d0)) .and. compBes) then
	print *, "Wrong node or weight."
    endif
end subroutine

end program
