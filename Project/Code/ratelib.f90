! reaction rate library
module rates
      use interfaces
      implicit none

contains
! torch rates
! li7(t,n)   a(an,g)    be9(p,d)    be9(p,n)    b10(a,n)   b11(a,n)
! n14(p,a)   c11(p,g)   c12(a,n)    c13(a,n)    c13(p,n)   c14(a,g)
! c14(p,n)   c14(p,g)   o16(p,a)    n14(p,n)    n14(a,n)   n15(p,n)
! n15(a,n)   n15(a,g)   o14(a,g)    o17(a,g)    o17(a,n)   o18(a,g)
! o18(a,n)   ne20(p,a)  f18(p,g)    f19(p,g)    f19(p,n)   f19(a,p)
! na22(n,a)  ne20(p,g)  na23(p,a)   ne20(n,g)   ne21(p,g)  ne21(a,g)
! ne22(p,g)  ne22(a,g)  na22(n,p)   ne22(a,n)   na21(p,g)  mg24(p,a)
! ne21(a,n)  na22(p,g)  na23(p,g)   na23(p,n)   mg24(p,g)  al27(p,a)
! mg25(p,g)  mg25(a,p)  mg25(a,g)   mg25(a,n)   mg26(p,g)  mg26(a,g)
! mg26(a,n)  al25(p,g)  al26(p,g)   al27(a,n)   si27(p,g)  si28(p,g)
! si29(p,g)  si30(p,g)

! bigbang rates:
! n(e-nu)p   p(e-,nu)n  d(p,n)      d(n,g)      d(d,p)     d(d,n)
! t(p,n)     d(d,g)     t(p,g)      t(d,n)      t(t,2n)    he3(d,p)
! he3(t,d)   he3(t,np)  he4(np,g)   he4(d,g)    he4(t,n)   li6(p,he3)
! li6(n,g)   li7(d,n)   lit(t,2n)   li7(he3,np) li6(p,g)   li7(p,n)
! be7(d,p)   be7(t,np)  be7(3he,2p) li6(a,g)    li7(a,n)   be9(p,g)
! b10(p,a)   li7(a,g)   b11(p,a)    be7(a,g)    b11(p,n)   b8(a,p)
! b10(p,g)   c11(n,a)   be9(a,n)    b11(p,g)    b11(a,p)

! pp123 rates:
! p(p,e+nu)  p(n,g)     d(p,g)      he3(n,g)    he3+he3    he3(a,g)
! be7(e-,nu) be7(p,g)   li7(p,g)    li7(p,a)    b8(e+,nu)

! cno rates:
! c12(p,g)   n13(e-nu)  c13(p,g)    n14(p,g)    o15(e-nu)  n14(a,g)
! n15(p,g)   n15(p,a)   o16(p,g)    o17(p,a)    o17(p,g)   o18(p,a)
! o18(p,g)   f17(e-nu)  f18(e-nu)   f19(p,a)

! hot cno rates
! n13(p,g)   o14(e-nu)  o14(a,p)    o15(a,g)    f17(p,g)   ne18(e-nu)
! f18(p,a)   ne18(a,p)  ne19(p,g)   ne19(e-nu)  si26(a,p)

! alfa chain rates:
! a(aa,g)    c12(a,g)   c12+c12     c12+o16     o16+o16    o16(a,g)
! ne20(a,g)  ne20(a,g)  mg24(a,g)   mg24(a,p)   al27(p,g)  si28(a,g)
! si28(a,p)  p31(p,g)   s32(a,g)    s32(a,p)    cl35(p,g)  ar36(a,g)
! ar36(a,p)  k39(p,g)   ca40(a,g)   ca40(a,p)   sc43(p,g)  ti44(a,g)
! ti44(a,p)  v47(p,g)   cr48(a,g)   cr(a,p)     mn51(p,g)  fe52(a,g)
! fe52(a,p)  co55(p,g)

! photodisintegration rates:
! fe52(n,g) fe53(n,g)  fe54(p,g)    fe54(n,g)   fe54(a,p)   fe55(n,g)  






      subroutine tfactors(temp)
      include 'implno.dek'
      include 'tfactors.dek'

! sets various popular temperature factors into common block
! this routine must be called before any of the rates are called

! declare the pass
      double precision temp

! all these are in common block

      t9    = temp * 1.0d-9
      t92   = t9*t9
      t93   = t9*t92
      t94   = t9*t93
      t95   = t9*t94
      t96   = t9*t95

      t912  = sqrt(t9)
      t932  = t9*t912
      t952  = t9*t932
      t972  = t9*t952

      t913  = t9**oneth
      t923  = t913*t913
      t943  = t9*t913
      t953  = t9*t923
      t973  = t953*t923
      t9113 = t973*t943

      t914  = t9**(0.25d0)
      t934  = t914*t914*t914
      t954  = t9*t914
      t974  = t9*t934

      t915  = t9**onefif
      t935  = t915*t915*t915
      t945  = t915 * t935
      t965  = t9 * t915

      t916  = t9**onesix
      t976  = t9 * t916
      t9i76 = 1.0d0/t976

      t917  = t9**onesev
      t927  = t917*t917
      t947  = t927*t927

      t918  = sqrt(t914)
      t938  = t918*t918*t918
      t958  = t938*t918*t918

      t9i   = 1.0d0/t9
      t9i2  = t9i*t9i
      t9i3  = t9i2*t9i

      t9i12 = 1.0d0/t912
      t9i32 = t9i*t9i12
      t9i52 = t9i*t9i32
      t9i72 = t9i*t9i52

      t9i13 = 1.0d0/t913
      t9i23 = t9i13*t9i13
      t9i43 = t9i*t9i13
      t9i53 = t9i*t9i23

      t9i14 = 1.0d0/t914
      t9i34 = t9i14*t9i14*t9i14
      t9i54 = t9i*t9i14

      t9i15 = 1.0d0/t915
      t9i35 = t9i15*t9i15*t9i15
      t9i45 = t9i15 * t9i35
      t9i65 = t9i*t9i15

      t9i17 = 1.0d0/t917
      t9i27 = t9i17*t9i17
      t9i47 = t9i27*t9i27

      t9i18 = 1.0d0/t918
      t9i38 = t9i18*t9i18*t9i18
      t9i58 = t9i38*t9i18*t9i18

      return
      end





      subroutine rate_aan(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd

! he4(an,g)be9
      aa  = 1.0d0 + 0.344*t9
      bb  = t92 * aa
      dbb = 2.0d0 * t9 * aa + t92*0.344

	print *, aa,t92,bb,dbb

      cc  = 1.0d0/bb
      dcc = -cc*cc*dbb

      dd  = 2.59e-6 * exp(-1.062*t9i)
      ddd = dd*1.062*t9i2

      term    = cc * dd
      dtermdt = dcc*dd + cc*ddd

! rates

      fr    = den * den * term
      dfrdt = den * den * dtermdt * 1.0d-9
      dfrdd = 2.0d0 * den * term

      rev      = 5.84e19 * t93 * exp(-18.260*t9i)
      drevdt   = rev*(3.0d0*t9i + 18.260*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_be9pd(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/0.2704d0)


! be9(p,d)be8 =>2a
      aa  = 2.11e+11 * t9i23 * exp(-10.359*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*10.359*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0  + 0.04*t913 + 1.09*t923 + 0.307*t9 &
            + 3.21*t943 + 2.30*t953
      dbb = oneth*0.04*t9i23 + twoth*1.09*t9i13 + 0.307 &
            + fourth*3.21*t913 + fiveth*2.30*t923

      cc  = 5.79e+08 * t9i * exp(-3.046*t9i)
      dcc = cc*(-t9i + 3.046*t9i2)

      dd  = 8.50e+08 * t9i34 * exp(-5.800*t9i)
      ddd = dd*(-0.75d0*t9i + 5.800*t9i2)

      term    = aa*bb + cc + dd
      dtermdt = daa*bb + aa*dbb + dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.07e-11 * t9i32 *exp(-7.555*t9i)
      drevdt   = rev*(-1.5d0*t9i + 7.555*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end






      subroutine rate_be9pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz


! be9(p,n)b9
      aa  = 5.58e7*(1.0d0 + 0.042*t912 + 0.985*t9)
      daa = 5.58e7*(0.5d0*0.042*t9i12 + 0.985)

      zz  = exp(-21.473*t9i)
      dzz = zz*21.473*t9i2

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc  = 1.02e+09 * t9i32 * exp(-26.725*t9i)
      dcc = cc*(-1.5d0*t9i + 26.725*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term


      bb  = 0.998 * aa
      dbb = 0.998 * daa

      cc  = 0.998 * 1.02e+09 * t9i32 * exp(-5.252*t9i)
      dcc = cc*(-1.5d0*t9i + 5.252*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_b10an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,q1
      parameter        (q1 = 1.0d0/91.948921d0)


! b10(a,n)n13
      term    = 1.20e+13 * t9i23 * exp(-27.989*t9i13 - t92*q1)
      dtermdt = -twoth*term*t9i &
                + term*(oneth*27.989*t9i43 - 2.0d0*t9*q1)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.34 * exp(-12.287*t9i)
      drevdt   = rev*12.287*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_b11an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.0196d0)


! b11(a,n)n14
      aa  = 6.97e+12 * t9i23 * exp(-28.234*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*28.234*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.015*t913 + 8.115*t923 + 0.838*t9 &
            + 39.804*t943 + 10.456*t953
      dbb = oneth*0.015*t9i23 + twoth*8.115*t9i13 + 0.838 &
            + fourth*39.804*t913 + fiveth*10.456*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.79 * t9i32 * exp(-2.827*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.827*t9i2)

      ee   = 1.71e+03 * t9i32 * exp(-5.178*t9i)
      dee  = ee*(-1.5d0*t9i + 5.178*t9i2)

      ff   = 4.49e+06 * t935 * exp(-8.596*t9i)
      dff  = ff*(0.6d0*t9i + 8.596*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.67 * exp(-1.835*t9i)
      drevdt   = rev*1.835*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_n14pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! n14(p,a)b11
       aa     = 1.0d0 + 0.0478*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0478

       zz     = 1.0d0/bb
       cc     = aa + 7.56e-03*t953*zz
       dcc    = 0.0478 + (fiveth*7.56e-3*t923 - 7.56e-3*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix*t9a56*zz

       dd      = 2.63e+16 * t9a56 * t9i32 * exp(-31.883/t9a13)
       ddd     = dd*(dt9a56/t9a56 - 1.5d0*t9i &
                 + 31.883/t9a13**2 * dt9a13)

       term    = dd * exp(-33.915*t9i)
       dtermdt = term*(ddd/dd + 33.915*t9i2)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.272 * dd
      drevdt   = 0.272 * ddd

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_c11pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,q1
      parameter        (q1 = 1.0d0/2.647129d0)


! c11(p,g)n12
      aa  = 4.24e+04 * t9i23 * exp(-13.658*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*13.658*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0  + 0.031*t913 + 3.11*t923 + 0.665*t9 &
            + 4.61*t943 + 2.50*t953
      dbb = oneth*0.031*t9i23 + twoth*3.11*t9i13 + 0.665 &
            + fourth*4.61*t913 + fiveth*2.50*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 8.84e+03 * t9i32 * exp(-7.021*t9i)
      ddd  = dd*(-1.5d0*t9i + 7.021*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.33e+10 * t932 * exp(-6.975*t9i)
      drevdt   = rev*(1.5d0*t9i + 6.975*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! c12(a,n)o15
      aa  = 2.48e7 * (1.0d0 + 0.188*t912 + 0.015*t9)
      daa = 2.48e7 * (0.5d0*0.188*t9i12 + 0.015)

      bb  = exp(-98.661*t9i)
      dbb = bb*98.661*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 1.41 * aa
      drrdt = den * 1.41 * daa  * 1.0d-9
      drrdd = 1.41 * aa

      return
      end





      subroutine rate_c13an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,q1
      parameter        (q1 = 1.0d0/1.648656d0)


! c13(a,n)o16
      aa  = 6.77e+15 * t9i23 * exp(-32.329*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*32.329*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.013*t913 + 2.04*t923 + 0.184*t9
      dbb = oneth*0.013*t9i23 + twoth*2.04*t9i13 + 0.184

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.82e+05 * t9i32 * exp(-9.373*t9i)
      ddd  = dd*(-1.5d0*t9i + 9.373*t9i2)

      ee   = 1.41e+06 * t9i32 * exp(-11.873*t9i)
      dee  = ee*(-1.5d0*t9i + 11.873*t9i2)

      ff   = 2.0e+09 * t9i32 * exp(-20.409*t9i)
      dff  = ff*(-1.5d0*t9i + 20.409*t9i2)

      gg   = 2.92e+09 * t9i32 * exp(-29.283*t9i)
      dgg  = gg*(-1.5d0*t9i + 29.283*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.79e+00 * exp(-25.711*t9i)
      drevdt = rev*25.711*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_c13pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! c13(p,n)n13
      aa  = 1.88e+08*(1.0d0 - 0.167*t912 + 0.037*t9)
      daa = 1.88e+08*(0.037 - 0.5d0*0.167*t9i12)

      bb  = exp(-34.846*t9i)
      dbb = bb*34.846*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 0.998 * aa
      drrdt = den * 0.998 * daa * 1.0d-9
      drrdd = 0.998 * aa

      return
      end






      subroutine rate_c14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/7.086244d0)


! c14(a,g)o18
      aa  = 1.528e+09 * t9i23 * exp(-32.513*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*32.513*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.0128*t913 - 0.869*t923 - 0.0779*t9 &
            + 0.321*t943 + 0.0732*t953
      dbb = oneth*0.0128*t9i23 - twoth*0.869*t9i13 - 0.0779 &
            + fourth*0.321*t913 + fiveth*0.0732*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.375e+08 * t9i2 * exp(-32.513*t9i13)
      ddd  = dd*(-2.0d0*t9i + oneth*32.513*t9i43)

      ee   = 9.29e-08 * t9i32 * exp(-2.048*t9i)
      dee  = ee*(-1.5d0*t9i + 2.048*t9i2)

      ff   = 2.77e+03 * t9i45 * exp(-9.876*t9i)
      dff  = ff*(-0.8d0*t9i + 9.876*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.42e+10 * t932 * exp(-72.262*t9i)
      drevdt   = rev*(1.5d0*t9i + 72.262*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_c14pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz


! c14(p,n)n14
      aa  = 7.19e+05*(1.0d0 + 0.361*t912 + 0.502*t9)
      daa = 7.19e+05*(0.5d0*0.361*t9i12 + 0.502)

      zz  = exp(-7.263*t9i)
      dzz = zz*7.263*t9i2

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc  = 3.34e+08 * t9i12 * exp(-12.246*t9i)
      dcc = cc*(-0.5d0*t9i + 12.246*t9i2)

      term    = bb + cc
      dtermdt = dbb + dcc

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      cc  = 3.34e+08 * t9i12 * exp(-4.983*t9i)
      dcc = cc*(-0.5d0*t9i + 4.983*t9i2)

      rr    = den * 0.333 * (aa + cc)
      drrdt = den * 0.333 * (daa + dcc) * 1.0d-9
      drrdd = 0.333 * (aa + cc)

      return
      end






      subroutine rate_c14pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/32.729841d0)


! c14(p,g)n14
      aa  = 6.80e+06 * t9i23 * exp(-13.741*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*13.741*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.03*t913 + 0.503*t923 + 0.107*t9 &
            + 0.213*t943 + 0.115*t953
      dbb = oneth*0.03*t9i23 + twoth*0.503*t9i13 + 0.107 &
            + fourth*0.213*t913 + fiveth*0.115*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 5.36e+03 * t9i32 * exp(-3.811*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.811*t9i2)

      ee   = 9.82e+04 * t9i13 * exp(-4.739*t9i)
      dee  = ee*(-oneth*t9i + 4.739*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.00e+09 * t932 * exp(-118.452*t9i)
      drevdt   = rev*(1.5d0*t9i + 118.452*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_o16pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! o16(p,a)n13
       aa     = 1.0d0 + 0.0776*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0776

       zz     = 1.0d0/bb
       cc     = aa + 0.0264*t953*zz
       dcc    = 0.0776 + (fiveth*0.0264*t923 - 0.0264*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       dd      = 1.88e+18 * t9a56 * t9i32 * exp(-35.829/t9a13)
       ddd     = dd*(dt9a56/t9a56 - 1.5d0*t9i &
                 + 35.829/t9a13**2 * dt9a13)

       term    = dd * exp(-60.561*t9i)
       dtermdt = term*(ddd/dd + 60.561*t9i2)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.172 * dd
      drevdt   = 0.172 * ddd

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_n14pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! n14(p,n)o14
      aa  = 6.74e+07 * (1.0d0 + 0.658*t912 + 0.379*t9)
      daa = 6.74e+07 * (0.5d0*0.658*t9i12 + 0.379)

      bb  = exp(-68.762*t9i)
      dbb = bb*68.762*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * 2.99 * aa
      drrdt = den * 2.99 * daa * 1.0d-9
      drrdd = 2.99 * aa

      return
      end





      subroutine rate_n14an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,zz,dzz,q1
      parameter        (q1 = 1.0d0/7.828804d0)


! n14(a,n)f17
      aa  = 5.24e9*(1.0d0 - 1.15*t912 + 0.365*t9)
      daa = 5.24e9*(0.365 - 0.5d0*1.15*t9i12)

      zz  = exp(-t92*q1)
      dzz = -zz*2.0d0*t9*q1

      bb  = aa * zz
      dbb = daa*zz + aa*dzz

      cc   = 3.28e10 * t9i32 * exp(-1.5766e1*t9i)
      dcc  = cc*(-1.5d0*t9i + 1.5766e1*t9i2)

      term     = (bb + cc) * exp(-54.942*t9i)
      dtermdt  = term*((dbb+dcc)/(bb+cc) + 54.942*t9i2)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 1.48 * (bb + cc)
      dtermdt  = 1.48 * (dbb + dcc)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_n15pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,aa,daa,bb,dbb


! n15(p,n)o15
      t9a = min(t9,10.0d0)
      aa  = 3.51e+08 * (1.0d0 + 0.452*t912 - 0.191*t9a)
      if (t9a .eq. 10.0) then
       daa = 3.51e+08 * 0.5d0*0.452*t9i12
      else
       daa = 3.51e+08 * (0.5d0*0.452*t9i12 - 0.191)
      end if

      bb  = exp(-41.032*t9i)
      dbb = bb*41.032*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * aa
      dtermdt = 0.998 * daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_n15an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! n15(a,n)f18
      aa  = 3.14e8 * (1.0d0 - 0.641*t912 + 0.108*t9)
      daa = 3.14e8 * (0.108 - 0.5d0*0.641*t9i12)

      bb  = exp(-74.479*t9i)
      dbb = bb*74.479*t9i2

      term     = aa * bb
      dtermdt  = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 2.0d0 * aa
      dtermdt  = 2.0d0 * daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_n15ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.379456d0)


! n15(a,g)f19
      aa  = 2.54e+10 * t9i23 * exp(-36.211*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*36.211*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.012*t913 + 1.69*t923 + 0.136*t9 &
            + 1.91*t943 + 0.391*t953
      dbb = oneth*0.012*t9i23 + twoth*1.69*t9i13 + 0.136 &
            + fourth*1.91*t913 + fiveth*0.391*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 9.83e-03 * t9i32 * exp(-4.232*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.232*t9i2)

      ee   = 1.52e+03 * t9 * exp(-9.747*t9i)
      dee  = ee*(t9i + 9.747*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.54e+10 * t932 * exp(-46.578*t9i)
      drevdt = rev*(1.5d0*t9i + 46.578*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.514089d0)


! o14(a,g)ne18
      aa  = 9.47e+08 * t9i23 * exp(-39.388*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*39.388*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.011*t913 + 1.974*t923 + 0.146*t9 &
            + 3.036*t943 + 0.572*t953
      dbb = oneth*0.011*t9i23 + twoth*1.974*t9i13 + 0.146 &
            + fourth*3.036*t913 + fiveth*0.572*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.16e-01 * t9i32 * exp(-11.733*t9i)
      ddd  = dd*(-1.5d0*t9i + 11.733*t9i2)

      ee   = 3.39e+01 * t9i32 * exp(-22.609*t9i)
      dee  = ee*(-1.5d0*t9i + 22.609*t9i2)

      ff   = 9.10e-03 * t95 * exp(-12.159*t9i)
      dff  = ff*(5.0d0*t9i + 12.159*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 5.42e+10 * t932 * exp(-59.328*t9i)
      drevdt = rev*(1.5d0*t9i + 59.328*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       ft9a,dft9a,fpt9a,dfpt9a,gt9x,dgt9x,zz


! o17(a,g)ne21
       aa    = 1.0d0 + 0.1646*t9
       zz    = 1.0d0/aa
       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.1646)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.786/t9a
       daa    = -aa*zz
       bb     = aa**3.51
       dbb    = 3.51*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       aa     = t9a/1.084
       bb     = aa**1.69
       dbb    = 1.69*bb/aa * dt9a/1.084
       fpt9a  = exp(-bb)
       dfpt9a = -fpt9a*dbb

       aa     = oneth*exp(-10.106*t9i)
       daa    = aa*10.106*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = 1.73e17 * fpt9a*zz
       daa    = (1.73e17*dfpt9a - aa*dgt9x)*zz

       bb     = 3.50e15 * ft9a*zz
       dbb    = (3.50e15*dft9a - bb*dgt9x)*zz

       term    = (aa+bb) * t9a56 * t9i32 * exp(-39.914/t9a13)
       dtermdt = term*((daa+dbb)/(aa+bb) + dt9a56/t9a56 &
                 - 1.5d0*t9i + 39.914/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.63e+10 * t932 * exp(-85.305*t9i)
      drevdt   = rev*(1.5d0*t9i + 85.305*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,gt9x,dgt9x,zz


! o17(a,n)ne20
       aa     = 1.0d0 + 0.0268*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0268

       zz     = 1.0d0/bb
       cc     = aa + 0.0232*t953*zz
       dcc    = 0.0268 + (fiveth*0.0232*t923 - 0.0232*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       dd     = oneth*exp(-10.106*t9i)
       gt9x   = 1.0d0 + dd
       dgt9x  = dd*10.106*t9i2

       term      = 1.03e+18/gt9x * t9a56 * t9i32 * exp(-39.914/t9a13)
       dtermdt   = term*(-dgt9x/gt9x + dt9a56/t9a56 &
                   - 1.5d0*t9i + 39.914/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.86e+01 * exp(-6.852*t9i)
      drevdt   = rev*6.852*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_o18ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.117649d0)


! o18(a,g)ne22
! giessen et al 1994 nuc phys a 567, 146 for t9 less than 0.3
! cf88 otherwise

      if (t9.lt.0.3) then
       aa   = 1.066d-41 * t9i32 * exp(-5.507d-01*t9i)
       daa  = aa*(-1.5d0*t9i + 5.507d-1*t9i2)

       bb   = 1.852d-13 * t9i32 * exp(-2.070*t9i)
       dbb  = bb*(-1.5d0*t9i + 2.070*t9i2)

       cc   = 1.431d-02 * t9i32 * exp(-4.462*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.462*t9i2)

       dd   = 2.055d-04 * t9i32 * exp(-5.374*t9i)
       ddd  = dd*(-1.5d0*t9i + 5.374*t9i2)

       ee   = 5.332d+00 * t9i32 * exp(-6.285*t9i)
       dee  = ee*(-1.5d0*t9i + 6.285*t9i2)

       ff   = 1.457d+00 * t9i32 * exp(-7.121*t9i)
       dff  = ff*(-1.5d0*t9i + 7.121*t9i2)

       gg   = 3.121d-02 * t9i32 * exp(-7.292*t9i)
       dgg  = gg*(-1.5d0*t9i + 7.292*t9i2)

       hh   = 6.23d+03 * t9 * exp(-16.987*t9i)
       dhh  = hh*(t9i + 16.987*t9i2)

       term    = aa + bb + cc + dd + ee + ff + gg + hh
       dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh


      else
       aa  = 1.82d+12 * t9i23 * exp(-40.057*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*40.057*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.01*t913 + 0.988*t923 + 0.072*t9 &
             + 3.17*t943 + 0.586*t953
       dbb = oneth*0.01*t9i23 + twoth*0.988*t9i13 + 0.072 &
            + fourth*3.17*t913 + fiveth*0.586*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 7.54 * t9i32 * exp(-6.228*t9i)
       ddd  = dd*(-1.5d0*t9i + 6.228*t9i2)

       ee   = 34.8 * t9i32 * exp(-7.301*t9i)
       dee  = ee*(-1.5d0*t9i + 7.301*t9i2)

       ff   = 6.23d+03 * t9 * exp(-16.987*t9i)
       dff  = ff*(t9i + 16.987*t9i2)

       gg   = theta * 1.0d-11 * t9i32 * exp(-1.994*t9i)
       dgg  = gg*(-1.5d0*t9i + 1.994*t9i2)

       term    = cc + dd + ee + ff + gg
       dtermdt = dcc + ddd + dee + dff + dgg
      end if

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.85d+10 * t932 * exp(-112.208*t9i)
      drevdt   = rev*(1.5d0*t9i + 112.208*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_o18an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,bb,dbb,cc,dcc,dd, &
                       ee,dee,ff,dff,gg,dgg,hh,dhh,ft9a,dft9a,gt9,dgt9, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,gt9i,zz


! o18(a,n)ne21
      aa     = 1.0d0 + 0.0483*t9
      bb     = aa**twoth
      dbb    = twoth*bb/aa*0.0483

      zz     = 1.0d0/bb
      cc     = aa + 0.00569*t953*zz
      dcc    = 0.0483 + (fiveth*0.00569*t923 - 0.00569*t953*zz*dbb)*zz

      zz     = 1.0d0/cc
      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*dcc)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix * t9a56*zz

      dd     = 5.0d0 * exp(-23.002*t9i)
      gt9    = 1.0d0 + dd
      gt9i   = 1.0d0/gt9
      dgt9   = dd*23.002*t9i2

      ee     = 0.431/t9a
      dee    = -ee*zz
      ff     = ee**3.89
      dff    = 3.89*ff/ee*dee
      ft9a   = exp(-ff)
      dft9a  = -ft9a*dff

      gg     = 7.22e+17 * ft9a*gt9i * t9a56 * t9i32 * exp(-40.056/t9a13)
      dgg    = gg*(-dff - gt9i*dgt9 + dt9a56/t9a56 - 1.5d0*t9i &
               + 40.056/t9a13**2 *dt9a13)

      hh     = 150.31 / gt9 * exp(-8.045*t9i)
      dhh    = hh*(-gt9i*dgt9 + 8.045*t9i2)

      term     = gg + hh
      dtermdt  = dgg + dhh

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

! must protect the 8.045*t9i from overflow, so write it this way

      gg     = 7.22e+17*gt9i * t9a56 * t9i32 &
               * exp(-ff - 40.056/t9a13 + 8.045*t9i)
      dgg    = gg*(-gt9i*dgt9 + dt9a56/t9a56 - 1.5d0*t9i &
                - dff + 40.056/t9a13**2*dt9a13 - 8.045*t9i2)

      hh     = 150.31 * gt9i
      dhh    = -hh*gt9i*dgt9

      term    = 0.784 * (gg + hh)
      dtermdt = 0.784 * (dgg + dhh)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_ne20pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc,dd,ddd, &
                       ee,dee,ff,dff,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       zz,t9b


! ne20(p,a)f17
       aa     = 1.0d0 + 0.0612*t9
       bb     = aa**twoth
       dbb    = twoth*bb/aa*0.0612

       zz     = 1.0d0/bb
       cc     = aa + 0.013*t953*zz
       dcc    = 0.0612 + (fiveth*0.013*t923 - 0.013*t953*zz*dbb)*zz

       zz     = 1.0d0/cc
       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*dcc)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       t9b     = min(t9,10.0d0)
       dd      = 5.31 + 0.544*t9b - 0.0523*t9b*t9b
       ddd     = 0.544 - 2.0d0*0.0523*t9b
       if (t9b .eq. 10.0) ddd     = 0.0d0

       ee      = 3.25e19 * dd * t9a56 * t9i32 * exp(-43.176/t9a13)
       dee     = ee*(ddd/dd + dt9a56/t9a56 - 1.5d0*t9i &
                 +   43.176/t9a13**2 * dt9a13)

       ff      = exp(-47.969*t9i)
       dff     = ff*47.969*t9i2

       term    = ee * ff
       dtermdt = dee*ff + ee*dff


! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0537 * ee
      drevdt   = 0.0537 * dee

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_f18pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg


! f18(p,g)ne19
! wiescher and kettner, apj 263, 891 1982
      aa  = 1.658e+7 * t9i23 * exp(-18.06*t9i13)
      daa = aa*(-twoth*t9i + oneth*18.06*t9i43)

      bb  = 4.604 + 0.106*t913 + 0.053*t923 + 0.009*t9 &
            - 0.036*t943 - 0.015*t953
      dbb = oneth*0.106*t9i23 + twoth*0.053*t9i13 + 0.009 &
            - fourth*0.036*t913 - fiveth*0.015*t923

! for temps greater than about t9 = 20, bb goes negative
      if (bb .le. 0.0) then
       bb = 0.0d0
       dbb = 0.0d0
      end if

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 4.55e-14 * t9i32* exp(-0.302*t9i)
      ddd  = dd*(-1.5d0*t9i + 0.302*t9i2)

      ee   = 327.0 * t9i32 * exp(-3.84*t9i)
      dee  = ee*(-1.5d0*t9i + 3.84*t9i2)

      ff   = 1.32e+04 * t9i32 * exp(-5.22*t9i)
      dff  = ff*(-1.5d0*t9i + 5.22*t9i2)

      gg   = 93.0 * t9i32 * exp(-4.29*t9i)
      dgg  = gg*(-1.5d0*t9i + 4.29*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.73e+10 * t932 * exp(-74.396*t9i)
      drevdt   = rev*(1.5d0*t9i + 74.396*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_f19pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,zz,q1
      parameter        (q1 = 1.0d0/0.173056d0)


! f19(p,g)ne20
      aa  = 6.04e+07 * t9i23 * exp(-18.113*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*18.113*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.023*t913 + 2.06*t923 + 0.332*t9 &
            + 3.16*t943 + 1.30*t953
      dbb = oneth*0.023*t9i23 + twoth*2.06*t9i13 + 0.332 &
            + fourth*3.16*t913 + fiveth*1.30*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 6.32e+02 * t9i32 * exp(-3.752*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.752*t9i2)

      ee   = 7.56e+04 * t9i27 * exp(-5.722*t9i)
      dee  = ee*(-twosev*t9i + 5.722*t9i2)

      ff   = 7.0*exp(-16.44*t9i)
      dff  = ff*16.44*t9i2

      gg   = 4.0 * exp(-2.09*t9i)
      dgg  = gg*2.09*t9i2

      hh   = 1.0d0 + ff + gg
      dhh  = dff + dgg

      zz      = 1.0d0/hh
      term    = (cc + dd + ee)*zz
      dtermdt = (dcc + ddd + dee - term*dhh)*zz

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.7e+10 * t932 * exp(-149.093*t9i)
      drevdt   = rev*(1.5d0*t9i + 149.093*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_f19pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! f19(p,n)ne19
      aa  = 1.27e+08 * (1.0d0 - 0.147*t912 + 0.069*t9)
      daa = 1.27e+08 * (0.069 - 0.5d0*0.147*t9i12)

      bb  = exp(-46.659*t9i)
      dbb = bb*46.659*t9i2

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * aa
      dtermdt = 0.998 * daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_f19ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,q1
      parameter        (q1 = 1.0d0/0.405769)


! f19(a,p)ne22
      aa  = 4.50e+18 * t9i23 * exp(-43.467*t9i13  - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*43.467*t9i43 - 2.0d0*t9*q1)

      bb   = 7.98e+04 * t932 * exp(-12.760*t9i)
      dbb  = 1.5d0*bb*t9i + bb*12.760*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.36 * exp(-19.439*t9i)
      drevdt   = rev*19.439*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_na22na(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa, &
                       t9b,t9b2,t9b3


! na22(n,a)f19
      t9b  = min(t9,10.0d0)
      t9b2 = t9b*t9b
      t9b3 = t9b2*t9b
      aa  = 1.0d0 + 0.8955*t9b - 0.05645*t9b2 + 7.302e-04*t9b3
      daa = 0.8955 - 2.0d0*0.05645*t9b + 3.0d0*7.302e-4*t9b2
      if (t9b .eq. 10.0) daa = 0.0d0

      term     = 1.21e6 * exp(aa)
      dtermdt  = term*daa

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.10 * exp(-22.620*t9i)
      drevdt   = rev*22.620*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne20pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ff,gg,dgg,zz


! ne20(p,g)na21
      aa  = 9.55e+06 * exp(-19.447*t9i13)
      daa = aa*oneth*19.447*t9i43

      bb  = 1.0d0 + 0.0127*t9i23
      dbb = -twoth*0.0127*t9i53

      cc  = t92 * bb * bb
      dcc = 2.0d0*cc*t9i + 2.0d0*t92*bb*dbb

      zz  = 1.0d0/cc
      dd  = aa*zz
      ddd = (daa - dd*dcc)*zz

      aa  = 2.05e+08 * t9i23 * exp(-19.447*t9i13)
      daa = aa*(-twoth*t9i + oneth*19.447*t9i43)

      bb  = sqrt (t9/0.21)
      dbb = 0.5d0/(bb * 0.21)

      cc  = 2.67 * exp(-bb)
      dcc = -cc*dbb

      ff  = 1.0d0 + cc

      gg  = aa*ff
      dgg = daa*ff + aa*dcc


      aa  = 18.0 * t9i32 * exp(-4.242*t9i)
      daa = aa*(-1.5d0*t9i + 4.242*t9i2)

      bb  = 10.2 * t9i32 * exp(-4.607*t9i)
      dbb = bb*(-1.5d0*t9i + 4.607*t9i2)

      cc  = 3.6e+04 * t9i14 * exp(-11.249*t9i)
      dcc = cc*(-0.25d0*t9i + 11.249*t9i2)

      term    = dd + gg + aa + bb + cc
      dtermdt = ddd + dgg + daa + dbb + dcc


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.63e+09 * t932 * exp(-28.216*t9i)
      drevdt   = rev*(1.5d0*t9i + 28.216*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na23pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta,q1,q2
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.0169d0, &
                        q2    = 1.0d0/0.017161d0)


! na23(p,a)ne20
! el eid & champagne 1995
      if (t9 .le. 2.0) then
       aa  = 1.26d+10 * t9i23 * exp(-20.758*t9i13 - t92*q1)
       daa = -twoth*aa*t9i + aa*(oneth*20.758*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0  + 0.02*t913 - 13.8*t923 - 1.93*t9 &
             + 234.0*t943 + 83.6*t953
       dbb = oneth*0.02*t9i23 - twoth*13.8*t9i13 - 1.93 &
             + fourth*234.0*t913 + fiveth*83.6*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 4.38 * t9i32 * exp(-1.979*t9i)
       ddd  = -1.5d0*dd*t9i + dd*1.979*t9i2

       ee   = 6.50d+06 * (t9**(-1.366d0)) * exp(-6.490*t9i)
       dee  = -1.366d0*ee*t9i + ee*6.490*t9i2

       ff   = 1.19d+08 * (t9**(-1.055d0)) * exp(-11.411*t9i)
       dff  = -1.055d0*ff*t9i + ff*11.411*t9i2

       gg   = theta * 9.91d-14 * t9i32 * exp(-0.418*t9i)
       dgg  = -1.5d0*gg*t9i + gg*0.418*t9i2

       term    = cc + dd + ee + ff + gg
       dtermdt = dcc + ddd + dee + dff + dgg



! cf88 + one term from gorres, wiesher & rolfs 1989, apj 343, 365
      else
       aa  = 8.56d+09 * t9i23 * exp(-20.766*t9i13 - t92*q2)
       daa = -twoth*aa*t9i + aa*(oneth*20.766*t9i43 - 2.0d0*t9*q2)

       bb  = 1.0d0  + 0.02*t913 + 8.21*t923 + 1.15*t9 &
             + 44.36*t943 + 15.84*t953
       dbb = oneth*0.02*t9i23 + twoth*8.21*t9i13 + 1.15 &
             + fourth*44.36*t913 + fiveth*15.84*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = 4.02 * t9i32 * exp(-1.99*t9i)
       ddd  = -1.5d0*dd*t9i + dd*1.99*t9i2

       ee   = 1.18d+04 * t9i54 * exp(-3.148*t9i)
       dee  = -1.25d0*ee*t9i + ee*3.148*t9i2

       ff   = 8.59d+05 * t943 * exp(-4.375*t9i)
       dff  = fourth*ff*t9i + ff*4.375*t9i2

       gg   = theta * 3.06d-12 * t9i32 * exp(-0.447*t9i)
       dgg  = -1.5d0*gg*t9i + gg*0.447*t9i2

       hh   = theta * 0.820 * t9i32 * exp(-1.601*t9i)
       dhh  = -1.5d0*hh*t9i + hh*1.601*t9i2

       term    = cc + dd + ee + ff + gg + hh
       dtermdt = dcc + ddd + dee + dff + dgg + dhh
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.25 * exp(-27.606*t9i)
      drevdt   = rev*27.606*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_ne20ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


! ne20(n,g)ne21
! wm88 Apj 239, 943; fit over range of experimental data, constant otherwise

      if (t9 .lt. 5.8025d-2) then
       term    = 5.449d+03
       dtermdt = 0.0d0

      else if (t9 .gt. 1.1605) then
       term    = 6.977d+04
       dtermdt = 0.0d0

      else if (t9 .ge. 5.8025d-2 .and. t9 .le. 2.9012d-1) then
       term    = 4.7219d+3 + 2.5248d+4*t9 - 2.7448d+5*t92 &
                 + 9.2848d+5*t93
       dtermdt = 2.5248d+4 - 2.0d0*2.7448d+5*t9 &
                 + 3.0d0*9.2848d+5*t92

      else

       aa  = 1.802d+04 * (t9/0.348)**4.43
       daa = 4.43 * aa * t9i

       bb  = -5.931 * (t9-0.348) + 1.268 * (t9-0.348)**2
       dbb = -5.931 + 2.0d0*1.268*(t9 - 0.348)

       cc  = exp(bb)
       dcc = cc*dbb

       term    = aa * cc
       dtermdt = daa*cc + aa*dcc

      end if

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      =  4.650d+09 * t932 * exp(-78.46*t9i)
      drevdt   = rev*(1.5d0*t9i + 78.46*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne21pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,xx,dxx,zz, &
                       theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.003364d0)


! ne21(p,g)na22

! el eid & champagne 1995

      if (t9.le.2.0) then
       aa  = 3.4d+08 * t9i23 * exp(-19.41*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.41*t9i43)

       bb  = (16.7*t9 - 1.0)**2
       dbb = 2.0d0*(16.7*t9 - 1.0)*16.7

       cc  = 0.56 * exp(-bb)
       dcc = -cc*dbb

       dd  = 1.0d0 + cc
       ddd = dcc

       ee   = aa * dd
       dee  = daa*dd + aa*ddd

       ff   = 6.12 * t9i32 * exp(-1.403*t9i)
       dff  = ff*(-1.5d0*t9i + 1.403*t9i2)

       gg   = 1.35d+04 * t9i32 * exp(-3.008*t9i)
       dgg  = gg*(-1.5d0*t9i + 3.008*t9i2)

       aa   = t9**0.67
       daa  = 0.67*aa*t9i
       zz   = 1.0d0/aa

       hh   = 3.12d+06 * t9**(-0.72) * exp(-8.268*zz)
       dhh  = hh*(-0.72d0*t9i + 8.268*zz*zz*daa)

       xx   = theta * 1.1d-03 * t9i32 * exp(-1.114*t9i)
       dxx  = xx*(-1.5d0*t9i + 1.114*t9i2)

       term    = ee + ff + gg + hh + xx
       dtermdt = dee + dff + dgg + dhh + dxx


! cf88
      else

       aa  = theta * 2.95d+08 * t9i23 * exp(-19.462*t9i13 -t92*q1)
       daa = aa*(-twoth*t9i + oneth*19.462*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.021*t913 + 13.29*t923 + 1.99*t9 &
             + 124.1*t943 + 47.29*t953
       dbb = oneth*0.021*t9i23 + twoth*13.29*t9i13 + 1.99 &
             + fourth*124.1*t913 + fiveth*47.29*t923

       cc   = aa * bb
       dcc  = daa*bb + aa*dbb

       dd   = theta * 7.80d-01 * t9i32 * exp(-1.085*t9i)
       ddd  = dd*(-1.5d0*t9i + 1.085*t9i2)

       ee   = 4.37d+08 * t9i23 * exp(-19.462*t9i13)
       dee  = ee*(-twoth*t9i + oneth*19.462*t9i43)

       ff   = 5.85 * t9i32 * exp(-1.399*t9i)
       dff  = ff*(-1.5d0*t9i + 1.399*t9i2)

       gg   = 1.29d+04 * t9i32 * exp(-3.009*t9i)
       dgg  = gg*(-1.5d0*t9i + 3.009*t9i2)

       hh   = 3.15d+05 * t9i35 * exp(-5.763*t9i)
       dhh  = hh*(-0.6d0*t9i + 5.763*t9i2)

       term    = cc + dd + ee + ff + gg + hh
       dtermdt = dcc + ddd + dee + dff + dgg + dhh
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.06d+10 * t932 * exp(-78.194*t9i)
      drevdt   = rev*(1.5d0*t9i + 78.194*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ne21ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,gg,dgg,hh,dhh,zz, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! ne21(a,g)mg25
       aa    = 1.0d0 + 0.0537*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0537)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 8.72e-03*t9 - 6.87e-04*t92 + 2.15e-05*t93
       daa    = 8.72e-3 - 2.0d0*6.87e-4*t9 + 3.0d0*2.15e-5*t92

       bb     = 1.52e-04 * exp(-46.90*t9i13*aa)
       dbb    = bb*46.90*(oneth*t9i43*aa - t9i13*daa)

       cc     = 1.5*exp(-4.068*t9i)
       dcc    =  cc*4.068*t9i2

       gg     = 2.0 * exp(-20.258*t9i)
       dgg    = gg*20.258*t9i2

       hh     = 1.0d0 + cc + gg
       dhh    = dcc + dgg

       zz     = 1.0d0/hh
       dd     = bb*zz
       ddd    = (dbb - dd*dhh)*zz

       aa     = 4.94e+19 * t9a56 * t9i32 * exp(-46.89/t9a13)
       daa    = aa*(dt9a56/t9a56 - 1.5d0*t9i &
                    + 46.89/t9a13**2 * dt9a13)

       bb     =  2.66e+07 * t9i32 * exp(-22.049*t9i)
       dbb    = bb*(-1.5d0*t9i + 22.049*t9i2)

       cc     = aa + bb
       dcc    = daa + dbb

       term    = dd * cc
       dtermdt = ddd*cc + dd*dcc


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 4.06e+10 * t932 * exp(-114.676*t9i)
      drevdt = rev*(1.5d0*t9i + 114.676*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_ne21an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,zz, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! ne21(a,n)mg24
       aa    = 1.0d0 + 0.0537*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0537)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 4.94e+19 * t9a56 * t9i32 * exp(-46.89/t9a13)
       daa    = aa*(dt9a56/t9a56 - 1.5d0*t9i &
                    + 46.89/t9a13**2 * dt9a13)

       bb     =  2.66e+07 * t9i32 * exp(-22.049*t9i)
       dbb    = bb*(-1.5d0*t9i + 22.049*t9i2)

       cc     = 2.0d0*exp(-20.258*t9i)
       dcc    = cc*20.258*t9i2

       dd     = 1.5*exp(-4.068*t9i)
       ddd    = dd*4.068*t9i2

       ee     = 1.0d0 + cc + dd
       dee    = dcc + ddd

       zz      = 1.0d0/ee
       term    = (aa + bb)*zz
       dtermdt = (daa + dbb - term*dee)*zz


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 12.9 * exp(-29.606*t9i)
      drevdt   = rev*29.606*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_ne22pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,theta
      parameter        (theta = 0.1d0)


! ne22(p,g)na23

! el eid & champagne 1995

      if (t9.le.2.0) then
       aa  = 1.05d+09 * t9i23 * exp(-19.431*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.431*t9i43)

       bb  = 1.24d-09 * t9i32 * exp(-0.414*t9i)
       dbb = bb*(-1.5d0*t9i + 0.414*t9i2)

       cc  = 2.90d-02 * t9i32 * exp(-1.752*t9i)
       dcc = cc*(-1.5d0*t9i + 1.752*t9i2)

       dd  = 9.30d+04 * t9**(-1.174) * exp(-5.100*t9i)
       ddd = dd*(-1.174*t9i + 5.100*t9i2)

       ee   = 5.71d+05 * t9**(0.249) * exp(-7.117*t9i)
       dee  = ee*(0.249*t9i + 7.117*t9i2)

       ff   = theta * 3.25d-04 * t9i32 * exp(-0.789*t9i)
       dff  = ff*(-1.5d0*t9i + 0.789*t9i2)

       gg   = theta * 0.10 * t9i32 * exp(-1.161*t9i)
       dgg  = gg*(-1.5d0*t9i + 1.161*t9i2)

       term    = aa + bb + cc + dd + ee + ff + gg
       dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


! cf88
      else

       aa  = 1.15d+09 * t9i23 * exp(-19.475*t9i13)
       daa = aa*(-twoth*t9i + oneth*19.475*t9i43)

       bb  = 9.77d-12 * t9i32 * exp(-0.348*t9i)
       dbb = bb*(-1.5d0*t9i + 0.348*t9i2)

       cc   = 8.96d+03 * t9i32 * exp(-4.84*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.84*t9i2)

       dd   = 6.52d+04 * t9i32 * exp(-5.319*t9i)
       ddd  = dd*(-1.5d0*t9i + 5.319*t9i2)

       ee   = 7.97d+05 * t9i12 * exp(-7.418*t9i)
       dee  = ee*(-0.5d0*t9i + 7.418*t9i2)

       ff   = theta * 1.63d-01 * t9i32 * exp(-1.775*t9i)
       dff  = ff*(-1.5d0*t9i + 1.775*t9i2)

       term    = aa + bb + cc + dd + ee + ff
       dtermdt = daa + dbb + dcc + ddd + dee + dff

      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.67d+09 * t932 * exp(-102.048*t9i)
      drevdt   = rev*(1.5d0*t9i + 102.048*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne22ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,res1,dres1, &
                       ft9a,dft9a,fpt9a,dfpt9a,gt9x,dgt9x, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       rdmass,res2,zz
      parameter        (rdmass = 22.0d0*4.0d0/26.0d0, &
                        res2   = -11.604d0 * 22.0d0/26.0d0)


! ne22(a,g)mg26
! kappeler 1994 apj 437, 396

      if (t9 .lt. 1.25) then

       res1 = 1.54d-01*(t9*rdmass)**(-1.5)
       dres1 = -1.5d0 * res1 * t9i

       aa    = 1.7d-36 * res1 * exp(res2*t9i*0.097)
       daa   = aa/res1*dres1 - aa*res2*0.097*t9i2

       bb    = 1.5d-7 * res1 * exp(res2*t9i*0.400)
       dbb   = bb/res1*dres1 - bb * res2 * 0.400 * t9i2

       cc    = 0.5 * res1 * 3.7d-2 * exp(res2*t9i*0.633)
       dcc   = cc/res1*dres1 - cc*res2*0.633*t9i2

       dd    = res1 * 3.6d+1 * exp(res2*t9i*0.828)
       ddd   = dd/res1*dres1 - dd*res2*0.828*t9i2

       term    = aa + bb + cc + dd
       dtermdt = daa + dbb + dcc + ddd


! cf88
      else
       aa    = 1.0d0 + 0.0548*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0548)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.197/t9a
       daa    = -aa*zz
       bb     = aa**4.82
       dbb    = 4.82*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       aa     = t9a/0.249
       bb     = aa**2.31
       dbb    = 2.31*bb/aa * dt9a/0.249
       fpt9a  = exp(-bb)
       dfpt9a = -fpt9a*dbb

       aa     = 5.0d0*exp(-14.791*t9i)
       daa    = aa*14.791*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = 4.16e19 * fpt9a*zz
       daa    = (4.16e19*dfpt9a - aa*dgt9x)*zz

       bb     = 2.08e16 * ft9a*zz
       dbb    = (2.08e16*dft9a - bb*dgt9x)*zz

       term    = (aa+bb) * t9a56 * t9i32 * exp(-47.004/t9a13)
       dtermdt = term*((daa+dbb)/(aa+bb) + dt9a56/t9a56 &
                       - 1.5d0*t9i + 47.004/t9a13**2 * dt9a13)
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 6.15d+10 * t932 * exp(-123.151*t9i)
      drevdt = rev*(1.5d0*t9i + 123.151*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na22np(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! na22(n,p)ne22
      aa  = 1.0d0 - 3.037e-02*t9 + 8.380e-03*t92 - 7.101e-04*t93
      daa =  -3.037e-02 + 2.0d0*8.380e-03*t9 - 3.0d0*7.101e-04*t92

      term    = 1.24e+08 * exp(aa)
      dtermdt = term*daa

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev     = 7.01*exp(-42.059*t9i)
      drevdt  = rev*42.059*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne22an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,gg,dgg,ft9a,dft9a,gt9x,dgt9x, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,res1,res2, &
                       zz
      parameter        (res1 = 2.4731857150793075d-2, &
                        res2 = -9.8187694549560547d0)

! note: res1=1.54d-1*(88./26.)**(-1.5)   res2=-11.604*(22./26.)

! ne22(a,n)mg25
! kappeler 1994 apj 437, 396 ; wiescher suggest only 828 kev, ignore 633 kev

      if (t9 .lt. 0.6) then
       term = res1*1.64d+02 * t9i32 * exp(t9i*0.828*res2)
       dtermdt = -1.5d0*term*t9i - term*res2*0.828*t9i2

! cf88
      else
       aa    = 1.0d0 + 0.0548*t9
       zz    = 1.0d0/aa

       t9a   = t9*zz
       dt9a  = (1.0d0 - t9a*0.0548)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a56  = t9a**fivsix
       dt9a56 = fivsix * t9a56*zz

       aa     = 0.197/t9a
       daa    = -aa*zz
       bb     = aa**4.82
       dbb    = 4.82*bb/aa * daa
       ft9a   = exp(-bb)
       dft9a  = -ft9a*dbb

       gg     = bb
       dgg    = dbb

       aa     = 5.0d0*exp(-14.791*t9i)
       daa    = aa*14.791*t9i2
       gt9x   = 1.0d0 + aa
       dgt9x  = daa

       zz     = 1.0d0/gt9x
       aa     = ft9a*zz
       daa    = (dft9a - aa*dgt9x)*zz

       bb     = 4.16e+19 * t9a56 * t9i32 * exp(-47.004/t9a13)
       dbb    = bb*(dt9a56/t9a56 - 1.5d0*t9i &
                    + 47.004/t9a13**2 * dt9a13)

       cc     = aa*bb
       dcc    = daa*bb + aa*dbb

       dd      = 1.44e-04*zz * exp(-5.577*t9i)
       ddd     = -dd*zz*dgt9x + dd*5.577*t9i2

       term    = cc + dd
       dtermdt = dcc + ddd
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 7.833d-5
      drevdt = 0.0d0
      if (t9 .gt. 0.008) then
       rev    = 0.544 * exp(5.577*t9i)
       drevdt = -rev*5.577*t9i2
      end if

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_na21pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,q1
      parameter        (q1 = 1.0d0/0.133956d0)


! na21(p,g)mg22
      aa  = 1.41e+05 * t9i23 * exp(-20.739*t9i13 -  t92*q1)
      daa = aa*(-twoth*t9i + oneth*20.739*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.020*t913 + 4.741*t923 + 0.667*t9 &
            + 16.380*t943 + 5.858*t953
      dbb = oneth*0.020*t9i23 + twoth*4.741*t9i13 + 0.667 &
            + fourth*16.380*t913 + fiveth*5.858*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 6.72e+02 * t9i34 * exp(-2.436*t9i)
      ddd  = dd*(-0.75d0*t9i + 2.436*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.44e+10 * t932 * exp(-63.790*t9i)
      drevdt   = rev*(1.5d0*t9i + 63.790*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg24pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,gg,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! mg24(p,a)na21
       aa     = 1.0d0 + 0.127*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.127)*zz

       zz      = dt9a/t9a
       t9a13   = t9a**oneth
       dt9a13  = oneth*t9a13*zz

       t9a56   = t9a**fivsix
       dt9a56  = fivsix * t9a56*zz

       gg      = min(t9,12.0d0)
       aa      = 4.43 + 3.31*gg - 0.229*gg*gg
       daa     = 3.31 - 2.0d0*0.229*gg
       if (gg .eq. 12.0) daa = 0.0d0

       bb      = 1.81e21 * t9a56 * t9i32 * exp(-49.967/t9a13)
       dbb     = bb*(dt9a56/t9a56 - 1.5d0*t9i &
                     + 49.967/t9a13**2 * dt9a13)

       cc      = aa*bb
       dcc     = daa*bb + aa*dbb

       dd      = exp(-79.843*t9i)
       ddd     = dd*79.843*t9i2

       term    = cc * dd
       dtermdt = dcc*dd + cc*ddd

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0771 * cc
      drevdt   = 0.0771 * dcc

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end






      subroutine rate_na22pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! na22(p,g)mg23
      aa  = 9.63e-05 * t932 * exp(-0.517*t9i)
      daa = aa*(1.5d0*t9i + 0.517*t9i2)

      bb  = 2.51e+04 * t9 * exp(-2.013*t9i)
      dbb = bb*(t9i + 2.013*t9i2)

      term    = aa + bb
      dtermdt = daa + dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.27e+10 * t932 * exp(-87.933*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.933*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_na23pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,hhi,xx,dxx, &
                       theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.088209d0)


! na23(p,g)mg24

! el eid & champagne 1995
      if (t9 .le. 2.0) then
       aa  = 2.47d+09 * t9i23 * exp(-20.758*t9i13)
       daa = aa*(-twoth*t9i + oneth*20.758*t9i43)

       bb  = 9.19d+01 * t9i32 * exp(-2.789*t9i)
       dbb = bb*(-1.5d0*t9i + 2.789*t9i2)

       cc  = 1.72d+04 * t9i32 * exp(-3.433*t9i)
       dcc = cc*(-1.5d0*t9i + 3.433*t9i2)

       dd  = 3.44d+04 * t9**0.323 * exp(-5.219*t9i)
       ddd = dd*(0.323*t9i + 5.219*t9i2)

       ee   = theta * 2.34d-04 * t9i32 * exp(-1.590*t9i)
       dee  = ee*(-1.5d0*t9i + 1.590*t9i2)

       term    = aa + bb + cc + dd + ee
       dtermdt = daa + dbb + dcc + ddd + dee


! cf88 + gorres, wiesher & rolfs 1989, apj 343, 365
      else

       aa  = 2.93d+08 * t9i23 * exp(-20.766*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*20.766*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.02*t913 + 1.61*t923 + 0.226*t9 &
             + 4.94*t943 + 1.76*t953
       dbb = oneth*0.02*t9i23 + twoth*1.61*t9i13 + 0.226 &
            + fourth*4.94*t913 + fiveth*1.76*t923

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 9.34d+01 * t9i32 * exp(-2.789*t9i)
       dcc  = cc*(-1.5d0*t9i + 2.789*t9i2)

       dd   = 1.89d+04 * t9i32 * exp(-3.434*t9i)
       ddd  = dd*(-1.5d0*t9i + 3.434*t9i2)

       ee   = 5.1d+04 * t915 * exp(-5.51*t9i)
       dee  = ee*(0.2d0*t9i + 5.51*t9i2)

       ff   = theta * 0.820 * t9i32 * exp(-1.601*t9i)
       dff  = ff*(-1.5d0*t9i + 1.601*t9i2)

       gg   = 1.5 * exp(-5.105*t9i)
       dgg  = gg*5.105*t9i2

       hh   = 1.0d0 + gg
       hhi  = 1.0d0/hh

       term    = (xx + cc + dd + ee + ff) * hhi
       dtermdt = (dxx + dcc + ddd + dee + dff - term*dgg)*hhi

      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.49d+10 * t932 * exp(-135.665*t9i)
      drevdt   = rev*(1.5d0*t9i + 135.665*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_na23pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,bb,dbb,cc,dcc, &
                       t9a,dt9a,t9a32,dt9a32,zz


! na23(p,n)mg24
      aa  = 1.0d0 + 0.141*t9
      zz  = 1.0d0/aa

      t9a = t9*zz
      dt9a = (1.0d0 - t9a*0.141)*zz

      aa    = sqrt(t9a)
      t9a32 = t9a * aa
      dt9a32 = 1.5d0 * aa * dt9a

      bb   = 9.29d8 * (1.0d0 - 0.881d0 * t9a32 * t9i32)
      dbb  = -9.29d8 * 0.881d0 * t9i32*(dt9a32 - 1.5d0*t9a32*t9i)

      cc   = exp(-56.173*t9i)
      dcc  = cc*56.173*t9i2

      term    = bb * cc
      dtermdt = dbb*cc + bb*dcc

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term    = 0.998 * bb
      dtermdt = 0.998 * dbb

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end






      subroutine rate_mg24pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,ggi


! mg24(p,g)al25
      aa  = 5.60e+08 * t9i23 * exp(-22.019*t9i13)
      daa = aa*(-twoth*t9i + oneth*22.019*t9i43)

      bb  = 1.0d0 + 0.019*t913 - 0.173*t923 - 0.023*t9
      dbb = oneth*0.019*t9i23 - twoth*0.173*t9i13 - 0.023

! stop negative rates above t9 = 10
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd   = 1.48e+03 * t9i32 * exp(-2.484*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.484*t9i2)

      ee   = 4.00e+03 * exp(-4.180*t9i)
      dee  = ee*4.180*t9i2

      ff   = 5.0 * exp(-15.882*t9i)
      dff  = ff*15.882*t9i2

      gg   = 1.0d0 + ff
      ggi  = 1.0d0/gg

      term    = (cc + dd + ee) * ggi
      dtermdt = (dcc + ddd + dee - term*dff)*ggi


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.13e+09 * t932 * exp(-26.358*t9i)
      drevdt   = rev*(1.5d0*t9i + 26.358*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_al27pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,theta
      parameter        (theta = 0.1d0)


! al27(p,a)mg24
! champagne 1996
      aa  = 4.71d+05 * t9i23 * exp(-23.25*t9i13 - 3.57*t92)
      daa = -twoth*aa*t9i + aa*(oneth*23.25*t9i43 - 2.0d0*3.57*t9)

      bb  = 1.0d0 + 0.018*t913 - 7.29*t923 - 0.914*t9 &
            + 77.2*t943 + 24.6*t953
      dbb = oneth*0.018*t9i23 - twoth*7.29*t9i13 - 0.914 &
            + fourth*77.2*t913 + fiveth*24.6*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 2.23d+04 * (t9**3.989) * exp(-2.148 * t9**(-1.293))
      ddd = 3.989*dd*t9i + 1.293*dd*2.148*t9**(-2.293)

      ee   = 0.17 * 1.29d-09 * t9i32 * exp(-0.836*t9i)
      dee  = -1.5d0*ee*t9i + ee*0.836*t9i2

      ff   = theta * 2.73d-03 * t9i32 * exp(-2.269*t9i)
      dff  = -1.5d0*ff*t9i + ff*2.269*t9i2

      gg   = theta * 2.60d-02 * t9i32 * exp(-2.492*t9i)
      dgg  = -1.5d0*gg*t9i + gg*2.492*t9i2

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.81*exp(-18.572*t9i)
      drevdt   = rev*18.572*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_mg25pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,q1,q2
      parameter        (q1 = 1.0d0/0.0036d0, &
                        q2 = 1.0d0/169.0d0)


! mg25(p,g)al26
      aa  = 3.57e+09 * t9i23 * exp(-22.031*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*22.031*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.019*t913 + 7.669*t923 + 1.015*t9 &
            + 167.4*t943 + 56.35*t953
      dbb = oneth*0.019*t9i23 + twoth*7.669*t9i13 + 1.015 &
            + fourth*167.4*t913 + fiveth*56.35*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.07e-13 * t9i32 * exp(-0.435*t9i)
      ddd = dd*(-1.5d0*t9i + 0.435*t9i2)

      ee   = 1.94e-07 * t9i32 * exp(-0.673*t9i)
      dee  = ee*(-1.5d0*t9i + 0.673*t9i2)

      ff   = 3.15e-05 * t9**(-3.40)* exp(-1.342*t9i - t92*q2)
      dff  = ff*(-3.40d0*t9i + 1.342*t9i2 - 2.0d0*t9*q2)

      gg   = 1.77e+04 * t958 * exp(-3.049*t9i - t92*q2)
      dgg  = gg*(0.625*t9i + 3.049*t9i2 - 2.0d0*t9*q2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.03e+10 * t932 * exp(-73.183*t9i)
      drevdt   = rev*(1.5d0*t9i + 73.183*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_mg25ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc


! mg25(a,p)si28
      aa  = -23.271*t9i13 + 6.46*t9 - 2.39*t92 + 0.506*t93 &
            - 6.04e-2*t94 + 3.75e-3*t95 - 9.38e-5*t96

      daa = oneth*23.271*t9i43 + 6.46 - 2.0d0*2.39*t9 + 3.0d0*0.506*t92 &
            - 4.0d0*6.04e-2*t93 + 5.0d0*3.75e-3*t94 - 6.0d0*9.38e-5*t95

      bb  = 3.23e8 * t9i23 * exp(aa)
      dbb  = -twoth*bb*t9i + bb*daa

! dbb/bb
      cc   = -twoth*t9i + daa

      term    = bb * exp(-13.995*t9i)
      dtermdt = term*cc + term*13.995*t9i2

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.86 * bb
      drevdt   = 2.86 * dbb

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_mg25ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13, &
                       t9a56,dt9a56,zz


! mg25(a,g)si29
      aa    = 1.0d0 + 0.0630*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0630)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = oneth*10.0d0*exp(-13.180*t9i)
      daa    = aa*13.180*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      cc     = 3.59e+20 * bb * t9a56 * t9i32 * exp(-53.41/t9a13)
      dcc    = cc*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
                   + 53.41/t9a13**2 * dt9a13)

      dd     = 0.0156*t9 - 1.79e-03*t92 + 9.08e-05*t93
      ddd    = 0.0156 - 2.0d0*1.79e-03*t9 + 3.0d0*9.08e-05*t92

      ee     = 5.87e-04*exp(-53.42*t9i13*dd)
      dee    = ee*53.42*(oneth*t9i43*dd - t9i13*ddd)

      term    = cc * ee
      dtermdt = dcc*ee + cc*dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 1.90e+11 * t932 * exp(-129.128*t9i)
      drevdt = rev*(1.5d0*t9i + 129.128*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg25an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,zz, &
                       gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! mg25(a,n)si28
      aa    = 1.0d0 + 0.0630*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0630)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = oneth*10.0d0*exp(-13.180*t9i)
      daa    = aa*13.180*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      term    = 3.59e+20 * bb * t9a56 * t9i32 * exp(-53.41/t9a13)
      dtermdt = term*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
                + 53.41/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 20.0*exp(-30.792*t9i)
      drevdt   = rev*30.792*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_mg26pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,theta
      parameter        (theta = 0.1d0)


! mg26(p,g)al27
! champagne 1996
      aa  = 8.54d-12 * t9i32 * exp(-0.605*t9i)
      daa = aa*(-1.5d0*t9i + 0.605*t9i2)

      bb  = 2.75d-06 * t9i32 * exp(-1.219*t9i)
      dbb = bb*(-1.5d0*t9i + 1.219*t9i2)

      cc  = 1.30d-02 * t9i32 * exp(-1.728*t9i)
      dcc = cc*(-1.5d0*t9i + 1.728*t9i2)

      dd  = 8.06d+00 * t9i32 * exp(-2.537*t9i)
      ddd = dd*(-1.5d0*t9i + 2.537*t9i2)

      ee   = 1.45d+03 * t9i32 * exp(-3.266*t9i)
      dee  = ee*(-1.5d0*t9i + 3.266*t9i2)

      ff   = 4.03d+04 * t9i32 * exp(-3.784*t9i)
      dff  = ff*(-1.5d0*t9i + 3.784*t9i2)

      gg   = 8.82d+04 * t9**(-0.21) * exp(-4.194*t9i)
      dgg  = gg*(-0.21*t9i + 4.194*t9i2)

      hh   = theta * 1.93d-05 * t9i32 * exp(-1.044*t9i)
      dhh  = hh*(-1.5d0*t9i + 1.044*t9i2)

      term    = aa + bb + cc + dd + ee + ff + gg + hh
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.14d+09 * t932 * exp(-95.99*t9i)
      drevdt   = rev*(1.5d0*t9i + 95.99*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg26ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13, &
                       t9a56,dt9a56,zz


! mg26(a,g)si30
      aa    = 1.0d0 + 0.0628*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0628)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = 5.0d0*exp(-20.990*t9i)
      daa    = aa*20.990*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      cc     = 2.93e+20 * bb * t9a56 * t9i32 * exp(-53.505/t9a13)
      dcc    = cc*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
                   + 53.505/t9a13**2 * dt9a13)

      dd     = 0.0751*t9 - 0.0105*t92 + 5.57e-04*t93
      ddd    = 0.0751 - 2.0d0*0.0105*t9 + 3.0d0*5.57e-04*t92

      ee     = 4.55e-2 * exp(-53.51*t9i13*dd)
      dee    = ee*53.51*(oneth*t9i43*dd - t9i13*ddd)

      term    = cc * ee
      dtermdt = dcc*ee + cc*dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 6.38e+10 * t932 * exp(-123.52*t9i)
      drevdt = rev*(1.5d0*t9i + 123.52*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg26an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,zz, &
                       gt9x,dgt9x,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56


! mg26(a,n)si29
      aa    = 1.0d0 + 0.0628*t9
      zz    = 1.0d0/aa

      t9a   = t9*zz
      dt9a  = (1.0d0 - t9a*0.0628)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix * t9a56*zz

      aa     = 5.0d0*exp(-20.990*t9i)
      daa    = aa*20.990*t9i2
      gt9x   = 1.0d0 + aa
      dgt9x  = daa

      bb     = 1.0d0/gt9x
      dbb    = -bb*bb*dgt9x

      term   = 2.93e+20 * bb * t9a56 * t9i32 * exp(-53.505/t9a13)
      dtermdt= term*(dbb*gt9x + dt9a56/t9a56 - 1.5d0*t9i &
               + 53.505/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.68*exp(-0.401*t9i)
      drevdt   = rev*0.401*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_al25pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee

! al25(p,g)si26
! coc et al 1995 a&a 299, 479 , case b

      aa  = 8.98d+1 * t9i32 * exp(-4.874*t9i)
      daa = aa*(-1.5d0*t9i + 4.874*t9i2)

      bb  = 1.568d+3 * t9i32 * exp(-9.632*t9i)
      dbb = bb*(-1.5d0*t9i + 9.632*t9i2)

      cc  = 2.42d+8 * t9i23 * exp(-23.18*t9i13)
      dcc = cc*(-twoth*t9i + oneth*23.18*t9i43)

      dd  = 4.10d-02 * t9i32 * exp(-1.741*t9i)
      ddd = dd*(-1.5d0*t9i + 1.741*t9i2)

      ee  = 2.193d+3 * t9i32 * exp(-4.642*t9i)
      dee = ee*(-1.5d0*t9i + 4.642*t9i2)

      term    = aa + bb + cc + dd + ee
      dtermdt = daa + dbb + dcc + ddd + dee


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.117d+11 * t932 * exp(-64.048*t9i)
      drevdt   = rev*(1.5d0*t9i + 64.048*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_al26pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,theta
      parameter        (theta = 0.1d0)


! al26(p,g)si27
! coc et al 1995 a&a 299, 479

      aa  = 1.53d+9 * t9**(-1.75) * exp(-23.19*t9i13)
      daa = aa*(-1.75*t9i + oneth*23.19*t9i43)

      bb  = theta*8.7d-7 * t9i32 * exp(-0.7845*t9i)
      dbb = bb*(-1.5d0*t9i + 0.7845*t9i2)

      cc  = theta*1.00d-3 * t9i32 * exp(-1.075*t9i)
      dcc = cc*(-1.5d0*t9i + 1.075*t9i2)

      dd  = 9.00d+00 * t9i32 * exp(-2.186*t9i)
      ddd = dd*(-1.5d0*t9i + 2.186*t9i2)

      ee  = 5.05d+02 * t9i32 * exp(-3.209*t9i)
      dee = ee*(-1.5d0*t9i + 3.209*t9i2)

      ff  = 9.45d+03 * t9i * exp(-4.008*t9i)
      dff = ff*(-t9i + 4.008*t9i2)

      term    = aa + bb + cc + dd + ee + ff
      dtermdt = daa + dbb + dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.46d+10 * t932 * exp(-86.621*t9i)
      drevdt   = rev*(1.5d0*t9i + 86.621*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_al27an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! al27(a,n)p30

      aa    = 8.2e+04*exp(-30.588*t9i)
      daa   = aa*30.588*t9i2

      bb    = 5.21e+05 * t974 * exp(-33.554*t9i)
      dbb   = 1.75d0*bb*t9i + bb*33.554*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term


      aa  = 5.21e+05 * t974 * exp(-2.966*t9i)
      daa = aa*(1.75d0*t9i + 2.966*t9i2)

      rev      = 6.75d0 * (8.20e4 + aa)
      drevdt   = 6.75d0 * daa

      rr    = den * rev
      drrdt = den * drevdt * 1.0d-9
      drrdd = rev

      return
      end





      subroutine rate_si27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd


! si27(p,g)p28

      aa  = 1.64e+09 * t9i23 * exp(-24.439*t9i13)
      daa = aa*(-twoth*t9i + oneth*24.439*t9i43)

      bb  = 2.00e-08 * t9i32 * exp(-0.928*t9i)
      dbb = bb*(-1.5d0*t9i + 0.928*t9i2)

      cc  = 1.95e-02 * t9i32 * exp(-1.857*t9i)
      dcc = cc*(-1.5d0*t9i + 1.857*t9i2)

      dd  = 3.70e+02 * t9i47 * exp(-3.817*t9i)
      ddd = dd*(-foursev*t9i + 3.817*t9i2)

      term    = aa + bb + cc + dd
      dtermdt = daa + dbb + dcc + ddd


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.62e+10 * t932 * exp(-23.960*t9i)
      drevdt   = rev*(1.5d0*t9i + 23.960*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si28pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,xx,dxx,q1
      parameter        (q1 = 1.0d0/8.4681d0)



! si28(p,g)p29

! champagne et al 96

      if (t9.le.5.0) then

       aa  = 8.44d+08 * t9i23 * exp(-24.389*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*24.389*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.17*t913 + 0.113*t923 + 0.0135*t9 &
             + 0.194*t943 + 0.0591*t953
       dbb = oneth*0.17*t9i23 + twoth*0.113*t9i13 + 0.0135 &
            + fourth*0.194*t913 + fiveth*0.0591*t923

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 2.92d+02 * t9i32 * exp(-4.157*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.157*t9i2)

       dd   = 4.30d+05 * t9i32 * exp(-18.51*t9i)
       ddd  = dd*(-1.5d0*t9i + 18.51*t9i2)

       ee   = 6.05d+03 * t9i32 * exp(-18.17*t9i)
       dee  = ee*(-1.5d0*t9i + 18.17*t9i2)

       term    = xx + cc + dd + ee
       dtermdt = dxx + dcc + ddd + dee


! cf88
      else

       aa  = 1.64d+08 * t9i23 * exp(-24.449*t9i13 - t92*q1)
       daa = aa*(-twoth*t9i + oneth*24.449*t9i43 - 2.0d0*t9*q1)

       bb  = 1.0d0 + 0.017*t913 - 4.11*t923 - 0.491*t9 &
             + 5.22*t943 + 1.58*t953
       dbb = oneth*0.017*t9i23 - twoth*4.11*t9i13 - 0.491 &
            + fourth*5.22*t913 + fiveth*1.58*t923

       xx  = aa * bb
       dxx = daa*bb + aa*dbb

       cc   = 3.52d+02 * t9i32 * exp(-4.152*t9i)
       dcc  = cc*(-1.5d0*t9i + 4.152*t9i2)

       dd   = 6.3d+05 * t9i32 * exp(-18.505*t9i)
       ddd  = dd*(-1.5d0*t9i + 18.505*t9i2)

       ee   = 1.69d+03 * exp(-14.518*t9i)
       dee  = ee*14.518*t9i2

       term    = xx + cc + dd + ee
       dtermdt = dxx + dcc + ddd + dee

      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.46d+09 * t932 * exp(-31.879*t9i)
      drevdt   = rev*(1.5d0*t9i + 31.879*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si29pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,xx,dxx,q1
      parameter        (q1 = 1.0d0/0.065536d0)



! si29(p,g)p30

      aa  = 3.26e+09 * t9i23 * exp(-24.459*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*24.459*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.017*t913 + 4.27*t923 + 0.509*t9 &
            + 15.40*t943 + 4.67*t953
      dbb = oneth*0.017*t9i23 + twoth*4.27*t9i13 + 0.509 &
           + fourth*15.40*t913 + fiveth*4.67*t923

      xx  = aa * bb
      dxx = daa*bb + aa*dbb

      cc   = 2.98e+03 * t9i32 * exp(-3.667*t9i)
      dcc  = cc*(-1.5d0*t9i + 3.667*t9i2)

      dd   = 3.94e+04 * t9i32 * exp(-4.665*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.665*t9i2)

      ee   = 2.08e+04 * t912 * exp(-8.657*t9i)
      dee  = ee*(0.5d0*t9i + 8.657*t9i2)

      term    = xx + cc + dd + ee
      dtermdt = dxx + dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.26e+10 * t932 * exp(-65.002*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.002*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si30pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,xx,dxx,q1
      parameter        (q1 = 1.0d0/0.4489d0)



! si30(p,g)p31

      aa  = 4.25e8 * t9i23 * exp(-24.468*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*24.468*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.017*t913 + 0.150*t923 + 0.018*t9 &
            + 5.53*t943 + 1.68*t953
      dbb = oneth*0.017*t9i23 + twoth*0.150*t9i13 + 0.018 &
           + fourth*5.53*t913 + fiveth*1.68*t923

      xx  = aa * bb
      dxx = daa*bb + aa*dbb

      cc   = 1.86e4 * t9i32 * exp(-5.601*t9i)
      dcc  = cc*(-1.5d0*t9i + 5.601*t9i2)

      dd   = 3.15e5 * t9i32 * exp(-6.961*t9i)
      ddd  = dd*(-1.5d0*t9i + 6.961*t9i2)

      ee   = 2.75e5 * t9i12 * exp(-10.062*t9i)
      dee  = ee*(-0.5d0*t9i + 10.062*t9i2)

      term    = xx + cc + dd + ee
      dtermdt = dxx + dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.50e9 * t932 * exp(-84.673*t9i)
      drevdt   = rev*(1.5d0*t9i + 84.673*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end








      subroutine rate_weaknp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,cc,dcc, &
                       zm1,zm2,zm3,zm4,zm5, &
                       c1,c2
      parameter        (c1 = 1.0d0/5.93d0, &
                        c2 = 0.98d0/885.7d0)


! free decay of neutrons, n(e-nu)p and p(e-,nu)n
! fit formula from schramm and wagoner annual review 1977
! currently accepted best value for the neutron lifetime,
! 886.7 (+/- 1.9) seconds. P.R. Huffman et al., Nature, 6 January 2000.
! world average 885.7 +/- 0.8.

      zm1   = t9 * c1
      zm2   = zm1*zm1
      zm3   = zm1*zm2
      zm4   = zm1*zm3
      zm5   = zm1*zm4

      aa   = 27.512*zm5 + 36.492*zm4 + 11.108*zm3 &
             - 6.382*zm2 + 0.565*zm1 + 1.0d0
      daa  = (5.0d0*27.512*zm4 + 4.0d0*36.492*zm3 + 3.0d0*11.108*zm2 &
             - 2.0d0*6.382*zm1 + 0.565)*c1


! n=>p
      fr    = c2 * aa
      dfrdt = c2 * daa * 1.0d-9
      dfrdd = 0.0d0


      aa  = 27.617*zm5 + 34.181*zm4 + 18.059*zm3 &
            - 16.229*zm2 + 5.252*zm1
      daa = (5.0d0*27.617*zm4 + 4.0d0*34.181*zm3 + 3.0d0*18.059*zm2 &
            - 2.0d0*16.229*zm1 + 5.252)*c1

      bb = exp(-2.531d0/zm1)
      dbb = bb*2.531d0/zm2*c1

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

! p=>n
      rr    = c2 * cc
      drrdt = c2 * dcc * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_dpn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


! d(p,n)2p
       aa   = 3.35e7 * exp(-3.720*t9i13)
       daa  = aa*oneth*3.720d0*t9i43

       bb   = 1.0d0 + 0.784*t913 + 0.346*t923 + 0.690*t9
       dbb  = oneth*0.784*t9i23 + twoth*0.346*t9i13 + 0.690

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

! rate
      cc = exp(-25.815*t9i)
      dcc = cc*25.815*t9i2

      fr    = den * cc * term
      dfrdt = den * (dcc*term + cc*dtermdt) * 1.0d-9
      dfrdd = cc * term

      rev      =  4.24e-10 * t9i32
      drevdt   = -1.5d0*rev*t9i

      rr    = den**2 * rev * term
      drrdt = den**2 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0*den * rev * term

      return
      end




      subroutine rate_dng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,c1
      parameter        (c1 = 66.2d0*18.9d0)


! d(n,g)t
      term    = 66.2 * (1.0d0  + 18.9*t9)
      dtermdt = c1

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      =  1.63e+10 * t9i32 * exp(-72.62*t9i)
      drevdt   = rev*(-1.5d0*t9i + 72.62*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ddp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb



! d(d,p)t
      aa  = 4.13e8 * t9i23 * exp(-4.258*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*4.258*t9i43

      bb  = 1.0d0 + 0.098*t913 + 4.39e-2*t923 + 3.01e-2*t9 &
            + 0.543*t943 + 0.946*t953
      dbb = oneth*0.098*t9i23 + twoth*4.39e-2*t9i13 + 3.01e-2 &
            + fourth*0.543*t913 + fiveth*0.946*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.73 * exp(-46.798*t9i)
      drevdt   = rev*46.798*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ddn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb



! d(d,n)he3
      aa  = 3.88e8 * t9i23 * exp(-4.258*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*4.258*t9i43

      bb  = 1.0d0 + 0.098*t913 + 0.418*t923 + 0.287*t9 &
            + 0.638*t943 + 1.112*t953
      dbb = oneth*0.098*t9i23 + twoth*0.418*t9i13 + 0.287 &
            + fourth*0.638*t913 + fiveth*1.112*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.730 * exp(-37.935*t9i)
      drevdt   = rev*37.935*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_tpn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! t(p,n)he3
      term    = 7.07e8 * (1.0d0 - 0.15*t912 + 0.098*t9)
      dtermdt = 7.07e8 * (-0.5d0*0.15*t9i12 + 0.098)

      aa  = exp(-8.863*t9i)
      daa = aa*8.863*t9i2

! rate
      fr    = den * aa * term
      dfrdt = den * (daa*term + aa*dtermdt) * 1.0d-9
      dfrdd = aa * term

      rev      = 0.998
      drevdt   = 0.0d0

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_ddg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! d(d,g)he4
      aa  = 4.84e+01 * t9i23 * exp(-4.258*t9i13)
      daa = aa*(-twoth*t9i + oneth*4.258*t9i43)

      bb  = 1.0d0 + 0.098*t913 - 0.203*t923 - 0.139*t9 &
            + 0.106*t943 + 0.185*t953
      dbb = oneth*0.098*t9i23 - twoth*0.203*t9i13 - 0.139 &
            + fourth*0.106*t913 + fiveth*0.185*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.53e+10 * t932 * exp(-276.729*t9i)
      drevdt   = rev*(1.5d0*t9i + 276.729*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_tpg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! t(p,g)he4
      aa  = 2.20e+04 * t9i23 * exp(-3.869*t9i13)
      daa = aa*(-twoth*t9i + oneth*3.869*t9i43)

      bb  = 1. + 0.108*t913 + 1.68*t923 + 1.26*t9 &
            + 0.551*t943 + 1.06*t953
      dbb = oneth*0.108*t9i23 + twoth*1.68*t9i13 + 1.26 &
            + fourth*0.551*t913 + fiveth*1.06*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.61e+10 * t932 * exp(-229.932*t9i)
      drevdt   = rev*(1.5d0*t9i + 229.932*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_tdn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,q1
      parameter        (q1 = 1.0d0/0.0144d0)


! t(d,n)he4 ; the "dt" reaction
      aa  = 8.09e+10 * t9i23 * exp(-4.524*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*4.524*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.092*t913 + 1.80*t923 + 1.16*t9 &
            + 10.52*t943 + 17.24*t953
      dbb = oneth*0.092*t9i23 + twoth*1.80*t9i13 + 1.16 &
            + fourth*10.52*t913 + fiveth*17.24*t923

      cc  = 8.73e+08 * t9i23 * exp(-0.523*t9i)
      dcc = -twoth*cc*t9i + cc*0.523*t9i2

      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.54*exp(-204.117*t9i)
      drevdt   = rev*204.117*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_tt2n(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! t(t,2n)he4
      aa  = 1.67e+09 * t9i23 * exp(-4.872*t9i13)
      daa = aa*(-twoth*t9i + oneth*4.872*t9i43)

      bb  = 1.0d0 + 0.086*t913 - 0.455*t923 - 0.272*t9 &
            + 0.148*t943 + 0.225*t953
      dbb = oneth*0.086*t9i23 - twoth*0.455*t9i13 - 0.272 &
            + fourth*0.148*t913 + fiveth*0.225*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.38e-10 * t9i32 * exp(-131.504*t9i)
      drevdt   = rev*(-1.5d0*t9i + 131.504*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end





      subroutine rate_he3dp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,q1
      parameter        (q1 = 1.0d0/0.099225d0)


! he3(d,p)he4
      aa  = 5.86e+10 * t9i23 * exp(-7.181*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*7.181*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.058*t913 + 0.142*t923 + 0.0578*t9 &
            + 2.25*t943 + 2.32*t953
      dbb = oneth*0.058*t9i23 + twoth*0.142*t9i13 + 0.0578 &
            + fourth*2.25*t913 + fiveth*2.32*t923

      cc  = 4.36e+08 * t9i12 * exp(-1.72*t9i)
      dcc = -0.5d0*cc*t9i + cc*1.72*t9i2

      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.55*exp(-212.980*t9i)
      drevdt   = rev*212.980*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_he3td(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,t9a,dt9a, &
                       t9a13,dt9a13,t9a56,dt9a56,zz


! he3(t,d)he4
      aa       = 1.0d0 + 0.128*t9
      zz       = 1.0d0/aa

      t9a      = t9*zz
      dt9a     = (1.0d0 - t9a*0.128)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56    = t9a**fivsix
      dt9a56   = fivsix*t9a56*zz

      term     = 5.46e+09 * t9a56 * t9i32 * exp(-7.733/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 7.733/t9a13**2 * dt9a13)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.60*exp(-166.182*t9i)
      drevdt   = rev*166.182*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_he3tnp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,t9a,dt9a, &
                       t9a13,dt9a13,t9a56,dt9a56,zz


! he3(t,np)he4
      aa       = 1.0d0 + 0.115d0 * t9
      zz       = 1.0d0/aa

      t9a      = t9*zz
      dt9a     = (1.0d0 - t9a * 0.115d0)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56    = t9a**fivsix
      dt9a56   = fivsix*t9a56*zz

      term     = 7.71d9 * t9a56 * t9i32 * exp(-7.733d0/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 7.733d0/t9a13**2 * dt9a13)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.39d-10 * t9i32 * exp(-140.367*t9i)
      drevdt   = rev*(-1.5d0 * t9i + 140.367d0 * t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end




      subroutine rate_he4npg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! he4(np,g)li6
      aa  = 4.62e-6 * t9i2 * exp(-19.353*t9i)
      daa = aa*(-2.0d0*t9i + 19.353*t9i2)

      bb  = 1.0d0 + 0.075*t9
      dbb = 0.075

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.22e19 * t93 * exp(-42.933*t9i)
      drevdt   = rev*(3.0d0*t9i + 42.933*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he4dg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc


! he4(d,g)li6
      aa  = 3.01e1 * t9i23 * exp(-7.423*t9i13)
      daa = aa*(-twoth*t9i + oneth*7.423*t9i43)

      bb  = 1.0d0 + 0.056*t913 - 4.85*t923 + 8.85*t9 &
            - 0.585*t943 - 0.584*t953
      dbb = oneth*0.056*t9i23 - twoth*4.85*t9i13 + 8.850 &
            - fourth*0.585*t913 - fiveth*0.584*t923

! rate goes negative for t9 greater than about 15, so try this
      if (bb .le. 0.0) then
       bb = 0.0d0
       dbb = 0.0d0
      end if

      cc =  8.55e1 * t9i32 * exp(-8.228*t9i)
      dcc = cc*(-1.5d0*t9i + 8.228*t9i2)

      term    = aa * bb + cc
      dtermdt = daa*bb + aa*dbb + dcc


! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.53e10 * t932 * exp(-17.1180*t9i)
      drevdt   = rev*(1.5d0*t9i + 17.1180*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_he4tn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,t9a,dt9a,t9a32,dt9a32,zz


! he4(t,n)li6
      aa     = 1.0d0 + 49.180*t9
      zz     = 1.0d0/aa

      t9a   = t9*zz
      dt9a   = (1.0d0 - t9a*49.180)*zz

      t9a32  = t9a * sqrt(t9a)
      dt9a32 = 1.5d0*t9a32/t9a * dt9a

      aa     = 1.80e8 * exp(-55.4940*t9i)
      daa    = aa*55.4940*t9i2

      bb     = 1.0d0 - 0.2610 * t9a32 * t9i32
      dbb    = -0.2610*(-1.5d0*t9a32*t9i52 + dt9a32*t9i32)

      cc     = aa*bb
      dcc    = daa*bb + aa*dbb

      dd     = 2.72e9 * t9i32 * exp(-57.8840*t9i)
      ddd    = dd*(-1.5d0*t9i + 57.8840*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.72e9 * t9i32 * exp(-2.39*t9i)
      drevdt   = rev*(-1.5d0*t9i + 2.39*t9i2)

      term = 0.935*(1.80e8*bb + rev)
      dtermdt = 0.935*(1.80e8*dbb + drevdt)

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_li6phe3(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,q1
      parameter        (q1 = 1.0d0/30.25d0)


! li6(p,he3)he4
      aa  = 3.73e10 * t9i23 * exp(-8.413*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*8.413*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.050*t913 - 0.061*t923 - 0.0210*t9 &
            + 0.0060*t943 + 0.0050*t953
      dbb = oneth*0.050*t9i23 - twoth*0.061*t9i13 - 0.0210 &
            + fourth*0.0060*t913 + fiveth*0.0050*t923

      cc  = 1.33e10 * t9i32 * exp(-17.7630*t9i)
      dcc = -1.5d0*cc*t9i + cc*17.7630*t9i2

      dd  =  1.29e9 * t9i * exp(-21.82*t9i)
      ddd = -dd*t9i + dd*21.82*t9i2

      term    = aa * bb + cc + dd
      dtermdt = daa*bb + aa*dbb + dcc + ddd

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.070d0 * exp(-46.6310*t9i)
      drevdt   = rev*46.6310*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_li6ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li6(n,g)li7
! malaney-fowler 1989

      term    = 5.10e3
      dtermdt = 0.0d0

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.19e10 * t932 * exp(-84.17*t9i)
      drevdt   = rev*(1.5d0*t9i + 84.17*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he4tg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! he4(t,g)li7
      aa  = 8.67e5 * t9i23 * exp(-8.08*t9i13)
      daa = aa*(-twoth*t9i + oneth*8.08*t9i43)

      bb  = 1.0d0 + 0.052*t913 - 0.448*t923 - 0.165*t9 &
            + 0.144*t943 + 0.134*t953
      dbb = oneth*0.052*t9i23 - twoth*0.448*t9i13 - 0.165 &
            + fourth*0.144*t913 + fiveth*0.134*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.11e10 * t932 * exp(-28.64*t9i)
      drevdt   = rev*(1.5d0*t9i + 28.64*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_li7dn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(d,n)2a
      term    = 2.92e11 * t9i23 * exp(-10.259*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*10.259*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.95e-10 * t9i32 * exp(-175.476*t9i)
      drevdt   = rev*(-1.5d0*t9i + 175.476*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end




      subroutine rate_li7tn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(t,n)be9
! malaney and fowler (apjl, 345, l5, 1989)
      term    = 1.46d+11 * t9i23 * exp(-11.333*t9i13)
      dtermdt = -twoth*term*t9i + oneth*term*11.333*t9i43

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0d0
      drevdt   = 0.0d0

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_li7t2n(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(t,2n)2a
      term    = 8.81e11 * t9i23 * exp(-11.333*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*11.333*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.22e-19 * t9i3 * exp(-102.864*t9i)
      drevdt   = rev*(-3.0d0*t9i + 102.864*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end




      subroutine rate_li7he3np(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! li7(he3,np)2a
      term    = 1.11e13 * t9i23 * exp(-17.989*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*17.989*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.09e-20 * t9i3 * exp(-111.727*t9i)
      drevdt   = rev*(-3.0d0*t9i + 111.727*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end




      subroutine rate_li6pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb,cc,dcc, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz


! li6(p,g)be7
      if (t9 .gt. 10.0) then
       t9a  = 1.0d0
       dt9a = 0.0d0
      else
       aa   = 1.0d0 - 0.0969*t9

       bb   = aa**(-twoth)
       dbb  = twoth*bb/aa*0.0969

       cc   = aa + 0.0284*t953*bb
       dcc  = -0.0969 + 0.0284*(fiveth*t923*bb + t953*dbb)

       zz   = 1.0d0/cc
       t9a  = t9*zz
       dt9a = (1.0d0 - t9a*dcc)*zz
      end if

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      term    = 6.69e+05 * t9a56 * t9i32 * exp(-8.413/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 8.413/t9a13**2 * dt9a13)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.19e+10 * t932 * exp(-65.054*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.054*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0


      return
      end




      subroutine rate_li7pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! li7(p,n)be7
      aa  = 5.15e+09 * exp(-1.167*t913 - 19.081*t9i)
      daa = aa*(-oneth*1.167*t9i23 + 19.081*t9i2)

      bb  = 7.84e+09 * t9i32 * exp(-22.832*t9i)
      dbb = -1.5d0*bb*t9i + bb*22.832*t9i2

      term    = aa + bb
      dtermdt = daa + dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      aa  = 5.15e+09 * exp(-1.167*t913)
      daa = -aa*oneth*1.167*t9i23

      bb  = 0.998 * 7.84e+09 * t9i32 * exp(-3.751*t9i)
      dbb = -1.5d0*bb*t9i + bb*3.751*t9i2

      term    = aa + bb
      dtermdt = daa + dbb
      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end




      subroutine rate_li7ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! li7(n,g)li8
! apj 372, 1
      aa      = 4.26d+03 * t9i32 * exp(-2.576*t9i)
      daa     = aa*(-1.5d0*t9i + 2.576*t9i2)

      term    = 3.144d+03 + aa
      dtermdt = daa

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 1.2923d+10 * t932 * exp(-2.359d+01*t9i)
      drevdt = rev*(1.5d0*t9i + 2.359d+01*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_be7dp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! be7(d,p)2a
      term    = 1.07e12 * t9i23 * exp(-12.428*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*12.428*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.97e-10 * t9i32 * exp(-194.557*t9i)
      drevdt   = rev*(-1.5d0*t9i + 194.557*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end




      subroutine rate_be7tnp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! be7(t,np)2a
      term    = 2.91e12 * t9i23 * exp(-13.729*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*13.729*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.09e-20 * t9i3 * exp(-121.944*t9i)
      drevdt   = rev*(-3.0d0*t9i + 121.944*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end



      subroutine rate_be7he32p(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! be7(he3,2p)2a
      term    = 6.11e13 * t9i23 * exp(-21.793*t9i13)
      dtermdt = term*(-twoth*t9i + oneth*21.793*t9i43)

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.22e-19 * t9i3 * exp(-130.807*t9i)
      drevdt   = rev*(-3.0d0*t9i + 130.807*t9i2)

      rr    = den**3 * rev * term
      drrdt = den**3 * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 3.0d0 * den**2 * rev * term

      return
      end




      subroutine rate_be9pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.2704d0)


! be9(p,a)li6
      aa  = 2.11e11 * t9i23 * exp(-10.359*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*10.359*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.04*t913 + 1.09*t923 + 0.307*t9 &
            + 3.21*t943 + 2.30*t953
      dbb  = oneth*0.04*t9i23 + twoth*1.09*t9i13 + 0.307 &
             + fourth*3.21*t913 + fiveth*2.30*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 4.51e8 * t9i * exp(-3.046*t9i)
      ddd  = -dd*t9i + dd*3.046*t9i2

      ee   = 6.70e8 * t9i34 * exp(-5.160*t9i)
      dee  = -0.75d0*ee*t9i + ee*5.160*t9i2

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.18e-1 * exp(-24.674*t9i)
      drevdt   = rev*24.674*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li6ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/1.758276d0)


! li6(a,g)b10
      aa  = 4.06e6 * t9i23 * exp(-18.790*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*18.790*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.022*t913 + 1.54*t923 + 0.239*t9 &
             +  2.20*t943 + 0.869*t953
      dbb  = oneth*0.022*t9i23 + twoth*1.54*t9i13 + 0.239 &
             + fourth*2.20*t913 + fiveth*0.869*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1910.0 * t9i32 * exp(-3.484*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.484*t9i2)

      ee   = 1.01e4 * t9i * exp(-7.269*t9i)
      dee  = ee*(-t9i + 7.269*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.58e10 * t932 * exp(-51.753*t9i)
      drevdt   = rev*(1.5d0*t9i + 51.753*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_li7an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt

! li7(a,n)b10
      term    = 3.84e8 * exp(-32.382*t9i)
      dtermdt = term*32.382*t9i2

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.84e8*1.32
      drevdt   = 0.0d0

      rr    = den * rev
      drrdt = 0.0d0
      drrdd = rev

      return
      end





      subroutine rate_be9pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.715716d0)


! be9(p,g)b10
      aa  = 1.33e7 * t9i23 * exp(-10.359*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*10.359*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.04*t913 + 1.52*t923 + 0.428*t9 &
             + 2.15*t943 + 1.54*t953
      dbb  = oneth*0.04*t9i23 + twoth*1.52*t9i13 + 0.428 &
             + fourth*2.15*t913 + fiveth*1.54*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 9.64e4 * t9i32 * exp(-3.445*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.445*t9i2)

      ee   = 2.72e6 * t9i32 * exp(-10.62*t9i)
      dee  = ee*(-1.5d0*t9i + 10.62*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.73e9 * t932 * exp(-76.427*t9i)
      drevdt   = rev*(1.5d0*t9i + 76.427*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_b10pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/19.377604d0)

! b10(p,a)li7
      aa  = 1.26e11 * t9i23 * exp(-12.062*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*12.062*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.035*t913 - 0.498*t923 - 0.121*t9 &
             + 0.3*t943 + 0.184*t953
      dbb  = oneth*0.035*t9i23 - twoth*0.498*t9i13 - 0.121 &
             + fourth*0.3*t913 + fiveth*0.184*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.59e9 * t9i * exp(-12.260*t9i)
      ddd  = -dd*t9i + dd*12.260*t9i2

      term    = cc + dd
      dtermdt = dcc + ddd

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.54e-01 * exp(-13.301*t9i)
      drevdt   = rev*13.301*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_li7ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/17.598025d0)


! li7(a,g)b11
      aa  = 3.55e7 * t9i23 * exp(-19.161*t9i13 -t92*q1)
      daa = aa*(-twoth*t9i + oneth*19.161*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.022*t913 + 0.775*t923 + 0.118*t9 &
             + 0.884*t943 + 0.342*t953
      dbb  = oneth*0.022*t9i23 + twoth*0.775*t9i13 + 0.118 &
             + fourth*0.884*t913 + fiveth*0.342*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 3.33e2 * t9i32 * exp(-2.977*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.977*t9i2)

      ee   = 4.10e4 * t9i * exp(-6.227*t9i)
      dee  = ee*(-t9i + 6.227*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.02e10 * t932 * exp(-100.538*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.538*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_b11pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/2.702736d0)

! b11(p,a)be8=>2a
      aa  = 2.20e12 * t9i23 * exp(-12.095*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*12.095*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0  + 0.034*t913 + 0.14*t923 + 0.034*t9 &
             + 0.19*t943 + 0.116*t953
      dbb  = oneth*0.034*t9i23 + twoth*0.14*t9i13 + 0.034 &
             + fourth*0.19*t913 + fiveth*0.116*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb


      dd   = 4.03e6 * t9i32 * exp(-1.734*t9i)
      ddd  = dd*(-1.5d0*t9i + 1.734*t9i2)

      ee   = 6.73e9 * t9i32 * exp(-6.262*t9i)
      dee  = ee*(-1.5d0*t9i + 6.262*t9i2)

      ff   = 3.88e9*t9i * exp(-14.154*t9i)
      dff  = ff*(-t9i + 14.154*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd +dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.5e-10* t9i32 *exp(-100.753*t9i)
      drevdt   = rev*(-1.5d0*t9i + 100.753*t9i2)

      rr    = den * den * rev * term 
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term 

      return
      end




      subroutine rate_be7ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/22.743361d0)

! be7(a,g)c11
      aa  = 8.45e+07 * t9i23 * exp(-23.212*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*23.212*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.018*t913 + 0.488*t923 + 0.061*t9 &
             + 0.296*t943 + 0.095*t953
      dbb  = oneth*0.018*t9i23 + twoth*0.488*t9i13 + 0.061 &
             + fourth*0.296*t913 + fiveth*0.095*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.25e+04 * t9i32 * exp(-6.510*t9i)
      ddd  = dd*(-1.5d0*t9i + 6.510*t9i2)

      ee   = 1.29e+05 * t9i54 * exp(-10.039*t9i)
      dee  = ee*(-1.25d0*t9i + 10.039*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.02e+10 * t932 * exp(-87.539*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.539*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_b11pn(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! b11(p,n)c11
      aa  = 1.69e8*(1.0d0 - 0.048*t912 + 0.010*t9)
      daa = 1.69e8*(-0.5d0*0.048*t9i12 + 0.010)

      bb  = exp(-32.080*t9i)
      dbb = bb*32.080*t9i2

      term    = aa*bb
      dtermdt = daa*bb + aa*dbb

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      term     = 0.998*aa
      dtermdt  = 0.998*daa

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_b8ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh, &
                       pp,dpp,qq,dqq,drr


! b8(a,p)c11
      aa  = 1.67e-09 * t9i32 * exp(-1.079*t9i)
      daa = -1.5d0*aa*t9i + aa*1.079*t9i2

      bb  = 9.55e-08 * t9i32 * exp(-1.311*t9i)
      dbb = -1.5d0*bb*t9i + bb*1.311*t9i2

      cc  = 1.98e-01 * t9i32 * exp(-2.704*t9i)
      dcc = -1.5d0*cc*t9i + cc*2.704*t9i2

      dd  = 1.34e+00 * t9i32 * exp(-4.282*t9i)
      ddd = -1.5d0*dd*t9i + dd*4.282*t9i2

      ee  = 3.22e+04 * t9i32 * exp(-6.650*t9i)
      dee = -1.5d0*ee*t9i + ee*6.650*t9i2

      ff  = 2.33e+05 * t9i32 * exp(-8.123*t9i)
      dff = -1.5d0*ff*t9i + ff*8.123*t9i2

      gg  = 2.55e+06 * t9i32 * exp(-11.99*t9i)
      dgg = -1.5d0*gg*t9i + gg*11.99*t9i2

      hh  = 9.90e+06 * t9i32 * exp(-13.50*t9i)
      dhh = -1.5d0*hh*t9i + hh*13.50*t9i2

      pp  = 1.41e+06 * t9i32 * exp(-16.51*t9i)
      dpp = -1.5d0*pp*t9i + pp*16.51*t9i2

      qq  = 1.99e+07 * t9i32 * exp(-18.31*t9i)
      dqq = -1.5d0*qq*t9i + qq*18.31*t9i2

      rr  = 6.01e+07 * t9i32 * exp(-20.63*t9i)
      drr = -1.5d0*rr*t9i + rr*20.63*t9i2

      term    = aa + bb + cc + dd + ee + ff + gg + hh + pp + qq + rr
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg + dhh &
                + dpp + dqq + drr

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.101 * exp(-85.95*t9i)
      drevdt   = rev*85.95*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev*term

      return
      end





      subroutine rate_b10pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/19.377604d0)


! b10(p,g)c11
      aa  = 4.61e+05 * t9i23 * exp(-12.062*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*12.062*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.035*t913 + 0.426*t923 + 0.103*t9 &
             + 0.281*t943 + 0.173*t953
      dbb  = oneth*0.035*t9i23 + twoth*0.426*t9i13 + 0.103 &
             + fourth*0.281*t913 + fiveth*0.173*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.93e+05 * t9i32 * exp(-12.041*t9i)
      ddd  = dd*(-1.5d0*t9i + 12.041*t9i2)

      ee   = 1.14e+04 * t9i32 * exp(-16.164*t9i)
      dee  = ee*(-1.5d0*t9i + 16.164*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.03e+10 * t932 * exp(-100.840*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.840*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_c11na(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd


! c11(n,a)be8=>2a   get a real rate in here

!      fr    = den * 7.0e4
!      dfrdt = 0.0d0
!      dfrdd = 7.0e4

      fr    = 0.0d0
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_be9an(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh


! be9(a,n)c12
! Wrean 94 Phys Rev C (1994) vol 49, #2, 1205
      aa  = 6.476d+13 * t9i23 * exp(-23.8702*t9i13)
      daa = -twoth*aa*t9i + oneth*aa*23.8702*t9i43

      bb  = (1.0d0 - 0.3270*t913)
      dbb = -oneth*0.3270*t9i23

! rate goes negative for t9 greater than about 15, so try this
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 6.044d-3*t9i32*exp(-1.401*t9i)
      ddd = -1.5d0*dd*t9i + dd*1.401*t9i2

      ee  = 7.268*t9i32*exp(-2.063*t9i)
      dee = -1.5d0*ee*t9i + ee*2.063*t9i2

      ff  = 3.256d+4*t9i32*exp(-3.873*t9i)
      dff = -1.5d0*ff*t9i + ff*3.873*t9i2

      gg  = 1.946d+5*t9i32*exp(-4.966*t9i)
      dgg = -1.5d0*gg*t9i + gg*4.966*t9i2

      hh  = 1.838e9*t9i32*exp(-15.39*t9i)
      dhh = -1.5d0*hh*t9i + hh*15.39*t9i2

      term    = cc + dd + ee + ff + gg + hh
      dtermdt = dcc + ddd + dee + dff + dgg + dhh


! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 10.3d0 * exp(-66.160*t9i)
      drevdt   = rev*66.160d0*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term


      return
      end




      subroutine rate_b11pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/0.57121d0)


! b11(p,g)c12
      aa  = 4.62e+07 * t9i23 * exp(-12.095*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*12.095*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.035*t913 + 3.0*t923 + 0.723*t9 &
             + 9.91*t943 + 6.07*t953
      dbb  = oneth*0.035*t9i23 + twoth*3.0*t9i13 + 0.723 &
             + fourth*9.91*t913 + fiveth*6.07*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 7.89e+03 * t9i32 * exp(-1.733*t9i)
      ddd  = dd*(-1.5d0*t9i + 1.733*t9i2)

      ee   = 9.68e+04 * t9i15 * exp(-5.617*t9i)
      dee  = ee*(-0.2d0*t9i + 5.617*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.01e+10 * t932 * exp(-185.173*t9i)
      drevdt   = rev*(1.5d0*t9i + 185.173*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_b11ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.120409d0)


! b11(a,p)c14
      aa  = 5.37e+11 * t9i23 * exp(-28.234*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*28.234*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.015*t913 + 5.575*t923 + 0.576*t9 &
             + 15.888*t943 + 4.174*t953
      dbb  = oneth*0.015*t9i23 + twoth*5.575*t9i13 + 0.576 &
             + fourth*15.888*t913 + fiveth*4.174*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 5.44e-03 * t9i32 * exp(-2.827*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.827*t9i2)

      ee   = 3.36e+02 * t9i32 * exp(-5.178*t9i)
      dee  = ee*(-1.5d0*t9i + 5.178*t9i2)

      ff   = 5.32e+06 * t9i38 * exp(-11.617*t9i)
      dff  = ff*(-0.375*t9i + 11.617*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.10e+01 * exp(-9.098*t9i)
      drevdt   = rev*9.098*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_pp(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! p(p,e+nu)d
      if (t9 .le. 3.0) then
       aa   = 4.01d-15 * t9i23 * exp(-3.380d0*t9i13)
       daa  = aa*(-twoth*t9i + oneth*3.380d0*t9i43)

       bb   = 1.0d0 + 0.123d0*t913 + 1.09d0*t923 + 0.938d0*t9
       dbb  = oneth*0.123d0*t9i23 + twoth*1.09d0*t9i13 + 0.938d0

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

      else
       term    = 1.1581136d-15
       dtermdt = 0.0d0
      end if

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_pep(temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! p(e-p,nu)d
      if (t9 .le. 3.0) then

       aa   = 1.36e-20 * t9i76 * exp(-3.380*t9i13)
       daa  = aa*(-sevsix*t9i + oneth*3.380d0*t9i43)

       bb   = (1.0d0 - 0.729d0*t913 + 9.82d0*t923)
       dbb  = -oneth*0.729d0*t9i23 + twoth*9.82d0*t9i13

       term    = aa * bb
       dtermdt = daa * bb + aa * dbb

      else
       term    = 7.3824387e-21
       dtermdt = 0.0d0
      end if

! rate
      fr    = ye * den * den * term
      dfrdt = ye * den * den * dtermdt * 1.0d-9
      dfrdd = ye * 2.0d0 * den * term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_hep(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb


! he3(p,e+nu)he4
      if (t9 .le. 3.0) then

       aa   = 8.78e-13 * t9i23 * exp(-6.141d0*t9i13)
       daa  = aa*(-twoth*t9i + oneth*6.141d0*t9i43)

       term    = aa
       dtermdt = daa

      else
       term    = 5.9733434e-15
       dtermdt = 0.0d0
      end if

! rate
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end







      subroutine rate_png(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa


! p(n,g)d
! smith,kawano,malany 1992

       aa      = 1.0d0 - 0.8504*t912 + 0.4895*t9 &
                 - 0.09623*t932 + 8.471e-3*t92 &
                 - 2.80e-4*t952

       daa     =  -0.5d0*0.8504*t9i12 + 0.4895 &
                 - 1.5d0*0.09623*t912 + 2.0d0*8.471e-3*t9 &
                 - 2.5d0*2.80e-4*t932

       term    = 4.742e4 * aa
       dtermdt = 4.742e4 * daa


! wagoner,schramm 1977
!      aa      = 1.0d0 - 0.86*t912 + 0.429*t9
!      daa     =  -0.5d0*0.86*t9i12 + 0.429

!      term    = 4.4d4 * aa
!      dtermdt = 4.4d4 * daa



! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.71d+09 * t932 * exp(-25.82*t9i)
      drevdt   = rev*(1.5d0*t9i + 25.82*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_dpg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! d(p,g)he3
      aa      = 2.24d+03 * t9i23 * exp(-3.720*t9i13)
      daa     = aa*(-twoth*t9i + oneth*3.720*t9i43)

      bb      = 1.0d0 + 0.112*t913 + 3.38*t923 + 2.65*t9
      dbb     = oneth*0.112*t9i23 + twoth*3.38*t9i13 + 2.65

      term    = aa * bb
      dtermdt = daa * bb + aa * dbb


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.63d+10 * t932 * exp(-63.750*t9i)
      drevdt   = rev*(1.5d0*t9i + 63.750*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_he3ng(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt


! he3(n,g)he4
      term    = 6.62 * (1.0d0 + 905.0*t9)
      dtermdt = 5.9911d3

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.61d+10 * t932 * exp(-238.81*t9i)
      drevdt   = rev*(1.5d0*t9i + 238.81*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_he3he3(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! he3(he3,2p)he4
      aa   = 6.04d+10 * t9i23 * exp(-12.276*t9i13)
      daa  = aa*(-twoth*t9i + oneth*12.276*t9i43)

      bb   = 1.0d0 + 0.034*t913 - 0.522*t923 - 0.124*t9 &
             + 0.353*t943 + 0.213*t953
      dbb  = oneth*0.034*t9i23 - twoth*0.522*t9i13 - 0.124 &
             + fourth*0.353*t913 + fiveth*0.213*t923

      term    = aa * bb
      dtermdt = daa*bb + aa*dbb

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.39e-10 * t9i32 * exp(-149.230*t9i)
      drevdt   = rev*(-1.5d0*t9i + 149.230*t9i2)

      rr    = den * den * rev * term
      drrdt = den * den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 2.0d0 * den * rev * term

      return
      end





      subroutine rate_he3he4(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,t9a,dt9a, &
                       t9a13,dt9a13,t9a56,dt9a56,zz


! he3(he4,g)be7
      aa      = 1.0d0 + 0.0495*t9
      daa     = 0.0495

      zz      = 1.0d0/aa
      t9a     = t9*zz
      dt9a    = (1.0d0 - t9a*daa)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      term    = 5.61d+6 * t9a56 * t9i32 * exp(-12.826/t9a13)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 12.826/t9a13**2 * dt9a13)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.11e+10 * t932 * exp(-18.423*t9i)
      drevdt   = rev*(1.5d0*t9i + 18.423*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_be7em(temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,ye,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb

! be7(e-,nu+g)li7
      if (t9 .le. 3.0) then
       aa  = 0.0027 * t9i * exp(2.515e-3*t9i)
       daa = -aa*t9i - aa*2.515e-3*t9i2

       bb  = 1.0d0 - 0.537*t913 + 3.86*t923 + aa
       dbb = -oneth*0.537*t9i23 + twoth*3.86*t9i13 + daa

       term    = 1.34e-10 * t9i12 * bb
       dtermdt = -0.5d0*term*t9i + 1.34e-10*t9i12*dbb

      else
       term    = 0.0d0
       dtermdt = 0.0d0
      endif

! rates
      fr    = ye * den * term
      dfrdt = ye * den * dtermdt * 1.0d-9
      dfrdd = ye * term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_be7pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb


! be7(p,g)b8
      aa      = 3.11e+05 * t9i23 * exp(-10.262*t9i13)
      daa     = aa*(-twoth*t9i + oneth*10.262*t9i43)

      bb      = 2.53e+03 * t9i32 * exp(-7.306*t9i)
      dbb     = bb*(-1.5d0*t9i + 7.306*t9i2)

      term    = aa + bb
      dtermdt = daa + dbb


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.30e+10 * t932 * exp(-1.595*t9i)
      drevdt   = rev*(1.5d0*t9i + 1.595*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_li7pag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc,dd,ddd, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56,zz, &
                       term1,dterm1,term2,dterm2,rev1,drev1dt,rev2,drev2dt,q1
      parameter        (q1 = 1.0d0/2.876416d0)



! 7li(p,g)8be=>2a
      aa   = 1.56d5 * t9i23 * exp(-8.472d0 * t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*8.472*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.049d0*t913 + 2.498d0*t923 + 0.86d0*t9 &
             + 3.518d0*t943 + 3.08*t953
      dbb  = oneth*0.049*t9i23 + twoth*2.498*t9i13 + 0.86 &
             + fourth*3.518*t913 + fiveth*3.08*t923

      cc   = aa*bb
      dcc  = daa*bb + aa*dbb

      dd   =  1.55d6 * t9i32 * exp(-4.478d0 * t9i)
      ddd  = dd*(-1.5d0*t9i + 4.478*t9i2)

      term1  = cc + dd
      dterm1 = dcc + ddd

      rev1    = 6.55e+10 * t932 * exp(-200.225*t9i)
      drev1dt = rev1*(1.5d0*t9i + 200.225*t9i2)


! 7li(p,a)a
      aa     = 1.0d0 + 0.759*t9

      zz     = 1.0d0/aa
      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*0.759)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      aa     = 1.096e+09 * t9i23 * exp(-8.472*t9i13)
      daa    = aa*(-twoth*t9i + oneth*8.472*t9i43)

      bb     = -4.830e+08 * t9a56 * t9i32 * exp(-8.472/t9a13)
      dbb    = bb*(dt9a56/t9a56 - 1.5d0*t9i + 8.472/t9a13**2*dt9a13)

      cc     = 1.06e+10 * t9i32 * exp(-30.442*t9i)
      dcc    = cc*(-1.5d0*t9i + 30.442*t9i2)

      term2   = aa + bb + cc
      dterm2  = daa + dbb + dcc

      rev2    = 4.69 * exp(-201.291*t9i)
      drev2dt = rev2*201.291*t9i2


! sum the two forward rates (per f35 cf88)

      term    = term1 + term2
      dtermdt = dterm1 + dterm2


! per cf88 the reverse rate just the second one
      rev     = rev2
      drevdt  = drev2dt


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_b8ep(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 0.77d0, &
                        con      = lntwo/halflife)


! b8(e+,nu)be8 => 2a

      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/2.25d0)


! c12(p,g)13n
      aa   = 2.04e+07 * t9i23 * exp(-13.69*t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*13.69*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.03*t913 + 1.19*t923 + 0.254*t9 &
             + 2.06*t943 + 1.12*t953
      dbb  = oneth*0.03*t9i23 + twoth*1.19*t9i13 + 0.254 &
             + fourth*2.06*t913 + fiveth*1.12*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.08e+05 * t9i32 * exp(-4.925*t9i)
      ddd  = dd*(-1.5d0*t9i + 4.925*t9i2)

      ee   = 2.15e+05 * t9i32 * exp(-18.179*t9i)
      dee  = ee*(-1.5d0*t9i + 18.179*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 8.84e+09 * t932 * exp(-22.553*t9i)
      drevdt   = rev*(1.5d0*t9i + 22.553*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n13em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 597.9d0, &
                        con      = lntwo/halflife)

! n13(e-nu)c13
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_c13pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/4.0d0)


! c13(p,g)13n
      aa   = 8.01e+07 * t9i23 * exp(-13.717*t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*13.717*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.030*t913 + 0.958*t923 + 0.204*t9 &
             + 1.39*t943 + 0.753*t953
      dbb  = oneth*0.030*t9i23 + twoth*0.958*t9i13 + 0.204 &
             + fourth*1.39*t913 + fiveth*0.753*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.21e+06 * t9i65 * exp(-5.701*t9i)
      ddd  = dd*(-sixfif*t9i + 5.701*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.19e+10 * t932 * exp(-87.621*t9i)
      drevdt   = rev*(1.5d0*t9i + 87.621*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_n14pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,q1
      parameter        (q1 = 1.0d0/10.850436d0)


! n14(p,g)o15
      aa  = 4.90e+07 * t9i23 * exp(-15.228*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.228*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.027*t913 - 0.778*t923 - 0.149*t9 &
             + 0.261*t943 + 0.127*t953
      dbb  = oneth*0.027*t9i23 - twoth*0.778*t9i13 - 0.149 &
             + fourth*0.261*t913 + fiveth*0.127*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.37e+03 * t9i32 * exp(-3.011*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.011*t9i2)

      ee   = 2.19e+04 * exp(-12.530*t9i)
      dee  = ee*12.530*t9i2

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 2.70e+10 * t932 * exp(-84.678*t9i)
      drevdt = rev*(1.5d0*t9i + 84.678*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_o15em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 122.24d0, &
                        con      = lntwo/halflife)

! o15(e-nu)n15
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_n14ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.776161d0)


! n14(a,g)f18
      aa  = 7.78d+09 * t9i23 * exp(-36.031*t9i13- t92*q1)
      daa = aa*(-twoth*t9i + oneth*36.031*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.012*t913 + 1.45*t923 + 0.117*t9 &
             + 1.97*t943 + 0.406*t953
      dbb  = oneth*0.012*t9i23 + twoth*1.45*t9i13 + 0.117 &
             + fourth*1.97*t913 + fiveth*0.406*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 2.36d-10 * t9i32 * exp(-2.798*t9i)
      ddd  = dd*(-1.5d0*t9i + 2.798*t9i2)

      ee   = 2.03 * t9i32 * exp(-5.054*t9i)
      dee  = ee*(-1.5d0*t9i + 5.054*t9i2)

      ff   = 1.15d+04 * t9i23 * exp(-12.310*t9i)
      dff  = ff*(-twoth*t9i + 12.310*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.42e+10 * t932 * exp(-51.236*t9i)
      drevdt   = rev*(1.5d0*t9i + 51.236*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n15pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.2025d0)


! n15(p,g)o16
      aa  = 9.78e+08 * t9i23 * exp(-15.251*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.251*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0  + 0.027*t913 + 0.219*t923 + 0.042*t9 &
             + 6.83*t943 + 3.32*t953
      dbb  = oneth*0.027*t9i23 + twoth*0.219*t9i13 + 0.042 &
             + fourth*6.83*t913 + fiveth*3.32*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.11e+04*t9i32*exp(-3.328*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.328*t9i2)

      ee   = 1.49e+04*t9i32*exp(-4.665*t9i)
      dee  = ee*(-1.5d0*t9i + 4.665*t9i2)

      ff   = 3.8e+06*t9i32*exp(-11.048*t9i)
      dff  = ff*(-1.5d0*t9i + 11.048*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.62e+10 * t932 * exp(-140.734*t9i)
      drevdt   = rev*(1.5d0*t9i + 140.734*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_n15pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg, &
                       theta,q1
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.272484d0)


! n15(p,a)c12
      aa  = 1.08d+12*t9i23*exp(-15.251*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.251*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.027*t913 + 2.62*t923 + 0.501*t9 &
             + 5.36*t943 + 2.60*t953
      dbb  = oneth*0.027*t9i23 + twoth*2.62*t9i13 + 0.501 &
             + fourth*5.36*t913 + fiveth*2.60*t923

      cc   = aa * bb
      dcc  = daa*bb + aa*dbb

      dd   = 1.19d+08 * t9i32 * exp(-3.676*t9i)
      ddd  = dd*(-1.5d0*t9i + 3.676*t9i2)

      ee   = 5.41d+08 * t9i12 * exp(-8.926*t9i)
      dee  = ee*(-0.5d0*t9i + 8.926*t9i2)

      ff   = theta * 4.72d+08 * t9i32 * exp(-7.721*t9i)
      dff  = ff*(-1.5d0*t9i + 7.721*t9i2)

      gg   = theta * 2.20d+09 * t9i32 * exp(-11.418*t9i)
      dgg  = gg*(-1.5d0*t9i + 11.418*t9i2)

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.06d-01*exp(-57.625*t9i)
      drevdt   = rev*57.625*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o16pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,zz


! o16(p,g)f17
      aa  = exp(-0.728*t923)
      daa = -twoth*aa*0.728*t9i13

      bb  = 1.0d0 + 2.13 * (1.0d0 - aa)
      dbb = -2.13*daa

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923*dbb

      dd   = exp(-16.692*t9i13)
      ddd  = oneth*dd*16.692*t9i43

      zz   = 1.0d0/cc
      ee   = dd*zz
      dee  = (ddd - ee*dcc)*zz

      term    = 1.50d+08 * ee
      dtermdt = 1.50d+08 * dee


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.03e+09*t932*exp(-6.968*t9i)
      drevdt   = rev*(1.5d0*t9i + 6.968*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_o17pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,res1,dres1,res2,dres2,res3,dres3, &
                       res4,dres4,res5,dres5,res6,dres6,zz, &
                       theta,q1,q2
      parameter        (theta = 0.1d0, &
                        q1    = 1.0d0/0.319225d0, &
                        q2    = 1.0d0/0.0016d0)



! o17(p,a)n14
! rate from jeff blackmons thesis, includes terms from fowler 75,
! landre 1990 (a&a 240, 85), and new results
! use rev factor from cf88 rate

      aa  = 1.53d+07 * t9i23 * exp(-16.712*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*16.712*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.025*t913 + 5.39*t923 + 0.940*t9 &
             + 13.5*t943 + 5.98*t953
      dbb  = oneth*0.025*t9i23 + twoth*5.39*t9i13 + 0.940 &
             + fourth*13.5*t913 + fiveth*5.98*t923

      res1  = aa * bb
      dres1 = daa*bb + aa*dbb

      res2  = 2.92d+06 * t9 * exp(-4.247*t9i)
      dres2 = res2*(t9i + 4.247*t9i2)


      aa    = 0.479 * t923 + 0.00312
      daa   = twoth*0.479*t9i13

      bb    = aa*aa
      dbb   = 2.0d0 * aa * daa

      cc    =  1.78d+05 * t9i23 * exp(-16.669*t9i13)
      dcc   = cc*(-twoth*t9i + oneth*16.669*t9i43)

      zz    = 1.0d0/bb
      res3  = cc*zz
      dres3 = (dcc - res3*dbb)*zz

      res4  = 8.68d+10 * t9 * exp(-16.667*t9i13 - t92*q2)
      dres4 = res4*(t9i + oneth*16.667*t9i43 - 2.0d0*t9*q2)

      res5  = 9.22d-04 * t9i32 * exp(-0.767*t9i)
      dres5 = res5*(-1.5d0*t9i + 0.767*t9i2)

      res6  = theta * 98.0 * t9i32 * exp(-2.077*t9i)
      dres6 = res6*(-1.5d0*t9i + 2.077*t9i2)

      term    = res1 + res2 + res3 + res4 + res5 + res6
      dtermdt = dres1 + dres2 + dres3 + dres4 + dres5 + dres6

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.676 * exp(-13.825*t9i)
      drevdt   = rev*13.825*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end






      subroutine rate_o17pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg, &
                       t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       zz,theta
      parameter        (theta = 0.1d0)


! o17(p,g)18f
! from landre et al 1990 a&a 240, 85
      aa     = 1.0d0 + 2.69*t9
      zz     = 1.0d0/aa

      t9a    = t9*zz
      dt9a   = (1.0d0 - t9a*2.69)*zz

      zz     = dt9a/t9a
      t9a13  = t9a**oneth
      dt9a13 = oneth*t9a13*zz

      t9a56  = t9a**fivsix
      dt9a56 = fivsix*t9a56*zz

      aa  = 7.97d+07 * t9a56 * t9i32 * exp(-16.712/t9a13)
      daa = aa*(dt9a56/t9a56 - 1.5d0*t9i + 16.712/t9a13**2*dt9a13)

      bb  = 1.0d0  + 0.025*t913 - 0.051*t923 - 8.82d-3*t9
      dbb = oneth*0.025*t9i23 - twoth*0.051*t9i13 - 8.82d-3
      if (bb .le. 0.0) then
       bb  = 0.0d0
       dbb = 0.0d0
      end if

      cc  = 1.51d+08 * t9i23 * exp(-16.712*t9i13)
      dcc = cc*(-twoth*t9i + oneth*16.712*t9i43)

      dd  = bb*cc
      ddd = dbb*cc + bb*dcc

      ee  = 1.56d+5 * t9i * exp(-6.272*t9i)
      dee = ee*(-t9i + 6.272*t9i2)

      ff  = 2.0d0 * theta * 3.16d-05 * t9i32 * exp(-0.767*t9i)
      dff = ff*(-1.5d0*t9i + 0.767*t9i2)

      gg  = theta * 98.0 * t9i32 * exp(-2.077*t9i)
      dgg = gg*(-1.5d0*t9i + 2.077*t9i2)

      term    = aa + dd + ee + ff + gg
      dtermdt = daa + ddd + dee + dff + dgg

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.66d+10 * t932 * exp(-65.061*t9i)
      drevdt   = rev*(1.5d0*t9i + 65.061*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_o18pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,q1
      parameter        (q1 = 1.0d0/1.852321d0)


! o18(p,a)n15
      aa  = 3.63e+11 * t9i23 * exp(-16.729*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*16.729*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.025*t913 + 1.88*t923 + 0.327*t9 &
            + 4.66*t943 + 2.06*t953
      dbb = oneth*0.025*t9i23 + twoth*1.88*t9i13 + 0.327 &
            + fourth*4.66*t913 + fiveth*2.06*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 9.90e-14 * t9i32 * exp(-0.231*t9i)
      ddd = -1.5d0*dd*t9i + dd*0.231*t9i2

      ee  = 2.66e+04 * t9i32 * exp(-1.670*t9i)
      dee = -1.5d0*ee*t9i + ee*1.670*t9i2

      ff  = 2.41e+09 * t9i32 * exp(-7.638*t9i)
      dff = -1.5d0*ff*t9i + ff*7.638*t9i2

      gg  = 1.46e+09 * t9i * exp(-8.310*t9i)
      dgg = -gg*t9i + gg*8.310*t9i2

      term    = cc + dd + ee + ff + gg
      dtermdt = dcc + ddd + dee + dff + dgg


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.66e-01 * exp(-46.191*t9i)
      drevdt   = rev*46.191*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o18pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.019321d0)


! o18(p,g)19f
      aa  = 3.45e+08 * t9i23 * exp(-16.729*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*16.729*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.025*t913 + 2.26*t923 + 0.394*t9 &
            + 30.56*t943 + 13.55*t953
      dbb = oneth*0.025*t9i23 + twoth*2.26*t9i13 + 0.394 &
            + fourth*30.56*t913 + fiveth*13.55*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 1.25e-15 * t9i32 * exp(-0.231*t9i)
      ddd = dd*(-1.5d0*t9i + 0.231*t9i2)

      ee  = 1.64e+02 * t9i32 * exp(-1.670*t9i)
      dee = ee*(-1.5d0*t9i + 1.670*t9i2)

      ff  = 1.28e+04 * t912 * exp(-5.098*t9i)
      dff = ff*(0.5d0*t9i + 5.098*t9i2)

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 9.20e+09 * t932 * exp(-92.769*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.769*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_f17em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 64.49d0, &
                        con      = lntwo/halflife)

! f17(e-nu)o17
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_f18em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 6586.2d0, &
                        con      = lntwo/halflife)

! f18(e-nu)o18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_f19pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,q1
      parameter        (q1 = 1.0d0/0.714025d0)


! f19(p,a)o16
      aa  = 3.55d+11 * t9i23 * exp(-18.113*t9i13 - t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*18.113*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.023*t913 + 1.96*t923 + 0.316*t9 &
            + 2.86*t943 + 1.17*t953
      dbb = oneth*0.023*t9i23 + twoth*1.96*t9i13 + 0.316 &
            + fourth*2.86*t913 + fiveth*1.17*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.67d+06 * t9i32 * exp(-3.752*t9i)
      ddd = -1.5d0*dd*t9i + dd*3.752*t9i2

      ee  = 3.07d+08 * exp(-6.019*t9i)
      dee = ee*6.019*t9i2

      ff  = 4.0*exp(-2.090*t9i)
      dff = ff*2.090*t9i2

      gg  = 7.0*exp(-16.440*t9i)
      dgg = gg*16.440*t9i2

      hh  = 1.0d0 + ff + gg
      dhh = dff + dgg

      term    = (cc + dd + ee)/hh
      dtermdt = ((dcc + ddd + dee) - term*dhh)/hh


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.54e-01 * exp(-94.159*t9i)
      drevdt   = rev*94.159*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_n13pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,q1
      parameter        (q1 = 1.0d0/0.69288976d0)


! n13(p,g)o14
! Keiner et al 1993 Nucl Phys A552, 66
      aa  = -1.727d+7 * t9i23 * exp(-15.168*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*15.168*t9i43 -2.0d0*t9*q1)

      bb  = 1.0d0 + 0.027*t913 - 17.54*t923 - 3.373*t9 &
            + 0.0176*t943 + 0.766d-2*t953
      dbb = oneth*0.027*t9i23 - twoth*17.54*t9i13 - 3.373 &
            + fourth*0.0176*t913 + fiveth*0.766d-2*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 3.1d+05 * t9i32 * exp(-6.348*t9i)
      ddd = dd*(-1.5d0*t9i + 6.348*t9i2)

      term    = cc + dd
      dtermdt = dcc + ddd

! goes negative below about t7=1.5
! note cf88 rate stays positive
      if (term .lt. 0.0) then
       term    = 0.0d0
       dtermdt = 0.0d0
      end if


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.57d+10*t932*exp(-53.706*t9i)
      drevdt   = rev*(1.5d0*t9i + 53.706*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_o14em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 70.606d0, &
                        con      = lntwo/halflife)

! o14(e-nu)n14
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_o14ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,q1
      parameter        (q1 = 1.0d0/0.514089d0)


! o14(a,p)f17
      aa  = 1.68e+13 * t9i23 * exp(-39.388*t9i13- t92*q1)
      daa = -twoth*aa*t9i + aa*(oneth*39.388*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.011*t913 + 13.117*t923 + 0.971*t9 &
            + 85.295*t943 + 16.061*t953
      dbb = oneth*0.011*t9i23 + twoth*13.117*t9i13 + 0.971 &
            + fourth*85.295*t913 + fiveth*16.061*t923

      cc  = aa * bb
      dcc = daa*bb + aa*dbb

      dd  = 3.31e+04 * t9i32 * exp(-11.733*t9i)
      ddd = -1.5d0*dd*t9i + dd*11.733*t9i2

      ee  = 1.79e+07 * t9i32 * exp(-22.609*t9i)
      dee = -1.5d0*ee*t9i + ee*22.609*t9i2

      ff  = 9.00e+03 * t9113 * exp(-12.517*t9i)
      dff = elvnth*ff*t9i + ff*12.517*t9i2

      term    = cc + dd + ee + ff
      dtermdt = dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.93e-01*exp(-13.820*t9i)
      drevdt   = rev*13.820*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end




      subroutine rate_o15ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh, &
                       q1,q2,q3
      parameter        (q1 = 1.0d0/9.0d0, &
                        q2 = 1.0d0/3.751969d0, &
                        q3 = 1.0d0/64.0d0)


! o15(a,g)ne19

      aa  = 3.57d+11 * t9i23 * exp(-39.584d+0*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*39.584d0*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.011*t913 - 0.273*t923 - 0.020*t9
      dbb = oneth*0.011*t9i23 - twoth*0.273*t9i13 - 0.020

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 5.10d+10 * t9i23 * exp(-39.584d+0*t9i13 - t92*q2)
      ddd = dd*(-twoth*t9i + oneth*39.584*t9i43 - 2.0d0*t9*q2)

      ee  = 1.0d0 + 0.011*t913 + 1.59*t923 + 0.117*t9 &
            + 1.81*t943 + 0.338*t953
      dee = oneth*0.011*t9i23 + twoth*1.59*t9i13 + 0.117 &
            + fourth*1.81*t913 + fiveth*0.338*t923

      ff  = dd*ee
      dff = ddd*ee + dd*dee

      gg  = 3.95d-1 * t9i32 * exp(-5.849*t9i)
      dgg = gg*(-1.5d0*t9i + 5.849*t9i2)

      hh  = 1.90d+1 * t9**2.85 * exp(-7.356*t9i - t92*q3)
      dhh = hh*(2.85*t9i + 7.356*t9i2 - 2.0d0*t9*q3)


      term    = cc + ff + gg + hh
      dtermdt = dcc + dff + dgg + dhh


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.54e+10 * t932 * exp(-40.957*t9i)
      drevdt   = rev*(1.5d0*t9i + 40.957*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_f17pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee


! f17(p,g)ne18
! wiescher and kettner, ap. j., 263, 891 (1982)

      aa  = 1.66e+07 * t9i23 * exp(-18.03*t9i13)
      daa = aa*(-twoth*t9i + oneth*18.03*t9i43)

      bb  = 2.194 + 0.050*t913 - 0.376*t923 - 0.061*t9 &
            + 0.026*t943 + 0.011*t953
      dbb = oneth*0.050*t9i23 - twoth*0.376*t9i13 - 0.061 &
            + fourth*0.026*t913 + fiveth*0.011*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 839.0 * t9i32 * exp(-6.93*t9i)
      ddd = dd*(-1.5d0*t9i + 6.93*t9i2)

      ee  = 33.56 * t9i32 * exp(-7.75*t9i)
      dee = ee*(-1.5d0*t9i + 7.75*t9i2)

      term    = cc + dd + ee
      dtermdt = dcc + ddd + dee

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.087e+11 * t932 * exp(-45.501*t9i)
      drevdt   = rev*(1.5d0*t9i + 45.501*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne18em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 1.672d0, &
                        con      = lntwo/halflife)

! ne18(e-nu)f18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_f18pa(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff


! f18(p,a)o15
! wiescher and kettner, ap. j., 263, 891 (1982)

      aa  = 1.66e-10 * t9i32 * exp(-0.302*t9i)
      daa = aa*(-1.5d0*t9i + 0.302*t9i2)

      bb  = 1.56e+05 * t9i32 * exp(-3.84*t9i)
      dbb = bb*(-1.5d0*t9i + 3.84*t9i2)

      cc  = 1.36e+06 * t9i32 * exp(-5.22*t9i)
      dcc = cc*(-1.5d0*t9i + 5.22*t9i2)

      dd  = 8.1e-05 * t9i32 * exp(-1.05*t9i)
      ddd = dd*(-1.5d0*t9i + 1.05*t9i2)

      ee  = 8.9e-04 * t9i32 * exp(-1.51*t9i)
      dee = ee*(-1.5d0*t9i + 1.51*t9i2)

      ff  = 3.0e+05 * t9i32 * exp(-4.29*t9i)
      dff = ff*(-1.5d0*t9i + 4.29*t9i2)

      term    = aa + bb + cc + dd + ee + ff
      dtermdt = daa + dbb + dcc + ddd + dee + dff


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 4.93e-01 * exp(-33.433*t9i)
      drevdt   = rev*33.433*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_ne18ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,zz

      double precision z1,a1,ztot,ared,r,c1,c2,c3,c4
      parameter        (z1   = 10.0d0, &
                        a1   = 18.0d0, &
                        ztot = 2.0d0 * z1, &
                        ared = 4.0d0*a1/(4.0d0 + a1), &
                        r    = 5.1566081196876965d0, &
                        c1   = 4.9080044545315392d10, &
                        c2   = 4.9592784569936502d-2, &
                        c3   = 1.9288564401521285d1, &
                        c4   = 4.6477847042196437d1)

! note:
!      r    = 1.09 * a1**oneth + 2.3
!      c1   = 7.833e9 * 0.31 * ztot**fourth/(ared**fivsix)
!      c2   = 0.08617 * 0.1215 * sqrt(ared*r**3/ztot)
!      c3   = 2.0d0 * 0.52495 * sqrt(ared*r*ztot)
!      c4   = 4.2487 * (ztot**2*ared)**oneth


! ne18ap(a,p)na21
! was a call to aprate

      aa  = 1.0d0 + c2*t9
      zz  = c2/aa

      bb  = aa**fivsix
      dbb = fivsix*bb*zz

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923 * dbb

      dd = aa**oneth
      ddd = oneth*dd*zz

      ee  = t9i13 * dd
      dee = -oneth*ee*t9i + t9i13 * ddd

      zz      = 1.0d0/cc
      term    = c1*zz * exp(c3 - c4*ee)
      dtermdt = -term*(zz*dcc + c4*dee)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev    = 0.0d0
      drevdt = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne19pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,ff,dff,gg,dgg,q1
      parameter        (q1 = 1.0d0/1.304164d0)


! ne19(p,g)na20

      aa  = 1.71d+6 * t9i23 * exp(-19.431d0*t9i13)
      daa = aa*(-twoth*t9i + oneth*19.431*t9i43)

      bb  = 1.0d0 + 0.021*t913 + 0.130*t923 + 1.95d-2*t9 &
            + 3.86d-2*t943 + 1.47d-02*t953
      dbb = oneth*0.021*t9i23 + twoth*0.130*t9i13 + 1.95d-2 &
            + fourth*3.86d-2*t913 + fiveth*1.47d-2*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb


      dd  = 1.89d+5 * t9i23 * exp(-19.431d0*t9i13 - t92*q1)
      ddd = dd*(-twoth*t9i + oneth*19.431*t9i43 - 2.0d0*t9*q1)

      ee  = 1.0d0 + 0.021*t913 + 2.13*t923 + 0.320*t9 &
            + 2.80*t943 + 1.07*t953
      dee = oneth*0.021*t9i23 + twoth*2.13*t9i13 + 0.320 &
            + fourth*2.80*t913 + fiveth*1.07*t923

      ff  = dd*ee
      dff = ddd*ee + dd*dee

      gg  = 8.45d+3 * t9i54 * exp(-7.64d0*t9i)
      dgg = gg*(-fivfour*t9i + 7.64d0*t9i2)


      term    = cc + ff + gg
      dtermdt = dcc + dff + dgg


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.39e+09 * t932 * exp(-25.519*t9i)
      drevdt   = rev*(1.5d0*t9i + 25.519*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ne19em(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision lntwo,halflife,con
      parameter        (lntwo    = 0.6931471805599453d0, &
                        halflife = 17.3982d0, &
                        con      = lntwo/halflife)

! ne18(e-nu)f18
      fr    = con
      dfrdt = 0.0d0
      dfrdd = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_si26ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,bb,dbb, &
                       cc,dcc,dd,ddd,ee,dee,zz

      double precision z1,a1,ztot,ared,r,c1,c2,c3,c4
      parameter        (z1   = 14.0d0, &
                        a1   = 26.0d0, &
                        ztot = 2.0d0 * z1, &
                        ared = 4.0d0*a1/(4.0d0 + a1), &
                        r    = 5.5291207145640335d0, &
                        c1   = 7.3266779970543091d10, &
                        c2   = 4.7895369289991982d-02, &
                        c3   = 2.4322657793918662d1, &
                        c4   = 5.9292366232997814d1)

! note:
!      r    = 1.09 * a1**oneth + 2.3
!      c1   = 7.833e9 * 0.31 * ztot**fourth/(ared**fivsix)
!      c2   = 0.08617 * 0.1215 * sqrt(ared*r**3/ztot)
!      c3   = 2.0d0 * 0.52495 * sqrt(ared*r*ztot)
!      c4   = 4.2487 * (ztot**2*ared)**oneth



! si26ap(a,p)p29
! was a call to aprate

      aa  = 1.0d0 + c2*t9
      zz  = c2/aa

      bb  = aa**fivsix
      dbb = fivsix*bb*zz

      cc  = t923 * bb
      dcc = twoth*cc*t9i + t923 * dbb

      dd = aa**oneth
      ddd = oneth*dd*zz

      ee  = t9i13 * dd
      dee = -oneth*ee*t9i + t9i13 * ddd

      zz      = 1.0d0/cc
      term    = c1*zz * exp(c3 - c4*ee)
      dtermdt = -term*(zz*dcc + c4*dee)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.0d0
      drevdt   = 0.0d0

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end





      subroutine rate_c12ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg,hh,dhh,f1,df1,f2,df2, &
                       zz,q1
      parameter        (q1 = 1.0d0/12.222016d0)


! c12(a,g)o16
      aa   = 1.0d0 + 0.0489d0*t9i23
      daa  = -twoth*0.0489d0*t9i53

      bb   = t92*aa*aa
      dbb  = 2.0d0*(bb*t9i + t92*aa*daa)

      cc   = exp(-32.120d0*t9i13 - t92*q1)
      dcc  = cc * (oneth*32.120d0*t9i43 - 2.0d0*t9*q1)

      dd   = 1.0d0 + 0.2654d0*t9i23
      ddd  = -twoth*0.2654d0*t9i53

      ee   = t92*dd*dd
      dee  = 2.0d0*(ee*t9i + t92*dd*ddd)

      ff   = exp(-32.120d0*t9i13)
      dff  = ff * oneth*32.120d0*t9i43

      gg   = 1.25d3 * t9i32 * exp(-27.499*t9i)
      dgg  = gg*(-1.5d0*t9i + 27.499*t9i2)

      hh   = 1.43d-2 * t95 * exp(-15.541*t9i)
      dhh  = hh*(5.0d0*t9i + 15.541*t9i2)

      zz   = 1.0d0/bb
      f1   = cc*zz
      df1  = (dcc - f1*dbb)*zz

      zz   = 1.0d0/ee
      f2   = ff*zz
      df2  = (dff - f2*dee)*zz

      term    = 1.04d8*f1  + 1.76d8*f2 + gg + hh
      dtermdt = 1.04d8*df1 + 1.76d8*df2 + dgg + dhh


! 1.7 times cf88 value
      term     = 1.7d0 * term
      dtermdt  = 1.7d0 * dtermdt

      fr    = term * den
      dfrdt = dtermdt * den * 1.0d-9
      dfrdd = term

      rev    = 5.13d10 * t932 * exp(-83.111*t9i)
      drevdt = rev*(1.5d0*t9i + 83.111*t9i2)

      rr     = rev * term
      drrdt  = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd  = 0.0d0

      return
      end






      subroutine rate_tripalf(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,r2abe,dr2abedt,rbeac, &
                       drbeacdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,xx,dxx,yy,dyy,zz,dzz,uu,vv,f1,df1,rc28, &
                       q1,q2
      parameter        (rc28   = 0.1d0, &
                        q1     = 1.0d0/0.009604d0, &
                        q2     = 1.0d0/0.055225d0)



! triple alfa to c12
! this is a(a,g)be8
      aa    = 7.40d+05 * t9i32 * exp(-1.0663*t9i)
      daa   = aa*(-1.5d0*t9i  + 1.0663*t9i2)

      bb    = 4.164d+09 * t9i23 * exp(-13.49*t9i13 - t92*q1)
      dbb   = bb*(-twoth*t9i + oneth*13.49*t9i43 - 2.0d0*t9*q1)

      cc    = 1.0d0 + 0.031*t913 + 8.009*t923 + 1.732*t9 &
              + 49.883*t943 + 27.426*t953
      dcc   = oneth*0.031*t9i23 + twoth*8.009*t9i13 + 1.732 &
              + fourth*49.883*t913 + fiveth*27.426*t923

      r2abe    = aa + bb * cc
      dr2abedt = daa + dbb*cc + bb*dcc


! this is be8(a,g)c12
      dd    = 130.0d0 * t9i32 * exp(-3.3364*t9i)
      ddd   = dd*(-1.5d0*t9i + 3.3364*t9i2)

      ee    = 2.510d+07 * t9i23 * exp(-23.57*t9i13 - t92*q2)
      dee   = ee*(-twoth*t9i + oneth*23.57*t9i43 - 2.0d0*t9*q2)

      ff    = 1.0d0 + 0.018*t913 + 5.249*t923 + 0.650*t9 + &
              19.176*t943 + 6.034*t953
      dff   = oneth*0.018*t9i23 + twoth*5.249*t9i13 + 0.650 &
              + fourth*19.176*t913 + fiveth*6.034*t923

      rbeac    = dd + ee * ff
      drbeacdt = ddd + dee * ff + ee * dff


! a factor
      xx    = rc28 * 1.35d-07 * t9i32 * exp(-24.811*t9i)
      dxx   = xx*(-1.5d0*t9i + 24.811*t9i2)


! high temperature rate
      if (t9.gt.0.08) then
       term    = 2.90d-16 * r2abe * rbeac + xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac &
                 + 2.90d-16 * r2abe * drbeacdt &
                 + dxx


! low temperature rate
      else
       uu   = 0.8d0*exp(-(0.025*t9i)**3.263)
       yy   = 0.2d0 + uu
!       yy   = 0.01 + 0.2d0 + uu
       dyy  = uu * 3.263*(0.025*t9i)**2.263 * (0.025*t9i2)
       vv   = 4.0d0*exp(-(t9/0.025)**9.227)
       zz   = 1.0d0 + vv
       dzz  = vv * 9.227*(t9/0.025)**8.227 * 40.0d0
       aa   = 1.0d0/zz
       f1   = 0.01d0 + yy * aa
!       f1   = yy * aa
       df1  = (dyy - f1*dzz)*aa
       term = 2.90d-16 * r2abe * rbeac * f1 +  xx
       dtermdt =   2.90d-16 * dr2abedt * rbeac * f1 &
                 + 2.90d-16 * r2abe * drbeacdt * f1 &
                 + 2.90d-16 * r2abe * rbeac * df1 &
                 + dxx
      end if


! rates
!      term    = 1.2d0 * term
!      dtermdt = 1.2d0 * term

      fr    = term * den * den
      dfrdt = dtermdt * den * den * 1.0d-9
      dfrdd = 2.0d0 * term * den

      rev    = 2.00d+20*t93*exp(-84.424*t9i)
      drevdt = rev*(3.0d0*t9i + 84.424*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_c12c12(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a56,dt9a56, &
                       aa,zz


! c12 + c12 reaction
      aa      = 1.0d0 + 0.0396*t9
      zz      = 1.0d0/aa

      t9a     = t9*zz
      dt9a    = (1.0d0 -  t9a*0.0396)*zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      term    = 4.27d+26 * t9a56 * t9i32 * &
                exp(-84.165/t9a13 - 2.12d-03*t93)
      dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                      + 84.165/t9a13**2*dt9a13 - 6.36d-3*t92)

! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_c12c12npa(temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd

! locals
      double precision term,dtermdt,rev,drevdt,t9a,dt9a,t9a13,dt9a13, &
                       t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,dd,ddd,zz, &
                       b24n,db24n,b24p,db24p,b24a,db24a


! c12(c12,n)mg23
! c12(c12,p)na23
! c12(c12,a)ne20


      aa      = 1.0d0 + 0.0396*t9
      zz      = 1.0d0/aa

      t9a     = t9 * zz
      dt9a    = (1.0d0 -  t9a*0.0396) * zz

      zz      = dt9a/t9a
      t9a13   = t9a**oneth
      dt9a13  = oneth*t9a13*zz

      t9a56   = t9a**fivsix
      dt9a56  = fivsix*t9a56*zz

      aa = 4.27d+26 * t9a56 * t9i32 * &
           exp(-84.165/t9a13 - 2.12d-03*t93)

      daa = aa * (dt9a56/t9a56 - 1.5d0*t9i &
                + 84.165/t9a13**2*dt9a13 - 6.36d-3*t92)


! neutron branching from dayras switkowski and woosley 1976
      if (t9 .ge. 1.5) then

       bb    =  0.055 * exp(0.976 - 0.789*t9)
       dbb   = -bb*0.789

       b24n  = 0.055  - bb
       db24n = -dbb

      else

       bb    = 1.0d0 + 0.0789*t9 + 7.74*t92
       dbb   = 0.0789 + 2.0d0*7.74*t9

       cc    = 0.766*t9i3
       dcc   = -3.0d0*cc*t9i

       dd    = bb * cc
       ddd   = dbb*cc + bb*dcc

       b24n  = 0.859*exp(-dd)
       db24n = -b24n*ddd
      end if


! proton branching ratio
      if (t9.gt.3.) then
        b24p  = oneth*(1.0d0 - b24n)
        db24p = -oneth*db24n

        b24a  = 2.0d0 * b24p
        db24a = 2.0d0 * db24p

       else
        b24p  = 0.5d0*(1.0d0 - b24n)
        db24p = -0.5d0*db24n

        b24a  = b24p
        db24a = db24p

       end if

! rates

! c12(c12,n)mg23
      term    = aa * b24n
      dtermdt = daa*b24n + aa*db24n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.93 * exp(30.16100515d0*t9i)
       drevdt = -rev*30.16100515d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr1dd = rev*term


! c12(c12,p)na23
      term    = aa * b24p
      dtermdt = daa*b24p + aa*db24p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.93 * exp(-25.98325915d0*t9i)
       drevdt = rev*25.98325915d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


! c12(c12,a)ne20
      term    = aa * b24a
      dtermdt = daa*b24a + aa*db24a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 2.42 * exp(-53.576110995d0*t9i)
       drevdt = rev*53.576110995d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term

      return
      end






      subroutine rate_c12o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,t9a,dt9a,t9a13,dt9a13,t9a23,dt9a23, &
                       t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc,zz


! c12 + o16 reaction; see cf88 references 47-4

      if (t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.055)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0d0 * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa      = exp(-0.18*t9a*t9a)
       daa     = -aa * 0.36 * t9a * dt9a

       bb      = 1.06d-03*exp(2.562*t9a23)
       dbb     = bb * 2.562 * dt9a23

       cc      = aa + bb
       dcc     = daa + dbb

       zz      = 1.0d0/cc
       term    = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) * zz
       dtermdt = term*(dt9a56/t9a56 - 1.5d0*t9i &
                       + 106.594/t9a23*dt9a13 - zz*dcc)

      else
!       term    = 2.6288035d-29
       term    = 0.0d0
       dtermdt = 0.0d0
      endif


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end




      subroutine rate_c12o16npa(temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd

! locals
      double precision term,dtermdt,rev,drevdt,t9a,dt9a,t9a13,dt9a13, &
                       t9a23,dt9a23,t9a56,dt9a56,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,b27n,b27p,b24a,zz


! c12(o16,n)si27
! c12(o16,p)al27
! c12(o16,a)mg24


      if (t9.ge.0.5) then
       aa     = 1.0d0 + 0.055*t9
       zz     = 1.0d0/aa

       t9a    = t9*zz
       dt9a   = (1.0d0 - t9a*0.055)*zz

       zz     = dt9a/t9a
       t9a13  = t9a**oneth
       dt9a13 = oneth*t9a13*zz

       t9a23  = t9a13*t9a13
       dt9a23 = 2.0d0 * t9a13 * dt9a13

       t9a56  = t9a**fivsix
       dt9a56 = fivsix*t9a56*zz

       aa     = exp(-0.18*t9a*t9a)
       daa    = -aa * 0.36 * t9a * dt9a

       bb     = 1.06d-03*exp(2.562*t9a23)
       dbb    = bb * 2.562 * dt9a23

       cc     = aa + bb
       dcc    = daa + dbb

       zz     = 1.0d0/cc
       dd     = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) *zz
       ddd    = dd*(dt9a56/t9a56 - 1.5d0*t9i &
                 + 106.594/t9a23 * dt9a13 - dcc*zz)

      else
!       dd     = 2.6288035d-29
       dd     = 0.0d0
       ddd    = 0.0d0
      endif


! branching ratios from pwnsz data
        b27n = 0.1d0
        b27p = 0.5d0
        b24a = 0.4d0


! rates

! c12(o16,n)si27
      term    = dd * b27n
      dtermdt = ddd * b27n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 1.58d0 * exp(4.8972467d0*t9i)
       drevdt = -rev*4.8972467d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = den*(drevdt*term + rev*dtermdt)*1.0d-9
      drr1dd = rev*term


! c12(o16,p)al27
      term    = dd * b27p
      dtermdt = ddd * b27p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 1.58d0 * exp(-59.9970745d0*t9i)
       drevdt = rev*59.9970745d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


! c12(o16,a)mg24
      term    = dd * b24a
      dtermdt = ddd * b24a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 2.83d0 * exp(-78.5648345d0*t9i)
       drevdt = rev*78.5648345d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term


      return
      end








      subroutine rate_o16o16(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt


! o16 + o16
      term  = 7.10d36 * t9i23 * &
              exp(-135.93 * t9i13 - 0.629*t923 &
                   - 0.445*t943 + 0.0103*t9*t9)

      dtermdt = -twoth*term*t9i &
                + term * (oneth*135.93*t9i43 - twoth*0.629*t9i13 &
                          - fourth*0.445*t913 + 0.0206*t9)


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rr    = 0.0d0
      drrdt = 0.0d0
      drrdd = 0.0d0

      return
      end






      subroutine rate_o16o16npad(temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd, &
                      fr4,dfr4dt,dfr4dd,rr4,drr4dt,drr4dd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                      fr1,dfr1dt,dfr1dd,rr1,drr1dt,drr1dd, &
                      fr2,dfr2dt,dfr2dd,rr2,drr2dt,drr2dd, &
                      fr3,dfr3dt,dfr3dd,rr3,drr3dt,drr3dd, &
                      fr4,dfr4dt,dfr4dd,rr4,drr4dt,drr4dd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa, &
                       b32n,db32n,b32p,db32p,b32a,db32a,b32d,db32d, &
                       ezro,dezro,dlt,ddlt,xxt,dxxt,thrs,dthrs


! o16(o16,n)s31
! o16(o16,p)p31
! o16(o16,a)si28
! o16(o16,d)p30

      aa  = 7.10d36 * t9i23 * exp(-135.93 * t9i13 - 0.629*t923 &
                   - 0.445*t943 + 0.0103*t9*t9)

      daa = -twoth*aa*t9i + aa * (oneth*135.93*t9i43 - twoth*0.629*t9i13 &
                          - fourth*0.445*t913 + 0.0206*t9)


! branching ratios highly uncertain;  guessed using fcz 1975
! deuteron channel is endoergic. apply error function cut-off.
       ezro = 3.9*t923
       dezro = twoth*ezro*t9i

       dlt  = 1.34*t9**fivsix
       ddlt = fivsix*dlt*t9i

       xxt  = 2.0d0*(2.406 - ezro)/dlt
       dxxt = -(2.0d0*dezro + xxt*ddlt)/dlt

       call fowthrsh(xxt,thrs,dthrs)

       b32d  = 0.05d0*thrs
       db32d = 0.05d0*dthrs*dxxt

! shut down deutrium channel
!       b32d  = 0.0d0
!       db32d = 0.0d0

       b32n  = 0.1d0
       db32n = 0.0d0

       b32a  = 0.25d0
       db32a = 0.0d0

       b32p  = 1.0d0 - b32d - b32a - b32n
       db32p = -db32d



! rates

! o16(o16,n)s31
      term    = aa * b32n
      dtermdt = daa*b32n + aa*db32n
      fr1     = den * term
      dfr1dt  = den * dtermdt * 1.0d-9
      dfr1dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 5.92 * exp(-16.8038228d0*t9i)
       drevdt = rev*16.8038228d0*t9i2
      end if
      rr1    = den * rev * term
      drr1dt = den*(drevdt*term + rev*dtermdt) * 1.0d-9
      drr1dd = rev*term

! o16(o16,p)p31
      term    = aa * b32p
      dtermdt = daa*b32p + aa*db32p
      fr2     = den * term
      dfr2dt  = den * dtermdt * 1.0d-9
      dfr2dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 5.92*exp(-89.0788286d0*t9i)
       drevdt = rev*89.0788286d0*t9i2
      end if
      rr2    = den * rev * term
      drr2dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr2dd = rev*term


! o16(o16,a)si28
      term    = aa * b32a
      dtermdt = daa*b32a + aa*db32a
      fr3     = den * term
      dfr3dt  = den * dtermdt * 1.0d-9
      dfr3dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 3.46*exp(-111.3137212d0*t9i)
       drevdt = rev*111.3137212d0*t9i2
      end if
      rr3    = den * rev * term
      drr3dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr3dd = rev*term


! o16(o16,d)p30
      term    = aa * b32d
      dtermdt = daa*b32d + aa*db32d
      fr4     = den * term
      dfr4dt  = den * dtermdt * 1.0d-9
      dfr4dd  = term

      rev    = 0.0d0
      drevdt = 0.0d0
      if (t9 .gt. 0.1) then
       rev    = 0.984*exp(27.9908982d0*t9i)
       drevdt = -rev*27.9908982d0*t9i2
      end if
      rr4    = den * rev * term
      drr4dt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drr4dd = rev*term

      return
      end






      subroutine fowthrsh(x,thrs,dthrs)
      include 'implno.dek'

! fowler threshold fudge function.
! err func rational (abramowitz p.299)7.1.25 and its derivative

! declare
      double precision x,thrs,dthrs, &
                       ag,z,z2,t,t2,t3,tt,er,aa,daa,dt,dtt,der

      ag   = sign(1.0d0,x)
      z    = abs(x)
      z2   = z*z
      aa   = 1.0d0 + 0.47047d0*z

      t    = 1.0d0/aa
      dt   = -t*t*0.47047*ag

      t2   = t*t
      t3   = t2*t

      tt   = 0.3480242d0*t - 0.0958798d0*t2 + 0.7478556d0*t3
      dtt  = dt * (0.3480242d0 - 2.0d0*0.0958798d0*t &
                   + 3.0d0*0.7478556d0*t2)

      thrs  = 0.5d0
      dthrs = -5.6452433937900004d-1

      if (z .ne. 0) then
       aa   = exp(-z2)
       daa  = -2.0d0*aa*z*ag

       er   = 1.0d0 - tt * aa
       der  = -dtt*aa - tt*daa

       thrs  = 0.5d0 * (1.0d0 - ag*er)
       dthrs = -0.5d0*ag*der
      end if
      return
      end





      subroutine rate_o16ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                       cc,dcc,term2,dterm2,rev,drevdt,q1
      parameter        (q1 = 1.0d0/2.515396d0)



! o16(a,g)ne20
      term1   = 9.37d9 * t9i23 * exp(-39.757*t9i13 - t92*q1)
      dterm1  = term1*(-twoth*t9i + oneth*39.757*t9i43 - 2.0d0*t9*q1)

      aa      = 62.1 * t9i32 * exp(-10.297*t9i)
      daa     = aa*(-1.5d0*t9i + 10.297*t9i2)

      bb      = 538.0d0 * t9i32 * exp(-12.226*t9i)
      dbb     = bb*(-1.5d0*t9i + 12.226*t9i2)

      cc      = 13.0d0 * t92 * exp(-20.093*t9i)
      dcc     = cc*(2.0d0*t9i + 20.093*t9i2)

      term2   = aa + bb + cc
      dterm2  = daa + dbb + dcc

      term    = term1 + term2
      dtermdt = dterm1 + dterm2


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 5.65d+10*t932*exp(-54.937*t9i)
      drevdt   = rev*(1.5d0*t9i + 54.937*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_ne20ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,term1,dterm1,aa,daa,bb,dbb, &
                       term2,dterm2,term3,dterm3,rev,drevdt,zz,rc102,q1
      parameter        (rc102 = 0.1d0, &
                        q1    = 1.0d0/4.923961d0)



! ne20(a,g)mg24

      aa   = 4.11d+11 * t9i23 * exp(-46.766*t9i13 - t92*q1)
      daa  = aa*(-twoth*t9i + oneth*46.766*t9i43 - 2.0d0*t9*q1)

      bb   = 1.0d0 + 0.009*t913 + 0.882*t923 + 0.055*t9 &
             + 0.749*t943 + 0.119*t953
      dbb  = oneth*0.009*t9i23 + twoth*0.882*t9i13 + 0.055 &
             + fourth*0.749*t913 + fiveth*0.119*t923

      term1  = aa * bb
      dterm1 = daa * bb + aa * dbb


      aa   = 5.27d+03 * t9i32 * exp(-15.869*t9i)
      daa  = aa*(-1.5d0*t9i + 15.869*t9i2)

      bb   = 6.51d+03 * t912 * exp(-16.223*t9i)
      dbb  = bb*(0.5d0*t9i + 16.223*t9i2)

      term2  = aa + bb
      dterm2 = daa + dbb


      aa   = 42.1 * t9i32 * exp(-9.115*t9i)
      daa  = aa*(-1.5d0*t9i + 9.115*t9i2)

      bb   =  32.0 * t9i23 * exp(-9.383*t9i)
      dbb  = bb*(-twoth*t9i + 9.383*t9i2)

      term3  = rc102 * (aa + bb)
      dterm3 = rc102 * (daa + dbb)


      aa  = 5.0d0*exp(-18.960*t9i)
      daa = aa*18.960*t9i2

      bb  = 1.0d0 + aa
      dbb = daa

      zz      = 1.0d0/bb
      term    = (term1 + term2 + term3)*zz
      dtermdt = ((dterm1 + dterm2 + dterm3) - term*dbb)*zz


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.01d+10 * t932 * exp(-108.059*t9i)
      drevdt   = rev*(1.5d0*t9i + 108.059*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_mg24ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,hh,hhi,rev,drevdt,rc121
      parameter        (rc121 = 0.1d0)



! 24mg(a,g)28si

      aa    = 4.78d+01 * t9i32 * exp(-13.506*t9i)
      daa   = aa*(-1.5d0*t9i + 13.506*t9i2)

      bb    =  2.38d+03 * t9i32 * exp(-15.218*t9i)
      dbb   = bb*(-1.5d0*t9i + 15.218*t9i2)

      cc    = 2.47d+02 * t932 * exp(-15.147*t9i)
      dcc   = cc*(1.5d0*t9i + 15.147*t9i2)

      dd    = rc121 * 1.72d-09 * t9i32 * exp(-5.028*t9i)
      ddd   = dd*(-1.5d0*t9i + 5.028*t9i2)

      ee    = rc121* 1.25d-03 * t9i32 * exp(-7.929*t9i)
      dee   = ee*(-1.5d0*t9i + 7.929*t9i2)

      ff    = rc121 * 2.43d+01 * t9i * exp(-11.523*t9i)
      dff   = ff*(-t9i + 11.523*t9i2)

      gg    = 5.0d0*exp(-15.882*t9i)
      dgg   = gg*15.882*t9i2

      hh    = 1.0d0 + gg
      hhi   = 1.0d0/hh

      term    = (aa + bb + cc + dd + ee + ff) * hhi
      dtermdt = (daa + dbb + dcc + ddd + dee + dff - term*dgg) * hhi


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.27d+10 * t932 * exp(-115.862*t9i)
      drevdt   = rev*(1.5d0*t9i + 115.862*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_mg24ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,term1,dterm1,term2,dterm2, &
                       rev,drevdt,rc148,q1
      parameter        (rc148 = 0.1d0, &
                        q1    = 1.0d0/0.024649d0)



! 24mg(a,p)al27
      aa     = 1.10d+08 * t9i23 * exp(-23.261*t9i13 - t92*q1)
      daa    = -twoth*aa*t9i + aa*(23.261*t9i43 - 2.0d0*t9*q1)

      bb     =  1.0d0 + 0.018*t913 + 12.85*t923 + 1.61*t9 &
               + 89.87*t943 + 28.66*t953
      dbb    = oneth*0.018*t9i23 + twoth*12.85*t9i13 + 1.61 &
                + fourth*89.87*t913 + fiveth*28.66*t923

      term1  = aa * bb
      dterm1 = daa * bb + aa * dbb

      aa     = 129.0d0 * t9i32 * exp(-2.517*t9i)
      daa    = -1.5d0*aa*t9i + aa*2.517*t9i2

      bb     = 5660.0d0 * t972 * exp(-3.421*t9i)
      dbb    = 3.5d0*bb*t9i +  bb*3.421*t9i2

      cc     = rc148 * 3.89d-08 * t9i32 * exp(-0.853*t9i)
      dcc    = -1.5d0*cc*t9i + cc*0.853*t9i2

      dd     = rc148 * 8.18d-09 * t9i32 * exp(-1.001*t9i)
      ddd    = -1.5d0*dd*t9i + dd*1.001*t9i2

      term2  = aa + bb + cc + dd
      dterm2 = daa + dbb + dcc + ddd

      ee     = oneth*exp(-9.792*t9i)
      dee    = ee*9.792*t9i2

      ff     =  twoth * exp(-11.773*t9i)
      dff    = ff*11.773*t9i2

      gg     = 1.0d0 + ee + ff
      dgg    = dee + dff

      term    = (term1 + term2)/gg
      dtermdt = ((dterm1 + dterm2) - term*dgg)/gg


! the rates
      rev      = 1.81 * exp(-18.572*t9i)
      drevdt   = rev*18.572*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end






      subroutine rate_al27pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'

! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,bb,dbb,cc,dcc, &
                       dd,ddd,ee,dee,ff,dff,gg,dgg


! al27(p,g)si28
! champagne 1996

      aa  = 1.32d+09 * t9i23 * exp(-23.26*t9i13)
      daa = aa*(-twoth*t9i + oneth*23.26*t9i43)

      bb  = 3.22d-10 * t9i32 * exp(-0.836*t9i)*0.17
      dbb = bb*(-1.5d0*t9i + 0.836*t9i2)

      cc  = 1.74d+00 * t9i32 * exp(-2.269*t9i)
      dcc = cc*(-1.5d0*t9i + 2.269*t9i2)

      dd  = 9.92d+00 * t9i32 * exp(-2.492*t9i)
      ddd = dd*(-1.5d0*t9i + 2.492*t9i2)

      ee  = 4.29d+01 * t9i32 * exp(-3.273*t9i)
      dee = ee*(-1.5d0*t9i + 3.273*t9i2)

      ff  = 1.34d+02 * t9i32 * exp(-3.654*t9i)
      dff = ff*(-1.5d0*t9i + 3.654*t9i2)

      gg  = 1.77d+04 * (t9**0.53) * exp(-4.588*t9i)
      dgg = gg*(0.53*t9i + 4.588*t9i2)

      term    = aa + bb + cc + dd + ee + ff + gg
      dtermdt = daa + dbb + dcc + ddd + dee + dff + dgg


! rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.13d+11 * t932 * exp(-134.434*t9i)
      drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

      rr    = rev * term
      drrdt = (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end







      subroutine rate_al27pg_old(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,bb,dbb,cc,dcc,dd,ddd,ee,dee, &
                       ff,dff,gg,dgg,hh,dhh,xx,dxx,yy,dyy,zz,dzz,pp, &
                       rev,drevdt,rc147,q1
      parameter        (rc147 = 0.1d0, &
                        q1    = 1.0d0/0.024025d0)



! 27al(p,g)si28  cf88

      aa  = 1.67d+08 * t9i23 * exp(-23.261*t9i13 - t92*q1)
      daa = aa*(-twoth*t9i + oneth*23.261*t9i43 - 2.0d0*t9*q1)

      bb  = 1.0d0 + 0.018*t913 + 5.81*t923 + 0.728*t9 &
            + 27.31*t943 + 8.71*t953
      dbb = oneth*0.018*t9i23 + twoth*5.81*t9i13 + 0.728 &
            + fourth*27.31*t913 + fiveth*8.71*t923

      cc  = aa*bb
      dcc = daa*bb + aa*dbb

      dd  = 2.20d+00 * t9i32 * exp(-2.269*t9i)
      ddd = dd*(-1.5d0*t9i + 2.269*t9i2)

      ee  = 1.22d+01 * t9i32 * exp(-2.491*t9i)
      dee = ee*(-1.5d0*t9i + 2.491*t9i2)

      ff  =  1.50d+04 * t9 * exp(-4.112*t9i)
      dff = ff*(t9i + 4.112*t9i2)

      gg  = rc147 * 6.50d-10 * t9i32 * exp(-0.853*t9i)
      dgg = gg*(-1.5d0*t9i + 0.853*t9i2)

      hh  = rc147 * 1.63d-10 * t9i32 * exp(-1.001*t9i)
      dhh = hh*(-1.5d0*t9i + 1.001*t9i2)

      xx     = oneth*exp(-9.792*t9i)
      dxx    = xx*9.792*t9i2

      yy     =  twoth * exp(-11.773*t9i)
      dyy    = yy*11.773*t9i2

      zz     = 1.0d0 + xx + yy
      dzz    = dxx + dyy

      pp      = 1.0d0/zz
      term    = (cc + dd + ee + ff + gg + hh)*pp
      dtermdt = ((dcc + ddd + dee + dff + dgg + dhh) - term*dzz)*pp


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.13d+11*t932*exp(-134.434*t9i)
      drevdt   = rev*(1.5d0*t9i + 134.434*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_si28ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! si28(a,g)s32
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 6.340d-2*z + 2.541d-3*z2 - 2.900d-4*z3
      if (z .eq. 10.0) then
       daa = 0
      else
       daa   = 6.340d-2 + 2.0d0*2.541d-3*t9 - 3.0d0*2.900d-4*t92
      end if

      term    = 4.82d+22 * t9i23 * exp(-61.015 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 61.015*t9i13*(oneth*t9i*aa - daa))

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.461d+10 * t932 * exp(-80.643*t9i)
      drevdt   = rev*(1.5d0*t9i + 80.643*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_si28ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! si28(a,p)p31
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 2.798d-3*z + 2.763d-3*z2 - 2.341d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 2.798d-3 + 2.0d0*2.763d-3*t9 - 3.0d0*2.341d-4*t92
      end if

      term    = 4.16d+13 * t9i23 * exp(-25.631 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*25.631*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 0.5825d0 * exp(-22.224*t9i)
      drevdt   = rev*22.224*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt * term + rev * dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_p31pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! p31(p,g)s32
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.928d-1*z - 1.540d-2*z2 + 6.444d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.928d-1 - 2.0d0*1.540d-2*t9 + 3.0d0*6.444d-4*t92
      end if

      term    = 1.08d+16 * t9i23 * exp(-27.042 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 27.042*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 3.764d+10 * t932 * exp(-102.865*t9i)
      drevdt   = rev*(1.5d0*t9i + 102.865*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_s32ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! s32(a,g)ar36
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 4.913d-2*z + 4.637d-3*z2 - 4.067d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 4.913d-2 + 2.0d0*4.637d-3*t9 - 3.0d0*4.067d-4*t92
      end if

      term    = 1.16d+24 * t9i23 * exp(-66.690 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 66.690*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.616d+10 * t932 * exp(-77.080*t9i)
      drevdt   = rev*(1.5d0*t9i + 77.080*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_s32ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! s32(a,p)cl35
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.041d-1*z - 1.368d-2*z2 + 6.969d-4*z3
      if (z .eq. 10) then
       daa = 0.0d0
      else
       daa   = 1.041d-1 - 2.0d0*1.368d-2*t9 + 3.0d0*6.969d-4*t92
      end if

      term    = 1.27d+16 * t9i23 * exp(-31.044 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*31.044*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.144 * exp(-21.643*t9i)
      drevdt   = rev*21.643*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_cl35pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt


! cl35(p,g)ar36
      aa    = 1.0d0 + 1.761d-1*t9 - 1.322d-2*t92 + 5.245d-4*t93
      daa   = 1.761d-1 - 2.0d0*1.322d-2*t9 + 3.0d0*5.245d-4*t92


      term    =  4.48d+16 * t9i23 * exp(-29.483 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 29.483*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.568d+10*t932*exp(-98.722*t9i)
      drevdt   = rev*(1.5d0*t9i + 98.722*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ar36ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ar36(a,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.458d-1*z - 1.069d-2*z2 + 3.790d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.458d-1 - 2.0d0*1.069d-2*t9 + 3.0d0*3.790d-4*t92
      end if

      term    = 2.81d+30 * t9i23 * exp(-78.271 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 78.271*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.740d+10 * t932 * exp(-81.711*t9i)
      drevdt   = rev*(1.5d0*t9i + 81.711*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ar36ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ar36(a,p)k39
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 4.826d-3*z - 5.534d-3*z2 + 4.021d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 4.826d-3 - 2.0d0*5.534d-3*t9 + 3.0d0*4.021d-4*t92
      end if

      term    = 2.76d+13 * t9i23 * exp(-34.922 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*34.922*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.128*exp(-14.959*t9i)
      drevdt   = rev*14.959*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_k39pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! k39(p,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.622d-1*z - 1.119d-2*z2 + 3.910d-4*z3
      if (z .eq. 10) then
       daa = 0.0d0
      else
       daa   = 1.622d-1 - 2.0d0*1.119d-2*t9 + 3.0d0*3.910d-4*t92
      end if

      term    = 4.09d+16 * t9i23 * exp(-31.727 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 31.727*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.600d+10 * t932 * exp(-96.657*t9i)
      drevdt   = rev*(1.5d0*t9i + 96.657*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ca40ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ca40(a,g)ti44
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.650d-2*z + 5.973d-3*z2 - 3.889d-04*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.650d-2 + 2.0d0*5.973d-3*z - 3.0d0*3.889d-4*z2
      end if

      term    = 4.66d+24 * t9i23 * exp(-76.435 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 76.435*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.843d+10 * t932 * exp(-59.510*t9i)
      drevdt   = rev*(1.5d0*t9i + 59.510*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ca40ap(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ca40(a,p)sc43
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 - 1.206d-2*z + 7.753d-3*z2 - 5.071d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = -1.206d-2 + 2.0d0*7.753d-3*t9 - 3.0d0*5.071d-4*t92
      end if

      term    = 4.54d+14 * t9i23 * exp(-32.177 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*32.177*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 2.229 * exp(-40.966*t9i)
      drevdt   = rev*40.966*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_sc43pg(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! sc43(p,g)ca40
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.023d-1*z - 2.242d-3*z2 - 5.463d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.023d-1 - 2.0d0*2.242d-3*t9 - 3.0d0*5.463d-5*t92
      end if

      term    = 3.85d+16 * t9i23 * exp(-33.234 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 33.234*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.525d+11 * t932 * exp(-100.475*t9i)
      drevdt   = rev*(1.5d0*t9i + 100.475*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_ti44ag(temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den,fr,dfrdt,dfrdd,rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ti44(a,g)cr48
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.066d-1*z - 1.102d-2*z2 + 5.324d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.066d-1 - 2.0d0*1.102d-2*t9 + 3.0d0*5.324d-4*t92
      end if

      term    = 1.37d+26 * t9i23 * exp(-81.227 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 81.227*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 6.928d+10*t932*exp(-89.289*t9i)
      drevdt   = rev*(1.5d0*t9i + 89.289*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_ti44ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! ti44(a,p)v47
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 2.655d-2*z - 3.947d-3*z2 + 2.522d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 2.655d-2 - 2.0d0*3.947d-3*t9 + 3.0d0*2.522d-4*t92
      end if

      term    = 6.54d+20 * t9i23 * exp(-66.678 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*66.678*t9i13*(oneth*t9i*aa - daa)


! the rates
      rev      = 1.104 * exp(-4.723*t9i)
      drevdt   = rev*4.723*t9i2

      fr    = den * rev * term
      dfrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      dfrdd = rev * term

      rr    = den * term
      drrdt = den * dtermdt * 1.0d-9
      drrdd = term

      return
      end





      subroutine rate_v47pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! v47(p,g)cr48
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.979d-2*z - 2.269d-3*z2 - 6.662d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.979d-2 - 2.0d0*2.269d-3*t9 - 3.0d0*6.662d-5*t92
      end if

      term    = 2.05d+17 * t9i23 * exp(-35.568 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 35.568*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.649d+10*t932*exp(-93.999*t9i)
      drevdt   = rev*(1.5d0*t9i + 93.999*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_cr48ag(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! cr48(a,g)fe52
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 6.325d-2*z - 5.671d-3*z2 + 2.848d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 6.325d-2 - 2.0d0*5.671d-3*t9 + 3.0d0*2.848d-4*t92
      end if

      term    = 1.04d+23 * t9i23 * exp(-81.420 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 81.420*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.001d+10 * t932 * exp(-92.177*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.177*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_cr48ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! cr48(a,p)mn51
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.384d-2*z + 1.081d-3*z2 - 5.933d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.384d-2 + 2.0d0*1.081d-3*t9 - 3.0d0*5.933d-5*t92
      end if

      term    = 1.83d+26 * t9i23 * exp(-86.741 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*86.741*t9i13*(oneth*t9i*aa - daa)


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.6087*exp(-6.510*t9i)
      drevdt   = rev*6.510*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_mn51pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! mn51(p,g)fe52
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 8.922d-2*z - 1.256d-3*z2 - 9.453d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 8.922d-2 - 2.0d0*1.256d-3*t9 - 3.0d0*9.453d-5*t92
      end if

      term    = 3.77d+17 * t9i23 * exp(-37.516 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 37.516*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.150d+11*t932*exp(-85.667*t9i)
      drevdt   = rev*(1.5d0*t9i + 85.667*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end






      subroutine rate_fe52ag(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! fe52(a,g)ni56
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 7.846d-2*z - 7.430d-3*z2 + 3.723d-4*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 7.846d-2 - 2.0d0*7.430d-3*t9 + 3.0d0*3.723d-4*t92
      end if

      term    = 1.05d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 91.674*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 7.064d+10*t932*exp(-92.850*t9i)
      drevdt   = rev*(1.5d0*t9i + 92.850*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_fe52ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! fe52(a,p)co55
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 1.367d-2*z + 7.428d-4*z2 - 3.050d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 1.367d-2 + 2.0d0*7.428d-4*t9 - 3.0d0*3.050d-5*t92
      end if

      term    = 1.30d+27 * t9i23 * exp(-91.674 * t9i13 * aa)
      dtermdt = -twoth*term*t9i + term*91.674*t9i13*(oneth*t9i*aa - daa)


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 0.4597*exp(-9.470*t9i)
      drevdt   = rev*9.470*t9i2

      rr    = den * rev * term
      drrdt = den * (drevdt*term + rev*dtermdt) * 1.0d-9
      drrdd = rev * term

      return
      end





      subroutine rate_co55pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,aa,daa,rev,drevdt,z,z2,z3


! co55(p,g)ni56
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.894d-2*z - 3.131d-3*z2 - 2.160d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.894d-2 - 2.0d0*3.131d-3*t9 - 3.0d0*2.160d-5*t92
      end if

      term    = 1.21d+18 * t9i23 * exp(-39.604 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 39.604*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.537d+11*t932*exp(-83.382*t9i)
      drevdt   = rev*(1.5d0*t9i + 83.382*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end





      subroutine rate_fe52ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,tq2


! fe52(n,g)fe53
      tq2     = t9 - 0.348d0
      term    = 9.604d+05 * exp(-0.0626*tq2)
      dtermdt = -term*0.0626

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.43d+09 * t932 * exp(-123.951*t9i)
      drevdt   = rev*(1.5d0*t9i + 123.951*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end




      subroutine rate_fe53ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,tq1,tq10,dtq10,tq2


! fe53(n,g)fe54
      tq1   = t9/0.348
      tq10  = tq1**0.10
      dtq10 = 0.1d0*tq10/(0.348*tq1)
      tq2   = t9 - 0.348d0

      term    = 1.817d+06 * tq10 * exp(-0.06319*tq2)
      dtermdt = term/tq10*dtq10 - term*0.06319

! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 1.56d+11 * t932 * exp(-155.284*t9i)
      drevdt   = rev*(1.5d0*t9i + 155.284*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_fe54ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe54(n,g)fe55
      aa   =  2.307390d+01 - 7.931795d-02 * t9i + 7.535681d+00 * t9i13 &
            - 1.595025d+01 * t913 + 1.377715d+00 * t9 - 1.291479d-01 * t953 &
            + 6.707473d+00 * log(t9)

      daa  =  7.931795d-02 * t9i2 - oneth * 7.535681d+00 * t9i43 &
            - oneth * 1.595025d+01 *t9i23 + 1.377715d+00 - fiveth * 1.291479d-01 *t923 &
            + 6.707473d+00 * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 4.800293d+09 * t932 * exp(-1.078986d+02 * t9i)
      dbb = bb*(1.5d0*t9i + 1.078986d+02 * t9i2)

! reverse rate
      rr    = term*bb
      drrdt = dtermdt*bb + term*dbb*1.0d-9
      drrdd = 0.0d0

! forward rate
      dfrdd = term
      fr    = term*den
      dfrdt = dtermdt*den

      return
      end



      subroutine rate_fe54pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision term,dtermdt,rev,drevdt,aa,daa,z,z2,z3


! fe54(p,g)co55
      z     = min(t9,10.0d0)
      z2    = z*z
      z3    = z2*z
      aa    = 1.0d0 + 9.593d-2*z - 3.445d-3*z2 + 8.594d-5*z3
      if (z .eq. 10.0) then
       daa = 0.0d0
      else
       daa   = 9.593d-2 - 2.0d0*3.445d-3*t9 + 3.0d0*8.594d-5*t92
      end if

      term    = 4.51d+17 * t9i23 * exp(-38.483 * t9i13 * aa)
      dtermdt = term*(-twoth*t9i + 38.483*t9i13*(oneth*t9i*aa - daa))


! the rates
      fr    = den * term
      dfrdt = den * dtermdt * 1.0d-9
      dfrdd = term

      rev      = 2.400d+09 * t932 * exp(-58.605*t9i)
      drevdt   = rev*(1.5d0*t9i + 58.605*t9i2)

      rr    = rev * term
      drrdt = (drevdt * term + rev * dtermdt) * 1.0d-9
      drrdd = 0.0d0

      return
      end



      subroutine rate_fe54ap(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe54(a,p)co57
      aa   =  3.97474900d+01 - 6.06543100d+00 * t9i + 1.63239600d+02 * t9i13 &
            - 2.20457700d+02 * t913 + 8.63980400d+00 * t9 - 3.45841300d-01 * t953 &
            + 1.31464200d+02 * log(t9)

      daa  =  6.06543100d+00 * t9i2 - oneth * 1.63239600d+02 * t9i43 &
            - oneth * 2.20457700d+02 * t9i23 + 8.63980400d+00 - fiveth * 3.45841300d-01 * t923 &
            + 1.31464200d+02  * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 2.16896000d+00  * exp(-2.05631700d+01 * t9i)
      dbb = bb * 2.05631700d+01 * t9i2

! reverse rate
      drrdd = term
      rr    = term*den
      drrdt = dtermdt*den

! forward rate
      fr    = rr*bb
      dfrdt = drrdt*bb + rr*dbb*1.0d-9
      dfrdd = drrdd*bb

      return
      end







      subroutine rate_fe55ng(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe55(n,g)fe56
      aa   =  1.954115d+01 - 6.834029d-02 * t9i + 5.379859d+00 * t9i13 &
            - 8.758150d+00 * t913 + 5.285107d-01 * t9 - 4.973739d-02  * t953 &
            + 4.065564d+00  * log(t9)

      daa  =  6.834029d-02 * t9i2 - oneth * 5.379859d+00 * t9i43 &
            - oneth * 8.758150d+00 * t9i23 + 5.285107d-01 - fiveth * 4.973739d-02  *t923 &
            + 4.065564d+00  * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 7.684279d+10  * t932 * exp(-1.299472d+02  * t9i)
      dbb = bb*(1.5d0*t9i + 1.299472d+02 * t9i2)

! reverse rate
      rr    = term*bb
      drrdt = dtermdt*bb + term*dbb*1.0d-9
      drrdd = 0.0d0

! forward rate
      dfrdd = term
      fr    = term*den
      dfrdt = dtermdt*den

      return
      end




      subroutine rate_fe56pg(temp,den, &
                            fr,dfrdt,dfrdd, &
                            rr,drrdt,drrdd)
      include 'implno.dek'
      include 'tfactors.dek'


! declare the pass
      double precision temp,den, &
                       fr,dfrdt,dfrdd, &
                       rr,drrdt,drrdd

! locals
      double precision aa,daa,bb,dbb,term,dtermdt


! fe56(p,g)co57
   
      aa   =  1.755960d+02 - 7.018872d+00 * t9i + 2.800131d+02 * t9i13 &
            - 4.749343d+02 * t913 + 2.683860d+01 * t9 - 1.542324d+00  * t953 &
            + 2.315911d+02  * log(t9)

      daa  =  7.018872d+00 * t9i2 - oneth * 2.800131d+02 * t9i43 &
            - oneth * 4.749343d+02 * t9i23 + 2.683860d+01 - fiveth * 1.542324d+00  *t923 &
            + 2.315911d+02  * t9i

      if (aa .lt. 200.0) then
       term    = exp(aa)
       dtermdt = term*daa*1.0d-9
      else
       term    = exp(200.0d0)
       dtermdt = 0.0d0
      end if

      bb  = 2.402486d+09 * t932 * exp(-6.995192d+01 * t9i)
      dbb = bb*(1.5d0*t9i + 6.995192d+01 * t9i2)


! reverse rate
      rr    = term*bb
      drrdt = dtermdt*bb + term*dbb*1.0d-9
      drrdd = 0.0d0

! forward rate
      dfrdd = term
      fr    = term*den
      dfrdt = dtermdt*den

      return
      end
end module rates