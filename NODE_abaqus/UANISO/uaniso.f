C-------------------------------------------------------------
      subroutine uanisohyper_inv (ainv, ua, zeta, nfibers, ninv,
     $     ui1, ui2, ui3, temp, noel, cmname, incmpflag, ihybflag,
     $     numstatev, statev, numfieldv, fieldv, fieldvinc,
     $     numprops, props)
C
      include 'aba_param.inc'
C
      character*80 cmname
      dimension ua(2), ainv(ninv), ui1(ninv),
     $     ui2(ninv*(ninv+1)/2), ui3(ninv*(ninv+1)/2),
     $     statev(numstatev), fieldv(numfieldv),
     $     fieldvinc(numfieldv), props(numprops)
C
C
C
      call UANISOHYPER_INVHGO(ainv, ua, zeta, nfibers, ninv,
     $     ui1, ui2, ui3, temp, noel, cmname, incmpflag, ihybflag,
     $     numstatev, statev, numfieldv, fieldv, fieldvinc,
     $     numprops, props)

C
C
C
      return
      end
c------------------------------------------------------------------
c
c     HGO model
c
      subroutine uanisohyper_invhgo (ainv, ua, zeta, nfibers, ninv,
     $     ui1, ui2, ui3, temp, noel, cmname, incmpflag, ihybflag,
     $     numstatev, statev, numfieldv, fieldv, fieldvinc,
     $     numprops, props)
C
      include 'aba_param.inc'
C
      character*80 cmname
      dimension ua(2), ainv(ninv), ui1(ninv),
     $     ui2(ninv*(ninv+1)/2), ui3(ninv*(ninv+1)/2),
     $     statev(numstatev), fieldv(numfieldv),
     $     fieldvinc(numfieldv), props(numprops)
C
c     ainv: invariants
c     ua  : energies ua(1): utot, ua(2); udev
c     ui1 : dUdI
c     ui2 : d2U/dIdJ
c     ui3 : d3U/dIdJdJ, not used for regular elements
C
      parameter ( half = 0.5d0,
     *            zero = 0.d0, 
     *            one  = 1.d0, 
     *            two  = 2.d0, 
     *            three= 3.d0, 
     *            four = 4.d0, 
     *            five = 5.d0, 
     *            six  = 6.d0,
c
     *            index_I1 = 1,
     *            index_J  = 3,
     *            asmall   = 2.d-16  )
C
C     HGO model
C
C       C10 = props(1)
C       rk1 = props(3)
C       rk2 = props(4)
C       rkp = props(5)
      C10 = 0.012902496702913772
      rk1 = 0.01724173170395558
      rk2 = 14.00442692847235
      rkp = 0.0
c
      ua(2) = zero
      om3kp = one - three * rkp
      do k1 = 1, nfibers
**         index_i4 = 4 + k1*(k1-1) + 2*(k1-1)
         index_i4 = indxInv4(k1,k1)
         E_alpha1 = rkp  * (ainv(index_i1) - three) 
     *           + om3kp * (ainv(index_i4) - one  )
         E_alpha = max(E_alpha1, zero)
         ht4a    = half + sign(half,E_alpha1 + asmall)
         aux     = exp(rk2*E_alpha*E_alpha)
c energy
         ua(2) = ua(2) +  aux - one
c ui1
         ui1(index_i1) = ui1(index_i1) + aux * E_alpha
         ui1(index_i4) = rk1 * om3kp * aux * E_alpha
c ui2
         aux2 = ht4a + two * rk2 * E_alpha * E_alpha
         ui2(indx(index_I1,index_I1)) = ui2(indx(index_I1,index_I1))
     *                                + aux * aux2
         ui2(indx(index_I1,index_i4)) = rk1*rkp*om3kp * aux * aux2
         ui2(indx(index_i4,index_i4)) = rk1*om3kp*om3kp*aux * aux2
      end do
c
c     deviatoric energy
c
      ua(2) = ua(2) * rk1 / (two * rk2)
      ua(2) = ua(2) + C10 * (ainv(index_i1) - three)
c
c     compute derivatives
c
      ui1(index_i1) = rk1 * rkp * ui1(index_i1) + C10
      ui2(indx(index_I1,index_I1))= ui2(indx(index_I1,index_I1)) 
     *                            * rk1 * rkp * rkp
c     
c     compressible case
      if(props(2).gt.zero) then
         Dinv = one / props(2)
         det = ainv(index_J)
         ua(1) = ua(2) + Dinv *((det*det - one)/two - log(det))
         ui1(index_J) = Dinv * (det - one/det)
         ui2(indx(index_J,index_J))= Dinv * (one + one / det / det)
         if (hybflag.eq.1) then
           ui3(indx(index_J,index_J))= - Dinv * two / (det*det*det)
         end if
      end if
c
      return
      end
C-------------------------------------------------------------
C     Function to map index from Square to Triangular storage 
C 		 of symmetric matrix
C
      integer function indx( i, j )
      include 'aba_param.inc'
      ii = min(i,j)
      jj = max(i,j)
      indx = ii + jj*(jj-1)/2
      return
      end
C-------------------------------------------------------------
C
C     Function to generate enumeration of scalar
C     Pseudo-Invariants of type 4

      integer function indxInv4( i, j )
      include 'aba_param.inc'
      ii = min(i,j)
      jj = max(i,j)
      indxInv4 = 4 + jj*(jj-1) + 2*(ii-1)
      return
      end
C-------------------------------------------------------------
C
C     Function to generate enumeration of scalar
C     Pseudo-Invariants of type 5
C
      integer function indxInv5( i, j )
      include 'aba_param.inc'
      ii = min(i,j)
      jj = max(i,j)
      indxInv5 = 5 + jj*(jj-1) + 2*(ii-1)
      return
      end
C----------------------------------------------------------------------
      subroutine aprd(A,B,C,n,m,k)
c
      include 'aba_param.inc'
c
      parameter (zero = 0.d0)
      dimension A(n,m),B(m,k),C(n,k)
c
      do k1=1,n
         do k2=1,k
            C(k1,k2) = zero
            do k3=1,m
               C(k1,k2)=C(k1,k2)+A(k1,k3)*B(k3,k2)
            end do
         end do
      end do
c
      return
      end
C----------------------------------------------------------------------
      subroutine aTprd(A,B,C,n,m,k)
c
      include 'aba_param.inc'
c
      parameter (zero = 0.d0)
      dimension A(n,m),B(m,k),C(n,k)
c
      do k1=1,n
         do k2=1,k
            C(k1,k2) = zero
            do k3=1,m
               C(k1,k2)=C(k1,k2)+A(k1,k3)*B(k2,k3)
            end do
         end do
      end do
c
      return
      end
