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
      call NODE(ninv, ainv, numprops, props, ua, ui1, ui2, noel)
C
C
C
      return
      end

C------------------------------------------------------------------
C     NODE based material model
      subroutine NODE(ninv, ainv, nprops, props, ua, ui1, ui2, noel)

C      implicit none
      integer ninv, nprops, noel
      real*8 ua(2), ainv(ninv), ui1(ninv), ui2(ninv*(ninv+1)/2),
     $      props(nprops)

      integer nlayers, n_input, weight_count, bias_count, ind
      integer n_neuronsperlayer(props(1)), i, j, k, l
      real*8 weightsI1(props(3)), weightsI2(props(3))
      real*8 weightsIv(props(3)), weightsIw(props(3)) 
      real*8 weightsJ1(props(3)), weightsJ2(props(3))
      real*8 weightsJ3(props(3)), weightsJ4(props(3))  
      real*8 weightsJ5(props(3)), weightsJ6(props(3))
      integer activtypes(props(1)-1)
      real*8 w1, w2, w3, w4, w5, w6, theta, Psi1_bias, Psi2_bias, Kvol, 
     *       I1, I2, Iv, Iw
      real*8 Psi1, Psi2, Psiv, Psiw, Phi1, Phi2, Phi3, Phi4, Phi5, Phi6
      real*8 Psi11, Psi12, Psi1v, Psi1w, Psi22, Psi2v, Psi2w, Psiww,
     *      Psivw, Psivv
      real*8 dPsi1, dPsi2, dPsiv, dPsiw, dPhi1, dPhi2, dPhi3, dPhi4,
     *      dPhi5, dPhi6
      real*8 biases(props(4)), dt, output_vector(props(4)+props(2))
      real*8 output_grad((props(2))*(props(2)+props(4)))

C     Only needed for the closed form HGO
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


      nlayers   = props(1) 
      n_input = props(2)
      weight_count = props(3)
      bias_count = props(4)

      do i = 1,nlayers 
        n_neuronsperlayer(i)    = props(4+i)
      end do

      ind = 1
      do i=1,nlayers-1
        do j=1,n_neuronsperlayer(i)
          do k=1,n_neuronsperlayer(i+1)
            weightsI1(ind) = props(4+nlayers+ind)
            weightsI2(ind) = props(4+nlayers+ind+weight_count*1)
            weightsIv(ind) = props(4+nlayers+ind+weight_count*2)
            weightsIw(ind) = props(4+nlayers+ind+weight_count*3)
            weightsJ1(ind) = props(4+nlayers+ind+weight_count*4)
            weightsJ2(ind) = props(4+nlayers+ind+weight_count*5)
            weightsJ3(ind) = props(4+nlayers+ind+weight_count*6)
            weightsJ4(ind) = props(4+nlayers+ind+weight_count*7)
            weightsJ5(ind) = props(4+nlayers+ind+weight_count*8)
            weightsJ6(ind) = props(4+nlayers+ind+weight_count*9)
            ind = ind + 1
          end do
        end do
      end do

      do i=1,nlayers-1
        activtypes(i)=int(props(4+nlayers + weight_count*10+i))
      end do

      ind = 4+nlayers + weight_count*10 + 3
      w1 = sigmoid(props(ind+1))
      w2 = sigmoid(props(ind+2))
      w3 = sigmoid(props(ind+3))
      w4 = sigmoid(props(ind+4))
      w5 = sigmoid(props(ind+5))
      w6 = sigmoid(props(ind+6))

      ind = 4+nlayers + weight_count*10 + 3 + 6
      theta = props(ind+1)
      Psi1_bias = props(ind+2)
      Psi2_bias = props(ind+3)
      Kvol  = props(ind+4)

      I1 = ainv(1)
      I2 = ainv(2)
      Iv = ainv(4)
      Iw = ainv(8)

      Psi1 = I1-3
      Psi2 = I2-3
      Psiv = Iv-1
      Psiw = Iw-1
      Phi1 = w1*Psi1+(1-w1)*Psi2
      Phi2 = w2*Psi1+(1-w2)*Psiv
      Phi3 = w3*Psi1+(1-w3)*Psiw
      Phi4 = w4*Psi2+(1-w4)*Psiv
      Phi5 = w5*Psi2+(1-w5)*Psiw
      Phi6 = w6*Psiv+(1-w6)*Psiw
      
C       Psi2 = 0.
C       Psiw = 0.
C       Phi1 = 0.
C       Phi3 = 0.
C       Phi4 = 0.
C       Phi5 = 0.
C       Phi6 = 0.

      dPsi1 = 1
      dPsi2 = 1
      dPsiv = 1
      dPsiw = 1
      dPhi1 = 1
      dPhi2 = 1
      dPhi3 = 1
      dPhi4 = 1
      dPhi5 = 1
      dPhi6 = 1


C     Call the NNs in a loop and perform forward Euler integration
      biases = 0
      do i=1, 10
        dt = 1.
        call NN(Psi1, weightsI1, biases, weight_count, bias_count, 
     #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
        Psi1  = Psi1 + output_vector(bias_count+1)*dt
        dPsi1 = dPsi1*(1+output_grad((0+bias_count)*n_input+1)*dt)

        call NN(Psi2, weightsI2, biases, weight_count, bias_count, 
     #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
        Psi2  = Psi2 + output_vector(bias_count+1)*dt
        dPsi2 = dPsi2*(1+output_grad((0+bias_count)*n_input+1)*dt)

        call NN(Psiv, weightsIv, biases, weight_count, bias_count, 
     #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
        Psiv  = Psiv + output_vector(bias_count+1)*dt
        dPsiv = dPsiv*(1+output_grad((0+bias_count)*n_input+1)*dt)

        call NN(Psiw, weightsIw, biases, weight_count, bias_count, 
     #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
        Psiw  = Psiw + output_vector(bias_count+1)*dt
        dPsiw = dPsiw*(1+output_grad((0+bias_count)*n_input+1)*dt)

        call NN(Phi1, weightsJ1, biases, weight_count, bias_count, 
     #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
        Phi1  = Phi1 + output_vector(bias_count+1)*dt
        dPhi1 = dPhi1*(1+output_grad((0+bias_count)*n_input+1)*dt)

        call NN(Phi2, weightsJ2, biases, weight_count, bias_count, 
     #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
        Phi2  = Phi2 + output_vector(bias_count+1)*dt
        dPhi2 = dPhi2*(1+output_grad((0+bias_count)*n_input+1)*dt)

        call NN(Phi3, weightsJ3, biases, weight_count, bias_count, 
     #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
        Phi3  = Phi3 + output_vector(bias_count+1)*dt
        dPhi3 = dPhi3*(1+output_grad((0+bias_count)*n_input+1)*dt)

        call NN(Phi4, weightsJ4, biases, weight_count, bias_count, 
     #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
        Phi4  = Phi4 + output_vector(bias_count+1)*dt
        dPhi4 = dPhi4*(1+output_grad((0+bias_count)*n_input+1)*dt)

        call NN(Phi5, weightsJ5, biases, weight_count, bias_count, 
     #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
        Phi5  = Phi5 + output_vector(bias_count+1)*dt
        dPhi5 = dPhi5*(1+output_grad((0+bias_count)*n_input+1)*dt)

        call NN(Phi6, weightsJ6, biases, weight_count, bias_count, 
     #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
        Phi6  = Phi6 + output_vector(bias_count+1)*dt
        dPhi6 = dPhi6*(1+output_grad((0+bias_count)*n_input+1)*dt)
      end do

      if (Phi1.lt.0) then
        Phi1 = 0
        dPhi1 = 0
      end if 
      if (Phi2.lt.0) then
        Phi2 = 0
        dPhi2 = 0
      end if 
      if (Phi3.lt.0) then
        Phi3 = 0
        dPhi3 = 0
      end if 
      if (Phi4.lt.0) then
        Phi4 = 0
        dPhi4 = 0
      end if 
      if (Phi5.lt.0) then
        Phi5 = 0
        dPhi5 = 0
      end if 
      if (Phi6.lt.0) then
        Phi6 = 0
        dPhi6 = 0
      end if 
      if (Psiv.lt.0) then
        Psiv = 0
        dPsiv = 0
      end if 
      if (Psiw.lt.0) then
        Psiw = 0
        dPsiw = 0
      end if 

      Psi1 = Psi1+    w1*Phi1+    w2*Phi2+    w3*Phi3+exp(Psi1_bias)
      Psi2 = Psi2+(1-w1)*Phi1+    w4*Phi4+    w5*Phi5+exp(Psi2_bias)
      Psiv = Psiv+(1-w2)*Phi2+(1-w4)*Phi4+    w6*Phi6
      Psiw = Psiw+(1-w3)*Phi3+(1-w5)*Phi5+(1-w6)*Phi6
      
C       Psi2 = 0.
      
      Psi11 = dPsi1 +    w1**2*dPhi1 +     w2**2*dPhi2 +     w3**2*dPhi3
      Psi22 = dPsi2 +(1-w1)**2*dPhi1 +     w4**2*dPhi4 +     w5**2*dPhi5
      Psivv = dPsiv +(1-w2)**2*dPhi2 + (1-w4)**2*dPhi4 +     w6**2*dPhi6
      Psiww = dPsiw +(1-w3)**2*dPhi3 + (1-w5)**2*dPhi5 + (1-w6)**2*dPhi6
      Psi12 = w1*(1-w1)*dPhi1
      Psi1v = w2*(1-w2)*dPhi2
      Psi1w = w3*(1-w3)*dPhi3
      Psi2v = w4*(1-w4)*dPhi4
      Psi2w = w5*(1-w5)*dPhi5
      Psivw = w6*(1-w6)*dPhi6

      ui1(1) = Psi1
      ui1(2) = Psi2
      ui1(4) = Psiv
      ui1(8) = Psiw

      ui2(indx(1,1)) = Psi11
      ui2(indx(1,2)) = Psi12
      ui2(indx(1,4)) = Psi1v
      ui2(indx(1,8)) = Psi1w
      ui2(indx(2,2)) = Psi22
      ui2(indx(2,4)) = Psi2v
      ui2(indx(2,8)) = Psi2w
      ui2(indx(4,4)) = Psivv
      ui2(indx(4,8)) = Psivw
      ui2(indx(8,8)) = Psiww

C ---------------------------------------------------------------
C     Analytical HGO model for comparison
      C10 = 1.02356332e-02/2
      rk1 = 5.13664702e-01
      rk2 = 5.91491834e+01
      rkp = 2.74447648e-01
c
      ua(2) = zero
      ui1 = zero
      ui2 = zero
      om3kp = one - three * rkp
      nfibers = 1
      do k1 = 1, nfibers
         index_i4 = indxInv4(k1,k1)
         E_alpha1 = rkp  * (ainv(index_i1) - three) 
     *           + om3kp * (ainv(index_i4) - one  )

C The correct killswitch:
C          if (ainv(index_i4).lt.one) then
C             E_alpha = zero
C             ui1(index_i1) = zero
C             ui1(index_i4) = zero
C             ui2(indx(index_I1,index_I1)) = zero
C             ui2(indx(index_I1,index_i4)) = zero
C             ui2(indx(index_i4,index_i4)) = zero
C          else
C             E_alpha = E_alpha1
C             aux     = exp(rk2*E_alpha*E_alpha)
C             ui1(index_i1) = ui1(index_i1) + aux * E_alpha
C             ui1(index_i4) = rk1 * om3kp * aux * E_alpha
C             aux2 = ht4a + two * rk2 * E_alpha * E_alpha
C             ui2(indx(index_I1,index_I1)) = ui2(indx(index_I1,index_I1))
C      *                                + aux * aux2
C             ui2(indx(index_I1,index_i4)) = rk1*rkp*om3kp * aux * aux2
C          end if

C C The incorrect killswitch:
         E_alpha = max(E_alpha1, zero)
         ht4a    = half + sign(half,E_alpha1 + asmall)
         aux     = exp(rk2*E_alpha*E_alpha)
         ui1(index_i1) = ui1(index_i1) + aux * E_alpha
         ui1(index_i4) = rk1 * om3kp * aux * E_alpha
         aux2 = ht4a + two * rk2 * E_alpha * E_alpha
            ui2(indx(index_I1,index_I1)) = ui2(indx(index_I1,index_I1))
     *                                + aux * aux2
         ui2(indx(index_I1,index_i4)) = rk1*rkp*om3kp * aux * aux2
         ui2(indx(index_i4,index_i4)) = rk1*om3kp*om3kp*aux * aux2
      end do

      ui1(index_i1) = rk1 * rkp * ui1(index_i1) + C10
      ui2(indx(index_I1,index_I1))= ui2(indx(index_I1,index_I1)) 
     *                            * rk1 * rkp * rkp
C ---------------------------------------------------------------

C Neo Hookean material
C       ui1 = 0
C       ui2 = 0
C       ui1(1) = 0.05

C     Debugging
      if (noel.eq.500) then
        print*, '---------------------------------------------'
        print*, 'I1, I2, I4v, I4w'
        print*, I1, I2, Iv, Iw
        print*, 'First derivatives with NODE'
        print*, Psi1, Psi2, Psiv, Psiw
        print*, 'First derivatives with GOH'
        print*, ui1(1), ui1(2), ui1(4), ui1(8)

        print*, 'Second derivatives with NODE'
        print*, 'Psi11, 12, 1v, 1w, 22, 2v, 2w, vv, vw, ww'
        print*, Psi11, Psi12, Psi1v, Psi1w, Psi22, Psi2v, Psi2w, Psivv,
     *          Psivw, Psiww
        print*, 'Second derivatives with GOH'
        print*, ui2(indx(1,1)), ui2(indx(1,2)), ui2(indx(1,4)),
     *          ui2(indx(1,8)), ui2(indx(2,2)), ui2(indx(2,4)),
     *          ui2(indx(2,8)), ui2(indx(4,4)), ui2(indx(4,8)),
     *          ui2(indx(8,8))
      end if

      return
      contains
C-------------------------------------------------------------
C     Function to map index from Square to Triangular storage 
C            of symmetric matrix
C
      integer function indx( i, j )
      implicit none
      integer i, j, ii, jj
C      include 'aba_param.inc'
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
      real function sigmoid(x)
      implicit none
      real*8 x
      sigmoid = 1./(1.+exp(-x))
      return
      end
C-------------------------------------------------------------
C     Subroutine that evaluates the neural network in every step
C
      subroutine NN(input, weights, biases, weight_count, 
     #              bias_count, output_vector, output_grad, nlayers,
     #              n_neuronsperlayer, activtypes)
        implicit none
        integer n_input
        parameter (n_input = 1)
        integer weight_count, bias_count
        real*8 input
        real*8 weights(weight_count), biases(bias_count)
        real*8 output_vector(n_input+bias_count)
        real*8 output_grad(n_input*(n_input+bias_count))
        real*8 activout, gradactivout
        integer activtypes(nlayers-1), flag
        integer i, j, k, l, io1, iw1, ig1, ib1, io2, iw2, ig2, ib2
        integer nlayers, n_neuronsperlayer(nlayers)

        output_vector = 0
        output_grad = 0

        do i=1,n_input
          output_vector(i) = input
          output_grad((i-1)*n_input+i) = 1
        end do

        io1 = 0
        iw1 = 0
        ig1 = 0
        ib1 = 0
        do i=1,nlayers-1
c...    Beginning and end of the chunk in output vector to be used as input
          io2 = io1 + n_neuronsperlayer(i) 
c...    Beginning and end of the chunk in weight array defining matrix 
          iw2 = iw1 + n_neuronsperlayer(i)*n_neuronsperlayer(i+1) 
c...    Beginning and end of the chunk for grad outputs 
          ig2 = ig1 + n_neuronsperlayer(i)*n_input
c...    do the matrix vector product and store in output chunk 
          do k = 1,n_neuronsperlayer(i+1)
            do j = 1,n_neuronsperlayer(i) 
c...        Matrix*vector + bias 
              output_vector(io2+k) = output_vector(io2+k) + weights(iw1+(j-1)* 
     #                n_neuronsperlayer(i+1)+k)*output_vector(io1+j)
c...        Matrix*Matrix for jacobian 
              do l = 1,n_input
                output_grad(ig2+(k-1)*n_input+l) = output_grad(ig2+(k-1)*n_input+l ) +
     #                weights(iw1+(j-1)*n_neuronsperlayer(i+1)+k)
     #                *output_grad(ig1+(j-1)*n_input+l)
              end do
            end do
            activout = output_vector(io2+k) + biases(ib1+k)
            gradactivout = activout
            call activation(activout,activtypes(i))
            output_vector(io2+k) = activout
            call grad_activation(gradactivout,activtypes(i))
            do l = 1,n_input
              output_grad(ig2+(k-1)*n_input+l) = gradactivout*output_grad(ig2+(k-1)*n_input+l)
            end do
          end do
          io1 = io2 
          iw1 = iw2
          ig1 = ig2
          ib1 = ib1 + n_neuronsperlayer(i+1)
        end do

        io2 = n_neuronsperlayer(i)
      end


      subroutine activation(value, typea)
      implicit none

      real*8 value
      integer typea
c...          | 0       =>    linear
c...          | 1       =>    ReLU
c...  typea = | 2       =>    sigmoid
c...          | 3       =>    tanh
c...          | 4       =>    not used yet
      if (typea==0) then
        value = value
      else if (typea==1) then
        if (value<0) then
          value = 0
        end if
      else if (typea==2) then
        value = 1/(1+exp(-value))
      else if (typea==3) then
        value = tanh(value)
      end if

      return
      end
c...  ------------------------------------------------------------------

      subroutine grad_activation(value, typea)
      implicit none

      real*8 value
      real*8 sigmoid
      integer typea
c...          | 0       =>    linear
c...          | 1       =>    ReLU
c...  typea = | 2       =>    sigmoid
c...          | 3       =>    tanh
c...          | 4       =>    not used yet
      if (typea==0) then
        value = 1
      else if (typea==1) then
        if (value<1e-9) then
          value = 0
        else
          value = 1
        end if
      else if (typea==2) then
        sigmoid =1/(1+exp(-value))
        value = sigmoid*(1-sigmoid)
      else if (typea==3) then
        value = 1 - tanh(value)**2
      end if

      return
      end


      end