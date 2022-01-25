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

      implicit none
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
      real*8 w1, w2, w3, w4, w5, w6, theta, Kvol, I1, I2, Iv, Iw
      real*8 Psi1, Psi2, Psiv, Psiw, Phi1, Phi2, Phi3, Phi4, Phi5, Phi6
      real*8 Psi11, Psi12, Psi1v, Psi1w, Psi22, Psi2v, Psi2w, Psiww,
     *      Psivw, Psivv
      real*8 dPsi1, dPsi2, dPsiv, dPsiw, dPhi1, dPhi2, dPhi3, dPhi4,
     *      dPhi5, dPhi6
      real*8 biases(props(4)), dt, output_vector(props(4)+props(2))
      real*8 output_grad((props(2))*(props(2)+props(4)))

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
      w1 = props(ind+1)
      w2 = props(ind+2)
      w3 = props(ind+3)
      w4 = props(ind+4)
      w5 = props(ind+5)
      w6 = props(ind+6)

      ind = 4+nlayers + weight_count*10 + 3 + 6
      theta = props(ind+1)
      Kvol  = props(ind+2)

      I1 = ainv(1)
      I2 = ainv(2)
      Iv = ainv(4)
      Iw = ainv(8)

      Psi1 = I1-3
      Psi2 = I2-3
      Psiv = Iv-1
      Psiw = Iw-1
      Phi1 = Psi1+Psi2
      Phi2 = Psi1+Psiv
      Phi3 = Psi1+Psiw
      Phi4 = Psi2+Psiv
      Phi5 = Psi2+Psiw
      Phi6 = Psiv+Psiw

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

C     Temporary. Only for neoHookean
C       Psi2 = 0
      Psiv = 0
      Psiw = 0
      Phi1 = 0
      Phi2 = 0
      Phi3 = 0
      Phi4 = 0
      Phi5 = 0
      Phi6 = 0
C       dPsi2 = 0
      dPsiv = 0
      dPsiw = 0
      dPhi1 = 0
      dPhi2 = 0
      dPhi3 = 0
      dPhi4 = 0
      dPhi5 = 0
      dPhi6 = 0

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

C         call NN(Psiv, weightsIv, biases, weight_count, bias_count, 
C      #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
C         Psiv  = Psiv + output_vector(bias_count+1)*dt
C         dPsiv = dPsiv*(1+output_grad((0+bias_count)*n_input+1)*dt)

C         call NN(Psiw, weightsIw, biases, weight_count, bias_count, 
C      #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
C         Psiw  = Psiw + output_vector(bias_count+1)*dt
C         dPsiw = dPsiw*(1+output_grad((0+bias_count)*n_input+1)*dt)

C         call NN(Phi1, weightsJ1, biases, weight_count, bias_count, 
C      #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
C         Phi1  = Phi1 + output_vector(bias_count+1)*dt
C         dPhi1 = dPhi1*(1+output_grad((0+bias_count)*n_input+1)*dt)

C         call NN(Phi2, weightsJ2, biases, weight_count, bias_count, 
C      #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
C         Phi2  = Phi2 + output_vector(bias_count+1)*dt
C         dPhi2 = dPhi2*(1+output_grad((0+bias_count)*n_input+1)*dt)

C         call NN(Phi3, weightsJ3, biases, weight_count, bias_count, 
C      #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
C         Phi3  = Phi3 + output_vector(bias_count+1)*dt
C         dPhi3 = dPhi3*(1+output_grad((0+bias_count)*n_input+1)*dt)

C         call NN(Phi4, weightsJ4, biases, weight_count, bias_count, 
C      #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
C         Phi4  = Phi4 + output_vector(bias_count+1)*dt
C         dPhi4 = dPhi4*(1+output_grad((0+bias_count)*n_input+1)*dt)

C         call NN(Phi5, weightsJ5, biases, weight_count, bias_count, 
C      #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
C         Phi5  = Phi5 + output_vector(bias_count+1)*dt
C         dPhi5 = dPhi5*(1+output_grad((0+bias_count)*n_input+1)*dt)

C         call NN(Phi6, weightsJ6, biases, weight_count, bias_count, 
C      #          output_vector, output_grad, nlayers, n_neuronsperlayer, activtypes)
C         Phi6  = Phi6 + output_vector(bias_count+1)*dt
C         dPhi6 = dPhi6*(1+output_grad((0+bias_count)*n_input+1)*dt)
      end do

C       if (Phi1.lt.0) then
C         Phi1 = 0
C         dPhi1 = 0
C       end if 
C       if (Phi2.lt.0) then
C         Phi2 = 0
C         dPhi2 = 0
C       end if 
C       if (Phi3.lt.0) then
C         Phi3 = 0
C         dPhi3 = 0
C       end if 
C       if (Phi4.lt.0) then
C         Phi4 = 0
C         dPhi4 = 0
C       end if 
C       if (Phi5.lt.0) then
C         Phi5 = 0
C         dPhi5 = 0
C       end if 
C       if (Phi6.lt.0) then
C         Phi6 = 0
C         dPhi6 = 0
C       end if 
C       if (Psiv.lt.0) then
C         Psiv = 0
C         dPsiv = 0
C       end if 
C       if (Psiw.lt.0) then
C         Psiw = 0
C         dPsiw = 0
C       end if 

C       Psi1 = Psi1 + w1*Phi1 + w2*Phi2 + w3*Phi3
C       Psi2 = Psi2 + w1*Phi1 + w4*Phi4 + w5*Phi5
C       Psiv = Psiv + w2*Phi2 + w4*Phi4 + w6*Phi6
C       Psiw = Psiw + w3*Phi3 + w5*Phi5 + w6*Phi6

      print*, 'Good till 234'
      Psi11 = dPsi1 + w1*dPhi1 + w2*dPhi2 + w3*dPhi3
      Psi12 = w1*dPhi1
      Psi1v = w2*dPhi2
      Psi1w = w3*dPhi3
      Psi22 = dPsi2 + w1*dPhi1 + w4*dPhi4 + w5*dPhi5
C      Psi21 = w1*dPhi1
      Psi2v = w4*dPhi4
      Psi2w = w5*dPhi5
      Psivv = dPsiv + w2*dPhi2 + w4*dPhi4 + w6*dPhi6
C      Psiv1 = w2*dPhi2
C      Psiv2 = w4*dPhi4
      Psivw = w6*dPhi6
      Psiww = dPsiw + w3*dPhi3 + w5*dPhi5 + w6*dPhi6
C      Psiw1 = w3*dPhi3
C      Psiw2 = w5*dPhi5
C      Psiwv = w6*dPhi6

      ui1(1) = Psi1
      ui1(2) = Psi2
      ui1(4) = Psiv
      ui1(8) = Psiw

      ui2(indx(1,1)) = Psi11
      ui2(indx(1,2)) = Psi12
      ui2(indx(1,4)) = Psi1v
      ui2(indx(1,8)) = Psi1w
      ui2(indx(2,4)) = Psi2v
      ui2(indx(2,8)) = Psi2w
      ui2(indx(4,8)) = Psivw

C     Debugging
      if (noel.eq.500) then
        print*, '---------------------------------------------'
        print*, 'I1, I2, I4v, I4w'
        print*, I1, I2, Iv, Iw
        print*, 'Psi1, Psi2, Psiv, Psiw'
        print*, Psi1, Psi2, Psiv, Psiw
        print*, 'Psi11, Psi12, Psi1v, Psi1w, Psi22, Psi2v, Psi2w'
        print*, Psi11, Psi12, Psi1v, Psi1w, Psi22, Psi2v, Psi2w
        print*, 'Psivv, Psivw, Psiww'
        print*, Psivv, Psivw, Psiww
      end if

      return
      contains
C-------------------------------------------------------------
C     Function to map index from Square to Triangular storage 
C            of symmetric matrix
C
      integer function indx( i, j )
      integer i, j, ii, jj
C      include 'aba_param.inc'
      ii = min(i,j)
      jj = max(i,j)
      indx = ii + jj*(jj-1)/2
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