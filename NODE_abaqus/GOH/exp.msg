1

   Abaqus 2020                                  Date 15-Dec-2021   Time 20:28:35
   For use by PURDUE UNIVERSITY under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     S T A T I C   A N A L Y S I S


                                                                                          

     AUTOMATIC TIME CONTROL WITH -
          A SUGGESTED INITIAL TIME INCREMENT OF                 1.00    
          AND A TOTAL TIME PERIOD OF                            10.0    
          THE MINIMUM TIME INCREMENT ALLOWED IS                1.000E-04
          THE MAXIMUM TIME INCREMENT ALLOWED IS                 1.00    

     LINEAR EQUATION SOLVER TYPE         DIRECT SPARSE

 CONVERGENCE TOLERANCE PARAMETERS FOR FORCE    
     CRITERION FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR DISP.    CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INITIAL VALUE OF TIME AVERAGE FORCE                               1.000E-02
     AVERAGE FORCE     IS TIME AVERAGE FORCE    
     ALTERNATE CRIT. FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO FORCE     RELATIVE TO TIME AVRG. FORCE         1.000E-05
     CRITERION FOR RESIDUAL FORCE     WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR DISP.    CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL FORCE     FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO FORCE     REL. TO TIME AVRG. MAX. FORCE        1.000E-05
     CRITERION FOR ZERO DISP.    RELATIVE TO CHARACTERISTIC LENGTH     1.000E-08

 CONVERGENCE TOLERANCE PARAMETERS FOR MOMENT   
     CRITERION FOR RESIDUAL MOMENT    FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR ROTATION CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INITIAL VALUE OF TIME AVERAGE MOMENT                              1.000E-02
     AVERAGE MOMENT    IS TIME AVERAGE MOMENT   
     ALTERNATE CRIT. FOR RESIDUAL MOMENT    FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO MOMENT    RELATIVE TO TIME AVRG. MOMENT        1.000E-05
     CRITERION FOR RESIDUAL MOMENT    WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR ROTATION CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL MOMENT    FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO MOMENT    REL. TO TIME AVRG. MAX. MOMENT       1.000E-05

 CONVERGENCE TOLERANCE PARAMETERS FOR CAV. VOL.
     CRITERION FOR RESIDUAL CAV. VOL. FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR H.PRESS. CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INITIAL VALUE OF TIME AVERAGE CAV. VOL.                           1.000E-02
     AVERAGE CAV. VOL. IS TIME AVERAGE CAV. VOL.
     ALTERNATE CRIT. FOR RESIDUAL CAV. VOL. FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO CAV. VOL. RELATIVE TO TIME AVRG. CAV. VOL.     1.000E-05
     CRITERION FOR RESIDUAL CAV. VOL. WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR H.PRESS. CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL CAV. VOL. FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO CAV. VOL. REL. TO TIME AVRG. MAX. CAV. VOL.    1.000E-05

     VOLUMETRIC STRAIN COMPATIBILITY TOLERANCE FOR HYBRID SOLIDS       1.000E-05
     AXIAL STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS             1.000E-05
     TRANS. SHEAR STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS      1.000E-05
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P>P0          5.000E-03
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P=0.0         0.100    
     CONTACT FORCE ERROR TOLERANCE FOR CONVERT SDI=YES                 1.00    
     DISPLACEMENT COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS           1.000E-05
     ROTATION COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS               1.000E-05

 EQUILIBRIUM WILL BE CHECKED FOR SEVERE DISCONTINUITY ITERATIONS

 TIME INCREMENTATION CONTROL PARAMETERS:
     FIRST EQUILIBRIUM ITERATION FOR CONSECUTIVE DIVERGENCE CHECK              4
     EQUILIBRIUM ITERATION AT WHICH LOG. CONVERGENCE RATE CHECK BEGINS         8
     EQUILIBRIUM ITERATION AFTER WHICH ALTERNATE RESIDUAL IS USED              9
     MAXIMUM EQUILIBRIUM ITERATIONS ALLOWED                                   16
     EQUILIBRIUM ITERATION COUNT FOR CUT-BACK IN NEXT INCREMENT               10
     MAXIMUM EQUILIB. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      4
     MAXIMUM ITERATIONS FOR SEVERE DISCONTINUITIES                            50
     MAXIMUM ATTEMPTS ALLOWED IN AN INCREMENT                                  5
     MAXIMUM DISCON. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      50
     MAXIMUM CONTACT AUGMENTATIONS FOR *SURFACE BEHAVIOR,AUGMENTED LAGRANGE   50
     CUT-BACK FACTOR AFTER DIVERGENCE                                 0.2500    
     CUT-BACK FACTOR FOR TOO SLOW CONVERGENCE                         0.5000    
     CUT-BACK FACTOR AFTER TOO MANY EQUILIBRIUM ITERATIONS            0.7500    
     CUT-BACK FACTOR AFTER TOO MANY SEVERE DISCONTINUITY ITERATIONS   0.2500    
     CUT-BACK FACTOR AFTER PROBLEMS IN ELEMENT ASSEMBLY               0.2500    
     INCREASE FACTOR AFTER TWO INCREMENTS THAT CONVERGE QUICKLY        1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED                       1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DYNAMICS)            1.250    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DIFFUSION)           2.000    
     MINIMUM TIME INCREMENT RATIO FOR EXTRAPOLATION TO OCCUR          0.1000    
     MAX. RATIO OF TIME INCREMENT TO STABILITY LIMIT                   1.000    
     FRACTION OF STABILITY LIMIT FOR NEW TIME INCREMENT               0.9500    
     TIME INCREMENT INCREASE FACTOR BEFORE A TIME POINT                1.000    
     GLOBAL STABILIZATION CONTROL IS NOT USED

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                     100

          LARGE DISPLACEMENT THEORY WILL BE USED

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH     0.237    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.24000    
       ELAPSED SYSTEM TIME (SEC)    =  6.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.30000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 203.201E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 284.957E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.89329E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 3.38079E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 271.375E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 6 RATIO = 435.047E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.460E-05   TIME AVG. FORCE       4.460E-05
 LARGEST RESIDUAL FORCE             3.834E-03   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.357       AT NODE       2463   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        0.357       AT NODE       2463   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.079E-02   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.824E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.824E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.61035E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 2.97061E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 908.708E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.520E-06   TIME AVG. FORCE       1.520E-06
 LARGEST RESIDUAL FORCE             3.260E-05   AT NODE       4229   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.         0.343       AT NODE       2463   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -1.458E-02   AT NODE       4157   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            8.288E-05   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.822E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.471E-05   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     7.048E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.17000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          1
 

 ***WARNING: THE SYSTEM MATRIX HAS 2 NEGATIVE EIGENVALUES.

     NEGATIVE EIGENVALUES MEAN THAT THE SYSTEM MATRIX IS NOT 
     POSITIVE DEFINITE:
     FOR EXAMPLE, A BIFURCATION (BUCKLING) LOAD MAY HAVE BEEN EXCEEDED.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.28062E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.07714E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.69991E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.588E-06   TIME AVG. FORCE       3.588E-06
 LARGEST RESIDUAL FORCE            -5.029E-04   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.632       AT NODE       4146   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        0.647       AT NODE       4230   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.120E-03   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.251E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.429E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     5.493E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 2.20587E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.17837E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 31.7879E+15 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.009E-06   TIME AVG. FORCE       1.009E-06
 LARGEST RESIDUAL FORCE             2.071E-06   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.755       AT NODE        221   DOF  3
   INSTANCE: ES_WALLS-1                                                                      
 LARGEST CORRECTION TO DISP.       -0.225       AT NODE       2463   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.268E-05   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.140E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.104E-04   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     5.363E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 832.422E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 3.84978E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.42738E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      9.904E-07   TIME AVG. FORCE       9.904E-07
 LARGEST RESIDUAL FORCE            -2.740E-06   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.728       AT NODE        221   DOF  3
   INSTANCE: ES_WALLS-1                                                                      
 LARGEST CORRECTION TO DISP.       -9.472E-02   AT NODE       4150   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.030E-05   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.296E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.560E-04   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     5.261E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.81995E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 3.04865E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 7.76424E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      9.877E-07   TIME AVG. FORCE       9.877E-07
 LARGEST RESIDUAL FORCE             2.229E-09   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.727       AT NODE        221   DOF  3
   INSTANCE: ES_WALLS-1                                                                      
 LARGEST CORRECTION TO DISP.       -2.563E-03   AT NODE        986   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.533E-08   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.292E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.682E-06   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     5.262E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.100    
 STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         1.00    


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 5.21836E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 757.412E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.51308E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.498E-06   TIME AVG. FORCE       1.743E-06
 LARGEST RESIDUAL FORCE            -9.191E-04   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.732       AT NODE        221   DOF  3
   INSTANCE: ES_WALLS-1                                                                      
 LARGEST CORRECTION TO DISP.        9.454E-02   AT NODE        586   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -6.094E-03   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.459E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.670E-04   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.004E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          1
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.03462E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.81505E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 3.52803E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.075E-06   TIME AVG. FORCE       3.031E-06
 LARGEST RESIDUAL FORCE             8.797E-05   AT NODE       1079   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.        -0.607       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        0.137       AT NODE       4227   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.646E-04   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.119E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     4.876E-04   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     8.234E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.17000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.02126E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.80217E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 892.291E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.781E-06   TIME AVG. FORCE       1.384E-06
 LARGEST RESIDUAL FORCE             2.726E-06   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.602       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.878E-03   AT NODE       4145   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.453E-05   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.162E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.641E-05   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     8.038E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.00399E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 2.18829E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.33913E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.733E-06   TIME AVG. FORCE       1.360E-06
 LARGEST RESIDUAL FORCE            -1.212E-06   AT NODE       4436   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.        -0.608       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -1.542E-02   AT NODE        769   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.015E-06   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.092E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.083E-05   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     7.877E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 593.204E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 7.12731E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 977.991E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.717E-06   TIME AVG. FORCE       1.352E-06
 LARGEST RESIDUAL FORCE            -1.318E-09   AT NODE       4118   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.        -0.608       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        6.449E-04   AT NODE       4165   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            7.381E-11   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.095E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.470E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     7.872E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.200    
 STEP TIME COMPLETED        2.00    ,  TOTAL TIME COMPLETED         2.00    


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          1
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.41973E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.02102E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 55.4630E+15 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.552E-06   TIME AVG. FORCE       1.752E-06
 LARGEST RESIDUAL FORCE            -1.056E-03   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.70       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        -2.69       AT NODE        807   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -6.668E-03   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.554E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     4.590E-04   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.694E-21   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.18000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.41228E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 32.4821E+15 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 556.469E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.452E-06   TIME AVG. FORCE       1.719E-06
 LARGEST RESIDUAL FORCE            -4.993E-06   AT NODE       4168   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.         -2.70       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        0.258       AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.504E-05   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.479E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     9.810E-05   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     9.559E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.89852E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 2.57730E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.72410E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.366E-06   TIME AVG. FORCE       1.690E-06
 LARGEST RESIDUAL FORCE             1.661E-08   AT NODE        863   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.         -2.70       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.972E-03   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -6.763E-08   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.493E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.380E-05   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     9.538E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  3.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 4.18792E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 678.531E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.08158E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.366E-06   TIME AVG. FORCE       1.690E-06
 LARGEST RESIDUAL FORCE            -2.724E-11   AT NODE       4198   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.         -2.70       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.097E-04   AT NODE       3879   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            6.354E-11   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.492E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.304E-07   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     9.538E-36   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.300    
 STEP TIME COMPLETED        3.00    ,  TOTAL TIME COMPLETED         3.00    


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 458.973E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.05172E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 4.93365E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.064E-06   TIME AVG. FORCE       2.034E-06
 LARGEST RESIDUAL FORCE            -1.664E-03   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.89       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        0.270       AT NODE        448   DOF  2
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.028E-02   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      5.068E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     5.760E-04   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.    -3.388E-21   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 329.577E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 6.91497E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.24529E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.989E-06   TIME AVG. FORCE       2.015E-06
 LARGEST RESIDUAL FORCE            -1.585E-06   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.89       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -4.995E-03   AT NODE        711   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.047E-05   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      5.075E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.548E-06   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.080E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.36091E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.86080E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 690.132E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.980E-06   TIME AVG. FORCE       2.013E-06
 LARGEST RESIDUAL FORCE             5.688E-10   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.89       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -7.932E-04   AT NODE       4116   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -4.032E-09   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      5.078E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.522E-06   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.080E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.400    
 STEP TIME COMPLETED        4.00    ,  TOTAL TIME COMPLETED         4.00    


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 2.58092E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 547.791E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.23259E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.660E-06   TIME AVG. FORCE       2.342E-06
 LARGEST RESIDUAL FORCE            -3.657E-03   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.75       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        0.298       AT NODE       3558   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.296E-02   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      5.758E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.805E-04   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.187E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 654.860E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 980.737E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.61630E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.571E-06   TIME AVG. FORCE       2.324E-06
 LARGEST RESIDUAL FORCE            -2.418E-06   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.76       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -1.268E-02   AT NODE        775   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.078E-05   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      5.837E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -8.362E-05   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.183E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.17000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.88875E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 341.851E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 875.687E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.570E-06   TIME AVG. FORCE       2.324E-06
 LARGEST RESIDUAL FORCE             3.475E-09   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.76       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        3.854E-04   AT NODE       3203   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.479E-08   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      5.838E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.149E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.183E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.500    
 STEP TIME COMPLETED        5.00    ,  TOTAL TIME COMPLETED         5.00    


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.62161E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.09299E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.64390E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.376E-06   TIME AVG. FORCE       2.666E-06
 LARGEST RESIDUAL FORCE            -8.909E-03   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.59       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -0.650       AT NODE        223   DOF  2
   INSTANCE: ES_WALLS-1                                                                      
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            5.844E-02   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      6.934E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.097E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     6.776E-21   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          1
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 6.89604E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.13826E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.72485E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.138E-06   TIME AVG. FORCE       2.626E-06
 LARGEST RESIDUAL FORCE             5.444E-07   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.55       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        0.169       AT NODE        449   DOF  1
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.833E-05   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.241E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.391E-04   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.272E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 440.867E+15  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 65.6254E+12 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 987.678E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.139E-06   TIME AVG. FORCE       2.627E-06
 LARGEST RESIDUAL FORCE             2.307E-09   AT NODE       3102   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.         -2.55       AT NODE       3393   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        3.279E-03   AT NODE       3143   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.872E-08   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.245E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.796E-06   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.272E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.600    
 STEP TIME COMPLETED        6.00    ,  TOTAL TIME COMPLETED         6.00    


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          1
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 2.06149E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 217.077E+15  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 403.812E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.451E-06   TIME AVG. FORCE       3.030E-06
 LARGEST RESIDUAL FORCE            -2.324E-02   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -5.27       AT NODE       4287   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.         2.91       AT NODE        449   DOF  1
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            0.157       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      6.795E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.062E-03   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.    -6.776E-21   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 2.19968E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 363.043E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.68007E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.794E-06   TIME AVG. FORCE       2.936E-06
 LARGEST RESIDUAL FORCE             9.099E-04   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.11       AT NODE       4223   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.         2.26       AT NODE       1702   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -4.300E-03   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.079E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.610E-03   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.349E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 2.16753E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.69192E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 967.412E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.687E-06   TIME AVG. FORCE       2.921E-06
 LARGEST RESIDUAL FORCE            -1.249E-06   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.12       AT NODE       4223   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -8.631E-03   AT NODE        289   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            8.938E-06   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.077E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.859E-05   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.368E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 882.345E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 2.01344E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.76537E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 34.3425E+09 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 117.183E+09  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.687E-06   TIME AVG. FORCE       2.921E-06
 LARGEST RESIDUAL FORCE            -1.743E-09   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.12       AT NODE       4223   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.984E-04   AT NODE       1939   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.586E-08   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.077E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.080E-07   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.368E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.700    
 STEP TIME COMPLETED        7.00    ,  TOTAL TIME COMPLETED         7.00    


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 2.69634E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.13330E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.68229E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.473E-06   TIME AVG. FORCE       3.240E-06
 LARGEST RESIDUAL FORCE            -9.706E-03   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.08       AT NODE       4222   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        4.490E-02   AT NODE        447   DOF  2
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            6.638E-02   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.096E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.115E-04   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.454E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          1
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 3.28034E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 696.691E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 22.8034E+15 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.221E-06   TIME AVG. FORCE       3.208E-06
 LARGEST RESIDUAL FORCE             5.196E-06   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.03       AT NODE       4222   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        5.032E-02   AT NODE        448   DOF  2
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.144E-05   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.117E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.774E-04   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.454E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.02276E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.09727E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.76799E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.222E-06   TIME AVG. FORCE       3.209E-06
 LARGEST RESIDUAL FORCE             1.433E-05   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -5.98       AT NODE       4253   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.757E-02   AT NODE        449   DOF  1
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -4.139E-04   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.080E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     5.689E-04   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.455E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 5.57574E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 2.01037E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 697.986E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.222E-06   TIME AVG. FORCE       3.209E-06
 LARGEST RESIDUAL FORCE            -4.537E-09   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -5.98       AT NODE       4253   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.060E-03   AT NODE       3394   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.410E-08   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.080E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.693E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.455E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.800    
 STEP TIME COMPLETED        8.00    ,  TOTAL TIME COMPLETED         8.00    


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.62767E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 5.04765E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 5.10386E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.140E-06   TIME AVG. FORCE       3.534E-06
 LARGEST RESIDUAL FORCE            -1.386E-02   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.54       AT NODE       4252   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -0.557       AT NODE          8   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            0.104       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.098E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.483E-04   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.534E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 5.38211E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 3.25347E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.05990E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.743E-06   TIME AVG. FORCE       3.490E-06
 LARGEST RESIDUAL FORCE             2.353E-05   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.61       AT NODE       4252   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -7.450E-02   AT NODE       3396   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.486E-04   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.126E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.053E-04   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.535E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 698.668E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.65889E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.73269E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.744E-06   TIME AVG. FORCE       3.490E-06
 LARGEST RESIDUAL FORCE            -4.569E-08   AT NODE          4   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.60       AT NODE       4252   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.340E-02   AT NODE       3925   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.572E-07   AT NODE          4   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.127E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.354E-05   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.535E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 3.64729E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 814.188E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.20014E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 33.7496E+09 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 66.9492E+09 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.744E-06   TIME AVG. FORCE       3.490E-06
 LARGEST RESIDUAL FORCE            -1.436E-10   AT NODE       3185   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.         -6.60       AT NODE       4252   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -4.168E-04   AT NODE       1938   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            3.715E-10   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.127E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.361E-08   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.534E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.900    
 STEP TIME COMPLETED        9.00    ,  TOTAL TIME COMPLETED         9.00    


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          1
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 799.222E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 740.118E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.66784E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.641E-06   TIME AVG. FORCE       3.805E-06
 LARGEST RESIDUAL FORCE            -1.381E-02   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.61       AT NODE       4251   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -2.168E-02   AT NODE       3913   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            0.103       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.144E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.040E-04   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.608E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 494.919E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.72241E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 810.456E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.254E-06   TIME AVG. FORCE       3.767E-06
 LARGEST RESIDUAL FORCE             2.516E-05   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.65       AT NODE       4220   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -3.952E-02   AT NODE       3396   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.433E-04   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.172E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.106E-04   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.609E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.86957E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.13478E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.44498E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.256E-06   TIME AVG. FORCE       3.767E-06
 LARGEST RESIDUAL FORCE            -6.161E-08   AT NODE          4   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -6.65       AT NODE       4220   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -3.960E-03   AT NODE       1700   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -3.522E-07   AT NODE          4   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.173E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.701E-05   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.609E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 2.84016E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 888.449E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 87.9956E+15 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 19.3539E+09 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 48.9842E+09 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.256E-06   TIME AVG. FORCE       3.767E-06
 LARGEST RESIDUAL FORCE            -1.322E-10   AT NODE       3185   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.         -6.65       AT NODE       4220   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        2.664E-03   AT NODE       1903   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -3.177E-10   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.173E-02   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.951E-08   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      1.000E-05   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.609E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED   1.00    
 STEP TIME COMPLETED        10.0    ,  TOTAL TIME COMPLETED         10.0    
1

   Abaqus 2020                                  Date 15-Dec-2021   Time 20:29:12
   For use by PURDUE UNIVERSITY under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    2     INCREMENT     1     STEP TIME    0.00    


                        S T E P       2     S T A T I C   A N A L Y S I S


                                                                                          

     AUTOMATIC TIME CONTROL WITH -
          A SUGGESTED INITIAL TIME INCREMENT OF                 1.00    
          AND A TOTAL TIME PERIOD OF                            10.0    
          THE MINIMUM TIME INCREMENT ALLOWED IS                1.000E-04
          THE MAXIMUM TIME INCREMENT ALLOWED IS                 1.00    

     LINEAR EQUATION SOLVER TYPE         DIRECT SPARSE

 CONVERGENCE TOLERANCE PARAMETERS FOR FORCE    
     CRITERION FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR DISP.    CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INIT. VALUE OF TIME AVG. FORCE     IS TIME AVG. FORCE     IN PREVIOUS STEP
     AVERAGE FORCE     IS TIME AVERAGE FORCE    
     ALTERNATE CRIT. FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO FORCE     RELATIVE TO TIME AVRG. FORCE         1.000E-05
     CRITERION FOR RESIDUAL FORCE     WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR DISP.    CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL FORCE     FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO FORCE     REL. TO TIME AVRG. MAX. FORCE        1.000E-05
     CRITERION FOR ZERO DISP.    RELATIVE TO CHARACTERISTIC LENGTH     1.000E-08

 CONVERGENCE TOLERANCE PARAMETERS FOR MOMENT   
     CRITERION FOR RESIDUAL MOMENT    FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR ROTATION CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INIT. VALUE OF TIME AVG. MOMENT    IS TIME AVG. MOMENT    IN PREVIOUS STEP
     AVERAGE MOMENT    IS TIME AVERAGE MOMENT   
     ALTERNATE CRIT. FOR RESIDUAL MOMENT    FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO MOMENT    RELATIVE TO TIME AVRG. MOMENT        1.000E-05
     CRITERION FOR RESIDUAL MOMENT    WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR ROTATION CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL MOMENT    FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO MOMENT    REL. TO TIME AVRG. MAX. MOMENT       1.000E-05

 CONVERGENCE TOLERANCE PARAMETERS FOR CAV. VOL.
     CRITERION FOR RESIDUAL CAV. VOL. FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR H.PRESS. CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INIT. VALUE OF TIME AVG. CAV. VOL. IS TIME AVG. CAV. VOL. IN PREVIOUS STEP
     AVERAGE CAV. VOL. IS TIME AVERAGE CAV. VOL.
     ALTERNATE CRIT. FOR RESIDUAL CAV. VOL. FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO CAV. VOL. RELATIVE TO TIME AVRG. CAV. VOL.     1.000E-05
     CRITERION FOR RESIDUAL CAV. VOL. WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR H.PRESS. CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL CAV. VOL. FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO CAV. VOL. REL. TO TIME AVRG. MAX. CAV. VOL.    1.000E-05

     VOLUMETRIC STRAIN COMPATIBILITY TOLERANCE FOR HYBRID SOLIDS       1.000E-05
     AXIAL STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS             1.000E-05
     TRANS. SHEAR STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS      1.000E-05
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P>P0          5.000E-03
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P=0.0         0.100    
     CONTACT FORCE ERROR TOLERANCE FOR CONVERT SDI=YES                 1.00    
     DISPLACEMENT COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS           1.000E-05
     ROTATION COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS               1.000E-05

 EQUILIBRIUM WILL BE CHECKED FOR SEVERE DISCONTINUITY ITERATIONS

 TIME INCREMENTATION CONTROL PARAMETERS:
     FIRST EQUILIBRIUM ITERATION FOR CONSECUTIVE DIVERGENCE CHECK              4
     EQUILIBRIUM ITERATION AT WHICH LOG. CONVERGENCE RATE CHECK BEGINS         8
     EQUILIBRIUM ITERATION AFTER WHICH ALTERNATE RESIDUAL IS USED              9
     MAXIMUM EQUILIBRIUM ITERATIONS ALLOWED                                   16
     EQUILIBRIUM ITERATION COUNT FOR CUT-BACK IN NEXT INCREMENT               10
     MAXIMUM EQUILIB. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      4
     MAXIMUM ITERATIONS FOR SEVERE DISCONTINUITIES                            50
     MAXIMUM ATTEMPTS ALLOWED IN AN INCREMENT                                  5
     MAXIMUM DISCON. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      50
     MAXIMUM CONTACT AUGMENTATIONS FOR *SURFACE BEHAVIOR,AUGMENTED LAGRANGE   50
     CUT-BACK FACTOR AFTER DIVERGENCE                                 0.2500    
     CUT-BACK FACTOR FOR TOO SLOW CONVERGENCE                         0.5000    
     CUT-BACK FACTOR AFTER TOO MANY EQUILIBRIUM ITERATIONS            0.7500    
     CUT-BACK FACTOR AFTER TOO MANY SEVERE DISCONTINUITY ITERATIONS   0.2500    
     CUT-BACK FACTOR AFTER PROBLEMS IN ELEMENT ASSEMBLY               0.2500    
     INCREASE FACTOR AFTER TWO INCREMENTS THAT CONVERGE QUICKLY        1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED                       1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DYNAMICS)            1.250    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DIFFUSION)           2.000    
     MINIMUM TIME INCREMENT RATIO FOR EXTRAPOLATION TO OCCUR          0.1000    
     MAX. RATIO OF TIME INCREMENT TO STABILITY LIMIT                   1.000    
     FRACTION OF STABILITY LIMIT FOR NEW TIME INCREMENT               0.9500    
     TIME INCREMENT INCREASE FACTOR BEFORE A TIME POINT                1.000    
     GLOBAL STABILIZATION CONTROL IS NOT USED

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                     100

          LARGE DISPLACEMENT THEORY WILL BE USED

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH     0.237    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  3.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 633.095E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.43157E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 3.12005E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 5.21221E+12.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 133.566E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 6 RATIO = 796.048E+09  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.540E-06   TIME AVG. FORCE       6.540E-06
 LARGEST RESIDUAL FORCE             2.151E-05   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         1.196E-02   AT NODE        769   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.196E-02   AT NODE        769   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.082E-04   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -5.085E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.085E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 581.624E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 400.263E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.59247E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.506E-06   TIME AVG. FORCE       6.506E-06
 LARGEST RESIDUAL FORCE             4.292E-09   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         1.167E-02   AT NODE        769   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -2.837E-04   AT NODE        769   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          DISP.    CORRECTION TOO LARGE COMPARED TO DISP.    INCREMENT

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.142E-08   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -5.243E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -7.384E-07   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.644E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          1
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 438.474E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 611.984E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 2.60181E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 65.3409E+09 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 62.2946E+09 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.506E-06   TIME AVG. FORCE       6.506E-06
 LARGEST RESIDUAL FORCE             2.683E-14   AT NODE        367   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.         1.169E-02   AT NODE        769   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        2.038E-04   AT NODE       4497   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.169E-13   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -5.243E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.398E-09   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.644E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.100    
 STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         11.0    


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.59870E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 897.852E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.91767E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.758E-06   TIME AVG. FORCE       6.632E-06
 LARGEST RESIDUAL FORCE             5.949E-09   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         1.116E-02   AT NODE        769   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.319E-04   AT NODE        769   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION   -7.542E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -6.458E-08   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -5.356E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.842E-07   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.    -1.355E-20   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.200    
 STEP TIME COMPLETED        2.00    ,  TOTAL TIME COMPLETED         12.0    


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.17000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 3.63648E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 17.1088E+15 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 4.94789E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.003E-06   TIME AVG. FORCE       6.756E-06
 LARGEST RESIDUAL FORCE             1.905E-09   AT NODE        756   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.         1.059E-02   AT NODE        769   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.838E-03   AT NODE       3818   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION   -3.270E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -6.575E-08   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -5.510E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.825E-07   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.    -1.355E-20   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.300    
 STEP TIME COMPLETED        3.00    ,  TOTAL TIME COMPLETED         13.0    


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 2.60769E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.49447E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.24025E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.248E-06   TIME AVG. FORCE       6.879E-06
 LARGEST RESIDUAL FORCE            -2.795E-09   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.042E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -2.436E-03   AT NODE       3818   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION    1.802E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -9.991E-08   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -5.648E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.816E-07   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.355E-20   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

 ITERATION SUMMARY FOR THE INCREMENT:   1 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  1 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.400    
 STEP TIME COMPLETED        4.00    ,  TOTAL TIME COMPLETED         14.0    


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 38.9931E+15 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 270.107E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 3.19418E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.489E-06   TIME AVG. FORCE       7.001E-06
 LARGEST RESIDUAL FORCE            -6.338E-09   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.337E-02   AT NODE       4466   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -3.098E-03   AT NODE       3818   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION    5.457E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.318E-07   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -5.795E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.492E-07   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.772E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.07720E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 5.08863E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.41142E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 24.9891E+09 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 25.0597E+09 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.489E-06   TIME AVG. FORCE       7.001E-06
 LARGEST RESIDUAL FORCE             1.049E-13   AT NODE        725   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.        -1.225E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.136E-03   AT NODE       4497   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION   -1.881E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            5.554E-13   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -5.788E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.508E-09   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.772E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.500    
 STEP TIME COMPLETED        5.00    ,  TOTAL TIME COMPLETED         15.0    


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 3.80371E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 14.1245E+15 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 223.272E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.729E-06   TIME AVG. FORCE       7.122E-06
 LARGEST RESIDUAL FORCE            -9.781E-09   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.161E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        3.214E-03   AT NODE       4497   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION   -9.151E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.619E-07   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -5.909E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.581E-07   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.    -1.355E-20   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.15000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 275.610E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 950.149E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 4.11435E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 29.9208E+09 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 29.3632E+09 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.729E-06   TIME AVG. FORCE       7.122E-06
 LARGEST RESIDUAL FORCE            -7.585E-14   AT NODE       4113   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.        -1.207E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -4.515E-04   AT NODE       4114   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION   -3.502E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.673E-13   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -5.908E-06   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.469E-08   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.801E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.600    
 STEP TIME COMPLETED        6.00    ,  TOTAL TIME COMPLETED         16.0    


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 2.15262E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.29646E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.50465E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.968E-06   TIME AVG. FORCE       7.243E-06
 LARGEST RESIDUAL FORCE            -1.304E-08   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.214E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.299E-03   AT NODE       4497   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION   -5.161E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.893E-07   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -6.082E-06   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.240E-07   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.830E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.66339E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 996.296E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 3.21994E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 32.9105E+09 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 33.4605E+09 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.968E-06   TIME AVG. FORCE       7.243E-06
 LARGEST RESIDUAL FORCE             6.227E-14   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.291E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        2.002E-03   AT NODE       4529   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -4.925E-13   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -6.082E-06   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.223E-09   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.830E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.700    
 STEP TIME COMPLETED        7.00    ,  TOTAL TIME COMPLETED         17.0    


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 430.338E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 668.495E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.51973E+15.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.205E-06   TIME AVG. FORCE       7.363E-06
 LARGEST RESIDUAL FORCE            -1.632E-08   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.235E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        5.596E-04   AT NODE       4228   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION   -2.902E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.166E-07   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -6.372E-06   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.901E-07   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.858E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 13.4789E+15 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 810.207E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.70059E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 36.3667E+09 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 37.0130E+09 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.205E-06   TIME AVG. FORCE       7.363E-06
 LARGEST RESIDUAL FORCE             1.069E-13   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.277E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -2.569E-03   AT NODE       3785   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION    1.682E-08
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.680E-13   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -6.370E-06   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.174E-09   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.858E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.800    
 STEP TIME COMPLETED        8.00    ,  TOTAL TIME COMPLETED         18.0    


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.17000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.66415E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 4.16484E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 339.421E+12  .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.439E-06   TIME AVG. FORCE       7.483E-06
 LARGEST RESIDUAL FORCE            -1.862E-08   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.275E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        3.098E-04   AT NODE       4497   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION   -1.911E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.380E-07   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -6.648E-06   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.779E-07   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.711E-20   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.76828E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 3.33989E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 4.62787E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 39.4020E+09 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 41.9876E+09 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.438E-06   TIME AVG. FORCE       7.483E-06
 LARGEST RESIDUAL FORCE             9.209E-13   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.093E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        2.236E-03   AT NODE       4531   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION   -1.106E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -5.819E-12   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -6.637E-06   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.122E-08   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.885E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.900    
 STEP TIME COMPLETED        9.00    ,  TOTAL TIME COMPLETED         19.0    


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.17000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 1.00046E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 535.182E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 10.8230E+15 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.672E-06   TIME AVG. FORCE       7.602E-06
 LARGEST RESIDUAL FORCE            -2.169E-08   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.163E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        9.222E-04   AT NODE       1820   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION   -6.922E-07
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.620E-07   AT NODE          4   DOF  4
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -6.890E-06   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.533E-07   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.912E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 1 RATIO = 696.286E+12  .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 2 RATIO = 1.10733E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 3 RATIO = 1.73543E+15.
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 4 RATIO = 45.4916E+09 .
 

 ***WARNING: SOLVER PROBLEM. NUMERICAL SINGULARITY WHEN PROCESSING NODE 
             ASSEMBLY.4 D.O.F. 5 RATIO = 49.7884E+09 .

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.672E-06   TIME AVG. FORCE       7.602E-06
 LARGEST RESIDUAL FORCE            -1.120E-13   AT NODE          4   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -1.183E-02   AT NODE       4197   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        2.155E-04   AT NODE       4529   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 ESTIMATE OF DISP.    CORRECTION    1.113E-09
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            7.980E-13   AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -6.887E-06   AT NODE          4   DOF  6
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -7.368E-09   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.912E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED   1.00    
 STEP TIME COMPLETED        10.0    ,  TOTAL TIME COMPLETED         20.0    
1

   Abaqus 2020                                  Date 15-Dec-2021   Time 20:29:32
   For use by PURDUE UNIVERSITY under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    3     INCREMENT     1     STEP TIME    0.00    


                        S T E P       3     S T A T I C   A N A L Y S I S


                                                                                          

     AUTOMATIC TIME CONTROL WITH -
          A SUGGESTED INITIAL TIME INCREMENT OF                 1.00    
          AND A TOTAL TIME PERIOD OF                            10.0    
          THE MINIMUM TIME INCREMENT ALLOWED IS                1.000E-04
          THE MAXIMUM TIME INCREMENT ALLOWED IS                 1.00    

     LINEAR EQUATION SOLVER TYPE         DIRECT SPARSE

 CONVERGENCE TOLERANCE PARAMETERS FOR FORCE    
     CRITERION FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR DISP.    CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INIT. VALUE OF TIME AVG. FORCE     IS TIME AVG. FORCE     IN PREVIOUS STEP
     AVERAGE FORCE     IS TIME AVERAGE FORCE    
     ALTERNATE CRIT. FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO FORCE     RELATIVE TO TIME AVRG. FORCE         1.000E-05
     CRITERION FOR RESIDUAL FORCE     WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR DISP.    CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL FORCE     FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO FORCE     REL. TO TIME AVRG. MAX. FORCE        1.000E-05
     CRITERION FOR ZERO DISP.    RELATIVE TO CHARACTERISTIC LENGTH     1.000E-08

 CONVERGENCE TOLERANCE PARAMETERS FOR MOMENT   
     CRITERION FOR RESIDUAL MOMENT    FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR ROTATION CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INIT. VALUE OF TIME AVG. MOMENT    IS TIME AVG. MOMENT    IN PREVIOUS STEP
     AVERAGE MOMENT    IS TIME AVERAGE MOMENT   
     ALTERNATE CRIT. FOR RESIDUAL MOMENT    FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO MOMENT    RELATIVE TO TIME AVRG. MOMENT        1.000E-05
     CRITERION FOR RESIDUAL MOMENT    WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR ROTATION CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL MOMENT    FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO MOMENT    REL. TO TIME AVRG. MAX. MOMENT       1.000E-05

 CONVERGENCE TOLERANCE PARAMETERS FOR CAV. VOL.
     CRITERION FOR RESIDUAL CAV. VOL. FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR H.PRESS. CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INIT. VALUE OF TIME AVG. CAV. VOL. IS TIME AVG. CAV. VOL. IN PREVIOUS STEP
     AVERAGE CAV. VOL. IS TIME AVERAGE CAV. VOL.
     ALTERNATE CRIT. FOR RESIDUAL CAV. VOL. FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO CAV. VOL. RELATIVE TO TIME AVRG. CAV. VOL.     1.000E-05
     CRITERION FOR RESIDUAL CAV. VOL. WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR H.PRESS. CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL CAV. VOL. FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO CAV. VOL. REL. TO TIME AVRG. MAX. CAV. VOL.    1.000E-05

     VOLUMETRIC STRAIN COMPATIBILITY TOLERANCE FOR HYBRID SOLIDS       1.000E-05
     AXIAL STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS             1.000E-05
     TRANS. SHEAR STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS      1.000E-05
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P>P0          5.000E-03
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P=0.0         0.100    
     CONTACT FORCE ERROR TOLERANCE FOR CONVERT SDI=YES                 1.00    
     DISPLACEMENT COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS           1.000E-05
     ROTATION COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS               1.000E-05

 EQUILIBRIUM WILL BE CHECKED FOR SEVERE DISCONTINUITY ITERATIONS

 TIME INCREMENTATION CONTROL PARAMETERS:
     FIRST EQUILIBRIUM ITERATION FOR CONSECUTIVE DIVERGENCE CHECK              4
     EQUILIBRIUM ITERATION AT WHICH LOG. CONVERGENCE RATE CHECK BEGINS         8
     EQUILIBRIUM ITERATION AFTER WHICH ALTERNATE RESIDUAL IS USED              9
     MAXIMUM EQUILIBRIUM ITERATIONS ALLOWED                                   16
     EQUILIBRIUM ITERATION COUNT FOR CUT-BACK IN NEXT INCREMENT               10
     MAXIMUM EQUILIB. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      4
     MAXIMUM ITERATIONS FOR SEVERE DISCONTINUITIES                            50
     MAXIMUM ATTEMPTS ALLOWED IN AN INCREMENT                                  5
     MAXIMUM DISCON. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      50
     MAXIMUM CONTACT AUGMENTATIONS FOR *SURFACE BEHAVIOR,AUGMENTED LAGRANGE   50
     CUT-BACK FACTOR AFTER DIVERGENCE                                 0.2500    
     CUT-BACK FACTOR FOR TOO SLOW CONVERGENCE                         0.5000    
     CUT-BACK FACTOR AFTER TOO MANY EQUILIBRIUM ITERATIONS            0.7500    
     CUT-BACK FACTOR AFTER TOO MANY SEVERE DISCONTINUITY ITERATIONS   0.2500    
     CUT-BACK FACTOR AFTER PROBLEMS IN ELEMENT ASSEMBLY               0.2500    
     INCREASE FACTOR AFTER TWO INCREMENTS THAT CONVERGE QUICKLY        1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED                       1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DYNAMICS)            1.250    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DIFFUSION)           2.000    
     MINIMUM TIME INCREMENT RATIO FOR EXTRAPOLATION TO OCCUR          0.1000    
     MAX. RATIO OF TIME INCREMENT TO STABILITY LIMIT                   1.000    
     FRACTION OF STABILITY LIMIT FOR NEW TIME INCREMENT               0.9500    
     TIME INCREMENT INCREASE FACTOR BEFORE A TIME POINT                1.000    
     GLOBAL STABILIZATION CONTROL IS NOT USED

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                     100

          LARGE DISPLACEMENT THEORY WILL BE USED

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH     0.237    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  3.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.19000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.088E-06   TIME AVG. FORCE       9.088E-06
 LARGEST RESIDUAL FORCE             2.054E-05   AT NODE       1825   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4213   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.         3.24       AT NODE       4213   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.902E-06   TIME AVG. FORCE       8.902E-06
 LARGEST RESIDUAL FORCE             2.052E-08   AT NODE       1826   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4213   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -3.102E-04   AT NODE       1702   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.923E-44   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.946E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.100    
 STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         21.0    


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.185E-06   TIME AVG. FORCE       9.044E-06
 LARGEST RESIDUAL FORCE             4.021E-05   AT NODE       1825   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4213   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.570E-04   AT NODE       4157   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.673E-43   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.979E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      9.131E-06   TIME AVG. FORCE       9.017E-06
 LARGEST RESIDUAL FORCE            -8.196E-08   AT NODE        292   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4214   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.944E-04   AT NODE       1944   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.763E-42   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.979E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      9.132E-06   TIME AVG. FORCE       9.017E-06
 LARGEST RESIDUAL FORCE             1.475E-09   AT NODE       3204   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4214   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        2.062E-05   AT NODE        351   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.673E-48   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     1.979E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.200    
 STEP TIME COMPLETED        2.00    ,  TOTAL TIME COMPLETED         22.0    


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.414E-06   TIME AVG. FORCE       9.150E-06
 LARGEST RESIDUAL FORCE             4.000E-05   AT NODE       1825   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4214   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.391E-04   AT NODE       4157   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.536E-43   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.012E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      9.360E-06   TIME AVG. FORCE       9.131E-06
 LARGEST RESIDUAL FORCE             7.896E-08   AT NODE       3208   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4214   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.946E-04   AT NODE       1944   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.760E-42   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.013E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      9.361E-06   TIME AVG. FORCE       9.132E-06
 LARGEST RESIDUAL FORCE             1.555E-09   AT NODE       3204   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4214   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        2.049E-05   AT NODE        351   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.394E-48   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.013E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.300    
 STEP TIME COMPLETED        3.00    ,  TOTAL TIME COMPLETED         23.0    


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.639E-06   TIME AVG. FORCE       9.259E-06
 LARGEST RESIDUAL FORCE             3.979E-05   AT NODE       1825   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4214   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.289E-04   AT NODE       4157   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.939E-18   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.045E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      9.586E-06   TIME AVG. FORCE       9.245E-06
 LARGEST RESIDUAL FORCE             7.953E-08   AT NODE       3208   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4214   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.946E-04   AT NODE       1944   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.757E-42   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.046E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.17000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      9.586E-06   TIME AVG. FORCE       9.245E-06
 LARGEST RESIDUAL FORCE             8.841E-10   AT NODE       3206   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4214   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        2.038E-05   AT NODE        351   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.203E-48   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.046E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.400    
 STEP TIME COMPLETED        4.00    ,  TOTAL TIME COMPLETED         24.0    


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      9.864E-06   TIME AVG. FORCE       9.369E-06
 LARGEST RESIDUAL FORCE             3.958E-05   AT NODE       1825   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4215   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.205E-04   AT NODE        449   DOF  2
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.939E-18   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.    -2.711E-20   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      9.811E-06   TIME AVG. FORCE       9.358E-06
 LARGEST RESIDUAL FORCE             8.823E-08   AT NODE       1705   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4215   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.945E-04   AT NODE       1944   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.753E-42   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.079E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      9.811E-06   TIME AVG. FORCE       9.358E-06
 LARGEST RESIDUAL FORCE            -1.058E-09   AT NODE       3205   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4215   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.948E-05   AT NODE        351   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.373E-48   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.079E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.500    
 STEP TIME COMPLETED        5.00    ,  TOTAL TIME COMPLETED         25.0    


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.009E-05   TIME AVG. FORCE       9.480E-06
 LARGEST RESIDUAL FORCE             3.937E-05   AT NODE       1825   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4215   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.203E-04   AT NODE        449   DOF  2
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.337E-19   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.711E-20   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.003E-05   TIME AVG. FORCE       9.471E-06
 LARGEST RESIDUAL FORCE             8.071E-08   AT NODE       3208   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4215   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.945E-04   AT NODE       1944   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.750E-42   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.112E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.003E-05   TIME AVG. FORCE       9.471E-06
 LARGEST RESIDUAL FORCE            -1.429E-09   AT NODE       3205   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4215   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.963E-05   AT NODE        351   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.147E-48   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.112E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.600    
 STEP TIME COMPLETED        6.00    ,  TOTAL TIME COMPLETED         26.0    


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.17000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.031E-05   TIME AVG. FORCE       9.591E-06
 LARGEST RESIDUAL FORCE             3.916E-05   AT NODE       1825   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4215   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.201E-04   AT NODE        449   DOF  2
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.168E-19   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.    -2.711E-20   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.026E-05   TIME AVG. FORCE       9.583E-06
 LARGEST RESIDUAL FORCE             8.133E-08   AT NODE       3208   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4216   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.945E-04   AT NODE       1944   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.747E-42   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.145E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.026E-05   TIME AVG. FORCE       9.583E-06
 LARGEST RESIDUAL FORCE             9.135E-10   AT NODE       3206   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4216   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.979E-05   AT NODE        351   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.000E-48   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.145E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.700    
 STEP TIME COMPLETED        7.00    ,  TOTAL TIME COMPLETED         27.0    


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.053E-05   TIME AVG. FORCE       9.701E-06
 LARGEST RESIDUAL FORCE             3.895E-05   AT NODE       1825   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4216   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.198E-04   AT NODE        449   DOF  2
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.735E-18   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.178E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.15000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.048E-05   TIME AVG. FORCE       9.695E-06
 LARGEST RESIDUAL FORCE             8.185E-08   AT NODE       3208   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4216   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.943E-04   AT NODE       1944   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.744E-42   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.178E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.048E-05   TIME AVG. FORCE       9.695E-06
 LARGEST RESIDUAL FORCE             9.184E-10   AT NODE       3206   DOF  1
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4216   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.979E-05   AT NODE        351   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.011E-48   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.178E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.800    
 STEP TIME COMPLETED        8.00    ,  TOTAL TIME COMPLETED         28.0    


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.075E-05   TIME AVG. FORCE       9.812E-06
 LARGEST RESIDUAL FORCE             3.873E-05   AT NODE       1825   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4216   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.195E-04   AT NODE        449   DOF  2
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.084E-19   AT NODE          4   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.711E-20   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.17000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.070E-05   TIME AVG. FORCE       9.807E-06
 LARGEST RESIDUAL FORCE             8.236E-08   AT NODE       3208   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4216   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.941E-04   AT NODE       1944   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.741E-42   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.211E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.070E-05   TIME AVG. FORCE       9.807E-06
 LARGEST RESIDUAL FORCE            -9.371E-10   AT NODE       3206   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4216   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.957E-05   AT NODE        351   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.038E-48   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.211E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED  0.900    
 STEP TIME COMPLETED        9.00    ,  TOTAL TIME COMPLETED         29.0    


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.18000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.128E-05   TIME AVG. FORCE       9.954E-06
 LARGEST RESIDUAL FORCE             3.852E-05   AT NODE       1825   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4216   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -6.191E-04   AT NODE        449   DOF  2
   INSTANCE: ES_FACE-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     8.674E-19   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.    -2.711E-20   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM RESPONSE WAS LINEAR IN THIS INCREMENT

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.16000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.17000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.123E-05   TIME AVG. FORCE       9.949E-06
 LARGEST RESIDUAL FORCE             8.285E-08   AT NODE       3208   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4217   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.       -5.939E-04   AT NODE       1944   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.738E-42   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.244E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

  NUMBER OF EQUATIONS =     13054     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  6.82E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.15000    
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.16000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.123E-05   TIME AVG. FORCE       9.949E-06
 LARGEST RESIDUAL FORCE            -9.698E-10   AT NODE       3206   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST INCREMENT OF DISP.          3.24       AT NODE       4217   DOF  2
   INSTANCE: SE_SKIN-1                                                                       
 LARGEST CORRECTION TO DISP.        1.928E-05   AT NODE        351   DOF  3
   INSTANCE: SE_SKIN-1                                                                       
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             0.00       AT NODE          4   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION     -7.562E-03   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.148E-48   AT NODE          4   DOF  4
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 AVERAGE CAV. VOL.                  1.333E-02   TIME AVG. CAV. VOL.   1.333E-02
 ALL CAV. VOL. RESIDUALS ARE ZERO
 LARGEST INCREMENT OF H.PRESS.      5.000E-06   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO H.PRESS.     2.244E-35   AT NODE          3   DOF  8
   INSTANCE: {assembly node}                                                                 
          THE CAV. VOL. EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED   1.00    ,  FRACTION OF STEP COMPLETED   1.00    
 STEP TIME COMPLETED        10.0    ,  TOTAL TIME COMPLETED         30.0    


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF         30  INCREMENTS
                       0  CUTBACKS IN AUTOMATIC INCREMENTATION
                      87  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                      87  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                      87  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                    1398  WARNING MESSAGES DURING USER INPUT PROCESSING
                     221  WARNING MESSAGES DURING ANALYSIS
                     200  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                      21  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   73.550    
       SYSTEM TIME (SEC)    =   3.2200    
       TOTAL CPU TIME (SEC) =   76.770    
       WALLCLOCK TIME (SEC) =         84
