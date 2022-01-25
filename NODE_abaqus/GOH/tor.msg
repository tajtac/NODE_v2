1

   Abaqus 2020                                  Date 15-Dec-2021   Time 20:26:56
   For use by PURDUE UNIVERSITY under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     S T A T I C   A N A L Y S I S


                                                                                          

     AUTOMATIC TIME CONTROL WITH -
          A SUGGESTED INITIAL TIME INCREMENT OF                0.100    
          AND A TOTAL TIME PERIOD OF                            10.0    
          THE MINIMUM TIME INCREMENT ALLOWED IS                1.000E-04
          THE MAXIMUM TIME INCREMENT ALLOWED IS                 10.0    

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

     CHARACTERISTIC ELEMENT LENGTH      3.07    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.100    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  9.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  2.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  0.11000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.159E-04   TIME AVG. FORCE       1.159E-04
 LARGEST RESIDUAL FORCE            -2.989E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.343       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.343       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      3.072E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            0.612       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.000E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.000E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.789E-05   TIME AVG. FORCE       3.789E-05
 LARGEST RESIDUAL FORCE            -7.260E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.310       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.256E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      3.072E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -3.919E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.000E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.377E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.799E-05   TIME AVG. FORCE       3.799E-05
 LARGEST RESIDUAL FORCE            -5.811E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.299       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.144E-02   AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      3.072E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            2.597E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.000E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.818E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.801E-05   TIME AVG. FORCE       3.801E-05
 LARGEST RESIDUAL FORCE             1.558E-08   AT NODE       1079   DOF  1
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.        -0.299       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.424E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      3.072E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -8.995E-08   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.000E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.173E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.100    ,  FRACTION OF STEP COMPLETED  1.000E-02
 STEP TIME COMPLETED       0.100    ,  TOTAL TIME COMPLETED        0.100    


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.100    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.761E-05   TIME AVG. FORCE       6.281E-05
 LARGEST RESIDUAL FORCE             5.045E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.299       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.289E-03   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            0.679       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.000E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.677E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.639E-05   TIME AVG. FORCE       5.720E-05
 LARGEST RESIDUAL FORCE             5.760E-05   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.299       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.099E-03   AT NODE         74   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.266E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.000E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.585E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      7.639E-05   TIME AVG. FORCE       5.720E-05
 LARGEST RESIDUAL FORCE             1.835E-07   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.299       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.346E-04   AT NODE        113   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.100E-07   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.000E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.622E-07   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          TIME INCREMENT MAY NOW INCREASE TO    0.150    

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.100    ,  FRACTION OF STEP COMPLETED  2.000E-02
 STEP TIME COMPLETED       0.200    ,  TOTAL TIME COMPLETED        0.200    


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.150    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.561E-04   TIME AVG. FORCE       9.017E-05
 LARGEST RESIDUAL FORCE             9.459E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.450       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.162E-03   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             1.28       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.500E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.906E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.347E-04   TIME AVG. FORCE       8.302E-05
 LARGEST RESIDUAL FORCE            -1.239E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.450       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.947E-03   AT NODE         74   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.202E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.500E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     4.926E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.347E-04   TIME AVG. FORCE       8.302E-05
 LARGEST RESIDUAL FORCE             5.325E-07   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.450       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.122E-04   AT NODE         74   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -3.990E-06   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.500E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.604E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.347E-04   TIME AVG. FORCE       8.302E-05
 LARGEST RESIDUAL FORCE            -1.449E-12   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.450       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.564E-06   AT NODE        111   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.053E-08   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.500E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -7.227E-08   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THERE IS ZERO MOMENT    EVERYWHERE
          TIME INCREMENT MAY NOW INCREASE TO    0.225    

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.150    ,  FRACTION OF STEP COMPLETED  3.500E-02
 STEP TIME COMPLETED       0.350    ,  TOTAL TIME COMPLETED        0.350    


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.225    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.758E-04   TIME AVG. FORCE       1.312E-04
 LARGEST RESIDUAL FORCE             0.213       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.677       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.381E-02   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             2.89       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      2.250E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.114E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.237E-04   TIME AVG. FORCE       1.182E-04
 LARGEST RESIDUAL FORCE             4.427E-04   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.679       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.851E-03   AT NODE         74   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.530E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      2.250E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.296E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.237E-04   TIME AVG. FORCE       1.182E-04
 LARGEST RESIDUAL FORCE             1.573E-06   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.678       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.632E-03   AT NODE         73   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.877E-06   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      2.250E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.737E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.237E-04   TIME AVG. FORCE       1.182E-04
 LARGEST RESIDUAL FORCE            -1.102E-10   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.678       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.517E-05   AT NODE        106   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.312E-07   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      2.250E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.425E-07   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          TIME INCREMENT MAY NOW INCREASE TO    0.338    

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.225    ,  FRACTION OF STEP COMPLETED  5.750E-02
 STEP TIME COMPLETED       0.575    ,  TOTAL TIME COMPLETED        0.575    


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.338    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.862E-04   TIME AVG. FORCE       1.918E-04
 LARGEST RESIDUAL FORCE             0.478       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.083E-02   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             6.53       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.375E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.597E-04   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.611E-04   TIME AVG. FORCE       1.668E-04
 LARGEST RESIDUAL FORCE             1.514E-03   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.03       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.057E-02   AT NODE         92   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.239E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.375E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     5.837E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.611E-04   TIME AVG. FORCE       1.668E-04
 LARGEST RESIDUAL FORCE            -8.298E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.155E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -4.247E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.375E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.016E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.611E-04   TIME AVG. FORCE       1.668E-04
 LARGEST RESIDUAL FORCE             2.116E-07   AT NODE        126   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST INCREMENT OF DISP.         -1.02       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.505E-04   AT NODE         88   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.380E-06   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.375E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.136E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          TIME INCREMENT MAY NOW INCREASE TO    0.506    

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.338    ,  FRACTION OF STEP COMPLETED  9.125E-02
 STEP TIME COMPLETED       0.913    ,  TOTAL TIME COMPLETED        0.913    


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.506    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.721E-04   TIME AVG. FORCE       2.843E-04
 LARGEST RESIDUAL FORCE              1.07       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.54       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.851E-02   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             14.8       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      5.063E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.392E-04   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.

     NEGATIVE EIGENVALUES MEAN THAT THE SYSTEM MATRIX IS NOT 
     POSITIVE DEFINITE:
     FOR EXAMPLE, A BIFURCATION (BUCKLING) LOAD MAY HAVE BEEN EXCEEDED.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.772E-04   TIME AVG. FORCE       2.352E-04
 LARGEST RESIDUAL FORCE             2.648E-03   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.49       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.486E-02   AT NODE         92   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.777E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      5.063E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.956E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.765E-04   TIME AVG. FORCE       2.350E-04
 LARGEST RESIDUAL FORCE            -2.482E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.54       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.359E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -3.069E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      5.062E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.290E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.770E-04   TIME AVG. FORCE       2.351E-04
 LARGEST RESIDUAL FORCE            -6.015E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.54       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.274E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.583E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      5.062E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.836E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          TIME INCREMENT MAY NOW INCREASE TO    0.759    

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.506    ,  FRACTION OF STEP COMPLETED  0.142    
 STEP TIME COMPLETED        1.42    ,  TOTAL TIME COMPLETED         1.42    


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.759    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.613E-03   TIME AVG. FORCE       4.320E-04
 LARGEST RESIDUAL FORCE              2.41       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.31       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.153       AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             33.5       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.595E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.272E-03   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      9.263E-04   TIME AVG. FORCE       3.339E-04
 LARGEST RESIDUAL FORCE             7.146E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.31       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.139       AT NODE         74   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            6.876E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.379E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  5.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  6.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      9.250E-04   TIME AVG. FORCE       3.337E-04
 LARGEST RESIDUAL FORCE            -4.251E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.32       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.198E-02   AT NODE         90   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            9.016E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.592E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.252E-03   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      9.256E-04   TIME AVG. FORCE       3.338E-04
 LARGEST RESIDUAL FORCE            -8.903E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.31       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.359E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -7.542E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.859E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      9.256E-04   TIME AVG. FORCE       3.338E-04
 LARGEST RESIDUAL FORCE            -1.012E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.31       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.611E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.528E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.510E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.759    ,  FRACTION OF STEP COMPLETED  0.218    
 STEP TIME COMPLETED        2.18    ,  TOTAL TIME COMPLETED         2.18    


  INCREMENT     8 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.759    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.117E-03   TIME AVG. FORCE       5.567E-04
 LARGEST RESIDUAL FORCE              2.88       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.30       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.183       AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             40.2       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.595E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.471E-03   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.310E-03   TIME AVG. FORCE       4.557E-04
 LARGEST RESIDUAL FORCE             8.092E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.33       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.178       AT NODE        109   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.515E-02   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.595E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.055E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.308E-03   TIME AVG. FORCE       4.555E-04
 LARGEST RESIDUAL FORCE            -9.518E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.33       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.104       AT NODE         90   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.602E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.591E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.566E-03   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.309E-03   TIME AVG. FORCE       4.557E-04
 LARGEST RESIDUAL FORCE            -3.863E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.31       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.854E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.232E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.018E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.309E-03   TIME AVG. FORCE       4.557E-04
 LARGEST RESIDUAL FORCE            -1.325E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.30       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.819E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -7.859E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.443E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.309E-03   TIME AVG. FORCE       4.557E-04
 LARGEST RESIDUAL FORCE            -1.365E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.30       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.872E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            8.742E-06   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.082E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.759    ,  FRACTION OF STEP COMPLETED  0.294    
 STEP TIME COMPLETED        2.94    ,  TOTAL TIME COMPLETED         2.94    


  INCREMENT     9 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.759    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.510E-03   TIME AVG. FORCE       6.840E-04
 LARGEST RESIDUAL FORCE              2.88       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.28       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.181       AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             40.1       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.595E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.254E-03   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.734E-03   TIME AVG. FORCE       5.978E-04
 LARGEST RESIDUAL FORCE             7.939E-03   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.31       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.180       AT NODE        109   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.986E-02   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.596E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.954E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.732E-03   TIME AVG. FORCE       5.976E-04
 LARGEST RESIDUAL FORCE            -1.462E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.32       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.108       AT NODE         90   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -3.997E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.590E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     4.130E-03   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.734E-03   TIME AVG. FORCE       5.978E-04
 LARGEST RESIDUAL FORCE            -7.028E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.28       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.976E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -3.405E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.384E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.735E-03   TIME AVG. FORCE       5.978E-04
 LARGEST RESIDUAL FORCE            -5.467E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.28       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.985E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.139E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.100E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      1.735E-03   TIME AVG. FORCE       5.978E-04
 LARGEST RESIDUAL FORCE            -8.102E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.28       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.737E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.850E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.479E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.759    ,  FRACTION OF STEP COMPLETED  0.370    
 STEP TIME COMPLETED        3.70    ,  TOTAL TIME COMPLETED         3.70    


  INCREMENT    10 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.759    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.948E-03   TIME AVG. FORCE       8.329E-04
 LARGEST RESIDUAL FORCE              2.87       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.24       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.177       AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             39.8       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.595E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.279E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.207E-03   TIME AVG. FORCE       7.588E-04
 LARGEST RESIDUAL FORCE             9.980E-03   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.27       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.183       AT NODE        109   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.114E-02   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.596E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.830E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.206E-03   TIME AVG. FORCE       7.586E-04
 LARGEST RESIDUAL FORCE            -2.239E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.30       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.106       AT NODE         90   DOF  1
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.788E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.590E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     4.936E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.208E-03   TIME AVG. FORCE       7.589E-04
 LARGEST RESIDUAL FORCE            -1.178E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.25       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.332E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -5.437E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.796E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.209E-03   TIME AVG. FORCE       7.589E-04
 LARGEST RESIDUAL FORCE            -2.114E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.25       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.104E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -5.477E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.163E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      2.209E-03   TIME AVG. FORCE       7.589E-04
 LARGEST RESIDUAL FORCE            -3.767E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.25       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.073E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            7.928E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.182E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      2.209E-03   TIME AVG. FORCE       7.589E-04
 LARGEST RESIDUAL FORCE            -7.333E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.25       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.732E-05   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.034E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     4.154E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.759    ,  FRACTION OF STEP COMPLETED  0.446    
 STEP TIME COMPLETED        4.46    ,  TOTAL TIME COMPLETED         4.46    


  INCREMENT    11 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.759    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.442E-03   TIME AVG. FORCE       1.003E-03
 LARGEST RESIDUAL FORCE              2.86       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.19       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.173       AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             39.5       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.595E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.456E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.738E-03   TIME AVG. FORCE       9.388E-04
 LARGEST RESIDUAL FORCE             1.005E-02   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.22       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.185       AT NODE        109   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.833E-02   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.597E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.273E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.735E-03   TIME AVG. FORCE       9.386E-04
 LARGEST RESIDUAL FORCE            -3.065E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.27       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.108       AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -3.378E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.589E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     5.754E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.739E-03   TIME AVG. FORCE       9.389E-04
 LARGEST RESIDUAL FORCE            -1.886E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.20       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.949E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -7.574E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.277E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.740E-03   TIME AVG. FORCE       9.390E-04
 LARGEST RESIDUAL FORCE            -6.478E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.19       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.339E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.208E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.749E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      2.740E-03   TIME AVG. FORCE       9.390E-04
 LARGEST RESIDUAL FORCE            -1.512E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.20       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.228E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.011E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.375E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      2.740E-03   TIME AVG. FORCE       9.390E-04
 LARGEST RESIDUAL FORCE            -5.812E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.20       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.253E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            3.627E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.157E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.759    ,  FRACTION OF STEP COMPLETED  0.522    
 STEP TIME COMPLETED        5.22    ,  TOTAL TIME COMPLETED         5.22    


  INCREMENT    12 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.759    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.001E-03   TIME AVG. FORCE       1.194E-03
 LARGEST RESIDUAL FORCE              2.85       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.13       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.168       AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             38.9       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.580E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.335E-03   TIME AVG. FORCE       1.139E-03
 LARGEST RESIDUAL FORCE             9.917E-03   AT NODE          1   DOF  2
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.15       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.185       AT NODE        109   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.933E-02   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.597E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.939E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.331E-03   TIME AVG. FORCE       1.138E-03
 LARGEST RESIDUAL FORCE            -4.246E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.22       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.120       AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -5.182E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.588E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.956E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.336E-03   TIME AVG. FORCE       1.139E-03
 LARGEST RESIDUAL FORCE            -2.737E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.14       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.739E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.029E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.754E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.337E-03   TIME AVG. FORCE       1.139E-03
 LARGEST RESIDUAL FORCE            -2.005E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.13       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.384E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.643E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.586E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      3.337E-03   TIME AVG. FORCE       1.139E-03
 LARGEST RESIDUAL FORCE            -4.937E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.13       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.251E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.528E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.148E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      3.337E-03   TIME AVG. FORCE       1.139E-03
 LARGEST RESIDUAL FORCE            -4.324E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.13       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.676E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.224E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.147E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      3.338E-03   TIME AVG. FORCE       1.139E-03
 LARGEST RESIDUAL FORCE            -1.131E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.13       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.041E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.137E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.482E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.759    ,  FRACTION OF STEP COMPLETED  0.598    
 STEP TIME COMPLETED        5.98    ,  TOTAL TIME COMPLETED         5.98    


  INCREMENT    13 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.759    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.647E-03   TIME AVG. FORCE       1.409E-03
 LARGEST RESIDUAL FORCE              2.84       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.05       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.162       AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             38.1       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.660E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.011E-03   TIME AVG. FORCE       1.360E-03
 LARGEST RESIDUAL FORCE            -1.344E-02   AT NODE          1   DOF  3
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.07       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.182       AT NODE         74   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            6.482E-02   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.598E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.388E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.005E-03   TIME AVG. FORCE       1.359E-03
 LARGEST RESIDUAL FORCE            -5.578E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.17       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.133       AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -7.444E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.588E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.999E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.012E-03   TIME AVG. FORCE       1.360E-03
 LARGEST RESIDUAL FORCE            -3.694E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.07       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.105       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.315E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.593E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.207E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.014E-03   TIME AVG. FORCE       1.360E-03
 LARGEST RESIDUAL FORCE            -5.294E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.05       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.376E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -5.226E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.064E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      4.015E-03   TIME AVG. FORCE       1.360E-03
 LARGEST RESIDUAL FORCE            -1.421E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.05       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.532E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            9.490E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.941E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  5.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      4.015E-03   TIME AVG. FORCE       1.360E-03
 LARGEST RESIDUAL FORCE            -2.585E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.06       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.732E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            3.649E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.639E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      4.015E-03   TIME AVG. FORCE       1.360E-03
 LARGEST RESIDUAL FORCE            -7.331E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.06       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.435E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -6.649E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.390E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      4.015E-03   TIME AVG. FORCE       1.360E-03
 LARGEST RESIDUAL FORCE            -1.329E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.06       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.242E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.625E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.478E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   9 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  9 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.759    ,  FRACTION OF STEP COMPLETED  0.673    
 STEP TIME COMPLETED        6.73    ,  TOTAL TIME COMPLETED         6.73    


  INCREMENT    14 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.759    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.392E-03   TIME AVG. FORCE       1.648E-03
 LARGEST RESIDUAL FORCE              2.82       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.96       AT NODE         90   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.154       AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             36.9       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.633E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.781E-03   TIME AVG. FORCE       1.604E-03
 LARGEST RESIDUAL FORCE            -1.864E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.96       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.186       AT NODE         74   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.459E-02   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.598E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.134E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.772E-03   TIME AVG. FORCE       1.604E-03
 LARGEST RESIDUAL FORCE            -7.587E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.10       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.160       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            8.393E-02   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.587E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     9.474E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.783E-03   TIME AVG. FORCE       1.605E-03
 LARGEST RESIDUAL FORCE            -4.452E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.99       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.122       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.567E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.593E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.547E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.785E-03   TIME AVG. FORCE       1.605E-03
 LARGEST RESIDUAL FORCE            -1.408E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.95       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.167E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.015E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.736E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      4.785E-03   TIME AVG. FORCE       1.605E-03
 LARGEST RESIDUAL FORCE            -3.454E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.96       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.242E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.787E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.967E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      4.785E-03   TIME AVG. FORCE       1.605E-03
 LARGEST RESIDUAL FORCE            -1.498E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.97       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.498E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.072E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.837E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      4.785E-03   TIME AVG. FORCE       1.605E-03
 LARGEST RESIDUAL FORCE            -3.825E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.97       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.316E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.810E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.108E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      4.786E-03   TIME AVG. FORCE       1.605E-03
 LARGEST RESIDUAL FORCE            -1.659E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.96       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.734E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.130E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.934E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION    10


 AVERAGE FORCE                      4.786E-03   TIME AVG. FORCE       1.605E-03
 LARGEST RESIDUAL FORCE            -4.250E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.96       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.387E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED BASED ON ALTERNATIVE
          RESIDUAL CRITERION

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.915E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.277E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:  10 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND 10 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.759    ,  FRACTION OF STEP COMPLETED  0.749    
 STEP TIME COMPLETED        7.49    ,  TOTAL TIME COMPLETED         7.49    


  INCREMENT    15 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.759    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.251E-03   TIME AVG. FORCE       1.914E-03
 LARGEST RESIDUAL FORCE              2.80       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.88       AT NODE         90   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.153       AT NODE         72   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            -37.5       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.692E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 1 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.662E-03   TIME AVG. FORCE       1.875E-03
 LARGEST RESIDUAL FORCE            -2.849E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.84       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.192       AT NODE         74   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            6.729E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.599E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -7.199E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.647E-03   TIME AVG. FORCE       1.874E-03
 LARGEST RESIDUAL FORCE            -0.101       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -2.02       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -0.217       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            0.124       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.587E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.114E-02   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.664E-03   TIME AVG. FORCE       1.875E-03
 LARGEST RESIDUAL FORCE            -4.778E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.90       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.136       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.849E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.592E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.698E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      5.667E-03   TIME AVG. FORCE       1.876E-03
 LARGEST RESIDUAL FORCE            -3.539E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.83       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.072E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.859E-02   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.758E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      5.668E-03   TIME AVG. FORCE       1.876E-03
 LARGEST RESIDUAL FORCE            -7.559E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.85       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.910E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            3.020E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     4.130E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      5.668E-03   TIME AVG. FORCE       1.876E-03
 LARGEST RESIDUAL FORCE            -8.047E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.86       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.115E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            3.038E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     4.256E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      5.668E-03   TIME AVG. FORCE       1.876E-03
 LARGEST RESIDUAL FORCE            -1.752E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.86       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.933E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -4.381E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.215E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      5.668E-03   TIME AVG. FORCE       1.876E-03
 LARGEST RESIDUAL FORCE            -1.876E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.86       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.702E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -4.600E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      7.594E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.500E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO
 

 ***NOTE: THE RATE OF CONVERGENCE IS VERY SLOW. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    15 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  0.380    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.357E-03   TIME AVG. FORCE       1.855E-03
 LARGEST RESIDUAL FORCE              1.06       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.944       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.771E-02   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             13.8       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.797E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     5.871E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.196E-03   TIME AVG. FORCE       1.844E-03
 LARGEST RESIDUAL FORCE            -0.774       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.92       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        -1.09       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            0.956       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      4.889E-02   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.894E-02   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.195E-03   TIME AVG. FORCE       1.844E-03
 LARGEST RESIDUAL FORCE            -5.188E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         -1.68       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.272       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -0.198       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.966E-02   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -9.183E-03   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.134E-03   TIME AVG. FORCE       1.840E-03
 LARGEST RESIDUAL FORCE            -0.358       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.935       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        0.846       AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -7.657E-02   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.793E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.901E-02   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  5.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      5.211E-03   TIME AVG. FORCE       1.845E-03
 LARGEST RESIDUAL FORCE            -5.552E-03   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.917       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.868E-02   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.488E-02   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.798E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.642E-03   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      5.212E-03   TIME AVG. FORCE       1.845E-03
 LARGEST RESIDUAL FORCE            -1.513E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.931       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.550E-02   AT NODE        343   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            3.805E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.797E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.355E-04   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      5.211E-03   TIME AVG. FORCE       1.845E-03
 LARGEST RESIDUAL FORCE            -1.640E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.947       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.783E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.200E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.797E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.342E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      5.212E-03   TIME AVG. FORCE       1.845E-03
 LARGEST RESIDUAL FORCE            -1.653E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.945       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.943E-03   AT NODE        109   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -4.351E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.797E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.846E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      5.212E-03   TIME AVG. FORCE       1.845E-03
 LARGEST RESIDUAL FORCE            -2.649E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.943       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.268E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -5.170E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.797E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -8.060E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO
 

 ***NOTE: THE RATE OF CONVERGENCE IS VERY SLOW. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    15 STARTS. ATTEMPT NUMBER  3, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.034E-03   TIME AVG. FORCE       1.833E-03
 LARGEST RESIDUAL FORCE             0.442       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.476       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.509E-02   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             5.86       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.792E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.995E-03   TIME AVG. FORCE       1.831E-03
 LARGEST RESIDUAL FORCE             3.667E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.501       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.527E-02   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            6.442E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.008E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.995E-03   TIME AVG. FORCE       1.831E-03
 LARGEST RESIDUAL FORCE            -3.673E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.474       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.004E-02   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -5.743E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -9.509E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.995E-03   TIME AVG. FORCE       1.831E-03
 LARGEST RESIDUAL FORCE            -1.558E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.472       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.587E-03   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            9.877E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.570E-04   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.995E-03   TIME AVG. FORCE       1.831E-03
 LARGEST RESIDUAL FORCE            -4.549E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.475       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.468E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            6.795E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.098E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      4.995E-03   TIME AVG. FORCE       1.831E-03
 LARGEST RESIDUAL FORCE            -2.089E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.475       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.989E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.110E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.816E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      4.995E-03   TIME AVG. FORCE       1.831E-03
 LARGEST RESIDUAL FORCE            -6.100E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.475       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.016E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -7.792E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.272E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      4.995E-03   TIME AVG. FORCE       1.831E-03
 LARGEST RESIDUAL FORCE            -2.799E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.475       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.461E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.290E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.102E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.768    
 STEP TIME COMPLETED        7.68    ,  TOTAL TIME COMPLETED         7.68    


  INCREMENT    16 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.222E-03   TIME AVG. FORCE       2.043E-03
 LARGEST RESIDUAL FORCE             0.177       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.469       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.038E-02   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT             2.35       AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.483E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.212E-03   TIME AVG. FORCE       2.042E-03
 LARGEST RESIDUAL FORCE             9.416E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.472       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.667E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.882E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.466E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.212E-03   TIME AVG. FORCE       2.042E-03
 LARGEST RESIDUAL FORCE            -1.177E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.468       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.229E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.113E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.737E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.212E-03   TIME AVG. FORCE       2.042E-03
 LARGEST RESIDUAL FORCE            -2.976E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.468       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.128E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.658E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.564E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      5.212E-03   TIME AVG. FORCE       2.042E-03
 LARGEST RESIDUAL FORCE            -1.895E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.468       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.667E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.424E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.209E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      5.212E-03   TIME AVG. FORCE       2.042E-03
 LARGEST RESIDUAL FORCE            -4.812E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.468       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.249E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.097E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.260E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.787    
 STEP TIME COMPLETED        7.87    ,  TOTAL TIME COMPLETED         7.87    


  INCREMENT    17 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.447E-03   TIME AVG. FORCE       2.242E-03
 LARGEST RESIDUAL FORCE             0.176       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.463       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.030E-02   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            -2.37       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.524E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.436E-03   TIME AVG. FORCE       2.242E-03
 LARGEST RESIDUAL FORCE             1.242E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.465       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.728E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.836E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.571E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.436E-03   TIME AVG. FORCE       2.242E-03
 LARGEST RESIDUAL FORCE            -1.323E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.460       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.464E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.219E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.815E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.436E-03   TIME AVG. FORCE       2.242E-03
 LARGEST RESIDUAL FORCE            -3.371E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.461       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.457E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.888E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.781E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      5.436E-03   TIME AVG. FORCE       2.242E-03
 LARGEST RESIDUAL FORCE            -2.509E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.461       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.628E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.713E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.532E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      5.436E-03   TIME AVG. FORCE       2.242E-03
 LARGEST RESIDUAL FORCE            -6.560E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.461       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.217E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.619E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.880E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.806    
 STEP TIME COMPLETED        8.06    ,  TOTAL TIME COMPLETED         8.06    


  INCREMENT    18 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.678E-03   TIME AVG. FORCE       2.433E-03
 LARGEST RESIDUAL FORCE             0.176       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.459       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.020E-02   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            -2.43       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.621E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.668E-03   TIME AVG. FORCE       2.432E-03
 LARGEST RESIDUAL FORCE             1.072E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.        -0.458       AT NODE         18   DOF  3
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.740E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.828E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.695E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.668E-03   TIME AVG. FORCE       2.432E-03
 LARGEST RESIDUAL FORCE            -1.418E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.457       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.675E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.327E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.889E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.668E-03   TIME AVG. FORCE       2.432E-03
 LARGEST RESIDUAL FORCE            -3.811E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.457       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.797E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.142E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.007E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      5.668E-03   TIME AVG. FORCE       2.432E-03
 LARGEST RESIDUAL FORCE            -3.290E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.457       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.693E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.051E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.889E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      5.668E-03   TIME AVG. FORCE       2.432E-03
 LARGEST RESIDUAL FORCE            -8.911E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.457       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.334E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -3.255E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     4.599E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.825    
 STEP TIME COMPLETED        8.25    ,  TOTAL TIME COMPLETED         8.25    


  INCREMENT    19 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.918E-03   TIME AVG. FORCE       2.615E-03
 LARGEST RESIDUAL FORCE             0.175       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.464       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.001E-02   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            -2.50       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.578E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      5.908E-03   TIME AVG. FORCE       2.615E-03
 LARGEST RESIDUAL FORCE             1.004E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.464       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.696E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.827E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.421E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      5.908E-03   TIME AVG. FORCE       2.615E-03
 LARGEST RESIDUAL FORCE            -1.517E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.463       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.824E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.432E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.943E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      5.908E-03   TIME AVG. FORCE       2.615E-03
 LARGEST RESIDUAL FORCE            -4.134E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.462       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.068E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.386E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.193E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      5.908E-03   TIME AVG. FORCE       2.615E-03
 LARGEST RESIDUAL FORCE            -4.217E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.463       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.779E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.424E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.256E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      5.908E-03   TIME AVG. FORCE       2.615E-03
 LARGEST RESIDUAL FORCE            -1.161E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.463       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.493E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -3.972E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     5.352E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.844    
 STEP TIME COMPLETED        8.44    ,  TOTAL TIME COMPLETED         8.44    


  INCREMENT    20 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.167E-03   TIME AVG. FORCE       2.793E-03
 LARGEST RESIDUAL FORCE             0.175       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.470       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.001E-02   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            -2.55       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.684E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.156E-03   TIME AVG. FORCE       2.792E-03
 LARGEST RESIDUAL FORCE             8.921E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.469       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.704E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.004E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.613E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.156E-03   TIME AVG. FORCE       2.792E-03
 LARGEST RESIDUAL FORCE            -1.686E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.468       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.104E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.573E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.034E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  5.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.156E-03   TIME AVG. FORCE       2.792E-03
 LARGEST RESIDUAL FORCE            -4.775E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.468       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.492E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.723E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.474E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.156E-03   TIME AVG. FORCE       2.792E-03
 LARGEST RESIDUAL FORCE            -5.601E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.468       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.122E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.918E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.737E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  5.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      6.156E-03   TIME AVG. FORCE       2.792E-03
 LARGEST RESIDUAL FORCE            -1.612E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.468       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.009E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -4.966E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.382E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.863    
 STEP TIME COMPLETED        8.63    ,  TOTAL TIME COMPLETED         8.63    


  INCREMENT    21 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.424E-03   TIME AVG. FORCE       2.965E-03
 LARGEST RESIDUAL FORCE             0.174       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.475       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.923E-03   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            -2.61       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.745E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.413E-03   TIME AVG. FORCE       2.964E-03
 LARGEST RESIDUAL FORCE             1.334E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.475       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.940E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.806E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     8.161E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.413E-03   TIME AVG. FORCE       2.964E-03
 LARGEST RESIDUAL FORCE            -1.871E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.474       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.460E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.746E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.154E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  5.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.413E-03   TIME AVG. FORCE       2.964E-03
 LARGEST RESIDUAL FORCE            -5.810E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.473       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.075E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            3.167E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.852E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.413E-03   TIME AVG. FORCE       2.964E-03
 LARGEST RESIDUAL FORCE            -7.606E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.473       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.305E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            3.552E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     4.339E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      6.413E-03   TIME AVG. FORCE       2.964E-03
 LARGEST RESIDUAL FORCE            -2.356E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.473       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.224E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -6.325E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.757E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      6.413E-03   TIME AVG. FORCE       2.964E-03
 LARGEST RESIDUAL FORCE            -3.086E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.473       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.627E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -7.133E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -8.740E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      6.413E-03   TIME AVG. FORCE       2.964E-03
 LARGEST RESIDUAL FORCE            -9.546E-10   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.473       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.464E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.276E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.562E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.882    
 STEP TIME COMPLETED        8.82    ,  TOTAL TIME COMPLETED         8.82    


  INCREMENT    22 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.690E-03   TIME AVG. FORCE       3.134E-03
 LARGEST RESIDUAL FORCE             0.174       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.480       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.821E-03   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            -2.67       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.877E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.680E-03   TIME AVG. FORCE       3.133E-03
 LARGEST RESIDUAL FORCE             1.152E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.480       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.286E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.208E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     9.131E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.680E-03   TIME AVG. FORCE       3.133E-03
 LARGEST RESIDUAL FORCE            -2.130E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.479       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.902E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.936E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.284E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.680E-03   TIME AVG. FORCE       3.133E-03
 LARGEST RESIDUAL FORCE            -7.714E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.478       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.942E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            3.782E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.386E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.680E-03   TIME AVG. FORCE       3.133E-03
 LARGEST RESIDUAL FORCE            -1.031E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.479       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.527E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.330E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     5.045E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      6.680E-03   TIME AVG. FORCE       3.133E-03
 LARGEST RESIDUAL FORCE            -3.762E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.479       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.533E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -8.275E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     9.686E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      6.680E-03   TIME AVG. FORCE       3.133E-03
 LARGEST RESIDUAL FORCE            -5.029E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.479       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.373E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -9.529E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.114E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      6.680E-03   TIME AVG. FORCE       3.133E-03
 LARGEST RESIDUAL FORCE            -1.834E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.479       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.385E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            1.831E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.139E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.901    
 STEP TIME COMPLETED        9.01    ,  TOTAL TIME COMPLETED         9.01    


  INCREMENT    23 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.967E-03   TIME AVG. FORCE       3.300E-03
 LARGEST RESIDUAL FORCE             0.173       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.485       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.759E-03   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            -2.74       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.912E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      6.956E-03   TIME AVG. FORCE       3.300E-03
 LARGEST RESIDUAL FORCE             1.273E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.485       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.051E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.146E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     8.918E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      6.956E-03   TIME AVG. FORCE       3.300E-03
 LARGEST RESIDUAL FORCE            -2.368E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.484       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.217E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -2.126E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.396E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      6.956E-03   TIME AVG. FORCE       3.300E-03
 LARGEST RESIDUAL FORCE            -8.328E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.483       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.338E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            4.219E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.667E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      6.956E-03   TIME AVG. FORCE       3.300E-03
 LARGEST RESIDUAL FORCE            -1.374E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.483       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.745E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            5.205E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     5.786E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  5.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      6.956E-03   TIME AVG. FORCE       3.300E-03
 LARGEST RESIDUAL FORCE            -4.858E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.484       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.772E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.009E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.127E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      6.956E-03   TIME AVG. FORCE       3.300E-03
 LARGEST RESIDUAL FORCE            -8.016E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.484       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        4.214E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT           -1.252E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.398E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      6.956E-03   TIME AVG. FORCE       3.300E-03
 LARGEST RESIDUAL FORCE            -2.833E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.483       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.279E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.000E-02
 LARGEST RESIDUAL MOMENT            2.442E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.722E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.920    
 STEP TIME COMPLETED        9.20    ,  TOTAL TIME COMPLETED         9.20    


  INCREMENT    24 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.253E-03   TIME AVG. FORCE       3.464E-03
 LARGEST RESIDUAL FORCE             0.173       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.490       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.709E-03   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.014E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            -2.80       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.957E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.243E-03   TIME AVG. FORCE       3.464E-03
 LARGEST RESIDUAL FORCE             1.257E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.490       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.090E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.014E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.271E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     9.124E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      7.243E-03   TIME AVG. FORCE       3.464E-03
 LARGEST RESIDUAL FORCE            -2.628E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.488       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.583E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.014E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.338E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.516E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.242E-03   TIME AVG. FORCE       3.464E-03
 LARGEST RESIDUAL FORCE            -9.721E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.488       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.969E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.014E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            4.816E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.083E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.242E-03   TIME AVG. FORCE       3.464E-03
 LARGEST RESIDUAL FORCE            -1.833E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.488       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.006E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.014E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            6.270E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.651E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  5.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      7.242E-03   TIME AVG. FORCE       3.464E-03
 LARGEST RESIDUAL FORCE            -6.793E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.488       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.106E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.014E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -1.259E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.344E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      7.242E-03   TIME AVG. FORCE       3.464E-03
 LARGEST RESIDUAL FORCE            -1.281E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.488       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.304E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.014E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -1.650E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -1.758E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      7.242E-03   TIME AVG. FORCE       3.464E-03
 LARGEST RESIDUAL FORCE            -4.745E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.488       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.567E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.014E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            3.338E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.551E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.939    
 STEP TIME COMPLETED        9.39    ,  TOTAL TIME COMPLETED         9.39    


  INCREMENT    25 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.549E-03   TIME AVG. FORCE       3.627E-03
 LARGEST RESIDUAL FORCE             0.172       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.495       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.627E-03   AT NODE         18   DOF  2
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.064E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            -2.86       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.974E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.539E-03   TIME AVG. FORCE       3.627E-03
 LARGEST RESIDUAL FORCE             1.272E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.495       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.842E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.064E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.440E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     8.804E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      7.539E-03   TIME AVG. FORCE       3.627E-03
 LARGEST RESIDUAL FORCE            -2.942E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.493       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        7.925E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.064E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.567E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.639E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.539E-03   TIME AVG. FORCE       3.627E-03
 LARGEST RESIDUAL FORCE            -1.050E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.492       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.375E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.064E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            5.346E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.384E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.539E-03   TIME AVG. FORCE       3.627E-03
 LARGEST RESIDUAL FORCE            -2.458E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.493       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.298E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.064E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            7.553E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.646E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      7.539E-03   TIME AVG. FORCE       3.627E-03
 LARGEST RESIDUAL FORCE            -8.816E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.493       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.426E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.064E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -1.531E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.560E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      7.539E-03   TIME AVG. FORCE       3.627E-03
 LARGEST RESIDUAL FORCE            -2.064E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.493       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        6.660E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.064E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.177E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.216E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      7.539E-03   TIME AVG. FORCE       3.627E-03
 LARGEST RESIDUAL FORCE            -7.398E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.493       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.029E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.064E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            4.450E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -4.519E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.190    ,  FRACTION OF STEP COMPLETED  0.958    
 STEP TIME COMPLETED        9.58    ,  TOTAL TIME COMPLETED         9.58    


  INCREMENT    26 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.190    

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.857E-03   TIME AVG. FORCE       3.790E-03
 LARGEST RESIDUAL FORCE             0.172       AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.499       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.896E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            -2.92       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.184E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.789E-03
 LARGEST RESIDUAL FORCE             1.232E-04   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.499       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.245E-03   AT NODE        127   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.384E-03   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.011E-04   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.789E-03
 LARGEST RESIDUAL FORCE            -3.483E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.497       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.628E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.906E-03   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.855E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.789E-03
 LARGEST RESIDUAL FORCE            -1.440E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.497       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.767E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            6.505E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.249E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.846E-03   TIME AVG. FORCE       3.789E-03
 LARGEST RESIDUAL FORCE            -3.482E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.497       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.748E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            9.404E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     9.083E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.789E-03
 LARGEST RESIDUAL FORCE            -1.458E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.498       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        3.107E-04   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.042E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.988E-05   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.789E-03
 LARGEST RESIDUAL FORCE            -3.527E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.497       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        8.744E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.974E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.891E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.789E-03
 LARGEST RESIDUAL FORCE            -1.475E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.497       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -9.886E-05   AT NODE          1   DOF  2
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            6.525E-05   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.326E-06   AT NODE          1   DOF  6
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.789E-03
 LARGEST RESIDUAL FORCE            -3.572E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.497       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.783E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            9.484E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      1.898E-02   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     9.199E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO
 

 ***NOTE: THE RATE OF CONVERGENCE IS VERY SLOW. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT    26 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  9.492E-02

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.695E-03   TIME AVG. FORCE       3.783E-03
 LARGEST RESIDUAL FORCE             6.455E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.249       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.694E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            -1.09       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.158E-04   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.692E-03   TIME AVG. FORCE       3.783E-03
 LARGEST RESIDUAL FORCE             2.081E-05   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.249       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.808E-03   AT NODE         89   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -4.325E-04   AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.309E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      7.692E-03   TIME AVG. FORCE       3.783E-03
 LARGEST RESIDUAL FORCE            -3.380E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.248       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.327E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -8.499E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -8.463E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.692E-03   TIME AVG. FORCE       3.783E-03
 LARGEST RESIDUAL FORCE            -1.213E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.248       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.578E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            1.659E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.647E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.692E-03   TIME AVG. FORCE       3.783E-03
 LARGEST RESIDUAL FORCE            -3.124E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.248       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.069E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            2.590E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.569E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      7.692E-03   TIME AVG. FORCE       3.783E-03
 LARGEST RESIDUAL FORCE            -1.118E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.248       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.301E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.114E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -5.034E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.000E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  9.492E-02,  FRACTION OF STEP COMPLETED  0.968    
 STEP TIME COMPLETED        9.68    ,  TOTAL TIME COMPLETED         9.68    


  INCREMENT    27 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.492E-02

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      7.849E-03   TIME AVG. FORCE       3.934E-03
 LARGEST RESIDUAL FORCE             4.298E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.250       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -2.915E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.162E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -0.740       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     6.926E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.934E-03
 LARGEST RESIDUAL FORCE             9.206E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.249       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.314E-03   AT NODE         89   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.162E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -5.282E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.114E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.934E-03
 LARGEST RESIDUAL FORCE            -1.395E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.249       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.363E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.162E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -5.378E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.225E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.934E-03
 LARGEST RESIDUAL FORCE            -1.195E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.249       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -7.061E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.162E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            1.871E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.813E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.934E-03
 LARGEST RESIDUAL FORCE            -1.419E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.249       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -4.340E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.162E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            1.714E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.663E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.934E-03
 LARGEST RESIDUAL FORCE            -1.211E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.249       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.248E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.162E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -5.946E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.773E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      7.847E-03   TIME AVG. FORCE       3.934E-03
 LARGEST RESIDUAL FORCE            -1.437E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.249       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.381E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.162E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -5.452E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.292E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   7 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  7 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  9.492E-02,  FRACTION OF STEP COMPLETED  0.977    
 STEP TIME COMPLETED        9.77    ,  TOTAL TIME COMPLETED         9.77    


  INCREMENT    28 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.492E-02

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.007E-03   TIME AVG. FORCE       4.079E-03
 LARGEST RESIDUAL FORCE             4.289E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.251       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.576E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.208E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -0.744       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     9.097E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.005E-03   TIME AVG. FORCE       4.079E-03
 LARGEST RESIDUAL FORCE             8.694E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.251       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.192E-03   AT NODE         89   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.208E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -4.492E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.298E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.005E-03   TIME AVG. FORCE       4.079E-03
 LARGEST RESIDUAL FORCE            -1.840E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.250       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.654E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.208E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -6.434E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.113E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      8.005E-03   TIME AVG. FORCE       4.079E-03
 LARGEST RESIDUAL FORCE            -1.029E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.250       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -6.852E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.208E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            1.728E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.637E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      8.005E-03   TIME AVG. FORCE       4.079E-03
 LARGEST RESIDUAL FORCE            -2.040E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.250       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.503E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.208E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            2.146E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.033E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      8.005E-03   TIME AVG. FORCE       4.079E-03
 LARGEST RESIDUAL FORCE            -1.139E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.250       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.280E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.208E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -5.741E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.446E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  9.492E-02,  FRACTION OF STEP COMPLETED  0.987    
 STEP TIME COMPLETED        9.87    ,  TOTAL TIME COMPLETED         9.87    


  INCREMENT    29 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  9.492E-02

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.167E-03   TIME AVG. FORCE       4.220E-03
 LARGEST RESIDUAL FORCE             4.283E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.252       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.022E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.253E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -0.756       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.219E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.165E-03   TIME AVG. FORCE       4.220E-03
 LARGEST RESIDUAL FORCE             8.807E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.252       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.444E-03   AT NODE         89   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.253E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -5.946E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.583E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.165E-03   TIME AVG. FORCE       4.220E-03
 LARGEST RESIDUAL FORCE            -1.672E-06   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.251       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.481E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.253E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -6.103E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -5.663E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      8.165E-03   TIME AVG. FORCE       4.220E-03
 LARGEST RESIDUAL FORCE            -1.678E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.251       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -8.274E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.253E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            2.329E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.155E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      8.165E-03   TIME AVG. FORCE       4.220E-03
 LARGEST RESIDUAL FORCE            -2.037E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.251       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -5.158E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.253E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            2.129E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     1.972E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      8.165E-03   TIME AVG. FORCE       4.220E-03
 LARGEST RESIDUAL FORCE            -2.037E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.251       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        2.882E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.253E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -8.098E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -7.509E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      8.165E-03   TIME AVG. FORCE       4.220E-03
 LARGEST RESIDUAL FORCE            -2.472E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.251       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.797E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.253E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -7.412E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -6.871E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      8.165E-03   TIME AVG. FORCE       4.220E-03
 LARGEST RESIDUAL FORCE            -2.472E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.251       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.004E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.253E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            2.823E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      9.492E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     2.616E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   8 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  8 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  9.492E-02,  FRACTION OF STEP COMPLETED  0.996    
 STEP TIME COMPLETED        9.96    ,  TOTAL TIME COMPLETED         9.96    


  INCREMENT    30 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  3.828E-02

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      8.232E-03   TIME AVG. FORCE       4.354E-03
 LARGEST RESIDUAL FORCE             1.211E-02   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.102       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.266E-03   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.296E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -0.214       AT NODE          1   DOF  5
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.828E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     3.512E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  3.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      8.231E-03   TIME AVG. FORCE       4.354E-03
 LARGEST RESIDUAL FORCE             5.611E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.102       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        9.104E-04   AT NODE         89   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.296E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.665E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.828E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.559E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      8.231E-03   TIME AVG. FORCE       4.354E-03
 LARGEST RESIDUAL FORCE            -2.188E-07   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.102       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        5.150E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.296E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -2.189E-04   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.828E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -2.011E-05   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  3.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      8.231E-03   TIME AVG. FORCE       4.354E-03
 LARGEST RESIDUAL FORCE            -2.897E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.102       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -3.332E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.296E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            9.771E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.828E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     8.966E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =  1.00000E-02
       ELAPSED TOTAL CPU TIME (SEC) =  5.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      8.231E-03   TIME AVG. FORCE       4.354E-03
 LARGEST RESIDUAL FORCE            -2.758E-08   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.102       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.       -1.828E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.296E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT            7.771E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.828E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION     7.134E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          THE MOMENT    IS ZERO EVERYWHERE BUT THE MOMENT    RESIDUAL OR
          THE ROTATION CORRECTION IS NON-ZERO

  NUMBER OF EQUATIONS =      3678     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  1.39E+08

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  4.00000E-02
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  4.00000E-02
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      8.231E-03   TIME AVG. FORCE       4.354E-03
 LARGEST RESIDUAL FORCE            -3.647E-09   AT NODE          1   DOF  1
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF DISP.         0.102       AT NODE         72   DOF  2
   INSTANCE: PART-1-1                                                                        
 LARGEST CORRECTION TO DISP.        1.183E-04   AT NODE          1   DOF  3
   INSTANCE: PART-1-1                                                                        
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 AVERAGE MOMENT                      0.00       TIME AVG. MOMENT      1.296E-02
                                                (FROM FORCE     FIELD)
 LARGEST RESIDUAL MOMENT           -3.464E-05   AT NODE          1   DOF  6
   INSTANCE:                                                                                 
 LARGEST INCREMENT OF ROTATION      3.828E-03   AT NODE          1   DOF  4
   INSTANCE: {assembly node}                                                                 
 LARGEST CORRECTION TO ROTATION    -3.182E-06   AT NODE          1   DOF  5
   INSTANCE: {assembly node}                                                                 
          ALL ROTATION CORRECTIONS ARE ZERO

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  3.828E-02,  FRACTION OF STEP COMPLETED   1.00    
 STEP TIME COMPLETED        10.0    ,  TOTAL TIME COMPLETED         10.0    


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF         30  INCREMENTS
                       3  CUTBACKS IN AUTOMATIC INCREMENTATION
                     219  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     219  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                     219  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       1  WARNING MESSAGES DURING USER INPUT PROCESSING
                      10  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                      10  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   39.610    
       SYSTEM TIME (SEC)    =   2.0800    
       TOTAL CPU TIME (SEC) =   41.690    
       WALLCLOCK TIME (SEC) =         43
