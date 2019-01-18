# Multivariable Optimization
#### Function used to demonstrate all the methods:
![alt text](https://i.imgur.com/2ZH0uA7.png "Exemplary function used for all methods.")
#### Additional assumptions:
![alt text](https://i.imgur.com/RMdclwj.png "Assumed interval for x, y")

![alt text](https://i.imgur.com/WTVDHlU.png "Assumed initial point")



## Gauss-Seidel
_**GaussSeidel**_ is a MATLAB implementation of Gauss-Seidel algorithm for finding the minimum of unimodal (within given interval) multivariable function.

This particular implementation assumes the function to be of two variables.
After providing (hardcode) the program with function formula, constraints for both variables, the initial point and the precision of searching, the program produces the following output.
#### Exemplary output: 

```
/////////////////////////////////////////////////////////////////////
//                                                                 //
//    OaDM Laboratory                                17.01.2019    //
//                                                                 //
//                   Multivariable Optimization                    //
//                      Gauss-Seidel method                        //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  Initial point x0: (0, 0)

  Iteration [1]
      a1 = 0.818
      Moving along X axis by [0.818].
      x1 = (0.818, 0)
      F(x1) = 0.112586

  Iteration [2]
      a2 = 0.2
      Moving along Y axis by [0.2].
      x2 = (0.818, 0.2)
      F(x2) = 0.0736213

  Iteration [3]
      a3 = 0
      Moving along X axis by [0].
      x3 = (0.818, 0.2)
      F(x3) = 0.0736213

  Iteration [4]
      a4 = 0
      Moving along Y axis by [0].
      x4 = (0.818, 0.2)
      F(x4) = 0.0736213


##############################################################

    Displacement in both directions is equal to 0. Algorithm stops here.
      Minimum found at (0.818, 0.2) with value: 0.0736213
```


## Powell
Similarly to the Gauss-Seidel algorithm, _**Powell**_ method was implemented here for minimization a function of two variables having one local minimum in the given interval.
As previously, after function formula, constraints for both variables, the initial point and the precision of searching, the program produces the following output.
#### Exemplary output:
```
/////////////////////////////////////////////////////////////////////
//                                                                 //
//    OaDM Laboratory                                17.01.2019    //
//                                                                 //
//                   Multivariable Optimization                    //
//                          Powell method                          //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  Initial point x0: (0, 0)

  Iteration [1]
    Gauss-Seidel part:
        Moving along X axis by [0.82].
        x1 = (0.82, 0)
        Moving along Y axis by [0.2].
        x2 = (0.82, 0.2)

    Powell part:
        Moving in direction [0.82, 0.2] by [2.22045e-16].
        x3 = (0.82, 0.2)
        F(x) = 0.0736264

  Iteration [2]
        Moving in direction [-0.82, 0] by [-2.22045e-16].
        x4 = (0.82, 0.2)
        F(x) = 0.0736264

  Iteration [3]
        Moving in direction [-0.82, -0.2] by [-2.22045e-16].
        x5 = (0.82, 0.2)
        F(x) = 0.0736264

  Iteration [4]
        Moving in direction [0, 0] by [-1.2].
        x6 = (0.82, 0.2)
        F(x) = 0.0736264


##############################################################

    Displacement in chosen direction is equal to 0. Algorithm stops here.
      Minimum found at (0.82, 0.2) with value: 0.0736264
```

## Steepest Descent
_**SteepestDescent**_ is yet another implementation of minimizing algorithms. This method uses gradient of the function to choose new direction in which the current point will be moved in the iteration.
#### Exemplary output:
```
/////////////////////////////////////////////////////////////////////
//                                                                 //
//    OaDM Laboratory                                17.01.2019    //
//                                                                 //
//                   Multivariable Optimization                    //
//                    Steepest Descent Method                      //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  Initial point x0: (0, 0)

  Iteration [1]
        Moving in direction [-0.973848, -0.379607] by [0.8].
        x1 = (0.779078, 0.303685)
        F(x) = 0.0859215

  Iteration [2]
        Moving in direction [-0.0836775, 0.204435] by [0.496315].
        x2 = (0.820608, 0.202222)
        F(x) = 0.0736345

  Iteration [3]
        Moving in direction [0.00600503, 0.00444303] by [0.477778].
        x3 = (0.817739, 0.200099)
        F(x) = 0.0736213

  Iteration [4]
        Moving in direction [-0.000232115, 0.000197491] by [0.479901].
        x4 = (0.817851, 0.200004)
        F(x) = 0.0736213

  Iteration [5]
        Moving in direction [9.94381e-06, 7.93865e-06] by [0.469996].
        x5 = (0.817846, 0.2)
        F(x) = 0.0736213

  Iteration [6]
        Moving in direction [-2.12122e-07, 4.76382e-07] by [0.49].
        x6 = (0.817846, 0.2)
        F(x) = 0.0736213

  Iteration [7]
        Moving in direction [1.37456e-08, 9.52787e-09] by [0.28].
        x7 = (0.817846, 0.2)
        F(x) = 0.0736213

  Iteration [8]
        Moving in direction [5.38199e-09, 4.19226e-09] by [-2.09613e-09].
        x8 = (0.817846, 0.2)
        F(x) = 0.0736213

  Iteration [9]
        Moving in direction [5.38199e-09, 4.19226e-09] by [-2.09613e-09].
        x9 = (0.817846, 0.2)
        F(x) = 0.0736213

  Iteration [10]
        Moving in direction [5.38199e-09, 4.19226e-09] by [-2.09613e-09].
        x10 = (0.817846, 0.2)
        F(x) = 0.0736213

  Iteration [11]
        Moving in direction [5.38199e-09, 4.19226e-09] by [-2.09613e-09].
        x11 = (0.817846, 0.2)
        F(x) = 0.0736213


##############################################################

    Stopping the calculations because of too little differences in F(x).
      Minimum found at (0.817846, 0.2) with value: 0.0736213
```

## Conjugate Gradient
_**ConjugateGradient**_ is the last of the four algorithms in this repository. The difference between the Conjugate Gradient method and the Steepest Descent method is that this algorithm instead of calculating the gradient with each iteration as new direction, "improves" the previous gradient. Details to be found in the _**LAB_INSTURCTION**_.
#### Exemplary output:
```
/////////////////////////////////////////////////////////////////////
//                                                                 //
//    OaDM Laboratory                                17.01.2019    //
//                                                                 //
//                   Multivariable Optimization                    //
//                   Conjugate Gradient Method                     //
//                                                                 //
/////////////////////////////////////////////////////////////////////

  Function: atan(y - 1/5)^2 + sin(x - 9/10)^2 + x^2/10
 
  Initial point x0: (0, 0)


  Iteration [1]
        Moving in direction [0.973848, 0.379607] by [0.8].
        x1 = (0.779078, 0.303685)
        F(x) = 0.0859215

  Iteration [2]
        Moving in direction [0.123712, -0.188829] by [0.476315].
        x2 = (0.838004, 0.213743)
        F(x) = 0.0742525

  Iteration [3]
        Moving in direction [-0.0420436, -0.0303527] by [0.466257].
        x3 = (0.818401, 0.199591)
        F(x) = 0.0736218

  Iteration [4]
        Moving in direction [-0.000761192, 0.00113814] by [0.480409].
        x4 = (0.818035, 0.200138)
        F(x) = 0.0736214

  Iteration [5]
        Moving in direction [-0.00040152, -0.000288948] by [0.469862].
        x5 = (0.817846, 0.200002)
        F(x) = 0.0736213

  Iteration [6]
        Moving in direction [1.81523e-06, -2.76712e-06] by [0.499998].
        x6 = (0.817847, 0.200001)
        F(x) = 0.0736213

  Iteration [7]
        Moving in direction [-2.25934e-06, -1.57192e-06] by [0.469999].
        x7 = (0.817846, 0.2)
        F(x) = 0.0736213

  Iteration [8]
        Moving in direction [3.07247e-08, -4.80733e-08] by [0.43].
        x8 = (0.817846, 0.2)
        F(x) = 0.0736213

  Iteration [9]
        Moving in direction [-9.8332e-09, -8.13588e-09] by [0.39].
        x9 = (0.817846, 0.2)
        F(x) = 0.0736213

  Iteration [10]
        Moving in direction [2.85602e-09, -5.42201e-09] by [-3.26148e-09].
        x10 = (0.817846, 0.2)
        F(x) = 0.0736213

  Iteration [11]
        Moving in direction [1.52539e-09, -6.52296e-09] by [-3.26148e-09].
        x11 = (0.817846, 0.2)
        F(x) = 0.0736213


##############################################################

    Stopping the calculations because of too little differences in F(x).
      Minimum found at (0.817846, 0.2) with value: 0.0736213
```


## Contour
_**Contour**_ plots the contour graph of the function values for (x, y).
#### Exemplary plot
![alt text](https://i.imgur.com/eCiTdEo.png "Exemplary plot in Task B")
