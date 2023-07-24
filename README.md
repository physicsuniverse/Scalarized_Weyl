# Solutions of 4th-Order Equations in Scalarized Black Hole Systems in Higher Derivative Gravity

This repository contains the code and results corresponding to the appendix of paper titled "Scalarized Black Holes: Extending From and Beyond the Schwarzschild Solution". In this paper, we explore the scalarized black hole systems within the framework of higher derivative gravity by directly solving the fourth-order equations of motion.

## Table of Contents
- [Solutions of 4th-Order Equations in Scalarized Black Hole Systems in Higher Derivative Gravity](#solutions-of-4th-order-equations-in-scalarized-black-hole-systems-in-higher-derivative-gravity)
  - [Table of Contents](#table-of-contents)
  - [Pre-installation](#pre-installation)
    - [File Description](#file-description)
  - [Introduction](#introduction)
  - [Technical Details](#technical-details)
  - [Citation](#citation)

## Pre-installation<a name="pre-installation"></a>

Before using the code in this repository, please ensure that you have the following software and packages installed:
- [Mathematica](https://www.wolfram.com/mathematica/): We utilize the free Wolfram Kernel to run the code. You can also run the code with Mathematica notebook, just copy the code inside the .ipynb file, and things will work very much like in the .ipynb files.
- [Jupyter Notebook](https://jupyter.org/): This acts as the front end for executing the code.
- [WolframLanguageForJupyter](https://github.com/WolframResearch/WolframLanguageForJupyter): Install this package to enable the integration of the Wolfram Language with Jupyter Notebook.

### File Description

- Folders:
  - `./eoms/`: This folder contains the equations of motion (EOMs) for both second and fourth-order problems. The EOMs are mathematical representations of the motion equations used to describe the behavior of physical systems.
  - `./packages/`: The packages folder houses the necessary software packages or libraries required to calculate and compute the equations of motion (EOMs). These packages provide the necessary tools and functions to perform numerical computations.
  - `./solutions/`: Within this folder, you can find the numerical solutions obtained for both the second-order and fourth-order problems.

- Files:
  - `Deduce_EOMs.ipynb`: This file contains the process and steps involved in deducing the equations of motion for both the second and fourth-order problems. It likely includes code, explanations, and computations performed to derive the EOMs.
  - `Solutions_Check.ipynb`: This file verifies and validates the accuracy and correctness of the obtained solutions for the fourth-order equations of motion.

## Introduction<a name="introduction"></a>

The recent breakthroughs in gravitational wave detection and black hole imaging have opened up new possibilities for exploring strong gravity. Scientists are particularly interested in investigating the existence of scalar fields that could leave distinct signatures on black holes. However, the no-hair theorem, which applies to both general relativity and other theories of gravity, poses a significant challenge as it prohibits stationary black hole solutions with scalar hairs.

To bypass this limitation, researchers propose violating certain assumptions of the no-hair theorem. One approach involves introducing a scalar field coupled to the Gauss-Bonnet invariant, resulting in scalarized black holes that extend beyond the solutions predicted by general relativity. These scalarized black holes offer a unique opportunity to detect scalar fields and serve as a means of distinguishing intriguing theoretical models.

In this study, we focus on a higher-derivative gravity theory coupled to a scalar field and explore its behavior. We argue against simplifying the analysis by reducing the governing equations to second-order form, instead choosing to directly solve the challenging fourth-order differential equations using numerical methods.

## Technical Details<a name="technical-details"></a>

We are presenting numerical methods for solving second and fourth-order equations of motion. The variables $f$, $h$, and $\varphi$ are treated as functions of the variable $z$. To discretize the axis $z$, we utilize the Gauss-Lobatto quadrature. This discretization allows us to represent $f$, $h$, and $\varphi$ as vectors $\boldsymbol{W}$.

We compute derivatives up to the fourth order using standard methods such as the ``NDSolve`FiniteDifferenceDerivative`` function in Mathematica. This ensures precise calculations necessary for numerical methods. For instance, by applying the following code:
```wl
NDSolve`FiniteDifferenceDerivative[Derivative[4], collocation, "DifferenceOrder" -> "Pseudospectral"]["DifferentiationMatrix"]
```
we obtain the differentiation matrix of the fourth-order derivative within the Pseudospectral method. Here, the `colloation` should be the Gauss-Lobatto collocation.

To obtain expansion coefficients ($a_i$), we employ the Fast Fourier Transform (FFT) method. These coefficients are then used to reconstruct the function using Chebyshev polynomials $T^i(x)$.

$$F = a_i T^i(z)$$

We introduce the residual $\mathcal R^I(z)$ of the equations of motion as a measure to evaluate the accuracy of the solutions. Given the nonlinearity of the equations, we utilize the Newton-Raphson iteration method, which linearizes the equations and updates them iteratively until convergence is achieved. The success of the iteration process relies on the initial guess for the solution.

The Pseudospectral method exhibits remarkable convergence capabilities, allowing it to converge easily even with random initial guesses. To obtain our results, we utilize thousands of random initial guesses, consolidate identical solutions, and finally obtain distinct solutions.

When considering second-order equations of motion (EOM), we encounter multiple solutions that heavily rely on the initial guesses provided. It appears that an infinite number of solutions could potentially exist within this context. This suggests a lack of constraints within the EOM, as they are unable to uniquely determine the solutions. In contrast, when exploring fourth-order EOMs and subjecting them to thousands of random seeds, we discover only four branches of solutions. This implies that the fourth-order equations are significantly more adept at resolving a solution given the same boundary conditions compared to their second-order counterparts.

Moreover, to validate the solutions obtained from the fourth-order equations, we can test them by reintroducing them into both the fourth-order and second-order equations. We provide examples to demonstrate that the solutions derived from the fourth-order equations satisfy the reduced second-order equation systems. In contrast, solutions obtained from the second-order equations often fail to satisfy the fourth-order equations. This suggests that the fourth-order equations identify the correct solution space compared to the second-order reduction.

## Citation<a name="citation"></a>

If you find this work or the code in this repository useful for your research, please consider citing our paper:

```
Author(s). "Title of the Paper", arXiv:xxxxxx.
```

We appreciate your interest and hope that this repository proves helpful for your study. If you have any questions or encounter any issues, please don't hesitate to reach out to us.