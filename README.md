# ZerosContour

A Python implementation of an algorithm for finding all zeros of a homolmorphic function in a given region, from the paper:

Dellnitz, Michael, Oliver Schütze, and Qinghua Zheng. "Locating all the zeros of an analytic function in one complex variable." Journal of Computational and Applied mathematics 138.2:325-333 (2002).

- The notebook `Test.ipynb` tests the code on examples of functions with known zeros.
- The notebook `SchrodEigs.ipynb` applies the code to find zeros of non-self-adjoint Schr&ouml;dinger operators.

This code was used in my paper 

Stepanenko A. "Spectral inclusion and pollution for a class of dissipative perturbations". Journal of Mathematical Physics.62(1):013501 (2021)

to create the following figure:

<img src="https://github.com/alexeistepa/ZerosContour/blob/main/fig2_dpi300.png?raw=true" width="650" height="650">

Idea of algorithm
-----------------
Given a region $\Omega = (\Re(b),\Re(b)+w)\times i (\Im(b),\Im(b)+h)$ (here $b \in \mathbb C, w,h \in \mathbb R$) in the complex plane we can calculate the number of zeros $N$ in $\Omega$ of a holomorphic function $f$ by the argument principle,
$$N =  \frac{1}{2 \pi i}\int_{\Omega} \frac{f'}{f}.$$
We split $\Omega$ into four quarters and apply the argument principle to determine wether each quarter contains at least one zero. This gives us an enclosure for the containing the zeros of $f$ in $\Omega$. 
To refine this enclosure, we take each quarter that contains at least one zero and again apply the same procedure to it, splittling it into 4 quarters and applying the argument principle. We iterate this process to obtain an enclosure consisting of a union of small boxes. We stop the iterations if either each box in the encosure only contains one zero and a minimum nimber of iterations has been reached, or if some threshold value for the number of iterations has been reached. In our implementation, we use the computer algebra package SymPy to evaluate $f'/f$, given $f$.
