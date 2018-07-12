# Project Title

A function that computes gradient and the Hessian of the distance between point and triangle in 3D

## Getting Started

The input parameter x[12] contains coordinates of points p0, p1, p2, p3, where (p1,p2,p3) define the triange and p0 is the point to which the distance is computed.
The return value is the squared distance.
Output array fd[12] contains first derivatives of the squared distance with respect to coordinates of p0, p1, p2, p3.
Similarly, sd[12][12] contains second derivatives.
Zeta2 and zeta3 are barycentric coordinates of the closest point on the triangle.

```
double pt(double(&x)[12], double(&fd)[12], double(&sd)[12][12], double &zeta2, double &zeta3);
```


## Please cite as

Gribanov, Igor; Taylor, Rocky; Sarracino, Robert. 2018. "The Gradient and the Hessian of the Distance between Point and Triangle in 3D." Algorithms 11, no. 7: 104.



