oKDE-Java
=========

This is a Java implementation of the oKDE algorithm proposed by Matej Kristan
([oKDE](http://www.vicos.si/Research/Multivariate_Online_Kernel_Density_Estimation)).
This project is part of my master's thesis that uses the oKDE algorithm to estimate
a model of human mobility and exploits the estimated model to predict future locations
of human individuals (*Location Prediction Based on Mobility Patterns in Location Histories*:
[git repository](https://github.com/joluet/prepos)).

*******************************

1. [oKDE Algorithm](#okde)
2. [Build Instructions](#build)
3. [Quickstart](#start)
4. [External Libraries Used](#ext_libs)
5. [License](#license)

*******************************


<a name="okde">
## oKDE Algorithm

Given a set of n-dimensional samples, this algorithm estimates 
the distribution of the samples using kernel density estimation.
The output model is a *Mixture of Gaussians*.
Moreover, the model is compressed and can be updated as new samples arrive.

The basic principle of the oKDE algorithm is summarized by this graphic:

![oKDE algorithm](oKDE.png "")

For more details see ([oKDE](http://www.vicos.si/Research/Multivariate_Online_Kernel_Density_Estimation)).



<a name="build">
## Build Instructions

Just execute ant in project root to compile the project:

**$ ant**

Afterwards, the packed jar file can be found in the *dist*-folder.



<a name="start">
## Quickstart

[Here](https://github.com/joluet/okde-java-example) you can find a simple example that uses oKDE-Java
to estimate a distribution of randomly generated samples. This example illustrates the basic usage
of oKDE-java.

To use oKDE-Java in another project just include the jar file (see above [how to build](#build)).



<a name="ext_libs">
## External Libraries Used

The following libraries are used in oKDE-Java:
 *  [EJML](https://code.google.com/p/efficient-java-matrix-library/),
 	a linear algebra library for manipulating dense matrices
 	License: http://www.apache.org/licenses/LICENSE-2.0


<a name="license">
## License

The MIT License (MIT)

Copyright (c) 2014 Jonas Luthke

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
