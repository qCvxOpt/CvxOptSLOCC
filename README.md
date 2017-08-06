# CvxOptSLOCC

This collection of *MATLAB* files aims to prove separability or membership in a certain SLOCC entanglement class for a given multiparticle quantum state, see `test_SEP_rhoPH.m` for an example. It is also able to perform the likelihood-ratio test, see `test_LRT_rhoPH.m` for a demonstration.

 * `test_SEP_rhoPH` Test for separability, using the 3x3 bound entangled states by P. Horodecki.
 * `test_LRT_rhoPH.m` Demonstration of the likelihood-ratio test (LRT), using the 3x3 bound entangled states by P. Horodecki.
 * `apg_Gilbert.m` The apg_Gilbert algorithm.
 * `dg_Gilbert.m` The dg_Gilbert algorithm.
 * `Gilbert_SLOCC.m` The Gilbert_SLOCC algorithm.
 * `qmt.m` Quantum measurement transform.
 * `WSLOCCMaxG.m` W-state SLOCC maximization (n-parties).
 * `RandomLocalUnitary.m` Random local unitary (uniform with respect to the Haar measure).
 * `ShuffleInit.m` Permutation of matrix generation (AB....YZ --> ZAB....Y).
 
 
 
Reference
----
If you are using the source codes for your research, please consider citing:
 * J. Shang and O. Guehne, Convex optimization over classes of multiparticle entanglement, *arXiv:1707.02958 [quant-ph] (2017)*.


License
----

MIT License

Copyright (c) 2017 Jiangwei Shang and Otfried Guehne

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
