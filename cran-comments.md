## Resubmission

This is a resubmission. In this version I have:

* Added small executable examples for Rd-files
* Removed code that previously modified user's home filespace (from tests) - code now writes to tempdir() instead when necessary  
* Maintainer asked for potential references describing methods in package, but unfortunately the package doesn't have any heavy statistical methodology that warrants a stand-alone paper/citation. The main rationale for the package is simplification of code by wrapping typical workflows in a series of functions. There is an analysis based manuscript being prepared that uses the package, but that will only include package as a method - it will not be a methodology focused manuscript       

## Test environments 

* ubuntu 18.04 (local), R 3.6.1
* ubuntu 16.04 (travis-ci), R 3.5.3
* ubuntu 16.04 (travis-ci), R 3.6.1
* Mac OS X 10.13.3 (travis-ci) R 3.5.3
* Mac OS X 10.13.3 (travis-ci) R 3.6.1
* Windows Server 2012 R2 x64 (travis-ci) R 3.6.1
* Windows 10 x64 (local) R 3.6.1 

## R CMD check results

── R CMD check results ────────────────────────────────────── LoopRig 0.1.1 ────
Duration: 1m 4.5s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

R CMD check succeeded
