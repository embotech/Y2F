# Changelog

## 0.1.6

- Fixed bug affecting parsing of problems with parameters influencing cost


## 0.1.5

- Changed method Y2F uses to store Simulink blocks:  All blocks are now stored in a single library that can be accessed from the Simulink Library Browser. If the user re-generates a solver, the Simulink block gets updated automatically.



## 0.1.4

*skipped*


## 0.1.3

- Added basic ADMM example
- Fixed bug that occurred in Matlab versions < 2014a when generating valid names for solvers


## 0.1.2

- Fixed a glitch that caused the build to fail when no solver libs are shipped
