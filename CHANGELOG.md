# Changelog

## 0.1.10
- Fix a compilation bug under Windows when using the new Azure servers

## 0.1.9

- Added `y2f_version` function that returns currently installed version of Y2F

## 0.1.8

- Fixed compilation for Visual Studio 2015 users

## 0.1.7

- Fixed bug in MEX file compilation on Windows machines

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
