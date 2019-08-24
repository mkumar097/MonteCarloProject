# MolSSI Summer School 2019
A project that was completed or the 2019 MolSSI Software Summer School
 which was held from July 22-31 at The Texas Advanced Computing Center. The 
 School’s focus was on best practices in software engineering, version
 control, continuous integration, data management, programming paradigms.
 
 FYI: This is built off of an algorithm that was provided by the MolSSI
 staff. I have coded and organized the algorithm using classes to help
  demonstrate my capabilities with object-oriented programming (OOP). I will
   be updating this small project over time to demonstrate what I have learned 
    about Python throughout my PhD.

<p align="center">
  <img src="./images/MM2_logo.png")>
</p>

## Molecular Mechanics Project 
[![Build Status](https://travis-ci.org/MolSSI-Education/mm_2019_sss_2.svg?branch=master)](https://travis-ci.org/MolSSI-Education/mm_2019_sss_2)
[![codecov](https://codecov.io/gh/MolSSI-Education/mm_2019_sss_2/branch/master/graph/badge.svg)](https://codecov.io/gh/jiayeguo/mm_2019_sss_2)

This code is written to run a simple Monte Carlo simulation of an ideal gas. 
The code works by either setting up a random box of ideal gas atoms or reading
in a file. (Format for the file is provided in the examples folder.)

The code is broken into three different classes:

1. SystemSetup
2. EnergyFunctions
3. MonteCarlo (This inherits attributes from the EnergyFunctions class.)