# Compartment Model for the Study of Insulin Circulation in the Blood

## Project Overview
This project aims to develop a mathematical model for studying the kinetics of insulin circulation in the blood. 
Diabetic patients often struggle to maintain stable blood sugar levels due to various factors affecting insulin dynamics. 
By creating a predictive model, we aim to assist in achieving more predictable blood sugar levels over time, improving diabetes management.

## Course and Degree
- Course: Numerical Analysis Part A
- Degree: MSc in Biomedical Engineering at the University of Bologna (UNIBO)

## Mathematical Model Description
The mathematical model consists of a system of differential equations representing the concentration of insulin in different compartments of the body (blood, kidneys, pancreas, and abdomen).
Various flow rate constants determine the movement of insulin between compartments and its removal by organs.
The model focuses solely on insulin dynamics and excludes factors like glucose for analytical simplicity.

## Solution Methods
Mpdel is solved in Matlab using following methods to assess performance.
The flux equations are solved using the following methods:
1. ODE45
2. Euler Explicit Method
3. Euler Implicit Method
4. Heun's Method
