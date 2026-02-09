#!/bin/sh

polynome2pi scan --sector "001|u_d" 
polynome2pi scan --sector "011|u_d_s" 
polynome2pi scan --sector "112|nucleon" 
polynome2pi scan --sector "H-Atom" 
polynome2pi scan --sector "222|c_tau" 

polynome2pi sensitivity --sector "011|u_d_s" 

polynome2pi sensitivity --sector "011|u_d_s" --particle u
polynome2pi sensitivity --sector "011|u_d_s" --particle d