#!/bin/sh

polynome2pi scan --sector "001|u_d" --no-open
polynome2pi scan --sector "011|u_d_s" --no-open
polynome2pi scan --sector "112|nucleon" --no-open
polynome2pi scan --sector "H-Atom" --no-open
polynome2pi scan --sector "222|c_tau" --no-open

polynome2pi sensitivity --sector "011|u_d_s"  --no-open

polynome2pi sensitivity --sector "011|u_d_s" --particle u --no-open
polynome2pi sensitivity --sector "011|u_d_s" --particle d --no-open