#!/bin/sh

polynome2pi sensitivity --sector minimal --eps 0.001
polynome2pi sensitivity --sector minimal --eps 0.01
polynome2pi sensitivity --sector minimal --eps 0.1

polynome2pi sensitivity --sector minimal --eps 0.001 --particle d
polynome2pi sensitivity --sector minimal --eps 0.01 --particle d
polynome2pi sensitivity --sector minimal --eps 0.1 --particle d

polynome2pi sensitivity --sector minimal --eps 0.001 --particle u
polynome2pi sensitivity --sector minimal --eps 0.01 --particle u
polynome2pi sensitivity --sector minimal --eps 0.1 --particle u

polynome2pi sensitivity --sector minimal --eps 0.001 --particle electron
polynome2pi sensitivity --sector minimal --eps 0.01 --particle electron
polynome2pi sensitivity --sector minimal --eps 0.1 --particle electron

polynome2pi sensitivity --sector light --eps 0.001
polynome2pi sensitivity --sector light --eps 0.01
polynome2pi sensitivity --sector light --eps 0.1

polynome2pi sensitivity --sector light --eps 0.001 --particle d
polynome2pi sensitivity --sector light --eps 0.01 --particle d
polynome2pi sensitivity --sector light --eps 0.1 --particle d

polynome2pi sensitivity --sector light --eps 0.001 --particle u
polynome2pi sensitivity --sector light --eps 0.01 --particle u
polynome2pi sensitivity --sector light --eps 0.1 --particle u

polynome2pi sensitivity --sector light --eps 0.001 --particle electron
polynome2pi sensitivity --sector light --eps 0.01 --particle electron
polynome2pi sensitivity --sector light --eps 0.1 --particle electron

polynome2pi sensitivity --sector light --eps 0.001 --particle s
polynome2pi sensitivity --sector light --eps 0.01 --particle s
polynome2pi sensitivity --sector light --eps 0.1 --particle s