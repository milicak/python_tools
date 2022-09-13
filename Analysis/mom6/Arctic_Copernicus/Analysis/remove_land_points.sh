#!/bin/bash

# After that run the code there are still blow ups you can use the code below
grep -Eo 'i=.{0,5}' Arctic_remove > Arctic_i_ind
grep -Eo 'j=.{0,5}' Arctic_remove > Arctic_j_ind
cut -c3- Arctic_i_ind > Arctic_i_ind16
cut -c3- Arctic_j_ind > Arctic_j_ind16
