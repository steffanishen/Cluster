read_dcd
========

[![Build Status](https://travis-ci.org/FHedin/read_dcd.svg?branch=master)](https://travis-ci.org/FHedin/read_dcd)

<a href="https://scan.coverity.com/projects/4505">
  <img alt="Coverity Scan Build Status"
         src="https://scan.coverity.com/projects/4505/badge.svg"/>
</a>

c++ class + main file example for reading a charmm dcd

Notes:
1. How does the code function:
   a) random_walk: create the placeholder psf and pdb, read monomer info from monomers, create a box and perform random walk to create polymer chains
   b) rule of thumb for monomers: 1) don't include more than 2 backbone atoms in each monomer;
				  2) the 1st and last monomers contain 2 backbone atoms

2. Notes:
   a) segid has to be an integer
How to use:
read_dcd input 
# Cluster
