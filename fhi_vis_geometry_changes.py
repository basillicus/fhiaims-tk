#!/usr/bin/env python3
# #!python
"""
Joins the initial and final geometry files of a geometry optimization and generates an extxyz file
"""

from ase.io.aims import read_aims, write_aims
from ase.io import read, write

outfile = 'geom_optim_changes.extxyz'
init_geom = read_aims('geometry.in')
final_geom = read_aims('geometry.in.next_step')

write(outfile, init_geom)
write(outfile, final_geom, append=True)
