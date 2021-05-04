# Unit test MC integrators, weighted/unweighted event generation etc.
#
#
# Run with: pytest ./tests/{filename} -s
#
# (c) 2017-2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os
import pathlib

from test_helpers import *


XS_UNIT = 1E6 # microbarns
y_unit  = '$\\mu b$'


def printout(data):
	"""
	Helper function
	"""
	for obs in data['h1'].keys():
		weights   = np.array(data['h1'][obs]['weights'])
		print(f"obs: {obs:15s} | sum[weights]: {np.sum(weights)} | fills = {data['h1'][obs]['fills']} | overflow = {data['h1'][obs]['overflow']} | underflow = {data['h1'][obs]['underflow']}")


def process(data, obs='M'):
	"""
	Helper function for reading in fast histogram json data
	"""

	print(data.keys())
	print(data['h1'].keys())
	print(data['h1'][obs].keys())

	obs       = 'M'
	binedges  = np.array(data['h1'][obs]['binedges'])
	weights   = np.array(data['h1'][obs]['weights'])
	fills     = data['h1'][obs]['fills']


	bincenter = np.mean(binedges,1)           # bincenter
	binwidth  = binedges[:,1] - binedges[:,0] # binwidth per bin

	# Differential cross section
	dsdx      = weights / fills / binwidth * XS_UNIT

	return bincenter, dsdx


def test_integrators():
	"""
	Test MC integrators
	"""
	
	pid = '[[211, -211], [211, -211], [211, -211]]'
	
	## Execute MC generator
	weightmode = ['true',  'false']
	integrator = ['VEGAS', 'FLAT']
	
	N = 1000
	
	for w in weightmode:
		for g in integrator:

			if g == 'FLAT' and w == 'true':
				continue; # FLAT and weighted has too high variance
			
			print(f'Generating events with weightmode <{w}> and integrator <{g}> ...')
			
			cmd = f"make -j4 && ./bin/gr -i ./input/test.json -p 'PP[CON]<C> -> pi+ pi-' -l false -h 0 -w {w} -n {N} -g {g} -o {g}_w_{w}"
			execute(cmd)
	
	cmd = f"python ./python/iceshot --hepmc3 VEGAS_w_false VEGAS_w_true FLAT_w_false --unit ub --pid '{pid}' --output testbench_integrators"
	execute(cmd, expect='[iceshot: done]')


def test_hfast():
	"""
	Test fast histograms
	"""
	
	obs = 'M'

	## Execute MC generator
	print('Computing ...')

	modes = ['VEGAS', 'FLAT']

	for mode in modes:
		cmd = f"make -j4 && ./bin/gr -i ./input/test.json -p 'PP[CON]<C> -> pi+ pi-' -w true -l false -h 1 -n 0 -o {mode} -g {mode}"
		execute(cmd)
	
	## Read fast histogram output
	fig,ax = plt.subplots()

	for mode in modes:

		with open(f'./output/{mode}.hfast') as f:
		  data = json.load(f)

		printout(data)
		bincenter, dsdx = process(data=data, obs=obs)

		plt.plot(bincenter, dsdx, label=mode)

		plt.ylabel(f'$d\\sigma/\\d {obs}$  [{y_unit}/unit]')
		plt.xlabel(f'${obs}$ [unit]')
		plt.legend()

	pathlib.Path('./testfigs/').mkdir(parents=True, exist_ok=True)
	fig.savefig(f"./testfigs/testbench_hfast.pdf", bbox_inches='tight')

