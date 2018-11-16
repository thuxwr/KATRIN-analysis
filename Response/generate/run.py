# Generate response function for various z values. Arguments: Z, core#.
import sys, subprocess, multiprocessing

# core from 0 to 20, Z from -5 to 5 with width 0.5
def run(core):
	Z = core * 0.5 - 5
	subprocess.call("./Generate %10.1f %d"%(Z, core), shell=True)

p = multiprocessing.Pool(21)
for i in range(21):
	p.apply_async(run, args=(i,))
p.close()
p.join()
