from distutils.core import setup, Extension

libs = ['S4', 'stdc++']
libs.extend( [lib[2::] for lib in '-lblas -llapack -lfftw3 -lpthread -lcholmod -lamd -lcolamd -lcamd -lccolamd '.split()])

S4module = Extension('S4',
	sources = [
		'S4/main_python.c'
	],
	libraries = libs,
	library_dirs = ['./build'],
	extra_link_args = [
		'./build/libS4.a'
	],
	extra_compile_args=['-std=gnu99'] 
)

setup(name = 'S4',
	version = '1.1',
	description = 'Stanford Stratified Structure Solver (S4): Fourier Modal Method',
	ext_modules = [S4module]
)
