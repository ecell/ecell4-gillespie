# Build script for Gillespie Solver (draft)
# vim: syntax=python 
top = '.'
out = 'build_with_waf'
APPNAME = 'gillespie_solver'
VERSION = '0.1.0'

# Header files which this module requires.
header_list = ['vector', 'map', 'numeric']

def options(opt):
	opt.add_option('--unit_test', action='store_true', default=False, help='unit test')
	opt.add_option('--enable_debug', action='store_true', default=False, help='debug')
	opt.load('compiler_cxx')

def configure(conf):
	conf.load('compiler_cxx')

	conf.check_cfg(package='gsl', uselib_store='gsl', atleat_version='1.13', args='--cflags --libs')
	conf.check_cfg(package='pficommon', uselib_store='pficommon', atleat_version='1.0.0', args='--cflags --libs')

	# Checking the existence of header files.
	for header in header_list:
		conf.check(header_name = header, features = 'c cprogram')

	# Save option flags.
	conf.env.unit_test = conf.options.unit_test
	conf.env.enable_debug =  conf.options.enable_debug

	conf.env.append_unique(
		'CXXFLAGS', 
		['-Wall', '-g']
		)
	

def build(bld):
	# always build libgillespie.so or .dylib(mac)
	bld.shlib(
		source = ['./GillespieSolver.cpp', './GillespieWorld.cpp', './serialize.cpp'],
		includes = ['.'],
		uselib = ['gsl', 'pficommon'],
		target = 'gillespie'
	)
	
	# make executable for testing.
	if bld.env.unit_test == True:
		bld.program(
			source='./test.cpp', 
			includes = ['.'],
			target = 'gillespie_unit',
			defines = ['UNITTEST'],
			use = 'gillespie',
		)

