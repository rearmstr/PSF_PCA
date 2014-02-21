import glob

tmv_dir='/data2/home/rarmst/soft/tmv0.71_icpc/'
ccfits_dir='/usr/global/CCfits/'
fitsio_dir='/usr/global/cfitsio/'
flags='-openmp -g'

tmv_incdir=tmv_dir+'include/'
tmv_libdir=tmv_dir+'lib/'

ccfits_incdir=ccfits_dir+'include/'
ccfits_libdir=ccfits_dir+'lib/'

fitsio_incdir=fitsio_dir+'include/'
fitsio_libdir=fitsio_dir+'lib/'

libs=[ 'tmv','tmv_symband','mkl_intel_lp64','mkl_core',
       'mkl_sequential', 'pthread','cfitsio','CCfits']

inc_path=['.',tmv_incdir,ccfits_incdir,fitsio_incdir]
lib_path=[tmv_libdir,ccfits_libdir,fitsio_libdir,
          '/opt/intel/composerxe-2011.4.191/mkl/lib/intel64/']
env=Environment()
env['CXX']='/opt/intel/composerxe-2011.4.191/bin/intel64/icpc'
env['CPPPATH']=inc_path
env['LIBPATH']=lib_path
env['LIBS']=libs
env['CCFLAGS']=Split(flags)
env['LINKFLAGS']=Split(flags)


files=glob.glob('*cc')
for file in files:
    name=file[:-3]
    concat_exp=env.Program(name,[Glob('*.cpp'),file])


