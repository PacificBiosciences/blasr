#!/usr/bin/env python
"""Configure the build.

- Create defines.mk
"""
import commands
import contextlib
import optparse
import os
import sys
import warnings

#DEFAULTCXXFLAG := -O3
#DEBUGCXXFLAG := -g -ggdb -fno-inline
#PROFILECXXFLAG := -Os -pg
#GCXXFLAG := -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -fno-omit-frame-pointer

ROOT = os.path.abspath(os.path.dirname(__file__))

def log(msg):
    sys.stderr.write(msg)
    sys.stderr.write('\n')

def shell(cmd):
    log('`%s`'%cmd)
    status, output = commands.getstatusoutput(cmd)
    if status:
        raise Exception('%d <- %r' %(status, cmd))
    log(output)
    return output

def system(cmd):
    log(cmd)
    status = os.system(cmd)
    if status:
        raise Exception('%d <- %r' %(status, cmd))
    return

def mkdirs(path):
    if not os.path.isdir(path):
        os.makedirs(path)

@contextlib.contextmanager
def cd(nwd):
    cwd = os.getcwd()
    log('cd %r -> %r' %(cwd, nwd))
    os.chdir(nwd)
    yield
    os.chdir(cwd)
    log('cd %r <- %r' %(cwd, nwd))

def update_content(fn, content):
    current_content = open(fn).read() if os.path.exists(fn) else None
    if content != current_content:
        log('writing to %r:' %fn)
        log('"""\n' + content + '"""\n----')
        open(fn, 'w').write(content)

def get_OS_STRING():
    G_BUILDOS_CMD = """bash -c 'set -e; set -o pipefail; id=$(lsb_release -si | tr "[:upper:]" "[:lower:]"); rel=$(lsb_release -sr); case $id in ubuntu) printf "$id-%04d\n" ${rel/./};; centos) echo "$id-${rel%%.*}";; *) echo "$id-$rel";; esac' 2>/dev/null"""
    return shell(G_BUILDOS_CMD)
def get_PREBUILT():
    cmd = 'cd ../../../../prebuilt.out 2>/dev/null && pwd || echo -n notfound'
    return shell(cmd)
def ifenvf(env, key, func):
    if key in env:
        return env[key]
    else:
        return func()
def setifenvf(envout, envin, key, func):
    envout[key] = ifenvf(envin, key, func)
def setifenv(envout, envin, key, val):
    envout[key] = envin.get(key, val)
def setenv(envout, key, val):
    envout[key] = val
def update_env_if(envout, envin, keys):
    for key in keys:
        if key in envin:
            envout[key] = envin[key]
def compose_defs_env(env):
    # We disallow env overrides for anything with a default from GNU make.
    nons = ['CXX', 'CC', 'AR'] # 'SHELL'?
    ovr    = ['%-20s ?= %s' %(k, v) for k,v in env.items() if k not in nons]
    nonovr = ['%-20s := %s' %(k, v) for k,v in env.items() if k in nons]
    return '\n'.join(ovr + nonovr + [''])
def compose_defines_pacbio(envin):
    """
    This is used by mobs via buildcntl.sh.
    """
    env = dict()
    setenv(env, 'SHELL', 'bash')
    #setifenvf(env, envin, 'OS_STRING', get_OS_STRING)
    #setifenvf(env, envin, 'PREBUILT', get_PREBUILT)
    nondefaults = set([
            'CXX',
            'BLASR_INC',
            'LIBPBDATA_INC', 'LIBPBDATA_LIB', 'LIBPBDATA_LIBFLAGS',
            'LIBPBIHDF_INC', 'LIBPBIHDF_LIB', 'LIBPBIHDF_LIBFLAGS',
            'LIBBLASR_INC', 'LIBBLASR_LIB', 'LIBBLASR_LIBFLAGS',
            'HDF5_INC', 'HDF5_LIB', 'HDF5_LIBFLAGS',
            'PBBAM_INC', 'PBBAM_LIB', 'PBBAM_LIBFLAGS',
            'HTSLIB_INC', 'HTSLIB_LIB', 'HTSLIB_LIBFLAGS',
            'BOOST_INC',
            'GCC_LIB',
            'ZLIB_LIB', 'ZLIB_LIBFLAGS',
            'SZLIB_LIB', 'SZLIB_LIBFLAGS',
            'PTHREAD_LIBFLAGS',
            'DL_LIBFLAGS',
            'RT_LIBFLAGS',
    ])
    update_env_if(env, envin, nondefaults)
    return compose_defs_env(env)


def configure_pacbio(envin, shared, build_dir):
    content1 = compose_defines_pacbio(envin)
    if not shared:
        content1 += 'LDFLAGS+=-static\n'
    update_content(os.path.join(build_dir, 'defines.mk'), content1)

def set_defs_submodule_defaults(env, nopbbam):
    subdir = os.path.join(ROOT, 'libcpp')
    defaults = {
        'LIBPBDATA_INC': os.path.join(subdir, 'pbdata'),
        'LIBBLASR_INC':  os.path.join(subdir, 'alignment'),
        #'LIBPBIHDF_INC': '' if nopbbam else os.path.join(subdir, 'hdf'),
        'LIBPBDATA_LIB': os.path.join(subdir, 'pbdata'),
        'LIBBLASR_LIB':  os.path.join(subdir, 'alignment'),
        #'LIBPBIHDF_LIB': '' if nopbbam else os.path.join(subdir, 'hdf'),
    }
    for k in defaults:
        if k not in env:
            env[k] = defaults[k]

def update_defaults_for_os(env):
    OS = shell('uname')
    if 'Darwin' in OS:
        #-lsz (for static builds?)
        env['RT_LIBFLAGS'] = ''

def set_defs_defaults(env, nopbbam, with_szlib):
    defaults = {
        'BLASR_INC': os.path.join(ROOT, 'include'),
        'LIBBLASR_INC':  os.path.join(ROOT, 'libcpp', 'alignment'),
        'LIBPBDATA_INC':  os.path.join(ROOT, 'libcpp', 'pbdata'),
        'LIBPBIHDF_INC':  os.path.join(ROOT, 'libcpp', 'hdf'),
        'LIBBLASR_LIB':  os.path.join(ROOT, 'libcpp', 'alignment'),
        'LIBPBDATA_LIB':  os.path.join(ROOT, 'libcpp', 'pbdata'),
        'LIBPBIHDF_LIB':  os.path.join(ROOT, 'libcpp', 'hdf'),
        'LIBBLASR_LIBFLAGS':  '-lblasr',
        'LIBPBDATA_LIBFLAGS': '-lpbdata',
        'LIBPBIHDF_LIBFLAGS': '-lpbihdf',
        'HDF5_LIBFLAGS': '-lhdf5_cpp -lhdf5',
        'RT_LIBFLAGS': '-lrt',
        'ZLIB_LIBFLAGS': '-lz',
        'PTHREAD_LIBFLAGS': '-lpthread',
        'DL_LIBFLAGS': '-ldl', # neeeded by HDF5 always
        'SHELL': 'bash -xe',
    }
    try:
        update_defaults_for_os(defaults)
    except Exception as e:
        warnings.warn(e)
    #setifenvf(defaults, env, 'OS_STRING', get_OS_STRING)
    #setifenvf(defaults, env, 'PREBUILT', get_PREBUILT)
    pbbam_defaults = {
        'PBBAM_LIBFLAGS': '-lpbbam',
        'HTSLIB_LIBFLAGS': '-lhts',
        'ZLIB_LIBFLAGS': '-lz',
        #'PTHREAD_LIBFLAGS': '-lpthread',
        #'DL_LIBFLAGS': '-ldl', # neeeded by HDF5 always
    }
    if not nopbbam:
        defaults.update(pbbam_defaults)
    szlib_defaults = {
        'SZLIB_LIBFLAGS': '-lsz',
        #'ZLIB_LIBFLAGS': '-lz', # probably needed, but provided elsewhere
    }
    if with_szlib:
        defaults.update(szlib_defaults)
    for k in defaults:
        if k not in env:
            env[k] = defaults[k]

def get_make_style_env(envin, args):
    envout = dict()
    for arg in args:
        if '=' in arg:
            k, v = arg.split('=')
            envout[k] = v
    envout.update(envin)
    return envout

def parse_args(args):
    parser = optparse.OptionParser()
    parser.add_option('--no-pbbam', action='store_true',
            help='Avoid compiling anything which would need pbbam.')
    parser.add_option('--with-szlib', action='store_true',
            help='If HDF5 was built with --with-szlib, then -lz is needed for static binaries.')
    parser.add_option('--submodules', action='store_true',
            help='Set variables to use our git-submodules, which must be pulled and built first. (Implies --no-pbbam.)')
    parser.add_option('--shared', action='store_true',
            help='Build for dynamic linking. (Non-static binaries.)')
    parser.add_option('--build-dir',
            help='Can be different from source directory, but only when *not* also building submodule.')
    return parser.parse_args(list(args))

def symlink_makefile(build_dir_root, src_dir_root, makefilename, relpath):
    src_dir = os.path.join(src_dir_root, relpath)
    build_dir = os.path.join(build_dir_root, relpath)
    src_name = os.path.join(src_dir, 'makefile')
    dst_name = os.path.join(build_dir, 'makefile')
    if os.path.lexists(dst_name):
        os.unlink(dst_name)
    print('%r <- %r' %(src_name, dst_name))
    mkdirs(build_dir)
    os.symlink(src_name, dst_name)

def symlink_makefiles(build_dir):
    symlink_makefile(build_dir, ROOT, 'makefile', '.')
    symlink_makefile(build_dir, ROOT, 'makefile', 'utils')
    symlink_makefile(build_dir, ROOT, 'makefile', 'extrautils')

def main(prog, *args):
    """We are still deciding what env-vars to use, if any.
    """
    # Set up an alias, until everything uses one consistently.
    conf, makevars = parse_args(args)
    if conf.build_dir is not None:
        symlink_makefiles(conf.build_dir)
    else:
        conf.build_dir = '.'
    conf.build_dir = os.path.abspath(conf.build_dir)
    envin = get_make_style_env(os.environ, makevars)
    if 'HDF5_INCLUDE' in envin and 'HDF5_INC' not in envin:
        envin['HDF5_INC'] = envin['HDF5_INCLUDE']
    if conf.submodules:
        set_defs_submodule_defaults(envin, conf.no_pbbam)
        conf.no_pbbam = True
    set_defs_defaults(envin, conf.no_pbbam, conf.with_szlib)
    configure_pacbio(envin, conf.shared, conf.build_dir)


if __name__=="__main__":
    main(*sys.argv)
