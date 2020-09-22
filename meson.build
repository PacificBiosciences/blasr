project(
  'blasr',
  'cpp',
  version : '5.3.3',
  default_options : [
    'buildtype=release',
    'warning_level=3',
    'cpp_std=c++17',
    'b_ndebug=false'],
  license : 'BSD-3',
  meson_version : '>= 0.52.0')

# TODO:
# try to reenable 'b_ndebug=true'
# this will require replacing a number of asserts
# with proper throw statements

############
# CXXFLAGS #
############

blasr_warning_flags = []
cpp = meson.get_compiler('cpp')
foreach cflag: [
  '-Wno-delete-non-virtual-dtor',
  '-Wno-non-virtual-dtor',
  '-Wno-unused-but-set-variable',
  '-Wno-unused-variable']
    if cpp.has_argument(cflag)
      blasr_warning_flags += cflag
    endif
endforeach

if cpp.get_id() == 'clang'
  foreach cflag: [
    '-Wno-overloaded-virtual']
      if cpp.has_argument(cflag)
        blasr_warning_flags += cflag
      endif
  endforeach
endif

################
# dependencies #
################

# threads
blasr_thread_dep = dependency('threads', required : true)

# boost
blasr_boost_dep = dependency('boost', include_type : 'system', required : true)

# pbbam
blasr_pbbam_dep = dependency('pbbam', fallback : ['pbbam', 'pbbam_dep'])

# libblasr
blasr_libblasr_dep = dependency('libblasr', fallback : ['blasr_libcpp', 'libblasr_dep'])

# zlib
blasr_zlib_dep = dependency('zlib', required : true, fallback : ['zlib', 'zlib_dep'])

# htslib
blasr_htslib_dep = dependency('htslib', required : true, version : '>=1.4', fallback : ['htslib', 'htslib_dep'])

blasr_deps = [
  blasr_thread_dep,
  blasr_boost_dep,
  blasr_pbbam_dep,
  blasr_libblasr_dep,
  blasr_zlib_dep,
  blasr_htslib_dep]

########################
# sources + executable #
########################

# replace version strings
blasr_ver_arr = meson.project_version().split('.')
blasr_major_version = blasr_ver_arr[0]
blasr_minor_version = blasr_ver_arr[1]
blasr_patch_version = blasr_ver_arr[2]

blasr_config = configuration_data()
blasr_config.set('BLASR_VERSION', meson.project_version())
blasr_config.set('BLASR_VERSION_MAJOR', blasr_major_version)
blasr_config.set('BLASR_VERSION_MINOR', blasr_minor_version)
blasr_config.set('BLASR_VERSION_PATCH', blasr_patch_version)

blasr_BlasrVersion_h = configure_file(
  input : 'BlasrVersion.h.in',
  output : 'BlasrVersion.h',
  configuration : blasr_config)

# replace git commit id
blasr_BlasrGitHash_h = vcs_tag(
  input : 'BlasrGitHash.h.in',
  output : 'BlasrGitHash.h',
  command : ['git', 'describe', '--always', '--dirty=*'],
  fallback : '',
  replace_string : '@BLASR_GIT_SHA1@')

blasr_sources = []
blasr_include_directories = [include_directories('.')]

subdir('extrautils')
subdir('utils')

blasr_static_impl = static_library(
  'blasr_impl',
  blasr_sources,
  install : false,
  dependencies : blasr_deps,
  cpp_args : [blasr_warning_flags, '-DUSE_PBBAM=1'])

blasr_main = executable(
  'blasr', [
    blasr_BlasrGitHash_h,
    blasr_BlasrVersion_h,
    files('Blasr.cpp')],
  install : true,
  dependencies : blasr_deps,
  link_with : blasr_static_impl,
  cpp_args : [blasr_warning_flags, '-DUSE_PBBAM=1', '-DCMAKE_BUILD=1'])

blasr_utils_sawriter = executable(
  'sawriter', files([
    'utils/SAWriter.cpp']),
  install : true,
  dependencies : blasr_deps,
  link_with : blasr_static_impl,
  cpp_args : [blasr_warning_flags, '-DUSE_PBBAM=1', '-DCMAKE_BUILD=1'])

blasr_utils_toAfg = executable(
  'toAfg', files([
    'utils/ToAfg.cpp']),
  install : false,
  dependencies : blasr_deps,
  link_with : blasr_static_impl,
  cpp_args : [blasr_warning_flags, '-DUSE_PBBAM=1', '-DCMAKE_BUILD=1'])

#########
# tests #
#########

if get_option('tests')
  blasr_clang_formatter = find_program('tools/check-formatting')
  test(
    'blasr formatting check',
    blasr_clang_formatter,
    args : [
      '--all'],
    workdir : meson.source_root())

  subdir('ctest')
endif
