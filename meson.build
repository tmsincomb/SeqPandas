project(
  'seqpandas',
  'c', 'cython',
  version: '0.0.2',
  license: 'MIT',
  meson_version: '>=1.2.0',
  default_options: [
    'buildtype=release',
    'b_ndebug=if-release',
    'c_std=c99',
  ],
)

# Python and NumPy dependencies
py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# NumPy dependency handling - robust approach for all platforms including Windows
# Windows often has issues with pkg-config, so we prioritize Python-based detection

# Always try to get NumPy include directory from Python first (most reliable on Windows)
incdir_numpy = run_command(py,
  ['-c', 'import numpy; print(numpy.get_include())'],
  check: false
).stdout().strip()

numpy_found = false
if incdir_numpy != ''
  inc_np = include_directories(incdir_numpy)
  numpy_found = true
  message('NumPy include directory found:', incdir_numpy)

  # Also try to get NumPy version for information
  numpy_version = run_command(py,
    ['-c', 'import numpy; print(numpy.__version__)'],
    check: false
  ).stdout().strip()

  if numpy_version != ''
    message('NumPy version:', numpy_version)
  endif
else
  # Fall back to dependency system (may not work on Windows)
  np_dep = dependency('numpy', required: false)

  if np_dep.found()
    inc_np = []  # Will be handled by np_dep
    numpy_found = true
    message('NumPy found via dependency system')
  else
    # NumPy not available - build will fail
    inc_np = []
    numpy_found = false
    error('NumPy not found - please install numpy before building')
  endif
endif

# Add subdirectories
subdir('seqpandas')
