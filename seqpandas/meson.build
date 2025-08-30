# Install Python sources
python_sources = [
  '__init__.py',
  'seqpandas.py',
  'vcf.py',
  'bed.py',
  'cli.py',
  'pileup.py',
  'alignment.py',
  'logo.py',
]

py.install_sources(
  python_sources,
  subdir: 'seqpandas'
)

# Add subdirectories
subdir('tools')
subdir('colors')

# Build Cython extension
# Always try to build it - the module has its own error handling
subdir('cython_numpy')
