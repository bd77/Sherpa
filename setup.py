"""Simple setup script."""

from setuptools import setup

def readme():
    """Read the README file."""
    return open('README.rst').read()

setup(name='sherpadf',
      version='0.1',
      description='SHERPA using dataframe approach',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Education :: Testing',
      ],
      keywords='sherpadf',
      url='',
      author='Denise Pernigotti',
      author_email='',
      license='MIT',
      packages=['sherpadf'],
      # uncomment this and add package names if you need to have certain
      # packages installed
      #install_requires=[
      #    'module_name1', 'module_name2'
      #],
      #entry_points={
      #  'console_scripts': ['prova-play=prova.cmd:main'],
      #  },
      zip_safe=False)
