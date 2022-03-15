from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Simulate forward time population genetics with Python'
LONG_DESCRIPTION = 'This package will allow you to simulate evolutionary processes using the concept of forward time population simulation'

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="simuPy", 
        version=VERSION,
        author="Mario Ruiz",
        author_email="<marioruizperez00@email.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        license='GPL',
         install_requires=[], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        #prob will include pandas
        
        keywords=['python', 'evolution'],
        # classifiers= [
        #     "Development Status :: 3 - Alpha",
        #     "Intended Audience :: Education",
        #     "Programming Language :: Python :: 2",
        #     "Programming Language :: Python :: 3",
        #     "Operating System :: MacOS :: MacOS X",
        #     "Operating System :: Microsoft :: Windows",
        # ]
)