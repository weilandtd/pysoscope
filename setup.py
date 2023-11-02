from setuptools import setup, find_packages

setup(
    name='pysoscope',
    version='0.0.1',
    url='https://github.com/weilandtd/pysoscope.git',
    author='Daniel R.Weilandt',
    author_email='daniel.r.weilandt@princeton.edu',
    description='A lightweight package to import MALDI data into python',
    packages=find_packages(),
    # TODO @DRW update the install requirements
    install_requires=['numpy >= 1.11.1', 
                      "h5py >= 2.6.0",
                      "scipy >= 0.17.1",
                      "scikit-learn >= 0.17.1",
                      "anndata >= 0.6.22",
                      "seaborn"],
)