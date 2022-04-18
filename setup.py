from setuptools import setup


__version__ = '1.1.0'


with open("README.md", "r") as readme_file:
    readme = readme_file.read()


setup(
    name = 'pystapler',
    version = __version__,
    author = 'Adam Moyer',
    author_email = 'atom.moyer@gmail.com',
    description = 'A Motif Hash Based Method for Matching Staples/Crosslinks into Peptides and Proteins',
    packages = ['stapler'],
    package_dir={'stapler' : 'stapler'},
    package_data={'stapler' : ['hash_tables']},
    install_requires = ['pyrosetta', 'xbin', 'numpy', 'getpy'],
    include_package_data=True,
    zip_safe = False,
    long_description=readme,
    long_description_content_type='text/markdown',
    url='https://github.com/atom-moyer/stapler',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: C++',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix'
    ],
)
