from setuptools import setup

setup(
    name="tcbp",
    version="0.1",
    author='Gabriele Bozzola',
    author_email='gabrielebozzola@arizona.edu',
    packages=['tcbp'],
    description='Solve the Newtonian two body problem with charges and radiation reaction.',
    install_requires=["numpy", "scipy"],
    license='LICENSE.txt',
)
