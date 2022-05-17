from setuptools import setup

setup(name='text_extraction',
      version='0.1',
      description='some text processing',
      url='https://github.com/MatthewMong/text_extraction',
      author='Matthew Mong',
      author_email='matthew.mong1999@gmail.com',
      license='MIT',
      packages=['flask', 'biopython'],
      install_requires=[
          'markdown',
      ],
      zip_safe=False)