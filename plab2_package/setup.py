#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

requirements = ['appnope', 'attrs', 'backcall', 'certifi', 'charset_normalizer', 'click', 'cycler', 'decorator','flask','pybase64','pytest-shutil','werkzeug','bytesbufio',
                'greenlet', 'httplib2', 'idna', 'iniconfig', 'IPython', 'jedi', 'jinja2', 'jsonpickle', 'kiwisolver', 'markupsafe',
                'matplotlib', 'matplotlib_inline', 'networkx', 'numpy', 'packaging', 'pandas', 'parso', 'pexpect', 'pickleshare',
                'pillow', 'pip', 'pluggy', 'prompt_toolkit', 'ptyprocess', 'py', 'pygments','pyparsing', 'pytest',
                'pytz', 'pyvis', 'requests', 'scipy', 'setuptools', 'six', 'sqlalchemy', 'sqlalchemy_utils', 'tabulate', 'toml', 'tqdm','traitlets',
                'urllib3','wcwidth', 'wheel', 'pymysql', 'cryptography']

test_requirements = ['pytest>=3', ]

setup(
    author="Danqi Wang",
    author_email='wdanqi@live.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    description="Package",
    entry_points={
        'console_scripts': [
            'plab2=plab2.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='plab2',
    name='plab2',
    packages=find_packages(include=['plab2', 'plab2.*']),
    test_suite='tests',
    tests_require=test_requirements,
    version='0.1.0',
    zip_safe=False,
)
