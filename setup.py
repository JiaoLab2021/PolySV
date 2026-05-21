from setuptools import setup, find_namespace_packages

# Read the contents of your README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="polysv",
    version="1.0.6",
    author="He Jiabao",
    author_email="hejiabao2001@gmail.com or hejiabao@webmail.hzau.edu.cn",
    description="A high-quality SV detection, genotyping and optimization pipeline.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JiaB-He/PolySV",
    packages=find_namespace_packages(include=['include', 'include.*']),
    py_modules=['polysv'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.9',
    entry_points={
        'console_scripts': [
            'polysv = polysv:main',
        ],
    },
)
