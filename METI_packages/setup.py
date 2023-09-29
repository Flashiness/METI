#!/usr/bin/env python
# coding: utf-8

# In[2]:


import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='METI',
    version='0.1',
    author="Jiahui Jiang",
    author_email="jjiang6@mdanderson.org",
    description="METI: Deep profiling of tumor ecosystems by integrating cell morphology and spatial transcriptomics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Flashiness/METI",
    packages=setuptools.find_packages(),
    install_requires=["python-igraph","torch","pandas","numpy","scipy","scanpy","anndata","louvain","scikit-learn", "numba", "TESLAforST"],
)

