<h1><center>METI Tutorial</center></h1>

## Outline
1. [Installation]
2. [Import modules]
3. [Read in data]

### 1. Installation
The installation should take a few minutes on a normal computer. 

Now you can install the current release of METI by the following three ways:
#### 1.1 PyPI: Directly install the package from PyPI.


```python
pip3 install METI
#Note: you need to make sure that the pip is for python3，or you can install METI by
python3 -m pip install METI==0.1
#If you do not have permission (when you get a permission denied error), you can install METI by
pip3 install --user METI==0.1
```
#### 1.2 Github
Download the package from Github and install it locally:


```python
git clone https://github.com/Flashiness/METI
cd METI/METI_package/
python3 setup.py install --user
```




