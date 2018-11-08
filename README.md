# Random Effect General(REG) `v1.0.0`(beta)

`reg` is a command line tool for performinig meta-analysis using
GWAS summary statistics. `reg` also estimates cumulative density
funcion of its statistic

## Getting Started

In order to download `reg`, you should clone this repository via the
commands
```
git clone http://github.com/cuelee/reg.git
cd reg
```

In order to run `reg.py` you will need the following python dependencies:
Python v.3.6.5 <
numpy v.1.15.2 <
scipy v.1.1.0 <

Once the above has completed, you can run:


```
./reg.py -h
```

In order to run `xtractmats.py` you will need the `ldsc` software on the same directory of reg.
This can be done by cloning `ldsc`'s repository via the commands
```  
git clone https://github.com/bulik/ldsc.git
```
Also, you will need the [Anaconda](https://store.continuum.io/cshop/anaconda/) Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with xtractsmats' dependencies:

```
conda env create --file environment.yml
source activate ldsc
```

Once the above has completed, you can run:

```
./xtractmats.py -h
```


to print a list of all command-line options. If these commands fail with
an error, then something as gone wrong during the installation process.

## Citation

If you use the software of the Random Effect General model, please cite
[Cue Hyunkyu Lee. et al. Cross-Disease meta-analysis method can identify
the pleiotropic effects.]

## License 

This project is licensed under GNU GPL v3.

## Authors
Cue Hyunkyu Lee ( Seoul National University )

