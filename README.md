# PLEIO 
- Pleiotropic Locus Exploration and Interpretation using Optimal test

`pleio` is a command line tool for performinig meta-analysis of multiple GWAS summary statistics. 

## Getting Started

In order to download `pleio`, you should clone this repository via the command
```
git clone https://github.com/cuelee/pleio.git
cd pleio
```

In order to run `pleio`, the following python dependencies must be installed in your system.

- python >= 3.7.4
- pandas
- numpy
- scipy


Once the above has completed, you can run the following command:

```
./pleio,py -h
```

## Updating PLEIO
You can update to the newest version of `PLEIO` using git. First, navigate to your pleio/ directory (e.g., cd pleio), then run
```
git pull
```
If `PLEIO` is up to date, you will see
```
Already up-to-date.
```
otherwise, you will see `git` output similar to 
```
remote: Enumerating objects: 9, done.
remote: Counting objects: 100% (9/9), done.
remote: Compressing objects: 100% (3/3), done.
remote: Total 6 (delta 4), reused 5 (delta 3), pack-reused 0
Unpacking objects: 100% (6/6), done.
From git://github.com/hanlab-SNU/pleio
   e065a06..14c3399  master     -> origin/master
Updating e065a06..14c3399
Fast-forward
 README.md       | 2 +-
 ldsc_preprocess | 2 +-
 2 files changed, 2 insertions(+), 2 deletions(-)
```
## For first-time users
An executable example file and run command can be found on the WIKI page of this Github repository. [Link](https://github.com/cuelee/pleio/wiki)

## Citation

If you use the software, please cite  
Lee, C. H., Shi, H., Pasaniuc, B., Eskin, E., & Han, B. (2021). PLEIO: a method to map and interpret pleiotropic loci with GWAS summary statistics. The American Journal of Human Genetics, 108(1), 36–48. [Link](https://doi.org/10.1016/j.ajhg.2020.11.017)

## Support

Issues with PLEIO? Email hl3565@cumc.columbia.edu

## License 

This project has no license currently.

## Authors

Cue Hyunkyu Lee ( Columbia University Irving Medical Center)
