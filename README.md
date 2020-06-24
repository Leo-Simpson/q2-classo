# q2-classo


## Commands to run an example


#### Install qiime2 and build the environment

https://docs.qiime2.org/2020.2/install/native/#install-qiime-2-within-a-conda-environment

#### Activate your qiime2 environment
```shell
source activate qiime2-2020.2
```

#### The plugin depends on the packages : 

- zarr 
- plotly
- c-lasso

#### To find an example of taxonomy : 

For example, in the qiime2 tutorial of the parkinson mouse : 

https://docs.qiime2.org/2020.2/tutorials/pd-mice/

at the section "Taxonomic classification", a file called taxonomy.qza can be downloaded, which is a FeatureData[Taxonomy]. 

One can use this taxonomy and build random data "with respect to this taxonomy". 


#### Push the python code into qiime
```shell
python setup.py install

qiime dev refresh-cache
```

#### Build random data 
```shell
qiime classo generate-data \
  --i-taxa taxonomy.qza \
  --o-x randomx.qza \
  --o-c randomc.qza
```

#### CLR transform
```shell
qiime classo features-clr \
  --i-features randomx.qza \
  --o-x xclr.qza
```

#### TAXA transform
```shell
qiime classo add-taxa \
  --i-features xclr.qza \
  --i-taxa taxonomy.qza \
  --i-c randomc.qza \
  --o-x xtaxa.qza \
  --o-ca ctaxa.qza
```


#### Apply regress to those random data
```shell
qiime classo regress \
  --i-features xtaxa.qza\
  --i-c ctaxa.qza\
  --m-y-file randomy.tsv\
  --m-y-column col\
  --o-result problem.qza
```

#### Make a visualizer of the solution
```shell
qiime classo summarize \
  --i-taxa taxonomy.qza \
  --i-problem problem.qza \
  --o-visualization problem.qzv
```

#### View our visualization file
```shell
qiime tools view problem.qzv
```
