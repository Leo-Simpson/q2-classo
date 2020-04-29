# q2-classo


## Commands to run an example




#### Activate your qiime2 environment
```shell
source activate qiime2-2020.2
```

PS : it is possible that after this, one might need to do as well :
```shell
conda install zarr
```




#### Push the python code into qiime
```shell
python setup.py install

qiime dev refresh-cache
```

#### Build random data 
```shell
qiime classo generate-data --o-x randomx.qza --o-c randomc.qza
```

#### Apply regress to those random data
```shell
qiime classo regress --i-features randomx.qza --i-c randomc.qza --m-y-file randomy.tsv --m-y-column col --o-result problem.qza
```

#### Make a vizualizer of the solution
```shell
qiime classo summarize \
  --i-problem problem.qza \
  --o-visualization problem.qzv
```

#### View our vizualization file
```shell
qiime tools view problem.qzv
```
