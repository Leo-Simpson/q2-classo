# q2-classo - a QIIME 2 plugin for constrained sparse regression and classification

currently developed by Léo Simpson, Evan Bolyen, Christian L. Müller


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
qiime classo transform-features \
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


## Commands for workflow on real data : Immune marker prediction in HIV patients


Bien, J., Yan, X., Simpson, L. and Müller, C. (2020). Tree-Aggregated Predictive Modeling of Microbiome Data : 

"
we consider soluble CD14 (sCD14) measurements in HIV patients as the variable to predict and learn an interpretable regression model from gut microbial amplicon data. sCD14 is a marker of microbial translocation and has been shown to be an independent predictor of mortality in HIV infection (Sandler et al., 2011). Following Rivera-Pinto et al. (2018), we analyze a HIV cohort of n = 151 patients where sCD14 levels (in pg/ml units) and fecal 16S rRNA amplicon data were measured.
"

We provide here a q2-classo workflow to study possible prediction of sCD14 using all available p = 539 bacterial and archaeal OTUs. 



One can do the following commands in the folder example/data_qiime : 

The workflow starts from a file ```table.qza```already provided , a taxonomic table ```taxonomy.qza``` and a metadata sample-metadata-complete

#### CLR transform
```shell
qiime classo transform-features \
	 --i-features table.qza \
	 --o-x xclr.qza
```

#### Aggregate data thanks to a taxonomic table
```shell
qiime classo add-taxa \
	--i-features xclr.qza  \
	--i-taxa taxonomy.qza \
	--o-x xtaxa.qza --o-ca ctaxa.qza

```


#### Split data into training and testing sets
```shell
qiime sample-classifier split-table \
	--i-table xtaxa.qza \
	--m-metadata-file sample-metadata-complete.tsv \
	--m-metadata-column sCD14  \
	--p-test-size 0.2 \
	--p-random-state 123 \
	--p-stratify False \
	--o-training-table xtraining \
	--o-test-table xtest
```


#### Apply classo to the training set to solve the linear regression problem
```shell
qiime classo regress  \
	 --i-features xtraining.qza \
	--i-c ctaxa.qza \
	--m-y-file sample-metadata-complete.tsv \
	--m-y-column sCD14  \
	--p-concomitant False \
	--p-stabsel-threshold 0.5 \
	--p-cv-seed 123456 \
	--p-cv-one-se False \
	--o-result problemtaxa
```


#### Compute the prediction on the testing set, for each model selection chosen
```shell
qiime classo predict \
	--i-features xtest.qza \
	--i-problem problemtaxa.qza \
	--o-predictions predictions.qza
```


#### Compute the visualisation of the problem solved
```shell
qiime classo summarize \
  --i-problem problemtaxa.qza \
  --i-taxa taxonomy.qza \
 --i-predictions predictions.qza \
  --o-visualization problemtaxa.qzv
```


#### Visualization 

```shell
qiime tools view problemtaxa.qzv
```

Alternatively, one can drag&drop the file problemtaxa.qzv on : https://view.qiime2.org
Thanks to this alternative, one can also track the workflow that the qiime2 artifact did. 







