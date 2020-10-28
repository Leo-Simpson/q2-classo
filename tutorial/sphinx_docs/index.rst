q2-classo : a QIIME 2 plugin for sparse log-contrast regression and classification tasks for compositional microbiome data
============================================================================================================================


**Title:** q2-classo : a QIIME 2 plugin for sparse log-contrast regression and classification tasks for compositional microbiome data

**Running Title:** q2-classo : a QIIME 2 plugin for sparse log-contrast regression and classification tasks for compositional microbiome data

**Authors:** Leo Simpson :sup:`1,#`, Evan Bolyen :sup:`2,#`, Christian Müller :sup:`3,#`

:sup:`1` Technische Universität München,
:sup:`2` Department of Statistics, Ludwig-Maximilians-Universität München,
:sup:`3` The	Pathogen	and	Microbiome	Institute,	Northern	Arizona	University,	Flagstaff,	AZ,	USA


Abstract
=========

QIIME 2 is a microbiome informatics platform which facilitates
comprehensive and fully reproducible microbiome data science, improving
accessibility to diverse users by adding multiple user interfaces.
The following tutorial describes how to use q2-classo,
a QIIME 2 plugin that incorporate several methods of variable selection
for microbiome data..... 

*Keywords:* Microbiome, QIIME 2, contraint lasso, regression, classification

Introduction
============


USING QIIME 2 WITH MICROBIOME DATA
====================================

Activating QIIME 2
------------


If at any point during the analysis the QIIME 2 conda environment is closed
or deactivated, QIIME 2 2020.20 can be reactivated by running the following
command:

.. command-block::
   :no-exec:

   conda activate qiime2-2020.8

To determine the currently active conda environment, run the following
command and look for the line that starts with “active environment”:

.. command-block::
   :no-exec:

   conda info


Upload data
------------

Taxonomic classification
------------------------

log-contrast processing
------------------------



.. command-block::

    qiime classo features-clr \
        --i-features randomx.qza \
        --o-x xclr.qza




Tree aggregation
-----------------

TAXA transform

.. command-block::

    qiime classo add-taxa \
        --i-features xclr.qza \
        --i-taxa taxonomy.qza \
        --i-c randomc.qza \
        --o-x xtaxa.qza \
        --o-ca ctaxa.qza




Regression tasks
-----------------

Split data into training and testing sets : 

.. command-block::

    qiime sample-classifier split-table \
        --i-table xtaxa.qza \
        --m-metadata-file sample-metadata-complete.tsv \
        --m-metadata-column sCD14  \
        --p-test-size 0.2 \
        --p-random-state 123 \
        --p-stratify False \
        --o-training-table xtraining \
        --o-test-table xtest


Apply classo to the training set to solve the linear regression problem : 

.. command-block::

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


Compute the prediction on the testing set, for each model selection chosen :

.. command-block::

    qiime classo predict \
        --i-features xtest.qza \
        --i-problem problemtaxa.qza \
        --o-predictions predictions.qza

Visualization
--------------

.. command-block::
    qiime classo summarize \
    --i-problem problemtaxa.qza \
    --i-taxa taxonomy.qza \
    --i-predictions predictions.qza \
    --o-visualization problemtaxa.qzv


.. command-block::

    qiime tools view problemtaxa.qzv



Alternatively, one can drag&drop the file problemtaxa.qzv on : https://view.qiime2.org
Thanks to this alternative, one can also track the workflow that the qiime2 artifact did. 


Commentary
==========



Acknowledgements
================


Literature cited
================



