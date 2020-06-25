# CRISPRed

This repository contains the notebooks with the scripts used to analyze genomic data and train machine learning models to predict CRISPR-Cas9 gene editing outcomes.

Some of the most interesting and interactive notebooks are:

- 2.5. EffModel.ipynb. Used to train and compare the gene editing efficiency models.
- 2.6. OutcomesModel.ipynb. To train and compare the gene editing outcome frequency models.
- 2.2. LabelGen.ipynb. Simple script to classify the outcomes in categories. The efficiency outcomes are used in 2.5.EffModel.ipynb, and the outcomes quantified are used in 2.3.OutcomesProfiling.ipynb.
- 2.3. OutcomesProfiling.ipynb. Characterization of the simulated data, and preparation of the labels for the 2.6. OutcomesModel notebook.
- 1.2. gRNALibDistribution.ipynb. Analysis of the gRNA library distribution in plasmids as a quality control.
- 1.4. CoverageAnalysis.ipynb. Used to assess the coverage of the genomic regions.

Description of the notebook contents (in catalan):

1. DescriptiveAnalysis/

1.1. SequencingDataProcessing.ipynb. Conté els scripts utilitzats per el processament de qualitat de les dades genòmiques i per l’alineament amb el genoma de referència.

1.2. gRNALibDistribution.ipynb. Conté els scripts utilitzats per el processament de les dades de target sequencing en el control de qualitat de la llibreria. Conté codi que es pot executar per visualitzar les figures generades.

1.3. CoordiantesToC3H.ipynb. Conté la informació del procés utilitzat per convertir les coordenades d’interès del genoma mm10 al genoma C3H.

1.4. CoverageAnalysis.ipynb. Conté els scripts utilitzats en l’anàlisi del coverage, es pot executar per generar les figures mostrades a la memòria (cal esperar uns 10 minuts).

2. Model/

2.1. DataSimulation.ipynb. Conté els scripts utilitzats per a simular les dades, combinant el model d’eficiència de Doench et.al. i el model Indelphi2 per predir les freqüències dels productes d’edició. Es pot executar un exemple per simular 5000 reads de 3 regions d’interès.

2.2. LabelGen.ipynb. Conté l’script utilitzat per analitzar les dades simulades i quantificar l’eficiència de cada regió i la freqüència de les categories de productes d’edició genètica segons la regió. Es mostra la quantificació de 3 regions d’interès com a exemple.

2.3. OutcomesProfiling.ipynb. Mostra l’ús de Seaborn i Matplotlib per representar visualment els resultats de l’anàlisi descriptiu de les dades simulades.

2.4. Featurization.ipynb. Conté l’script utilitzat per convertir cada regió d’interès en un seguit de descriptors de la seva seqüència i utilitzat per entrenar els models d’aprenentatge automàtic.

2.5. EffModel.ipynb. Conté el procés d’entrenament dels classificadors per predir l’eficiència d’edició, incloent la preparació de les dades, l’entrenament dels quatre models de classificació i la comparació entre ells i amb el model de Doench et.al. Es pot executar íntegrament per reproduir els resultats mostrats a la memòria.

2.6. OutcomesModel.ipynb. Conté la planificació i entrenament dels classificadors per predir la freqüència dels resultats d’edició. Es pot executar íntegrament per reproduir les figures de la memòria. 
