import sys
sys.path.append('/home/mexposit/bioinf/simdata/inDelphi-model/')
import inDelphi

import csv
import numpy as np
import pandas as pd

from math import exp

from progress.bar import Bar

# Doench 2014 parameters
doench_params = [
# pasted/typed table from PDF and converted to zero-based positions
(1,'G',-0.2753771),(2,'A',-0.3238875),(2,'C',0.17212887),(3,'C',-0.1006662),
(4,'C',-0.2018029),(4,'G',0.24595663),(5,'A',0.03644004),(5,'C',0.09837684),
(6,'C',-0.7411813),(6,'G',-0.3932644),(11,'A',-0.466099),(14,'A',0.08537695),
(14,'C',-0.013814),(15,'A',0.27262051),(15,'C',-0.1190226),(15,'T',-0.2859442),
(16,'A',0.09745459),(16,'G',-0.1755462),(17,'C',-0.3457955),(17,'G',-0.6780964),
(18,'A',0.22508903),(18,'C',-0.5077941),(19,'G',-0.4173736),(19,'T',-0.054307),
(20,'G',0.37989937),(20,'T',-0.0907126),(21,'C',0.05782332),(21,'T',-0.5305673),
(22,'T',-0.8770074),(23,'C',-0.8762358),(23,'G',0.27891626),(23,'T',-0.4031022),
(24,'A',-0.0773007),(24,'C',0.28793562),(24,'T',-0.2216372),(27,'G',-0.6890167),
(27,'T',0.11787758),(28,'C',-0.1604453),(29,'G',0.38634258),(1,'GT',-0.6257787),
(4,'GC',0.30004332),(5,'AA',-0.8348362),(5,'TA',0.76062777),(6,'GG',-0.4908167),
(11,'GG',-1.5169074),(11,'TA',0.7092612),(11,'TC',0.49629861),(11,'TT',-0.5868739),
(12,'GG',-0.3345637),(13,'GA',0.76384993),(13,'GC',-0.5370252),(16,'TG',-0.7981461),
(18,'GG',-0.6668087),(18,'TC',0.35318325),(19,'CC',0.74807209),(19,'TG',-0.3672668),
(20,'AC',0.56820913),(20,'CG',0.32907207),(20,'GA',-0.8364568),(20,'GG',-0.7822076),
(21,'TC',-1.029693),(22,'CG',0.85619782),(22,'CT',-0.4632077),(23,'AA',-0.5794924),
(23,'AG',0.64907554),(24,'AG',-0.0773007),(24,'CG',0.28793562),(24,'TG',-0.2216372),
(26,'GT',0.11787758),(28,'GG',-0.69774)]

def calc_doench_score(seq):
    """
    Addapted FROM CRISPOR!! This is the Doench, 2014 score. Not the 2016.
    Code reproduced following paper's methods section. Thanks to Daniel McPherson for fixing it.
    Input is a 30mer: 4bp 5', 20bp guide, 3bp PAM, 3bp 5'
    Output is a gRNA activity score, scaled from 0 (low efficiency) to 100 (high efficiency).
    """
    intercept = 0.59763615
    gc_high = -0.1665878
    gc_low = -0.2026259

    assert(len(seq)==30)
    score = intercept

    grna_seq = seq[4:24]
    gc_count = grna_seq.count("G") + grna_seq.count("C")
    if gc_count <= 10:
        gc_weight = gc_low
    if gc_count > 10:
        gc_weight = gc_high
    score += abs(10-gc_count)*gc_weight

    for pos, model_seq, weight in doench_params:
        sub_seq = seq[pos:pos+len(model_seq)]
        if sub_seq == model_seq:
            score += weight
    exp_score = int(100*(1.0/(1.0+exp(-score))))

    return exp_score

# input data, define cutsite and number of reads
# all gRNAs must be in the same orientation and position (if not, doench score would fail)
target_regions_file = 'C3H_targets.csv'
cutsite = 60
sim_reads = 5000

# #target_seqs = pd.read_csv('target_regions.csv', names=['gRNA_id', 'target_seq'])
target_seqs_data = dict(csv.reader(open(target_regions_file, 'r')))

sim_data = pd.DataFrame()
sim_info = pd.DataFrame()

bar = Bar('Simulating sequences:', max=len(target_seqs_data))

for gRNA_id, target_seq in target_seqs_data.items():

    # Calculate activity score using doench for that guide
    ## 30mer: 4bp 5', 20bp guide, 3bp PAM, 3bp 5'
    seq_for_doench = target_seq[ cutsite - 21 : cutsite + 9]
    doench_score = calc_doench_score(seq_for_doench)
    ## doench score is from 0 to 100, scale to get numb of edited reads
    n_edit_seqs = round(doench_score * sim_reads / 100) #round to have an integer number of reads

    # Calculate editing outcomes
    inDelphi.init_model(celltype = 'mESC')
    pred_df, stats = inDelphi.predict(target_seq, cutsite)
    pred_df = inDelphi.add_mhless_genotypes(pred_df,stats)
    #pred_df = inDelphi.add_genotype_column(pred_df,stats)
    ## adds gaps in the deletions, use add_genotype_column to avoid gaps, but sequences could be confused
    pred_df = inDelphi.add_genotype_column_wgaps(pred_df,stats)
    pred_frequency = np.array(pred_df["Predicted frequency"])
    # normalize probabilities to sum exactly 1
    pred_frequency /= pred_frequency.sum()

    # Simulate data
    ## first, create the edited reads
    edit_seqs = np.random.choice(pred_df["Genotype"], p=pred_frequency, size=(n_edit_seqs))
    ## add non edited reads up to the sim_reads objective
    wt_seqs = np.repeat(target_seq, sim_reads-n_edit_seqs)
    ## merge both and keep in dataframe
    sim_seqs = np.append(edit_seqs,wt_seqs)
    sim_data[gRNA_id] = sim_seqs

    # Keep useful information about simulation
    ## pred_df contains summary of outcomes, but rows are in different order for each simulation
    ## stats has many important parameters
    ## append doench score from 0 to 1, and store in dataframe
    stats['doench_eff_score'] = doench_score
    sim_info[gRNA_id] = pd.Series(stats)

    # Next progress bar
    bar.next()

bar.finish()
sim_data.to_csv('sim_data_grna_lib2.csv', index=None)
sim_info.to_csv('sim_info_grna_lib2.csv')
