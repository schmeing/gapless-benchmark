#!/usr/bin/env python3
import numpy as np
import os.path
import pandas as pd
import sys

def usage():
    print("calculate_filter_performance.py {folder}")

def main(argv):
    if 0 == len(argv):
        usage()
        sys.exit(1)
    
    if 1 != len(argv):
        print("Wrong number of arguments")
        usage()
        sys.exit(1)
    
    folder = argv[0]
    
    conparts = pd.read_csv(os.path.join(folder, "contig_parts.csv"))
    ground_truth = pd.read_csv(os.path.join(folder, "ground_truth.paf"), sep='\t', names=["qname","qlen","qstart","qend","strand","tname","tlen","tstart","tend","matches","aln_len","mapq"])
    filters_csv = pd.read_csv(os.path.join(folder, "filtered_bridges.csv"))
    alt_filters_csv = pd.read_csv(os.path.join(folder, "alt_filter_bridges.csv"))
    
    conparts.rename(columns={'name':'qname'}, inplace=True)
    truth = conparts.merge(ground_truth.loc[ground_truth['mapq'] >= 30, ['qname','qstart','qend','strand','tname','tstart','tend']], on='qname', how='inner')
    truth = truth[(np.minimum(truth['end'],truth['qend']) - np.maximum(truth['start'],truth['qstart'])) / (truth['end']-truth['start']) >= 0.8].copy()
    truth['ratio'] = (truth['tend']-truth['tstart'])/(truth['qend']-truth['qstart'])
    truth['tstart'] += np.round(np.where(truth['strand'] == '+', truth['start']-truth['qstart'], truth['qend']-truth['end'])*truth['ratio']).astype(int)
    truth['tend'] += np.round(np.where(truth['strand'] == '+', truth['end']-truth['qend'], truth['qstart']-truth['start'])*truth['ratio']).astype(int)
    truth = truth[['conpart','tname','tstart','tend','strand']].copy()
    unique_parts = truth.groupby('conpart').size().reset_index(name='count')
    unique_parts = unique_parts.loc[unique_parts['count'] == 1, ['conpart']].copy()
    truth = truth.merge(unique_parts, on='conpart', how='inner')
    
    correct = pd.concat([filters_csv, alt_filters_csv], ignore_index=True)
    correct = correct.merge(truth.rename(columns={'conpart':'from','tname':'from_name','tstart':'from_start','tend':'from_end','strand':'from_strand'}), on='from', how='inner')
    correct = correct.merge(truth.rename(columns={'conpart':'to','tname':'to_name','tstart':'to_start','tend':'to_end','strand':'to_strand'}), on='to', how='inner')
    correct['dist'] = np.where((correct['from_side'] == 'r') == (correct['from_strand'] == '+'), correct['to_start']-correct['from_end'], correct['from_start']-correct['to_end'])
    correct['correct'] = ( (correct['from_name'] == correct['to_name']) &
                           ((correct['from_side'] == correct['to_side']) != (correct['from_strand'] == correct['to_strand'])) &
                           (correct['dist']*np.where(correct['dist']>=0,0.5,1.5)-100 < correct['mean_dist']) & (correct['mean_dist'] < correct['dist']*np.where(correct['dist']>=0,1.5,0.5)+100) )
    
    fdr = correct.groupby(['fabs','fprem','frel','min_len'], dropna=False)['correct'].agg(['sum','size']).reset_index()
    fdr['fdr'] = 1.0 - fdr['sum']/fdr['size']
    fdr.drop(columns=['sum','size'], inplace=True)
    
    complete = correct.loc[correct['correct'], ['fabs','fprem','frel','min_len','from_name','from_start','from_end','to_start','to_end']].rename(columns={'from_name':'contig'})
    complete['con_start'] = np.minimum(complete['from_start'],complete['to_start'])
    complete['con_end'] = np.maximum(complete['from_end'],complete['to_end'])
    complete.drop(columns=['from_start','from_end','to_start','to_end'], inplace=True)
    
    bin_size = 100000
    truth['start_bin'] = truth['tstart'] // bin_size
    truth['nbins'] = truth['tend'] // bin_size - truth['start_bin'] + 1
    truth.reset_index(inplace=True, drop=True)
    truth = truth.loc[np.repeat(truth.index.values, truth['nbins'].values)].reset_index()
    truth['bin'] = truth['start_bin'] + truth.groupby('index',sort=False).cumcount()
    truth.drop(columns=['index','start_bin','nbins'], inplace=True)
    
    max_comp = complete[['contig','con_start','con_end','min_len']].drop_duplicates()
    max_comp['required'] = np.isnan(max_comp['min_len'])
    max_comp = max_comp.groupby(['contig','con_start','con_end'])['required'].max().reset_index()
    max_comp['id'] = np.arange(len(max_comp))
    max_comp['bin'] = max_comp['con_start'] // bin_size
    covered = max_comp.merge(truth[['tname','tstart','tend','bin']].rename(columns={'tname':'contig'}), on=['contig','bin'], how='inner')
    covered = covered[(covered['con_start'] >= covered['tstart']) & (covered['con_end'] <= covered['tend'])].copy()
    max_comp = max_comp[np.isin(max_comp['id'], np.unique(covered['id'])) == False].copy()
    
    max_comp.reset_index(drop=True, inplace=True)
    max_binned = max_comp.loc[np.repeat(max_comp.index.values, max_comp['con_end'].values // bin_size - max_comp['bin'].values + 1)].reset_index()
    max_binned['bin'] += max_binned.groupby('index',sort=False).cumcount()
    max_binned.drop(columns=['index','required'], inplace=True)
    
    covered = max_binned.rename(columns={'con_start':'start1','con_end':'end1','id':'id1'}).merge(max_binned.rename(columns={'con_start':'start2','con_end':'end2','id':'id2'}), on=['contig','bin'], how='inner')
    covered = covered[(covered['start1'] >= covered['start2']) & (covered['end1'] <= covered['end2'])].rename(columns={'id1':'included','id2':'id'})
    covered = covered[['id','included']].drop_duplicates()
    
    complete = complete.merge(max_comp.drop(columns=['bin','required']), on=['contig','con_start','con_end'], how='inner')
    complete.drop(columns=['contig','con_start','con_end'], inplace=True)
    complete = complete.merge(covered, on='id', how='inner')
    complete.drop(columns=['id'], inplace=True)
    complete.drop_duplicates(inplace=True)
    complete = complete[np.isin(complete['included'], max_comp.loc[max_comp['required'], 'id'].values)].copy()
    complete = complete.groupby(['fabs','fprem','frel','min_len'], dropna=False).size().reset_index(name='covered')
    complete['completeness'] = complete['covered']/sum(max_comp['required'])
    complete.drop(columns='covered', inplace=True)
    
    performance = fdr.merge(complete, on=['fabs','fprem','frel','min_len'], how='left')
    performance.to_csv(os.path.join(folder, "performance.csv"), index=False)

if __name__ == "__main__":
    main(sys.argv[1:])
