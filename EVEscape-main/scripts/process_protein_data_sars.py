import pandas as pd
import numpy as np
import scipy.stats
from sklearn.preprocessing import StandardScaler
from Bio.PDB import PDBParser
from weighted_contact_number import *
from seq_utils import *

##############################################
# General Paths
##############################################

# AA properties ()

aa_charge_hydro = '../data/aa_properties/dissimilarity_metrics.csv'

# SARS2 RBD Paths
##############################################

# Experimental data
rbd_replication = '../data/experiments/starr2020/Starr2020_rbd_bind_expr.csv'
rbd_escape = '../data/experiments/bloom_rbd_escape/escape_data_20220109.csv'
rbd_chan = '../data/experiments/chan2020/abf1738_processed_data_file_from_deep_mutagenesis_of_sars-cov-2_protein_s.xlsx'
rbd_studies_to_drop = ['2021_Greaney_B1351']
rbd_experiment_range = (331, 531)

# Models
rbd_eve_pre2020 = '../results/evol_indices/P0DTC2_321-541_sc0.5_cc0.3_b0.3_pre2020_evol_indices.csv'

# Structure data
rbd_structure_list = [{
    'name': '6VXX',
    'chains': ['A'],
    'trimer_chains': ['A', 'B', 'C'],
    'pdb_path': '../data/structures/6vxx.pdb'
}, {
    'name': '6VYB',
    'chains': ['B'],
    'trimer_chains': ['A', 'B', 'C'],
    'pdb_path': '../data/structures/6vyb.pdb'
}, {
    'name': '7BNN',
    'chains': ['B'],
    'trimer_chains': ['A', 'B', 'C'],
    'pdb_path': '../data/structures/7bnn.pdb'
#}, {
#    'name': 'Delta',
#    'chains': ['A'],
#    'trimer_chains': ['A', 'B', 'C'],
#    'pdb_path': '../data/structures/delta-af2-trunc.cif'
#}, {
#    'name': 'BA1',
#    'chains': ['A'],
#    'trimer_chains': ['A', 'B', 'C'],
#    'pdb_path': '../data/structures/ba1-af2-trunc.cif' 
}, {
    'name': '7CAB',
    'chains': ['A'],
    'trimer_chains': ['A', 'B', 'C'],
    'pdb_path': '../data/structures/7cab.pdb'
}]

rbd_target_seq_path = '../data/sequences/SPIKE_SARS2.fasta'

## RBD metadata save paths ##

bloom_ab_list = '../data/antibody_properties/Bloom_abs_to_use.txt'
xie_ab_list = '../data/antibody_properties/Xie_abs_to_use.txt'
rbd_ab_metadata = '../data/antibody_properties/rbd_antibody_metadata.csv'

##############################################
# SARS2 Spike Paths
##############################################

spike_target_seq_path = '../data/sequences/SPIKE_SARS2.fasta'
spike_experiment_range = (1, 1273)

BA1_target_seq_path = '../../BA1.fasta'
Delta_target_seq_path = '../../Delta.fasta'
Starr_target_seq_path = '../../starr_wt.fasta'
#spike_experiment_range = (1, 1273)

# Models
spike_eve_pre2020 = '../results/evol_indices/P0DTC2_sc0.5_cc0.3_b0.1_pre2020_evol_indices.csv'

##############################################
# Data Processing Functions
##############################################


def process_eve_smm(eve_path):
    '''
    Processes EVE single mutation matrix table
    '''
    eve = pd.read_csv(eve_path)
    eve = eve[1:]
    eve.columns = eve.columns.str.replace("_ensemble", "")
    eve['wt'] = eve.mutations.str[0]
    eve['mut'] = eve.mutations.str[-1]
    eve['i'] = eve.mutations.str[1:-1].astype(int)
    eve['evol_indices'] = -eve.evol_indices
    to_drop = ['protein_name', 'mutations']
    to_drop.extend([col for col in eve.columns if "semantic_change" in col])
    eve = eve.drop(columns=to_drop)
    return eve


def add_model_outputs(exps, eve_path):
    '''
    Merges EVE predictions on to experimental data table
    '''
    exps = exps.merge(process_eve_smm(eve_path),
                      on=['wt', 'mut', 'i'],
                      how='outer')
    return exps


def get_wcn(exps, pdb_path, trimer_chains, target_chains, map_table):
    '''
    Computes weighted contact number by alpha-carbon and sidechain 
    center of mass and merges on to experimental data table
    '''

    wcn = add_wcn_to_site_annotations(pdb_path, ''.join(trimer_chains))
    wcn = wcn.rename(columns={'pdb_position': 'i', 'pdb_aa': 'wt'})
    wcn['i'] = wcn.i.apply(lambda x: alphanumeric_index_to_numeric_index(x)
                           if (x != '') else x)
    wcn['i'] = wcn.i.replace('', np.nan)
    wcn = remap_struct_df_to_target_seq(wcn, target_chains, map_table)

    exps = exps.merge(wcn[['i', 'wcn_sc']], how='left', on='i')
    exps = exps.sort_values('i')
    exps['wcn_bfil'] = exps.wcn_sc.fillna(method='bfill')
    exps['wcn_ffil'] = exps.wcn_sc.fillna(method='ffill')
    exps['wcn_fill'] = (
        exps[['wcn_ffil', 'wcn_bfil']].sum(axis=1, min_count=2) / 2)
    exps = exps.drop(columns=['wcn_bfil', 'wcn_ffil'])
    return exps


def hydrophobicity_charge(exps, table):

    props = pd.read_csv(table, index_col=0)

    scale = StandardScaler()
    props['eisenberg_weiss_diff_std'] = scale.fit_transform(
        props['eisenberg_weiss_diff'].abs().values.reshape(-1, 1))
    props['charge_diff_std'] = scale.fit_transform(
        props['charge_diff'].abs().values.reshape(-1, 1))
    exps = exps.merge(props, how='left', on=['wt', 'mut'])

    exps['charge_ew-hydro'] = exps[[
        'eisenberg_weiss_diff_std', 'charge_diff_std'
    ]].sum(axis=1)
    exps = exps.drop(columns=['eisenberg_weiss_diff_std', 'charge_diff_std'])
    return exps


def norm_to_wt(df, prefvar):
    '''
    Normalize experimental variables to wildtype (for "prefs" style data) 
    '''
    newvar = 'norm_' + prefvar

    def grp_func(grp):

        ref = grp[grp['wt'] == grp['mut']][prefvar].mean()
        grp[newvar] = grp[prefvar] / ref
        return grp

    df[newvar] = df[prefvar]
    df = df.groupby(['i', 'wt']).apply(grp_func)
    return df

def rbd_metadata(escape_df, bloom_path, xie_path, metadata_path):
    escape = escape_df[['condition','condition_type',
                        'condition_subtype','condition_year',
                        'eliciting_virus','study',
                        'lab']].drop_duplicates()
    with open(xie_path, "w") as textfile:
        for element in escape[
                          (escape.lab=='Xie_XS')].condition.tolist():
            textfile.write(element + "\n")
    with open(bloom_path, "w") as textfile:
        for element in escape[
                          (escape.lab=='Bloom_JD')].condition.tolist():
            textfile.write(element + "\n")
    escape.to_csv(metadata_path)
    return(escape)


##############################################
# Summary workbook functions
##############################################


def load_spike(f = ""):

    # Make starting dataframe
    data = make_mut_table(f)

    # Read in and combine pre-2020 model data
    data = add_model_outputs(data, spike_eve_pre2020)

    # Get rid of wt data
    data = data[data.wt != data.mut]

    # Get mapping to PDB
    # Add WCNs to dataframe
    map_table_dict = {}

    for struct in rbd_structure_list:

        map_table = remap_pdb_seq_to_target_seq(struct['pdb_path'],
                                                struct['chains'],
                                                f)
        data = get_wcn(data, struct['pdb_path'], struct['trimer_chains'],
                       struct['chains'], map_table)
        data = data.rename(
            columns={
                'wcn_sc': 'wcn_sc_' + struct['name'],
                'wcn_fill': 'wcn_fill_' + struct['name']
            })

        map_table_dict[struct['name']] = map_table

    # Take min of structures
    data['wcn_fill'] = data[[
        col for col in data.columns if 'wcn_fill_' in col
    ]].min(axis=1)

    # Add aa properties to data
    data = hydrophobicity_charge(data, aa_charge_hydro)
    data = data.sort_values(['i', 'mut'])

    # Drop any rows not in experiment
    data = data[(data.i >= spike_experiment_range[0])
                & (data.i <= spike_experiment_range[1])]

    return data, map_table_dict



#spike, _ = load_spike(spike_target_seq_path)
#spike.to_csv('../spike_scores.csv', index=False)

#spike_BA1, _ = load_spike(BA1_target_seq_path)
#spike_BA1.to_csv('../BA1_scores.csv', index=False)

#spike_Delta, _ = load_spike(Delta_target_seq_path)
#spike_Delta.to_csv('../Delta_scores.csv', index=False)

spike_Starr, _ = load_spike(Starr_target_seq_path)
spike_Starr.to_csv('../Starr_scores.csv', index=False)
