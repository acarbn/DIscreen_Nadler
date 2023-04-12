
from Bio import SeqIO
from Bio.PDB import *
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from pyaaisc import Aaindex
import numpy as np
import pandas as pd
import math
from Bio.SubsMat import MatrixInfo as matlist
from Bio import pairwise2
from sklearn.metrics import auc
from sklearn.metrics import RocCurveDisplay, PrecisionRecallDisplay
import matplotlib.pyplot as plt


plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['lines.color'] = 'black'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['legend.frameon'] = False
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['lines.markeredgewidth'] = 1.5
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Cantarell'



def process_input(input):
    '''
    processing of the data corresponding to the initial (unsorted) libraries
    '''
    data_input = pd.read_csv(input, delimiter=',')
    data_input.insertion_site = data_input.insertion_site/3
    data_input['insertion_site'] = data_input['insertion_site']-1 #remove linker S after M
    data_input = data_input.loc[(data_input['in_frame_insertion'] == True) & (data_input['forward_insertion'] == True) & 
                                (data_input['linker_seq'].isnull())]
    data_input = data_input[data_input.insertion_site >= 0]
    return data_input

def process_enriched(input_enriched):
    '''
    processing of the final enriched libraries
    all positions that disappear completely will are set to -10
    '''
    data_cond = pd.read_csv(input_enriched, delimiter=',')
    data_cond.insertion_site = data_cond.insertion_site/3
    data_cond['insertion_site'] = data_cond['insertion_site']-1 #remove linker S after M
    data_cond = data_cond.loc[(data_cond['in_frame_insertion'] == True) & (data_cond['forward_insertion'] == True) & 
                            (data_cond['linker_seq'].isnull())]
    data_norm = pd.DataFrame([data_cond['insertion_site'].value_counts(), data_input['insertion_site'].value_counts(), 
    pd.RangeIndex(1, len(prot_dict[combination.split('_')[0]])+1).to_series()]).T
    data_norm.index = data_norm.index.astype(int)
    data_norm.columns = ['cond', 'input', 'position']
    data_norm['cond'] = data_norm['cond']/np.sum(data_norm['cond'])
    data_norm['input'] = data_norm['input']/np.sum(data_norm['input'])
    data_norm['norm'] = data_norm['cond']/data_norm['input']
    data_norm['position'] = data_norm.index.astype(str)
    data_norm.fillna(0, inplace=True)
    data_norm = data_norm[data_norm.index >= 0]
    data_norm = data_norm[data_norm.index < len(prot_dict[combination.split('_')[0]])+1]
    data_norm['log'] = np.log2(data_norm['norm'])
    data_norm.loc[data_norm.input == 0, 'log'] = np.nan
    data_norm.loc[data_norm.log == -np.inf, 'log'] = -10
    return data_norm

def get_Sec_struct(structure, name, wd):
    '''
    calculate secondary structure of each AA
    '''
    dssp = DSSP(structure[0], f'{wd}/{name}.pdb', dssp='mkdssp')
    sequence = ''
    Sec_structure = ''
    for z in range(len(dssp)):
        a_key = list(dssp.keys())[z]
        sequence += dssp[a_key][1]
        Sec_structure += dssp[a_key][2]
    Sec_structure = Sec_structure.replace('-', 'C')
    Sec_structure = Sec_structure.replace('I', 'C')
    Sec_structure = Sec_structure.replace('T', 'C')
    Sec_structure = Sec_structure.replace('S', 'C')
    Sec_structure = Sec_structure.replace('G', 'H')
    Sec_structure = Sec_structure.replace('E', 'B')
    Sec_structure = [AA for AA in Sec_structure]
    return Sec_structure, dssp

def get_surface_accessibility(dssp, name, window, property_dict):
    '''
    Extracts the surface accessibility  (which is stored in the same variable) for the given 'structure'
    Averages the surface accessibility around within the window
    '''
    surface = []
    for i in dssp.property_list:
        surface.append(i[3])
    surface_smooth = []
    if window == 0:
        for i in range(len(surface)):
            surface_smooth.append(np.mean(surface[i:i+1]))      
    else:
        for position in property_dict['prox_AAs']:
            temp_list = [surface[AA] for AA in position]
            surface_smooth.append(sum(temp_list)/len(temp_list))
    return surface_smooth

def get_pLDDT(structure, window, property_dict):
    '''
    Extracts the b-factor or pLDDT score (which is stored in the same variable) for the given 'structure'
    Averages the b-factor around within the window
    '''
    #get b factor
    pLDDT = []
    for chain in structure[0]: 
        for residue in chain:
            for atom in residue:
                if atom.full_id[4][0] == 'CA':
                    pLDDT.append(atom.bfactor)

    #smoothen b-factor
    pLDDT_smooth= []
    if window == 0:
        for i in range(len(pLDDT)):
            pLDDT_smooth.append(np.mean(pLDDT[i:i+1]))      
    else:
        for position in property_dict['prox_AAs']:
            temp_list = [pLDDT[AA] for AA in position]
            pLDDT_smooth.append(sum(temp_list)/len(temp_list))
    return pLDDT_smooth

# AA feature extraction
def import_features(feature_dict):
    '''
    import AA features from AAindex DB
    feature_dict: provide a dict with feature names and the respective aaindex IDs
    '''
    aaindex = Aaindex()
    for feature, ID in feature_dict.items():
        record = aaindex.get(ID)
        feature_dict[feature] = record.index_data
    return feature_dict

def apply_features(input_features, protein_seq, property_dict, window,  prox_mode=False):
    '''
    - map the given features to the primary sequence of a protein
    - if the proximity flag is set, the values are defined for the residues that are selected as proximal and the mean will be returned
    provide a name for the resulting df, as well as the input sequence and feature_dict
    '''
    if not prox_mode:
        window = round(window/2)
        for name, feat in input_features.items():
            feature_list = [feat[AA] for AA in protein_seq]
            property_dict[name] = []
            for i in range(len(feature_list)):
                if i < window-1:
                    property_dict[name].append(np.mean(feature_list[0:i+window+1]))
                elif i > len(feature_list)-window-1:
                    property_dict[name].append(np.mean(feature_list[i-window+1:len(feature_list)+1]))
                else:
                    property_dict[name].append(np.mean(feature_list[i-window+1:i+window+1]))
    else:
        for name, feat in input_features.items():
            property_dict[name] = []
            for position in property_dict['prox_AAs']:
                temp_list = [feat[protein_seq[AA]] for AA in position]
                property_dict[name].append(sum(temp_list)/len(temp_list))
    return property_dict

def find_proximal_AAs(property_dict, radius):
    '''
    defines insertion site as the middle between the c-alphas of the residues between which a domain is inserted.
    stores proximal residues of each position.
    Proximal residues are defined as the residues the c-alpha of which is located within the given radius (angstrom).
    '''
    residues = [r for r in property_dict['structure'].get_residues()]
    property_dict['prox_AAs'] = []
    for r1 in range(len(residues)):
        if r1+1 <= len(residues)-1:
            position = (residues[r1]["CA"].get_coord() + residues[r1+1]["CA"].get_coord())/2 # insertion position is distance between the 2 Calphas
        else:
            position = residues[r1]["CA"].get_coord()
        property_dict['prox_AAs'].append([])
        for r2 in range(len(residues)):
            residues[r1]["CA"].get_coord()            
            distance = np.linalg.norm(position-residues[r2]['CA'].get_coord())
            if distance < radius:
                property_dict['prox_AAs'][r1].append(r2)
    return property_dict

def process_alignment(alignment_file, uniprot_query):
    '''
    import and process alignment
    alignment_file: file with the MSA
    uniprot_query: Name of the reference protein within the file (usually the uniprot ID)
    '''
    fasta_alignment = SeqIO.parse(open(f'/work/projects/project01640/jm/domain_insertion/20220328_NGS_enriched-2/MSAs/{uniprot_query}_ref_aln.afa'),'fasta')
    alignment_dict = {}
    for i in fasta_alignment:
        alignment_dict[i.name] = str(i.seq)
    alignment_df = pd.DataFrame.from_dict([alignment_dict]).T
    alignment_df.columns = ['seq']
    alignment_df = alignment_df['seq'].str.split('', expand=True)
    alignment_df = alignment_df.iloc[:, 1:-1]
    alignment_counts = alignment_df.apply(pd.Series.value_counts)
    query = alignment_df.loc[uniprot_query]
    alignment_counts.drop('X', inplace=True)
    return alignment_df, query, alignment_counts

def KLD(alignment_counts, background, query):
    '''
    Calculates the Kullback-Leibler-Divergence for each position of the query protein
    '''
    entropy = []
    for idx, values in alignment_counts.T.iterrows():
        local_ent = 0
        for AA, i in values.items():
            if np.isnan(i) == False and AA != '-' and AA != 'B':
                local_ent += i/np.nansum(values)*math.log((i/np.nansum(values))/background[AA],10)
        entropy.append(local_ent)
    query_idx = np.where(query != '-')[0]
    entropy = np.array(entropy)
    entropy = entropy[query_idx]
    return entropy, query_idx

def insertion_stats(query_idx, primary_AA, property_dict):
    '''
    Calculates the insertion frequencies from pairwise alignments between the wildtype sequence and inputs
    The degree of relation to which the sequences are taken into account can be varied by changing the alugnments[0][2] cutoff
    '''
    matrix = matlist.blosum62
    insertions = [0]*len(query_idx)
    deletions = [0]*len(query_idx)
    ins_length = [[0]]*len(query_idx)
    for idx, variant in property_dict['alignment_df'].reset_index(drop=True).iterrows():
        seq = variant[variant != '-'].str.cat()
        alignments = pairwise2.align.globalds(primary_AA, seq, matrix,  -12, -3)
        position = -1
        if alignments[0][2] > 100:
            for i in range(len(alignments[0][0])):
                if alignments[0][0][i] != '-' and alignments[0][0][i-1] != '-':
                    position += 1
                elif alignments[0][0][i] != '-' and alignments[0][0][i-1] == '-':
                    position += 1
                    if ins_length[position-ins_len][0] == 0:
                        ins_length[position-ins_len] = [ins_len]
                    else:
                        ins_length[position-ins_len].append(ins_len)
                elif alignments[0][0][i] == '-' and alignments[0][0][i-1] != '-':
                    ins_len = 1
                    insertions[position] += 1
                elif alignments[0][0][i] == '-' and alignments[0][0][i-1] == '-':
                    ins_len += 1
                if alignments[0][1][i] == '-':
                    deletions[position] += 1
    insert_freq = [x/len(property_dict['alignment_df']) for x in insertions]
    deletion_freq = [x/len(property_dict['alignment_df']) for x in deletions]
    ave_ins_len = [np.mean(x) for x in ins_length]
    median_ins_len = [np.median(x) for x in ins_length]
    property_dict.update({'Insert_frequency':insert_freq, 'Deletion_frequency':deletion_freq, 'Mean_ins_len':ave_ins_len, 'Median_ins_len':median_ins_len})
    return property_dict

def one_hot_insertion(seq, alphabet = 'ACDEFGHIKLMNPQRSTVWY'):
    '''
    one-hot encoding of AA stings
    '''
    # mapping of AAs to integers
    AA_to_int = dict((c, i) for i, c in enumerate(alphabet))
    
    # integer encode input data
    integer_encoded = [AA_to_int[AA] for AA in seq]

    # one hot encode
    onehot_encoded = []
    for idx, value in enumerate(integer_encoded):
        letter = [0 for _ in range(len(alphabet))]
        
        letter[value] = 1
        if idx < len(integer_encoded)-1:
            letter[integer_encoded[idx+1]] = 1
        onehot_encoded.append(letter)
    return onehot_encoded, AA_to_int

class plot_prec_rec():
    '''
    Run classifier with cross-validation and plot ROC curves
    '''
    def __init__(self, classifier, data_split, X, y, title_text, groups):
        self.classifier = classifier
        self.data_split = data_split
        self.X = X
        self.y = y 
        self.groups = groups
        self.title_text = title_text
        self.average_precision = []
        self.mean_recall = np.linspace(0, 1, 100)
        self.colors = ['#332345', '#40498E', '#357BA3', '#38AAAC', '#79D6AE']

    def calculate(self, X, y, data_split, classifier):
        fig, self.ax = plt.subplots(figsize=(5,5))
        for i, (train, test) in enumerate(data_split.split(X, y, groups=self.groups)):
            classifier.fit(X[train], y.iloc[train])
            viz = PrecisionRecallDisplay.from_estimator(
                classifier,
                X[test],
                y.iloc[test],
                name="AP fold {}".format(i),
                alpha=0.8,
                lw=3,
                ax=self.ax,
                color=self.colors[i],
            )
            self.average_precision.append(viz.average_precision)
    
    def plot(self, title_text):
        mean_prec_tot = np.mean(self.average_precision)
        std_prec = np.std(self.average_precision)
        plt.plot([], [], ' ', label=r"Mean ave. prec. = %0.2f $\pm$ %0.2f" % (mean_prec_tot, std_prec))

        self.ax.set(
            xlim=[-.02, 1.05],
            ylim=[0, 1.05],
            title=f"Precision recall: {title_text}",
        )
        self.ax.legend(loc='center', bbox_to_anchor=(0.5, -.3), ncol=2)
        self.ax.set_aspect('equal')
        #self.ax.legend(loc="lower right")
        plt.savefig(f"/work/projects/project01640/jm/domain_insertion/dipSeq/batch_scripts_JMA/scripts_paper/figures/RP_{title_text}.svg")
        plt.show()

    def forward(self): 
        self.calculate(self.X, self.y, self.data_split, self.classifier)
        self.plot(self.title_text)

def plot_auroc(classifier, data_split, X, y, title_text):
    '''
    Run classifier with cross-validation and plot ROC curves
    '''
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    colors = ['#202020', '#353535', '#505050', '#656565', '#808080']

    fig, ax = plt.subplots(figsize=(5,5))
    for i, (train, test) in enumerate(data_split.split(X, y)):
        classifier.fit(X[train], y.iloc[train])
        viz = RocCurveDisplay.from_estimator(
            classifier,
            X[test],
            y.iloc[test],
            name="ROC fold {}".format(i),
            alpha=0.7,
            lw=3,
            ax=ax,
            color=colors[i]
        )
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="#87001D",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=3,
        alpha=1,
    )
    plt.fill_between(mean_fpr, mean_tpr, alpha=.2, color="#87001D")
    ax.set(ylim=(0, 1.01), xlim=(-.02, 1))
    ax.set_aspect('equal')
    ax.legend(loc="lower right")
    plt.savefig(f"/work/projects/project01640/jm/domain_insertion/dipSeq/batch_scripts_JMA/scripts_paper/figures/AUROC_{title_text}.svg")

def plot_auroc_custom(classifier, data_split, groups, X, y, title_text):
    # Run classifier with cross-validation and plot ROC curves
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    fig, ax = plt.subplots()
    for i, (train, test) in enumerate(data_split.split(X, y, groups= groups)):
        classifier.fit(X[train], y.iloc[train])
        viz = RocCurveDisplay.from_estimator(
            classifier,
            X[test],
            y.iloc[test],
            name="ROC fold {}".format(i),
            alpha=0.5,
            lw=2,
            ax=ax,
        )
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)

    ax.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(
        mean_fpr,
        mean_tpr,
        color="b",
        label=r"Mean ROC (AUC = %0.2f $\pm$ %0.2f)" % (mean_auc, std_auc),
        lw=2,
        alpha=0.8,
    )

    ax.set(
        xlim=[0, 1.05],
        ylim=[0, 1.05],
        title=f"ROC: {title_text}",
    )
    ax.legend(loc="lower right")
    plt.show()

def metrics_test(classifier, X_train, y_train, X_test, y_test, title_text):
    # Run classifier with cross-validation and plot ROC curves
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    classifier.fit(X_train, y_train)
    fig, ax = plt.subplots(1,2,figsize=(11,5))
    viz = RocCurveDisplay.from_estimator(classifier, X_test, y_test, name=f'ROC testset', alpha=1, lw=3, ax=ax[0], color='#87001D')
    ax[0].fill_between(viz.fpr, viz.tpr, alpha=.2, color="#87001D")
    ax[0].set(ylim=(0, 1.01), xlim=(-.02, 1))
    ax[0].set_aspect('equal')
    ax[0].legend(loc="lower right")
    
    viz = PrecisionRecallDisplay.from_estimator(classifier, X_test, y_test, name=f'AP fold', alpha=1, lw=3, ax=ax[1], color='#7582FF')
    ax[1].fill_between(viz.recall, viz.precision, step='post', alpha=.2, color='#7582FF')
    ax[1].set(ylim=(0, 1.01), xlim=(0, 1))
    ax[1].set_aspect('equal')
    ax[1].legend(loc="lower right")
    plt.savefig(f"/work/projects/project01640/jm/domain_insertion/dipSeq/batch_scripts_JMA/scripts_paper/figures/metrics_test_set{title_text}.svg")

def collect_training_data(radius, alignment, name, feat_dict_radius, feat_dict_sequential, primary_proteins, analysis_dict, wd):
    parser = PDBParser()
    protein = {}
    protein['structure'] = parser.get_structure(name, f'{wd}/{name}_AF.pdb')
    protein = find_proximal_AAs(protein, radius)
    protein['Sec_structure'], protein['all_dssp'] = get_Sec_struct(protein['structure'], f'{name}_AF', wd)
    protein['ASA'] = get_surface_accessibility(protein['all_dssp'], protein, window=radius, property_dict=protein)
    protein['pLDDT'] = get_pLDDT(protein['structure'], window=radius, property_dict=protein)
    seq = primary_proteins[name]
    # trim the uniprot sequence to the sequence of the structure
    trimmed_seq = [seq[r.id[1]-1] for r in protein['structure'].get_residues()]
    
    protein = apply_features(feat_dict_sequential, trimmed_seq, protein, window=radius, prox_mode=False)
    protein = apply_features(feat_dict_radius, trimmed_seq, protein, window=radius, prox_mode=True)
    del protein['structure'], protein['all_dssp']
    try:
        protein = pd.DataFrame.from_dict(protein)
    except:
        for i, j in protein.items():
            print(i, len(j))
    protein.index = [i for i in range(1, len(protein.index)+1)]

    #print output dataframes
    for i in analysis_dict.items():
        if i[0] == f'{name}_PDZ':
            #select samples from second enrichment and average them, if 2 replicates are available
            temp_dict = {}
            for idx, frame in i[1].items():
                if len(idx) == 2 and idx[1] == '2':
                    temp_dict[idx] = frame
            try: 
                out_df = (temp_dict['12']['log'] + temp_dict['22']['log'])/2
            except:
                try:
                    out_df = temp_dict['12']['log']
                except: 
                    out_df = temp_dict['22']['log']
            out_df = out_df.to_frame()
            complete_df = pd.concat([out_df, protein], axis=1)
            plot_df = complete_df[['log', 'Sec_structure', 'ASA', 'pLDDT']]   
            plot_df = pd.melt(plot_df, id_vars=['Sec_structure', 'ASA', 'pLDDT'], value_vars=out_df.columns, value_name='enrichment')
            plot_df['ASA'] =plot_df['ASA'].astype('float')
            plot_df['pLDDT'] =plot_df['pLDDT'].astype('float')
            protein = complete_df


    # append the alignment features to property dict
    protein = pd.concat([protein.reset_index(), alignment], axis=1)
    KLD_smooth=[]
    for position in protein['prox_AAs']:
        temp_list = [protein['KLD'][AA] for AA in position]
        KLD_smooth.append(sum(temp_list)/len(temp_list))
    protein['KLD'] = pd.Series(KLD_smooth)
    del protein['index']


    return protein
