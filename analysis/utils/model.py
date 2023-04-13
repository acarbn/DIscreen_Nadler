#import packages
from utils.plotting import *
from utils.processing import *
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import cross_validate
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import StratifiedKFold
from sklearn import ensemble


# set styles
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


class Classifier:
    '''
    trains a gradient boosting classifier
    protein_seq: the sequenc of the candidate protein, which will be one-hot encoded
    property_dict: a dict containing different protein properties that should be taken into account
    The data is scaled prior to the model fit
    '''
    
    def __init__(self, property_dict, n_cross_val, title_text=''):
        self.title_text = title_text
        self.scores = ['recall', 'roc_auc', 'average_precision']
        self.dataset = property_dict
        self.n_cross_val = n_cross_val
        self.min_max = MinMaxScaler()
        self.clf = ensemble.GradientBoostingClassifier(n_estimators=100, criterion='squared_error', learning_rate=0.1, max_depth=4, max_features='auto', loss='exponential', random_state=42) 
        self.cross_val = StratifiedKFold(n_splits=n_cross_val)

    def encode_AA_sequences(self, protein_seq):
        one_hot_enc, AA_to_int = one_hot_insertion(protein_seq)
        self.one_hot_enc = pd.DataFrame(one_hot_enc)
        self.one_hot_enc.columns = [*AA_to_int.keys()]
        self.dataset = pd.concat([self.one_hot_enc, self.dataset.reset_index()], axis=1)
        self.dataset.drop(columns=['index'], inplace=True)

    def encode_sec_structures(self):
        one_hot_sec, AA_to_int = one_hot_insertion(self.dataset['Sec_structure'], alphabet='CHB')
        one_hot_sec = pd.DataFrame(one_hot_sec)
        one_hot_sec.columns = ['Coil', 'Sheet', 'Helix']
        self.dataset = pd.concat([one_hot_sec.reset_index(), self.dataset.reset_index()], axis=1)
        self.dataset.drop(columns=['index', 'Sec_structure'], inplace=True)
      
    def prepare_datasets(self):
        self.dataset = self.dataset.dropna()
        self.enrichment = self.dataset['log'].copy()
        self.dataset= self.dataset.drop(columns=['log'])
        # self.dataset.drop(columns=['log'], inplace=True)
        
        # Train-test-split
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.dataset, self.enrichment.squeeze(), test_size=0.2, random_state=42)
        
        # Feature Scaling so that the numbers in different categories are similar
        self.sc = self.min_max
        self.X_train = self.sc.fit_transform(self.X_train)
        self.X_test = self.sc.transform(self.X_test)

    def grid_search(self, tuned_params, model):     
        try:
            self.prepare_datasets()
        except:
            pass    
        for score in self.scores:
            print("# Tuning hyper-parameters for %s" % score)
            clf = GridSearchCV(model, tuned_params, scoring=score)
            clf.fit(self.X_train, self.y_train)
            print("Best parameters set found on development set:")
            print(clf.best_params_)
            print("Grid scores on development set:")
            means = clf.cv_results_["mean_test_score"]
            stds = clf.cv_results_["std_test_score"]
            for mean, std, params in zip(means, stds, clf.cv_results_["params"]):
                print(f"{mean} (+/-{std}) for {params}")
        
    def split_and_train(self):
        self.prepare_datasets()  
        # Fit the model
        self.clf.fit(self.X_train, self.y_train)
    
        # n-fold cross-validation
        self.kf = self.cross_val
        scores = cross_validate(self.clf, self.X_train, self.y_train, cv=self.kf, scoring=self.scores, return_estimator=True) 
        return self.clf, self.X_test, self.y_test, self.X_train, self.y_train, scores, self.kf, self.sc
    
    def metrics(self):
        plot_auroc(self.clf, data_split=self.kf, X=self.X_train, y=self.y_train, title_text=self.title_text)
        plot_prec_rec(self.clf, data_split=self.kf, X=self.X_train, y=self.y_train, title_text=self.title_text, groups=None).forward()