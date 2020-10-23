from collections import OrderedDict
import matplotlib.pyplot as plt

TOTAL_LEN_GENOME = 3095677412

CHANNELS = ['A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'A[C>G]A', 'A[C>G]C',
            'A[C>G]G', 'A[C>G]T', 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T',
            'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'A[T>C]A', 'A[T>C]C',
            'A[T>C]G', 'A[T>C]T', 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T',
            'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'C[C>G]A', 'C[C>G]C',
            'C[C>G]G', 'C[C>G]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T',
            'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'C[T>C]A', 'C[T>C]C',
            'C[T>C]G', 'C[T>C]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T',
            'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'G[C>G]A', 'G[C>G]C',
            'G[C>G]G', 'G[C>G]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T',
            'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'G[T>C]A', 'G[T>C]C',
            'G[T>C]G', 'G[T>C]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T',
            'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'T[C>G]A', 'T[C>G]C',
            'T[C>G]G', 'T[C>G]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T',
            'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'T[T>C]A', 'T[T>C]C',
            'T[T>C]G', 'T[T>C]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T']

COLORS_SIGNATURES = {'SBS1':'#d35f5fff',
          'SBS2':'#b2abd2',
          'SBS5':'#5da394ff',
          'SBS6':'#dfc27d',
          'SBS9':'#3288bd',
          'SBS13':'#ffeda0',
          'SBS17a':'#678796',
          'SBS17b':'#bf812d',
          'SBS18':'#6c535dff',
          'SBS32':'#decd87ff',
          'SBS34':'#ff0000ff',
          'SBS36':'#66bd63',
          'SBS37':'#f1b6da',
          'unknown':'#e0e0e0',
          'SBS_hscp':'#ffd5d5ff',
          'SBSA_new':'#82a1e2',
          'SBSB_new':'#6a51a3'}

COLORS_AGES = {'0-9':"#e6e6e6ff",
                   '10-15':"#999999ff",
                   '16-20':"#666666ff",
                   '21-39':"#4d4d4dff",
                   '>40':"#333333ff"}

COLORS_AGES_TALL = {'0-9':"#ffffd9",
                   '10-15':"#c7e9b4",
                   '16-20':"#41b6c4",
                   '21-39':"#0570b0",
                   '>40':"#081d58"}

COLORS_SUBTYPES = {'TALL Adult':'#d50402',
                         "TALL Pediatric":'#ff8080',
                         "BALL Pediatric":'#ffb380ff',
                         'DUX4-ERG':'#8a6f91ff',
                         'Other':'#993404',
                         'Hypodiploid':'#f1b6da',
                         'Infant MLL-R':'#4393c3',
                         'Hyperdiploid':'#d38d5fff',
                         'Ph positive':'#a6bddb',
                         'Ph-like':'#9aaf6b',
                         'iAMP21':'#decd87ff'}

COLORS_COHORTS = {'ADULT TALL AECC PROJECT':'#d50402',
                'PEDIATRIC TALL WXS (Oshima et al., 2016; PNAS)':'#ff8080',
                'PEDIATRIC ALL (Li et al., 2019, Blood)':'#ffbfbfff'}

COLORS_IMMUNOPHENO = {'ETP':'#afdde9ff',
                          'Pre-T':'#e9c6afff',
                          'Cortical':'#e9ddafff',
                          'Madura':"#c8beb7ff"}

# We sequenced the patients at different moments in time so each batch had its own folder
# with the corresponding patient calls. That's why in many scripts this dictionary is imported.
# the paths added in this dictionary should not include de PAT* folder, just the parent level above it
PATS_DIRS = OrderedDict([('PAT1',
              ''),
             ('PAT2',
              ''),
             ('PAT3',
              ''),
             ('PAT4',
              ''),
             ('PAT5',
              ''),
             ('PAT6',
              ''),
             ('PAT7',
              ''),
             ('PAT8',
              ''),
             ('PAT9',
              ''),
             ('PAT10',
              ''),
             ('PAT11',
              ''),
             ('PAT12',
              ''),
             ('PAT13',
              ''),
             ('PAT14',
              ''),
             ('PAT15',
              ''),
             ('PAT16',
              ''),
             ('PAT17',
              ''),
             ('PAT18',
              ''),
             ('PAT19',
              '')])

# SET EQUAL PLOT PARAMETERS
commonFontsize = 12

def config_rcparams():
    plt.rcParams['font.family'] = ['sans-serif']
    plt.rcParams['font.sans-serif'] = ['arial']
    plt.rcParams['font.size'] = commonFontsize
    plt.rcParams['axes.labelsize'] = commonFontsize
    plt.rcParams['xtick.labelsize'] = commonFontsize
    plt.rcParams['ytick.labelsize'] = commonFontsize
    plt.rcParams['axes.titlesize'] = 14
    # plt.rcParams['text'] = commonFontsize
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'arial'
    plt.rcParams['mathtext.rm'] = 'arial'
    plt.rcParams['axes.edgecolor'] = 'black'
    plt.rcParams['ytick.color'] = 'black'
    plt.rcParams['axes.linewidth'] = 0.7
    plt.rcParams['xtick.major.width'] = 0.7
    plt.rcParams['ytick.major.width'] = 0.7
    plt.rcParams['xtick.major.size'] = 3
    plt.rcParams['ytick.major.size'] = 3
