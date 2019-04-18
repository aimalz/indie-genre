import collections
import sys

eps = sys.float_info.epsilon

OPop = collections.namedtuple('OPop', ['frac', 'type', 'pars'])

bconst = -0.003
sconst = 0.02
oconst = 0.1


def ingest(in_info):
    """
    Function reading in parameter file to define functions necessary for
    generation of posterior probability distributions

    Parameters
    ----------
    in_info: string or dict
        string containing path to plaintext input file or dict containing
        likelihood input parameters

    Returns
    -------
    in_dict: dict
        dict containing keys and values necessary for posterior probability
        distributions
    """
    if type(in_info) == str:
        with open(in_info) as infile:
            lines = (line.split(None) for line in infile)
            in_dict = {defn[0] : defn[1:] for defn in lines}
    else:
        in_dict = in_info
    return in_dict

def check_sim_params(params={}):
    """
    Checks simulation parameter dictionary for various keywords and sets to
    default values if not present

    Parameters
    ----------
    params: dict, optional
        dictionary containing initial key/value pairs for simulation of catalog

    Returns
    -------
    params: dict
        dictionary containing final key/value pairs for simulation of catalog
    """
    params = check_basic_setup(params)
    params = check_bias_params(params)
    params = check_variable_sigmas(params)
    params = check_catastrophic_outliers(params)
    return params
