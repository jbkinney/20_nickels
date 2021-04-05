import numpy as np
import pandas as pd
import pdb
import time

# For converting an alignment to a ct_df matrix
def x_to_ohe(x,
             alphabet,
             check_seqs=True,
             check_alphabet=True,
             ravel_seqs=True):
    """
    Convert a sequence array to a one-hot encoded matrix.

    Parameters
    ----------
    x: (np.ndarray)
        (N,) array of input sequences, each of length L

    alphabet: (np.ndarray)
        (C,) array describing the alphabet sequences are drawn from.

    check_seqs: (bool)
        Whether to validate the sequences

    check_alphabet: (bool)
        Whether to validate the alphabet

    ravel_seqs: (bool)
        Whether to return an (N, L*C) array, as opposed to an (N, L, C) array.

    Returns
    -------
    x_ohe: (np.ndarray)
        Array of one-hot encoded sequences, stored as np.int8 values.
    """
    # Get dimensions
    L = len(x[0])
    N = len(x)
    C = len(alphabet)

    # Shape sequences as array of int8s
    x_arr = np.frombuffer(bytes(''.join(x), 'utf-8'),
                          np.int8, N * L).reshape([N, L])

    # Create alphabet as array of int8s
    alphabet_arr = np.frombuffer(bytes(''.join(alphabet), 'utf-8'),
                                 np.int8, C)

    # Compute (N,L,C) grid of one-hot encoded values
    x_nlc = (x_arr[:, :, np.newaxis] ==
             alphabet_arr[np.newaxis, np.newaxis, :]).astype(np.int8)

    # Ravel if requested
    if ravel_seqs:
        x_ohe = x_nlc.reshape([N, L * C])
    else:
        x_ohe = x_nlc

    return x_ohe
def x_to_stats(x, weights=None, verbose=False):
    """
    Identify the consensus sequence from a sequence alignment.

    Parameters
    ----------
    x: (np.ndarray)
        List of sequences. Sequences must all be the same length.
        
    weights: (None, np.ndarray)
        Weights for each sequence. E.g., count values, or numerical y values.
        If None, a value of 1 will be assumed for each sequence.

    verbose: (bool)
        Whether to print computation time.

    Returns
    -------
    consensus_seq: (str)
        Consensus sequence.
    """
    # Start timer
    start_time = time.time()

    # Validate alphabet
    alphabet = np.array(list('ACGT'))

    # Check weights and set if not provided
    if weights is None:
        weights = np.ones(len(x))
    else:
        weights = weights.astype(float)
        assert len(weights) == len(x),\
              f"len(weights)={len(weights)} does not match len(x)={len(x)}" 

    # Do one-hot encoding of sequences
    t = time.time()
    x_nlc = x_to_ohe(x,
                     alphabet,
                     check_seqs=False,
                     check_alphabet=False,
                     ravel_seqs=False)
    #print(f'Time for x_to_ohe: {time.time()-t:.3f} sec.')
    N, L, C = x_nlc.shape

    # Dictionary to hold results
    stats = {}

    # Compute x_ohe
    stats['x_ohe'] = x_nlc.reshape([N, L*C]).astype(np.int8)

    # Multiply by weights
    x_nlc = x_nlc.astype(float) * weights[:, np.newaxis, np.newaxis]

    # Compute lc encoding of consensus sequence
    x_sum_lc = x_nlc.sum(axis=0)
    x_sum_lc = x_sum_lc.reshape([L, C])
    x_support_lc = (x_sum_lc != 0)

    # Set number of sequences
    stats['N'] = N

    # Set sequence length
    stats['L'] = L

    # Set number of characters
    stats['C'] = C

    # Set alphabet
    stats['alphabet'] = alphabet

    # Compute probability matrix
    p_lc = x_sum_lc / x_sum_lc.sum(axis=1)[:, np.newaxis]
    stats['probability_df'] = pd.DataFrame(index=range(L),
                                           columns=alphabet,
                                           data=p_lc)

    # Compute sparsity factor
    stats['sparsity_factor'] = (x_nlc != 0).sum().sum() / x_nlc.size

    # Compute the consensus sequence and corresponding matrix.
    # Adding noise prevents ties
    x_sum_lc += 1E-1 * np.random.rand(*x_sum_lc.shape)
    stats['consensus_seq'] = \
        ''.join([alphabet[np.argmax(x_sum_lc[l, :])] for l in range(L)])

    # Compute mask dict
    missing_dict = {}
    for l in range(L):
        if any(~x_support_lc[l, :]):
            try:
                missing_dict[l] = ''.join(alphabet[~x_support_lc[l, :]])
            except:
                pdb.set_trace()
    stats['missing_char_dict'] = missing_dict

    # Provide feedback if requested
    duration_time = time.time() - start_time
    if verbose:
        print(f'Stats computation time: {duration_time:.5f} sec.')

    return stats
def x_to_ct_df(x, ct, N_max=100000):
    
    # Remove zero rows
    x = np.array(x)[ct != 0]
    ct = np.array(ct)[ct != 0]
    
    # Randomly sample without replacement if N is set and is less then len(x)
    if N_max < len(x):
        ix = np.random.choice(a=len(x), size=N_max, replace=False)
        x = x[ix].copy()
        ct = ct[ix].copy()
    
    # Comptue stats
    stats = x_to_stats(x, weights=ct)
    
    # Compute ct_df and return
    ct_df = (stats['probability_df'] * stats['N']).astype(int)
    return ct_df