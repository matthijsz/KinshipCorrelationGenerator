try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import statsmodels.formula.api as smf
    from scipy import optimize
    from scipy.stats import rankdata
except ImportError:
    print('Failed to import required libraries. Please run the following:')
    print(' pip install --user pandas numpy matplotlib json scipy')

import json
import pickle
import time
import argparse
import os

upper_boundary = 110
lower_boundary = 0
check_cutoff_drops = False
seed = 1415926536
explore_plot = False
save_separate_data = False
parallel = False

__version__ = '1.0.0'

if parallel:
    import multiprocessing as mp
    from joblib import Parallel, delayed


def add_to_to_do(to_do_, cortype, val_1, val_2, fid):
    to_do_[cortype][0] += [val_1]
    to_do_[cortype][1] += [val_2]
    to_do_[cortype]['FamID'] += [fid]
    return to_do_


def combine_fam_ids(fid1, fid2, fid_mult):
    if fid1 >= fid2:
        return fid1*fid_mult+fid2
    else:
        return fid2*fid_mult+fid1


def _reformat_step_0(to_do, idxs, pedigree, all_IDs, parallel):
    for idx in idxs:
        if not parallel:
            print('Reformatting pedigree: Parentoffspring: {0}/{1}'.format(idx, len(pedigree)), end='\r')
        dat = pedigree.loc[idx, :]
        if not dat['Father'] == '':
            FID = dat['Father']
            if FID in all_IDs:
                if dat['Gender'] == 'M':
                    to_do = add_to_to_do(to_do, 'FatherSon', FID, dat['ID'], dat['!FamID'])
                elif dat['Gender'] == 'F':
                    to_do = add_to_to_do(to_do, 'FatherDaughter', FID, dat['ID'], dat['!FamID'])
        if not dat['Mother'] == '':
            MID = dat['Mother']
            if MID in all_IDs:
                if dat['Gender'] == 'M':
                    to_do = add_to_to_do(to_do, 'MotherSon', MID, dat['ID'], dat['!FamID'])
                elif dat['Gender'] == 'F':
                    to_do = add_to_to_do(to_do, 'MotherDaughter', MID, dat['ID'], dat['!FamID'])
    return to_do


def _reformat_step_1(to_do, code, nrs, mx, pedigree, parallel):
    for nr in nrs:
        if not parallel:
            print('Reformatting pedigree: {0}: {1}/{2}'.format(code, nr, mx), end='\r')
        dat = pedigree.loc[pedigree[code] == nr, :].copy()
        if (code == 'Twincode') & (len(dat) == 2):
            if dat['Gender'].tolist() == ['M', 'M']:
                to_do = add_to_to_do(to_do, 'MZM', dat['ID'].tolist()[0], dat['ID'].tolist()[1],
                                     dat['!FamID'].tolist()[0])
            elif dat['Gender'].tolist() == ['F', 'F']:
                to_do = add_to_to_do(to_do, 'MZF', dat['ID'].tolist()[0], dat['ID'].tolist()[1],
                                     dat['!FamID'].tolist()[0])
        elif (code == 'DZtwincode') & (len(dat) == 2):
            if dat['Gender'].tolist() == ['M', 'M']:
                to_do = add_to_to_do(to_do, 'DZM', dat['ID'].tolist()[0], dat['ID'].tolist()[1],
                                     dat['!FamID'].tolist()[0])
            elif dat['Gender'].tolist() == ['F', 'F']:
                to_do = add_to_to_do(to_do, 'DZF', dat['ID'].tolist()[0], dat['ID'].tolist()[1],
                                     dat['!FamID'].tolist()[0])
            elif dat['Gender'].tolist() == ['F', 'M']:
                to_do = add_to_to_do(to_do, 'DOS', dat['ID'].tolist()[1], dat['ID'].tolist()[0],
                                     dat['!FamID'].tolist()[0])
            elif dat['Gender'].tolist() == ['M', 'F']:
                to_do = add_to_to_do(to_do, 'DOS', dat['ID'].tolist()[0], dat['ID'].tolist()[1],
                                     dat['!FamID'].tolist()[0])
        elif code == 'SibHousehold2':
            tw_code_idxs = []
            for tw_code in list(set([x for x in dat['Twincode'].tolist() if x > -1])):
                tw_code_idxs += [list(dat.loc[dat['Twincode'] == tw_code, :].index)]
            for tw_pair_idxs in tw_code_idxs:
                dat.drop(axis=0, labels=np.random.choice(tw_pair_idxs, 1), inplace=True)
            if len(dat) > 1:
                females = dat.loc[dat['Gender'] == 'F', :].copy()
                males = dat.loc[dat['Gender'] == 'M', :].copy()
                if (len(females) > 0) & (len(males) > 0):
                    for female_fis in females['ID'].tolist():
                        for male_fis in males['ID'].tolist():
                            to_do = add_to_to_do(to_do, 'BrotherSister', male_fis, female_fis,
                                                 dat['!FamID'].tolist()[0])
                elif len(females) > 1:
                    for n1, fis_1 in enumerate(females['ID'].tolist()):
                        for n2, fis_2 in enumerate(females['ID'].tolist()):
                            if n2 > n1:
                                to_do = add_to_to_do(to_do, 'SisterSister', fis_1, fis_2, dat['!FamID'].tolist()[0])
                elif len(males) > 1:
                    for n1, fis_1 in enumerate(males['ID'].tolist()):
                        for n2, fis_2 in enumerate(males['ID'].tolist()):
                            if n2 > n1:
                                to_do = add_to_to_do(to_do, 'BrotherBrother', fis_1, fis_2, dat['!FamID'].tolist()[0])
        elif (code == 'SpouseHousehold3') & (len(dat) == 2):
            if ('M' in dat['Gender'].tolist()) and ('F' in dat['Gender'].tolist()):
                to_do = add_to_to_do(to_do, 'Spouse', dat.loc[dat['Gender'] == 'M', 'ID'].tolist()[0],
                                     dat.loc[dat['Gender'] == 'F', 'ID'].tolist()[0], dat['!FamID'].tolist()[0])
    return to_do


def _reformat_step_2(pedigree, loop, brothers, sisters, othersibs, parallel):
    fullsibs = {}
    gender_col = [n for n, x in enumerate(pedigree.columns) if x == 'Gender'][0]
    for n, fis in loop:
        if not parallel:
            print('Reformatting for extended pedigree: {0}/{1}'.format(n, len(pedigree) * 2), end='\r')
        fullsibs[fis] = {'sex': pedigree.iloc[int(n), gender_col]}
        if fullsibs[fis]['sex'] == 'M':
            fullsibs[fis]['brothers'] = [[y for y in x if y != fis][0] for x in brothers if fis in x]
            fullsibs[fis]['sisters'] = [[y for y in x if y != fis][0] for x in othersibs if fis in x]
        elif fullsibs[fis]['sex'] == 'F':
            fullsibs[fis]['sisters'] = [[y for y in x if y != fis][0] for x in sisters if fis in x]
            fullsibs[fis]['brothers'] = [[y for y in x if y != fis][0] for x in othersibs if fis in x]
    return fullsibs


def _reformat_step_3(to_do2, pedigree, loop, fullsibs, father_col, mother_col, famid_col, famid_mult, parallel):
    for n, fis in loop:
        if not parallel:
            print('Reformatting for extended pedigree: {0}/{1}'.format(n + len(pedigree), len(pedigree) * 2),
                  end='\r')
        aunts, uncles = [], []
        father = pedigree.iloc[n, father_col]
        if isinstance(father, str):
            if father in fullsibs.keys():
                aunts += fullsibs[father]['sisters']
                uncles += fullsibs[father]['brothers']
        mother = pedigree.iloc[n, mother_col]
        if isinstance(mother, str):
            if mother in fullsibs.keys():
                aunts += fullsibs[mother]['sisters']
                uncles += fullsibs[mother]['brothers']
        if fullsibs[fis]['sex'] == 'M':
            for x in aunts:
                fid1 = pedigree.iloc[n, famid_col]
                fid2 = int(pedigree.loc[pedigree['ID'] == x, '!FamID'])
                to_do2 = add_to_to_do(to_do2, 'AuntNephew', x, fis, combine_fam_ids(fid1, fid2, famid_mult))
            for x in uncles:
                fid1 = pedigree.iloc[n, famid_col]
                fid2 = int(pedigree.loc[pedigree['ID'] == x, '!FamID'])
                to_do2 = add_to_to_do(to_do2, 'UncleNephew', x, fis, combine_fam_ids(fid1, fid2, famid_mult))
        elif fullsibs[fis]['sex'] == 'F':
            for x in aunts:
                fid1 = pedigree.iloc[n, famid_col]
                fid2 = int(pedigree.loc[pedigree['ID'] == x, '!FamID'])
                to_do2 = add_to_to_do(to_do2, 'AuntNiece', x, fis, combine_fam_ids(fid1, fid2, famid_mult))
            for x in uncles:
                fid1 = pedigree.iloc[n, famid_col]
                fid2 = int(pedigree.loc[pedigree['ID'] == x, '!FamID'])
                to_do2 = add_to_to_do(to_do2, 'UncleNiece', x, fis, combine_fam_ids(fid1, fid2, famid_mult))
        for rel in aunts + uncles:
            rel_children_idx = None
            if rel in pedigree['Mother'].tolist():
                rel_children_idx = pedigree.loc[pedigree['Mother'] == rel].index
            if rel in pedigree['Father'].tolist():
                rel_children_idx = pedigree.loc[pedigree['Father'] == rel].index
            if rel_children_idx is not None:
                for idx in rel_children_idx:
                    fid1 = pedigree.iloc[n, famid_col]
                    fid2 = pedigree.iloc[idx, famid_col]
                    if pedigree.loc[idx, 'Gender'] == 'M':
                        if fullsibs[fis]['sex'] == 'M':
                            to_do2 = add_to_to_do(to_do2, 'NephewNephew', fis, pedigree.loc[idx, 'ID'],
                                                  combine_fam_ids(fid1, fid2, famid_mult))
                        if fullsibs[fis]['sex'] == 'F':
                            to_do2 = add_to_to_do(to_do2, 'NephewNiece', pedigree.loc[idx, 'ID'], fis,
                                                  combine_fam_ids(fid1, fid2, famid_mult))
                    if pedigree.loc[idx, 'Gender'] == 'F':
                        if fullsibs[fis]['sex'] == 'M':
                            to_do2 = add_to_to_do(to_do2, 'NephewNiece', fis, pedigree.loc[idx, 'ID'],
                                                  combine_fam_ids(fid1, fid2, famid_mult))
                        if fullsibs[fis]['sex'] == 'F':
                            to_do2 = add_to_to_do(to_do2, 'NieceNiece', fis, pedigree.loc[idx, 'ID'],
                                                  combine_fam_ids(fid1, fid2, famid_mult))
    return to_do2


def reformat_pedigree(pedigreefilepath='pedigree-hh-all.ped', make_extended=False, parallel=False):
    if parallel:
        print('When reformatting pedigree in parallel progress prints are disabled, but rest assured it still works fine.')
    pedigree = pd.read_csv(pedigreefilepath, low_memory=False)
    pedigree['SibHousehold2'] = pedigree['SibHousehold2'].astype(str)
    pedigree.loc[pedigree['SibHousehold2'] == 'x1', 'SibHousehold2'] = '100000'
    pedigree.loc[pedigree['SibHousehold2'] == 'x2', 'SibHousehold2'] = '100001'
    pedigree.loc[pedigree['SibHousehold2'] == 'x3', 'SibHousehold2'] = '100002'
    pedigree.loc[pedigree['SibHousehold2'] == 'x4', 'SibHousehold2'] = '100003'
    pedigree['SibHousehold2'] = pedigree['SibHousehold2'].astype(float)

    to_do = {'Spouse': {0: [], 1: [], 'FamID': []},
             'Spouse-parents': {0: [], 1: [], 'FamID': []},
             'Spouse-twins': {0: [], 1: [], 'FamID': []},
                 'MZM': {0: [], 1: [], 'FamID': []},
                 'MZF': {0: [], 1: [], 'FamID': []},
                 'DZM': {0: [], 1: [], 'FamID': []},
                 'DZF': {0: [], 1: [], 'FamID': []},
                 'DOS': {0: [], 1: [], 'FamID': []},
                 'MotherDaughter': {0: [], 1: [], 'FamID': []},
                 'MotherSon': {0: [], 1: [], 'FamID': []},
                 'FatherDaughter': {0: [], 1: [], 'FamID': []},
                 'FatherSon': {0: [], 1: [], 'FamID': []},
                 'BrotherBrother': {0: [], 1: [], 'FamID': []},
                 'BrotherSister': {0: [], 1: [], 'FamID': []},
                 'SisterSister': {0: [], 1: [], 'FamID': []}}
    for k in ['Father', 'Mother', 'ID']:
        pedigree[k] = pedigree[k].astype(str)
    all_IDs = pedigree['ID'].tolist()
    if not parallel:
        to_do = _reformat_step_0(to_do, pedigree.index, pedigree, all_IDs, parallel)
    else:
        chunks = [list(x) for x in np.array_split(list(pedigree.index), mp.cpu_count()-1)]
        results = Parallel(n_jobs=mp.cpu_count()-1)(delayed(_reformat_step_0)(to_do, chunk, pedigree, all_IDs, parallel) for chunk in chunks)
        for r in results:
            for k_cortype, v_cortype in r.items():
                for k, v in v_cortype.items():
                    to_do[k_cortype][k] += v
    print('Reformatting pedigree: Parentoffspring: Finished')

    codes = ['Twincode', 'DZtwincode', 'SibHousehold2', 'SpouseHousehold3']
    for code in codes:
        nrs = pedigree[code].unique().tolist()
        mx = nrs[-1]
        if not parallel:
            to_do = _reformat_step_1(to_do, code, nrs, mx, pedigree, parallel)
        else:
            chunks = [list(x) for x in np.array_split(nrs, mp.cpu_count()-1)]
            results = Parallel(n_jobs=mp.cpu_count()-1)(delayed(_reformat_step_1)(to_do, code, chunk, mx, pedigree, parallel) for chunk in chunks)
            for r in results:
                for k_cortype, v_cortype in r.items():
                    for k, v in v_cortype.items():
                        to_do[k_cortype][k] += v
        print('Reformatting pedigree: {}: Finished                      '.format(code))


    spouse_fis_0 = to_do['Spouse'][0]
    spouse_fis_1 = to_do['Spouse'][1]
    twin_fis = []
    for t in ['MZM', 'MZF', 'DZM', 'DZF', 'DOS']:
        twin_fis += to_do[t][0] + to_do[t][1]
    twin_fis = list(set(twin_fis))

    spouse_twins_0_idx = np.where(pd.Series(spouse_fis_0).isin(twin_fis))[0]
    spouse_twins_1_idx = np.where(pd.Series(spouse_fis_1).isin(twin_fis))[0]
    spouse_twins_idx = list(set(list(spouse_twins_1_idx) + list(spouse_twins_0_idx)))
    to_do['Spouse-twins'][0] = list(np.array(spouse_fis_0)[spouse_twins_idx])
    to_do['Spouse-twins'][1] = list(np.array(spouse_fis_1)[spouse_twins_idx])
    to_do['Spouse-twins']['FamID'] = list(np.array(to_do['Spouse']['FamID'])[spouse_twins_idx])

    parent_fis = []
    for p in ['Father', 'Mother']:
        for p_bond in ['Daughter', 'Son']:
            parent_fis = list(set(parent_fis + to_do[p+p_bond][0]))
    pedigree_twins = pedigree.loc[pedigree['ID'].isin(twin_fis), :].copy()
    parent_fis_include = list(np.where(pd.Series(parent_fis).isin(pedigree_twins['Mother'].tolist()))[0])
    parent_fis_include += list(np.where(pd.Series(parent_fis).isin(pedigree_twins['Father'].tolist()))[0])
    parent_fis = list(np.array(parent_fis)[list(set(parent_fis_include))])

    spouse_parents_0_idx = np.where(pd.Series(spouse_fis_0).isin(parent_fis))[0]
    spouse_parents_1_idx = np.where(pd.Series(spouse_fis_1).isin(parent_fis))[0]
    spouse_parents_idx = list(set(list(spouse_parents_1_idx) + list(spouse_parents_0_idx)))
    to_do['Spouse-parents'][0] = list(np.array(spouse_fis_0)[spouse_parents_idx])
    to_do['Spouse-parents'][1] = list(np.array(spouse_fis_1)[spouse_parents_idx])
    to_do['Spouse-parents']['FamID'] = list(np.array(to_do['Spouse']['FamID'])[spouse_parents_idx])
    if len(to_do['Spouse-twins'][0]) == 0:
        to_do = {k: v for k, v in to_do.items() if k != 'Spouse-twins'}
    with open('reformatted_pedigree.pickle', 'wb') as f:
        pickle.dump(to_do, f)
    if make_extended:
        to_do2 = {
            'UncleNephew': {0: [], 1: [], 'FamID': []},
            'UncleNiece': {0: [], 1: [], 'FamID': []},
            'AuntNiece': {0: [], 1: [], 'FamID': []},
            'AuntNephew': {0: [], 1: [], 'FamID': []},
            'NephewNephew': {0: [], 1: [], 'FamID': []},
            'NieceNiece': {0: [], 1: [], 'FamID': []},
            'NephewNiece': {0: [], 1: [], 'FamID': []}
        }
        sisters = pd.DataFrame.from_dict(to_do['SisterSister'], orient='columns')[[0, 1]].values.tolist()
        brothers = pd.DataFrame.from_dict(to_do['BrotherBrother'], orient='columns')[[0, 1]].values.tolist()
        othersibs = pd.DataFrame.from_dict(to_do['BrotherSister'], orient='columns')[[0, 1]].values.tolist()

        if not parallel:
            fullsibs = _reformat_step_2(pedigree, enumerate(pedigree['ID'].tolist()), brothers, sisters, othersibs, parallel)
        else:
            fullsibs = {}
            chunks = [list(x) for x in np.array_split(list(enumerate(pedigree['ID'].tolist())), mp.cpu_count() - 1)]
            results = Parallel(n_jobs=mp.cpu_count()-1)(delayed(_reformat_step_2)(pedigree, chunk, brothers, sisters, othersibs, parallel) for chunk in chunks)
            for x in results:
                for fis, d in x.items():
                    fullsibs[fis] = d
            print('Reformatting for extended pedigree: 50%')
        father_col = [n for n, x in enumerate(pedigree.columns) if x == 'Father'][0]
        mother_col = [n for n, x in enumerate(pedigree.columns) if x == 'Mother'][0]
        famid_col = [n for n, x in enumerate(pedigree.columns) if x == '!FamID'][0]
        famid_mult = 10**(len(str(pedigree['!FamID'].max())))
        if not parallel:
            to_do2 = _reformat_step_3(to_do2, pedigree, enumerate(pedigree['ID'].tolist()), fullsibs, father_col, mother_col, famid_col, famid_mult, parallel)
        else:
            chunks = [list(x) for x in np.array_split(list(enumerate(pedigree['ID'].tolist())), mp.cpu_count() - 1)]
            results = Parallel(n_jobs=mp.cpu_count()-1)(delayed(_reformat_step_3)(to_do2, pedigree, chunk, fullsibs, father_col, mother_col, famid_col, famid_mult, parallel) for chunk in chunks)
            for r in results:
                for k_cortype, v_cortype in r.items():
                    for k, v in v_cortype.items():
                        to_do2[k_cortype][k] += v
        for k, v in to_do2.items():
            to_do[k] = v
        with open('../reformatted_extended_pedigree.pickle', 'wb') as f:
            pickle.dump(to_do, f)


def get_ped(extended):
    if extended:
        with open('reformatted_extended_pedigree.pickle', 'rb') as f:
            to_do = pickle.load(f)
    elif os.path.isfile('reformatted_pedigree.pickle'):
        raise FileNotFoundError("An extended pedigree was requested but no reformatted extended pedigree was found")
    else:
        if os.path.isfile('reformatted_pedigree.pickle'):
            with open('reformatted_pedigree.pickle', 'rb') as f:
                to_do = pickle.load(f)
        elif os.path.isfile('reformatted_extended_pedigree.pickle'):
            with open('reformatted_extended_pedigree.pickle', 'rb') as f:
                to_do2 = pickle.load(f)
            ks = ['Spouse', 'Spouse-parents', 'Spouse-twins', 'MZM', 'MZF', 'DZM', 'DZF', 'DOS', 'MotherDaughter',
                  'MotherSon', 'FatherDaughter', 'FatherSon', 'BrotherBrother', 'BrotherSister', 'SisterSister']
            to_do = {k: to_do2[k] for k in ks if k in to_do2.keys()}
        else:
            raise FileNotFoundError("No (extended) reformatted pedigree file found.")
    return to_do


def make_familybased_selection(datafile, pedigreefile=None,
                               outfileprefix='familyselected',
                               upper_boundary=110, lower_boundary=0, surveycompletionfile=None,
                               check_cutoff_drops=False):
    print('longitudinal: Reading data...')
    if pedigreefile is None:
        to_do = get_ped(False)
        cor_types = list(to_do.keys())
        ped = []
        for cortype in cor_types:
            output_dict = dict(ID=[str(x) for x in to_do[cortype][0]] + [str(x) for x in to_do[cortype][1]],
                               FamID=[int(x) for x in to_do[cortype]['FamID']]+[int(x) for x in to_do[cortype]['FamID']])
            ped.append(pd.DataFrame.from_dict(output_dict, orient='columns'))
        ped = pd.concat(ped, axis=0, ignore_index=True).drop_duplicates(subset='ID')
        ped.rename(columns={'FamID': '!FamID'}, inplace=True)
    else:
        ped = pd.read_csv(pedigreefile, low_memory=False)
        ped['ID'] = ped['ID'].astype(str)
    data = pd.read_csv(datafile)
    data['FISNumber'] = data['FISNumber'].astype(str)
    data = data.merge(ped[['ID', '!FamID']], left_on='FISNumber', right_on='ID', how='left')
    data = data.loc[data['age'] < 998, :]
    data['formerge'] = data['FISNumber'] + data['Source']
    if surveycompletionfile is not None:
        invjr = pd.read_csv(surveycompletionfile)
        invjr['FISNumber'] = invjr['FISNumber'].astype(str)
        invjr = invjr.merge(ped[['ID', '!FamID']], left_on='FISNumber', right_on='ID', how='left')
        invjr.sort_values(by='!FamID', inplace=True)
        invjr.drop(axis=1, labels='ID', inplace=True)
        invjr['formerge'] = invjr['FISNumber'] + invjr['Source']
        invjr = invjr.merge(data, on='formerge')
        data.drop(axis=1, labels=['formerge'], inplace=True)
        invjr.drop(axis=1, labels=['FISNumber_y', 'Source_y', 'formerge'], inplace=True)
        invjr.rename(columns={'FISNumber_x': 'FISNumber', 'Source_x': 'Source'}, inplace=True)
    else:
        data.drop(axis=1, labels=['formerge'], inplace=True)
        invjr = data
        invjr['invjr'] = data['age']

    if check_cutoff_drops:
        f = open('Age-cutoff_drops.txt', 'w')
        f.write('# N dropped is the number of unique individuals dropped from the full dataset.\n')
        f.write('# Subjects that have data above and below the cutoff are already excluded from this number.\n')
        f.write('cutoff\tN_dropped\n.')
        for upper_cutoff in range(60, 81):
            drop_list = list(set(invjr.loc[invjr['age'] > upper_cutoff, 'FISNumber'].tolist()))
            not_drop_list = list(set(invjr.loc[invjr['age'] <= upper_cutoff, 'FISNumber'].tolist()))
            drop_list = [x for x in drop_list if x not in not_drop_list]
            f.write('{0}\t{1}\n'.format(upper_cutoff, len(drop_list)))
        f.close()

    invjr = invjr.loc[invjr['age'] >= lower_boundary, :]
    invjr = invjr.loc[invjr['age'] <= upper_boundary, :]

    all_families = list(set(invjr.loc[invjr['!FamID'].notnull(), '!FamID'].tolist()))

    to_use = {}
    for n_fam, fam in enumerate(all_families):
        if n_fam % 250 == 0:
            print('longitudinal: Working on family {0} of {1}'.format(n_fam, len(all_families)), end='\r')
        csub = invjr.loc[invjr['!FamID'] == fam, :]
        ccounts = csub['invjr'].value_counts()
        ccounts = ccounts.sort_values(0, ascending=False)
        ccounts = ccounts.reset_index()
        ccounts.sort_values('index', inplace=True)
        ccounts.sort_values('invjr', inplace=True)
        best_year = ccounts.loc[0, 'index']
        best_set = csub.loc[csub['invjr'] == best_year, ['FISNumber', 'Source']]
        best_set_fis = best_set['FISNumber'].tolist()
        best_set_source = best_set['Source'].tolist()
        for x in range(len(best_set_fis)):
            to_use[best_set_fis[x]] = best_set_source[x]
        csub = csub.loc[~(csub['FISNumber'].isin(best_set_fis)), :]
        while len(csub) > 0:
            csub['delta'] = abs(csub['invjr'] - best_year)
            csub.sort_values('delta', inplace=True)
            csub.reset_index(inplace=True, drop=True)
            next_fis = csub.loc[0, 'FISNumber']
            next_source = csub.loc[0, 'Source']
            to_use[next_fis] = next_source
            csub = csub.loc[~(csub['FISNumber'].isin([next_fis])), :]

    print('longitudinal: Optimal lists found.                        ')

    with open('Lists_used_per_subject.json', 'w') as f:
        json.dump(to_use, f)

    selected_data = []
    n_loops = len(to_use)
    a = 1
    for FIS, SRC in to_use.items():
        if a % 500 == 0:
            print('Selecting subject {0} of {1}         '.format(a, n_loops), end='\r')
        _ = data.loc[data['FISNumber'] == FIS, :]
        _ = _.loc[_['Source'] == SRC, :]
        selected_data.append(_)
        a += 1
    selected_data = pd.concat(selected_data, axis=0, ignore_index=True)
    if outfileprefix is not None:
        selected_data.to_csv(outfileprefix+'_familyselected.csv', index=False)
    return selected_data


class WeightedCorr:
    # Implementation as described in https://files.eric.ed.gov/fulltext/ED585538.pdf
    def __init__(self, xyw=None):
        self.x, self.y, self.w = (pd.to_numeric(xyw[i], errors='coerce').values for i in xyw.columns)

    def _wcov(self, x, y, ms):
        return np.sum(self.w * (x - ms[0]) * (y - ms[1])) / np.sum(self.w)

    def _pearson(self, x=None, y=None):
        x, y = (self.x, self.y) if ((x is None) and (y is None)) else (x, y)
        mx, my = (np.sum(i * self.w) / np.sum(self.w) for i in [x, y])
        return self._wcov(x, y, [mx, my]) / np.sqrt(self._wcov(x, x, [mx, mx]) * self._wcov(y, y, [my, my]))

    def _wrank(self, x):
        (unique, arr_inv, counts) = np.unique(rankdata(x), return_counts=True, return_inverse=True)
        a = np.bincount(arr_inv, self.w)
        return (np.cumsum(a) - a)[arr_inv] + ((counts + 1) / 2 * (a / counts))[arr_inv]

    def _spearman(self, x=None, y=None):
        x, y = (self.x, self.y) if ((x is None) and (y is None)) else (x, y)
        return self._pearson(self._wrank(x), self._wrank(y))

    def __call__(self, method='pearson'):
        if method not in ['pearson', 'spearman']:
            raise ValueError('method should be one of [\'pearson\', \'spearman\']')
        cor = {'pearson': self._pearson, 'spearman': self._spearman}[method]
        return cor()


def make_cor_table(datafile='Family_selected_data.csv', seed=1415926536, explore_plot=False,
                   outfileprefix=None, save_separate_data=False, use_repeated_families=False, method='pearson',
                   correction='', exclude='', use_extended=False, randomsample=False, raw_n=False):
    phenotype = pd.read_csv(datafile)
    phenotype['FISNumber'] = phenotype['FISNumber'].astype(str)
    phenotype['FISNumber'] = phenotype['FISNumber'].str.replace('\.0', '')
    to_do = get_ped(use_extended)
    phenotype.columns = [x.replace('.', '_') for x in phenotype.columns]
    variables = [x for x in phenotype.columns if x not in ['FISNumber', 'sex', 'age', 'Source', 'index']]
    if exclude != '':
        variables = [x for x in variables if x not in ','.split(exclude)]
    if correction != '':
        print('Corercting for the following: {}'.format(correction))
        for v in variables:
            phenotype[v] = pd.to_numeric(phenotype[v], errors='coerce')
            res = smf.ols('{} ~ {}'.format(v, correction), data=phenotype).fit()
            phenotype[v] = res.resid
    results, resultsN, resultsNtot = {}, {}, {}
    for cortype in list(to_do.keys()):
        output_dict = dict(ID_0=[str(x) for x in to_do[cortype][0]], ID_1=[str(x) for x in to_do[cortype][1]],
                           FamID=[int(x) for x in to_do[cortype]['FamID']])
        output_df = pd.DataFrame.from_dict(output_dict, orient='columns').drop_duplicates()
        output_mrg = output_df.merge(phenotype, left_on='ID_0', right_on='FISNumber', how='left').drop(axis=1, labels='ID_0')
        rename_dict0, rename_dict1 = {}, {}
        for x in list(output_mrg.columns[2:]):
            rename_dict0[x] = x + '_0'
            rename_dict1[x] = x + '_1'
        output_mrg.rename(columns=rename_dict0, inplace=True)
        final = output_mrg.merge(phenotype, left_on='ID_1', right_on='FISNumber', how='left').drop(axis=1, labels='ID_1')
        final = final.rename(columns=rename_dict1).sort_values('FamID', ascending=False).reset_index()
        final = final.dropna(subset=[x+'_0' for x in variables]+[x+'_1' for x in variables], how='all')
        if randomsample:
            final = final.sample(frac=1, random_state=seed).drop_duplicates(subset='FamID')
        len_fin = len(final)
        resultsNtot[cortype] = dict(total=len_fin)
        resultsN[cortype] = {}
        if save_separate_data:
            final.to_csv('{0}_{1}_data.csv'.format(outfileprefix, cortype), index=False)
        results[cortype] = {}
        print('cor_table: Getting correlations: {}              '.format(cortype), end='\r')
        for n, i in enumerate(variables):
            final[i + '_0'] = pd.to_numeric(final[i + '_0'], errors='coerce')
            final[i + '_1'] = pd.to_numeric(final[i + '_1'], errors='coerce')
            if len(final[[i+'_0', i+'_1']].dropna()) < 30:
                results[cortype][i] = 'NaN'
                resultsN[cortype][i] = -1
                resultsNtot[cortype][i] = -1
            else:
                cdat = final[[i+'_0', i+'_1']].dropna(subset=[i + '_0', i + '_1'])
                resultsNtot[cortype][i] = len(cdat)
                if randomsample or use_repeated_families:
                    results[cortype][i] = cdat.corr(method=method).iloc[1, 0]
                else:
                    if method == 'kendall':
                        raise NotImplementedError('A weighted kendall tau correlation is not implemented (yet).')
                    cdat = final[['FISNumber_0', 'FISNumber_1', i+'_0', i+'_1']].dropna(subset=[i + '_0', i + '_1'])
                    (unique, arr_inv, counts) = np.unique(cdat[['FISNumber_0', 'FISNumber_1']].values, return_counts=True, return_inverse=True)
                    cdat['weight'] = (.5 / counts[arr_inv].reshape((len(cdat), 2))).sum(axis=1)
                    results[cortype][i] = WeightedCorr(cdat[[i+'_0', i+'_1', 'weight']])(method)
                    resultsN[cortype][i] = cdat['weight'].sum()
                if explore_plot:
                    fig = plt.figure(figsize=(10, 7.5), dpi=80, facecolor='w', edgecolor='k')
                    fig.add_subplot(1, 3, 1)
                    ax1 = final[i + '_0'].plot.hist(ylim=(0, int(len(final) * (2 / 3))), bins=15)
                    ax1.set_title('0: m:{0}, sd:{1} \nrange:[{2},{3}]'.format(
                        round(final[i + '_0'].mean(), 2), round(final[i + '_0'].std(), 2),
                        round(final[i + '_0'].min(), 2), round(final[i + '_0'].max(), 2)
                    ))
                    ax_3 = fig.add_subplot(1, 3, 2)
                    ax3 = final.plot(x=i + '_0', y=i + '_1', kind='scatter', ax=ax_3)
                    ax3.set_title('N={2}\n{0}: r:{1}'.format(cortype, round(results[cortype][i], 2), len_fin))
                    fig.add_subplot(1, 3, 3)
                    ax2 = final[i + '_1'].plot.hist(ylim=(0, int(len(final) * (2 / 3))), bins=15)
                    ax2.set_title('1: m:{0}, sd:{1}, \nrange:[{2},{3}]'.format(
                        round(final[i + '_1'].mean(), 2), round(final[i + '_1'].std(), 2),
                        round(final[i + '_1'].min(), 2),
                        round(final[i + '_1'].max(), 2),
                    ))
                    fig.savefig('{0}_Explore_{1}_{2}.png'.format(outfileprefix, cortype, i))
                    plt.close()
    results_df = pd.DataFrame.from_dict(results, orient='index')
    if randomsample or use_repeated_families:
        resultsn_df = pd.DataFrame.from_dict(resultsNtot, orient='index')
    else:
        resultsn_df = pd.DataFrame.from_dict(resultsN, orient='index')
        if raw_n:
            pd.DataFrame.from_dict(resultsNtot, orient='index').to_csv(outfileprefix + '_Fam_raw_N.csv')
    if outfileprefix is not None:
        results_df.to_csv(outfileprefix+'_Fam_Correlations.csv')
        resultsn_df.to_csv(outfileprefix+'_Fam_N.csv')
    return results_df, resultsn_df


def printtime(start_t, prefix='Analysis'):
    runtime = time.time() - start_t
    hrs, mins = divmod(runtime, 3600)
    mins, secs = divmod(mins, 60)
    if hrs > 0:
        print('{} finished after {:3d}:{:02d}:{:05.2f}'.format(prefix, int(hrs), int(mins), secs))
    elif mins > 0:
        print('{} finished after {:02d}:{:05.2f}'.format(prefix, int(mins), secs))
    else:
        print('{} finished after {:05.2f} seconds'.format(prefix, secs))


def morehelp(mhdict, k):
    if 'all' in k:
        mh = list(mhdict.keys())
    else:
        mh = k
    for x in mh:
        print('More help on {}:'.format(x))
        print('\n'.join(['  ' + i for i in mhdict[x]]))
        print()


if __name__ == '__main__':
    print('+----------------------------------------------+\n|         Kinship correlation generator        |')
    print('+----------------------------------------------+\n|    by Matthijs van der Zee & Eco de Geus     |')
    print('+----------------------------------------------+\n')
    parser = argparse.ArgumentParser(description='Python scripts to generate weighted kinship correlation table.'
                                                 'This script is made for Python3 and requires the following libraries:'
                                                 'pandas, numpy, matplotlib, json')
    parser.add_argument('--data', help='Path to the datafile. This datafile should at least contain:[\'FISNumber\', '
                                       '\'age\', \'sex\'] columns and 1 phenotype column. If you are using longitudinal '
                                       'data and want to use the included optimal-familybased selection it should also '
                                       'contain [] variables.', default=None)
    parser.add_argument('--outprefix', help='Prefix to use for output files.', default=None)
    parser.add_argument('--pedigree', help='Path to the pedigree file if a reformatted pedigree does not exist, or you'
                                           'want to create a new reformatted pedigree.', default=None)
    parser.add_argument('--longitudinal', help='When using longitudinal data add this argument. This will apply find an'
                                               'optimal survey to use for each subject based on the surveys completed by'
                                               'other family members. A surveycompletion file is advised here, and the datafile'
                                               'should contain columns : [\'FISNumber\', \'age\', \'sex\', \'index\', \'Source\']',
                        action='store_true', default=False)
    parser.add_argument('--extended', help='Use extended pedigrees so the output will unclude Aunt-Nephew, Niece-Nephew, etc. correlations.',
                        action='store_true', default=False)
    parser.add_argument('--surveycompletion', help='Datafile with survey completion years. Should have columns:['
                                                   '\'FISNumber\', \'Source\', \'invjr\']'),
    parser.add_argument('--method', help='Correlation method, should be one of [\'pearson\', \'spearman\']. Default is pearson.',
                        default='pearson')
    parser.add_argument('--correct', help='Formula to use for linear regression correction. Defaults to no corrections.',
                        default='')
    parser.add_argument('--exclude', help='Variables for which no correlation should be calculated, to be used mainly with custom covariates.',
                        default='')
    parser.add_argument('--randomsample', help='Use only 1 pair per family instead of weighting for multiple occurences of the same invididual.',
                        action='store_true', default=False)
    parser.add_argument('--use_repeated_families', help='Add this argument to include all participants in larger '
                                                        'famileies, i.e. don\'t drop or weight for duplicate samples within correlations.',
                        action='store_true', default=False)
    parser.add_argument('--raw_n', help='Store an additional csv file with the raw N samples used, in addtion to the weighted N file.',
                        action='store_true', default=False)
    parser.add_argument('--morehelp', nargs="+", help='Print more help about this script (\'all\') or specific arguments')
    args = parser.parse_args()
    morehelpdict = {
        'this script': ['This script is used intended to be used to create a kinship correlations table.',
                        'This script is made for Python3 and requires pandas, numpy, matplotlib, json',
                        'It will generate at least two files: prefix_Fam_Correlations.csv and _Fam_N.csv,',
                        'for correlations and sample sizes respectively (weighted by default, see --morehelp weights).',
                        'A reformatted pedigree is required for this script to work efficiently.',
                        'You can generate this once and use the script without specifying pedigree as often as you wish after.',
                        'Note: This script is NOT fast, if your run the full monty and have a large dataset it might take a while.',
                        '   That being said, generating correlations if you already have a reformatted pedigree should finish in seconds.'
                        'In the output: Spouse-parents are parents of twins, and Spouse-twins are twins and twin spouses that are not parents.'],
        'weights': ['Correlations and sample size outputs are weighted by the number of occurences of an individual within the full set',
                    'of that correlation. Here\'s a crude explanation of the process for 1 correlation type (DOS) and 1 phenotype(Ph):',
                    ' - Create a dataframe of all DOS pairs with phenotypic data (columns are FID, ID1, ID2, Ph1, Ph2',
                    ' - Assign base weights of 0.5 to weight1 and weight2.',
                    ' - Devide weight1 for row i by the number of occurences of ID1[i] in ID1 and ID2',
                    ' - Devide weight2 for row i by the number of occurences of ID2[i] in ID1 and ID2',
                    ' - weight = weight1 + weight2',
                    'The sum of this weight will be used as N, and to calculate weighted correlations. ',
                    'I\'ve published the method I use fore weighted correlations here: https://github.com/matthijsz/weightedcorr'],
        'data': ['This is your input datafile, it should be a comma separated .csv file with the following columns:',
                 ' - FISNumber (case sensitive): Personal identifier',
                 ' - age (case sensitive): age of the participant at time of survey completion',
                 ' - sex (case sensitive): sex of the participant 1=male, 2=female',
                 ' - PHENOTYPE (name is irrelivent): 1 or more columns with your phenotypic data.',
                 'It should only contain complete data or highly similar missing data in all phenotypes',
                 'If you have longitudinal data and want to use my method of selection detailed under longitudinal',
                 'this file should have 1 row per survey per subject and, in addition the columns above contain columns:',
                 ' - index: A within-person index numbering the survey, first completed survey is 1, second 2, etc.',
                 ' - Source: Name of the survey, so ANTR9, or DHBQ12, or YNTR6, etc.'],
        'outprefix': ['The prefix that should be used for output files. Nothing fancy. Files will always be saved as .csv files.',
                      'Output suffixes will always we _Fam_Correlations.csv and _Fam_N.csv.'],
        'pedigree': [
            'Specify a pedigree file to generate a reformatted pedigree. This should also be comma seperated file',
            'and should have the following columns:',
            ' - !FamID: Family ID',
            ' - ID: Subject ID (identical to FISNumber),',
            ' - Father: Subject ID of Father',
            ' - Mother: Subject ID of Mother',
            ' - Gender: Sex of the participant (M for male F for female)',
            ' - Twincode: 1 for MZ twins, empty for everything else',
            ' - DZtwincode: 1 for DZ twins, empty for everything else',
            ' - TwinHousehold3: Unique numbers for each twin household',
            ' - SibHousehold2: Unique numbers for each sib household',
            ' - SpouseHousehold3: Unique numbers for each spouse household',
            'To generate a reformatted extended pedigree use --pedigree with --extended',
            'If you have an extended pedigree file, you do not need a smaller pedigree file.'],
        'longitudinal': ['This script can do longitudinal selection.',
                         'WARNING: If you are going to do this I HIGHLY recommend including surveycompletion!'
                         'It will check within each family which list has',
                         'been completed the most, and use that survey for all family members.',
                         'If no perfect survey can be found the survey completed closest in time to the best survey',
                         'will be picked for those that do not have data on the optimal survey.',
                         'a file called prefix_familyselected.csv will be saved with results of the selection.',
                         'This file will automatically be used as input for generating correlation tables.',
                         'This will also creata JSON file (List_used_per_subject.json) detailing which surveys were used.'],
        'extended': ['Use the extended pedigree file. This will add  the following correlations to the output:',
                     'AuntNephew, AuntNiece, UnclueNephew, UncleNiece, NephewNephew, NieceNiece, NieceNephew',
                     'Note this requires an reformatted extended pedigree file.',
                     'To generate a reformatted extended pedigree file use --extended with --pedigree.'],
        'surveycompletion': ['An additional datafile with years of survey completion for longitudinal data.',
                             'This file should have columns FISNumber, index, Source as described in morehelp data',
                             'WARNING: If this file is not specified, but longitudinal is used, age will be used as a proxy',
                             ' This is really not ideal!'],
        'method': ['Method to be used for calculation of correlations. Pearson is the default method.',
                   'Alternatively Spearman rank (spearman) correlations can be calculated.',
                   'Kendall Tau (kendal) corelations can be calculated as long as randomsample or use_repeated_families is enabled.'],
        'correct': ['Formula-style string of corrections to be applied to every phenotype. Linear regression will',
                    'be performed using this formula on each phenotype individually. Some examples:',
                    ' --correction age; this will correct for age-effects',
                    ' --correction age+age**2; this will correct for age and squared age effects',
                    ' --correction age+sex+sex*age; this will correct for age, sex and an age*sex interaction.',
                    'Note: make sure you DO NOT use spaces in the formula.',
                    'Custom covariates can be used, just enter the column name as a covariate in the formula.',
                    'Do note however that any \'.\' in column names should be replaced with \'_\'!'
                    'Additionally any covariates added here should also be added to the exclude argument.'],
        'exclude': ['Comma seperated list (WITHOUT spaces) of variables for which no correlation is requested.',
                    'This should mainly be used when adding custom covariates to the correction argument.',
                    'Any custom covariates should be added here as the regression var1~age+var1 may result in errors.'],
        'randomsample': ['The default of this script is to calculated weighted correlations, and return the sum of weights',
                         'rather then the true N by default. If you do not want to weigh the correlations, but instead',
                         'would like to select one pair per family per correlation, add this argument.'
                         'Setting this will override raw_n.'],
        'use_repeated_families': ['Adding this argument will prevent weighting by occurences of samples during',
                                  'calculation of kinship correlations. This will yield a larger sample ',
                                  'but there is now no longer any correction for nested families and thus ',
                                  'results will probably bebiased.',
                                  'Setting this overrides raw_n'],
        'raw_n': ['The default of this script is to calculated weighted correlations, and return the sum of weights',
                  'rather then the true N by default. If you want an additional csv file (_Fam_raw_N.csv) detailing',
                  'the raw total number of samples per phenotype, add this argument.',
                  'This argument will have no effect if either use_repeated_families or randomsample is used.'],
        'other arguments': ['In the this .py file you can change some additional options:',
                            ' - lower_boundary: minimum age for selecting participants in longitudinal',
                            ' - upper_obundary: maximum age for selecting participants in longitudinal',
                            ' - check_cutoff_drops: save an Age-cutoff_drops.txt file detailing subjects',
                            '                       dropped by the cutoffs',
                            ' - seed: seed for random selection of subjects from larger families',
                            ' - explore_plot: generate scatterplots for each correlation',
                            ' - save_separate_data: Save the data used for each correlation in a file',
                            '                        this will generate 1 file per correlation like prefix_MZM.csv',
                            ' - parallel: Generate the reformatted pedigree using multiple processing threads.']
    }
    if args.morehelp is not None:
        morehelp(morehelpdict, args.morehelp)
    else:
        if (args.data is None) and (args.pedigree is None):
            morehelp(morehelpdict, ['this script', 'data', 'outprefix'])
            print('No output was requested. Please use at least (data and outprefix) or (pedigree) arguments, see help above.')
            quit()
        if (args.data is not None) and (args.outprefix is None):
            morehelp(morehelpdict, ['outprefix'])
            print('Please specify an output prefix. See help above')
            quit()
        datafile = args.data
        if (args.pedigree is None) and (not os.path.isfile('reformatted_pedigree.pickle')) and (not os.path.isfile('reformatted_extended_pedigree.pickle')):
            morehelp(morehelpdict, ['pedigree'])
            print('No pedigree file specified and no reformatted pedigree file found.')
            print('Please generate a reformatted pedigree first. See help on pedigree above')
            quit()
        if args.pedigree is not None:
            if (args.extended) and (os.path.isfile('Resources/reformatted_extended_pedigree.pickle')):
                print('Generation of a reformatted extended pedigree file was requested but one already exists.')
                print('Regenerating this reformatted file takes a LONG time!')
                print('If you are sure you want to regenerate this file please remove the old one first.')
                quit()
            if (not args.extended) and (os.path.isfile('reformatted_pedigree.pickle')):
                print('Generation of a reformatted pedigree file was requested but one already exists.')
                print('Regenerating this reformatted file takes a LONG time!')
                print('If you are sure you want to regenerate this file please remove the old one first.')
                quit()
            start = time.time()
            if parallel:
                print('The parallel implementation is bugged or something and not actually faster as of now.')
                print('For now it\'s disabled by default, so no need to cancel this.')
                parallel = False
            reformat_pedigree(args.pedigree, make_extended=args.extended, parallel=parallel)
            printtime(start, 'Reformatting pedigree')
        if args.data is not None:
            if not os.path.isfile(args.data):
                print('File {} not found.'.format(args.data))
                quit()
            if args.longitudinal:
                start = time.time()
                make_familybased_selection(datafile=datafile, pedigreefile=args.pedigree,
                                           outfileprefix=args.outprefix, upper_boundary=upper_boundary,
                                           lower_boundary=lower_boundary, surveycompletionfile=args.surveycompletion,
                                           check_cutoff_drops=check_cutoff_drops)
                printtime(start, 'Longitudinal selection')
                datafile = args.outprefix + '_familyselected.csv'
            start = time.time()
            results, resultsN = make_cor_table(datafile=datafile, seed=seed, explore_plot=explore_plot,
                           save_separate_data=save_separate_data, outfileprefix=args.outprefix,
                           use_repeated_families=args.use_repeated_families, method=args.method, correction=args.correct,
                           use_extended=args.extended, exclude=args.exclude)
            printtime(start, 'Generating correlation table')
